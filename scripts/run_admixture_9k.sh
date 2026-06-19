#!/usr/bin/env bash
# ============================================================================
# run_admixture_9k.sh
#
# Reusable entrypoint for ADAMIXTURE-based population stratification on the
# 9k Chinese cohort. Implements the workflow proposed in:
#
#     SDSlog/logs/demography/2026-06-19_admixture_proposal_for_9k.md
#
# Contract:
#   - All server-local paths come from SDSpipeline/config/paths.env
#     (sourced automatically if SDS_PIPELINE_ROOT is unset).
#   - The Python environment is uv-managed (see companion
#     setup_admixture_env.sh). We do NOT depend on conda/mamba.
#   - Each step is idempotent via .done sentinels; safe to re-run.
#   - Outputs land under $SDS_ADMIXTURE_RUN_ROOT.
#
# Usage:
#   ./run_admixture_9k.sh {all|step1|step2|step3|step4|step5|step6 [K]}
# ============================================================================
set -euo pipefail

# ----------------------------------------------------------------------------
# Load paths.env (server-local config) if not already sourced.
# ----------------------------------------------------------------------------
if [[ -z "${SDS_PIPELINE_ROOT:-}" ]]; then
  SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  # SELF_DIR is SDSpipeline/scripts/
  PATHS_ENV="${SELF_DIR}/../config/paths.env"
  if [[ ! -f "$PATHS_ENV" ]]; then
    echo "FATAL: cannot find config/paths.env next to this script." >&2
    exit 2
  fi
  # shellcheck disable=SC1090
  source "$PATHS_ENV"
fi

# ----------------------------------------------------------------------------
# Required paths (fail fast with a clear message if any is missing).
# ----------------------------------------------------------------------------
: "${SDS_COHORT_9K_ROOT:?set SDS_COHORT_9K_ROOT in config/paths.env}"
: "${SDS_COHORT_9K_FREEZE:?set SDS_COHORT_9K_FREEZE in config/paths.env}"
: "${SDS_ADMIXTURE_RUN_ROOT:?set SDS_ADMIXTURE_RUN_ROOT in config/paths.env}"
: "${SDS_ADMIXTURE_ENV_PREFIX:?set SDS_ADMIXTURE_ENV_PREFIX in config/paths.env}"
: "${SDS_PLINK2_BIN:?set SDS_PLINK2_BIN in config/paths.env}"
: "${SDS_BCFTOOLS_BIN:?set SDS_BCFTOOLS_BIN in config/paths.env}"

for bin in "$SDS_PLINK2_BIN" "$SDS_BCFTOOLS_BIN"; do
  [[ -x "$bin" ]] || { echo "FATAL: not executable: $bin" >&2; exit 2; }
done

PLINK2="$SDS_PLINK2_BIN"
BCFTOOLS="$SDS_BCFTOOLS_BIN"
ADAMIXTURE_BIN="${SDS_ADMIXTURE_ENV_PREFIX}/bin/adamixture"
ADAMIXTURE_PROJECT_BIN="${SDS_ADMIXTURE_ENV_PREFIX}/bin/adamixture-project"
ADAMIXTURE_SUPERVISED_BIN="${SDS_ADMIXTURE_ENV_PREFIX}/bin/adamixture-supervised"

RUN_ROOT="$SDS_ADMIXTURE_RUN_ROOT"
WORK_DIR="${RUN_ROOT}/work"
INPUT_DIR="${RUN_ROOT}/input"
OUT_UNSUP="${RUN_ROOT}/unsupervised"
OUT_SUP="${RUN_ROOT}/supervised"
LOG_DIR="${RUN_ROOT}/logs"
mkdir -p "$WORK_DIR" "$INPUT_DIR" "$OUT_UNSUP" "$OUT_SUP" "$LOG_DIR"

# ----------------------------------------------------------------------------
# Plan-A default: 11 ready autosomes (chr6 excluded; rerun after chr6 is QC'd).
# Override via env, e.g.:  CHRS_READY=(1 4 5 6 7 14 17 18 19 20 21 22) ./run_admixture_9k.sh
# ----------------------------------------------------------------------------
if [[ -z "${CHRS_READY+x}" ]]; then
  CHRS_READY=(4 5 7 14 17 18 19 20 21 22)
fi

LD_WINDOW=${LD_WINDOW:-1000}
LD_STEP=${LD_STEP:-50}
LD_R2=${LD_R2:-0.1}

MIN_K=${MIN_K:-2}
MAX_K=${MAX_K:-10}
CV_FOLDS=${CV_FOLDS:-10}
THREADS=${THREADS:-32}

# ADAMIXTURE Adam-optimizer tuning (recommended by the upstream
# troubleshooting doc for large cohorts / elevated K).
EXTRA_ADAM_ARGS=(
  --patience_adam "${PATIENCE_ADAM:-5}"
  --lr_decay "${LR_DECAY:-0.85}"
  --lr "${LR:-0.0075}"
)

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
log() { printf '[%s] %s\n' "$(date +%F\ %T)" "$*"; }

need_done() {
  if [[ -f "$2" ]]; then log "skip $1 (sentinel $2 exists)"; return 0; fi
  return 1
}

resolve_qc_vcf() {
  local c=$1 d
  # Prefer snp_singleton_qc_chr<CHROM>_v* if present.
  for d in "${SDS_COHORT_9K_ROOT}"/snp_singleton_qc_chr${c}_v*; do
    if [[ -f "${d}/snp_qc_pass.vcf.gz" ]]; then
      echo "${d}/snp_qc_pass.vcf.gz"
      return 0
    fi
  done
  # Fallback: top-level all.chr<N>.QC.vcf.gz (e.g. chr1, chr6).
  local top="${SDS_COHORT_9K_ROOT}/all.chr${c}.QC.vcf.gz"
  if [[ -f "$top" ]]; then echo "$top"; return 0; fi
  echo "ERROR: no QC VCF for chr${c} under ${SDS_COHORT_9K_ROOT}" >&2
  return 1
}

freeze_sample_list() {
  # Emit one-sample-per-line list (no header) from cohort_freeze.tsv.
  awk -F'\t' 'NR>1 && $1!="" {print $1}' "$SDS_COHORT_9K_FREEZE"
}

# ----------------------------------------------------------------------------
# STEP 1 — per-chromosome biallelic-SNP pgen restricted to frozen samples
#   Uses plink2 (streaming, binary-packed) instead of bcftools because the
#   chr4 VCF (2.3GB compressed, 3.2M variants x 4393 samples) OOMed bcftools
#   at ~62GB resident. Cap RAM at PLINK_MEM_MB (default 32000).
# ----------------------------------------------------------------------------
PLINK_MEM_MB=${PLINK_MEM_MB:-32000}

step1_per_chr() {
  log "STEP 1: per-chromosome pgen (plink2, --memory ${PLINK_MEM_MB}MB)"
  local c vcf out_prefix samples_file
  samples_file="${WORK_DIR}/_freeze_samples.txt"
  freeze_sample_list > "$samples_file"

  for c in "${CHRS_READY[@]}"; do
    vcf=$(resolve_qc_vcf "$c")
    out_prefix="${WORK_DIR}/chr${c}.biallelic"
    if need_done "chr${c}" "${out_prefix}.done"; then continue; fi
    "$PLINK2" --vcf "$vcf" \
              --memory "$PLINK_MEM_MB" \
              --threads "$THREADS" \
              --keep "$samples_file" \
              --snps-only \
              --max-alleles 2 \
              --maf 0.01 \
              --make-pgen --out "$out_prefix"
    touch "${out_prefix}.done"
    log "  chr${c}: $(wc -l < "${out_prefix}.pvar") variants (incl. header)"
  done
}

# ----------------------------------------------------------------------------
# STEP 2 — pmerge-list across chromosomes into one genome-wide pgen
# ----------------------------------------------------------------------------
step2_concat() {
  log "STEP 2: pmerge-list across chromosomes (plink2)"
  local out_prefix="${WORK_DIR}/9k_raw"
  local list_file="${WORK_DIR}/chr_list_for_pmerge.txt"
  need_done "concat" "${out_prefix}.done" && return 0

  : > "$list_file"
  local c
  for c in "${CHRS_READY[@]}"; do
    echo "${WORK_DIR}/chr${c}.biallelic" >> "$list_file"
  done

  "$PLINK2" --memory "$PLINK_MEM_MB" \
            --threads "$THREADS" \
            --pmerge-list "$list_file" \
            --make-pgen --out "$out_prefix"
  touch "${out_prefix}.done"
  log "  merged variants: $(wc -l < "${out_prefix}.pvar") (incl. header)"
}

# ----------------------------------------------------------------------------
# STEP 3 — pgen + LD pruning
# ----------------------------------------------------------------------------
step3_pgen_and_prune() {
  log "STEP 3: pgen + LD prune"
  local pgen_prefix="${WORK_DIR}/9k_raw"
  local prune_in="${WORK_DIR}/9k.prune.in"
  local final_prefix="${INPUT_DIR}/9k_ldpruned"

  # pgen already produced by step2 (9k_raw.{pgen,pvar,psam} + .done sentinel).
  # Skip the legacy VCF-to-pgen branch.

  if need_done "prune" "${prune_in}.done"; then :; else
    "$PLINK2" --pfile "$pgen_prefix" \
              --indep-pairwise "$LD_WINDOW" "$LD_STEP" "$LD_R2" \
              --out "${WORK_DIR}/9k" \
              --threads "$THREADS" \
              --memory "$PLINK_MEM_MB"
    touch "${prune_in}.done"
  fi
  log "  LD-pruned SNPs: $(wc -l < "$prune_in")"

  if need_done "final" "${final_prefix}.done"; then :; else
    "$PLINK2" --pfile "$pgen_prefix" \
              --extract "$prune_in" \
              --make-pgen --out "$final_prefix" \
              --threads "$THREADS" \
              --memory "$PLINK_MEM_MB"
    touch "${final_prefix}.done"
  fi
}

# ----------------------------------------------------------------------------
# STEP 4 — build labels.txt aligned to psam order (Region + Superpopulation)
# ----------------------------------------------------------------------------
step4_labels() {
  log "STEP 4: build labels.txt (Region) and labels_superpop.txt"
  local out_region="${INPUT_DIR}/9k_ldpruned.labels.txt"
  local out_super="${INPUT_DIR}/9k_ldpruned.labels_superpop.txt"
  need_done "labels" "${out_region}.done" && return 0

  "$SDS_ADMIXTURE_ENV_PREFIX/bin/python" - "$SDS_COHORT_9K_FREEZE" \
      "${INPUT_DIR}/9k_ldpruned.psam" "$out_region" "$out_super" <<'PY'
import sys, csv
meta_path, psam_path, out_region, out_super = sys.argv[1:]

meta = {}
with open(meta_path) as f:
    rd = csv.DictReader(f, delimiter='\t')
    for row in rd:
        s = row['Sample']
        reg = (row.get('Region') or '-').strip()
        sup = (row.get('Superpopulation') or '-').strip()
        if reg in ('', '0', 'China', 'EastAsia', 'East Asian Ancestry'):
            reg = '-'
        meta[s] = (reg, sup)

with open(psam_path) as f:
    rd = csv.DictReader(f, delimiter='\t')
    samples = [row.get('IID') or row.get('#FID') for row in rd]

with open(out_region, 'w') as fo_r, open(out_super, 'w') as fo_s:
    for s in samples:
        reg, sup = meta.get(s, ('-', '-'))
        fo_r.write(reg + '\n')
        fo_s.write(sup + '\n')
PY
  touch "${out_region}.done"
  log "  label distribution (Region):"
  sort "$out_region" | uniq -c | sort -rn | sed 's/^/    /'
}

# ----------------------------------------------------------------------------
# STEP 5 — unsupervised K-sweep with CV
# ----------------------------------------------------------------------------
step5_unsupervised() {
  log "STEP 5: ADAMIXTURE unsupervised K=${MIN_K}..${MAX_K}, CV=${CV_FOLDS}-fold"
  local sentinel="${OUT_UNSUP}/sweep.done"
  need_done "unsupervised sweep" "$sentinel" && return 0

  [[ -x "$ADAMIXTURE_BIN" ]] || {
    echo "FATAL: $ADAMIXTURE_BIN not found. Run setup_admixture_env.sh first." >&2
    exit 2
  }

  "$ADAMIXTURE_BIN" \
    --min_k "$MIN_K" --max_k "$MAX_K" \
    --cv "$CV_FOLDS" \
    --data_path "${INPUT_DIR}/9k_ldpruned.pgen" \
    --save_dir "$OUT_UNSUP" \
    --name 9k \
    -t "$THREADS" \
    --labels "${INPUT_DIR}/9k_ldpruned.labels.txt" \
    --plot --plot_single \
    --max_iter 10000 --tol_adam 0.1 \
    "${EXTRA_ADAM_ARGS[@]}" \
    2>&1 | tee "${LOG_DIR}/adamixture_unsupervised.log"

  touch "$sentinel"
}

# ----------------------------------------------------------------------------
# STEP 6 — (optional) supervised run against a chosen K
# ----------------------------------------------------------------------------
step6_supervised() {
  local K=${1:-6}
  log "STEP 6: ADAMIXTURE supervised k=${K}"
  local sentinel="${OUT_SUP}/sup_k${K}.done"
  need_done "supervised k=${K}" "$sentinel" && return 0

  [[ -x "$ADAMIXTURE_SUPERVISED_BIN" ]] || {
    echo "FATAL: $ADAMIXTURE_SUPERVISED_BIN not found." >&2
    exit 2
  }

  "$ADAMIXTURE_SUPERVISED_BIN" \
    --data_path "${INPUT_DIR}/9k_ldpruned.pgen" \
    --labels "${INPUT_DIR}/9k_ldpruned.labels.txt" \
    --save_dir "$OUT_SUP" \
    --name 9k_sup \
    -k "$K" \
    -t "$THREADS" --plot \
    "${EXTRA_ADAM_ARGS[@]}" \
    2>&1 | tee "${LOG_DIR}/adamixture_supervised_k${K}.log"

  touch "$sentinel"
}

# ----------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------
case "${1:-all}" in
  step1) step1_per_chr ;;
  step2) step2_concat ;;
  step3) step3_pgen_and_prune ;;
  step4) step4_labels ;;
  step5) step5_unsupervised ;;
  step6) step6_supervised "${2:-6}" ;;
  all)
    step1_per_chr
    step2_concat
    step3_pgen_and_prune
    step4_labels
    step5_unsupervised
    ;;
  *)
    echo "usage: $0 {all|step1|step2|step3|step4|step5|step6 [K]}" >&2
    exit 2
    ;;
esac
