#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SDS_PYTHON="/data/home/grp-wangyf/intern/miniforge3/envs/sds/bin/python"

OLD_ROOT="$BASE_DIR/data/processed/sds_output_gravel_chb_ne100k_newinput_20260511/SCN"
MATCHED_ROOT="$BASE_DIR/data/processed/sds_output_gravel_chb_ne100k_scn5432matched_genomewide_run2_20260526/SCN"
OLDDEFAULT_ROOT="$BASE_DIR/data/processed/sds_output_scn_olddefault_genomewide_20260527/SCN"
OUT_DIR="$BASE_DIR/data/processed/scn_olddefault_genomewide_compare_20260528"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --old-root) OLD_ROOT="$2"; shift 2 ;;
        --matched-root) MATCHED_ROOT="$2"; shift 2 ;;
        --olddefault-root) OLDDEFAULT_ROOT="$2"; shift 2 ;;
        --out-dir) OUT_DIR="$2"; shift 2 ;;
        *) echo "Unknown parameter: $1" >&2; exit 1 ;;
    esac
done

mkdir -p "$OUT_DIR"

"$SDS_PYTHON" "$SCRIPT_DIR/compare_sds_diagnostics.py" \
  --track "old=$OLD_ROOT/diagnostics/SCN.gravel_chb_ne100k_manual" \
  --track "matched=$MATCHED_ROOT/diagnostics/SCN.gravel_chb_ne100k_submitter_array" \
  --track "olddefault=$OLDDEFAULT_ROOT/diagnostics/SCN.gravel_chb_ne100k_submitter_array" \
  --output-prefix "$OUT_DIR/SCN_old_vs_matched_vs_olddefault"

"$SDS_PYTHON" "$SCRIPT_DIR/compare_sds_frequency_bins.py" \
  --label-a old \
  --label-b olddefault \
  --bins-a "$OLD_ROOT/SCN.frequency_bins.tsv" \
  --bins-b "$OLDDEFAULT_ROOT/SCN.frequency_bins.tsv" \
  --output-tsv "$OUT_DIR/SCN_old_vs_olddefault.frequency_bins_compare.tsv"

"$SDS_PYTHON" "$SCRIPT_DIR/compare_sds_frequency_bins.py" \
  --label-a matched \
  --label-b olddefault \
  --bins-a "$MATCHED_ROOT/SCN.frequency_bins.tsv" \
  --bins-b "$OLDDEFAULT_ROOT/SCN.frequency_bins.tsv" \
  --output-tsv "$OUT_DIR/SCN_matched_vs_olddefault.frequency_bins_compare.tsv"

printf 'compare_dir\t%s\n' "$OUT_DIR"
