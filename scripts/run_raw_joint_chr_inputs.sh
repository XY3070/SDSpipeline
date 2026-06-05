#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$SCRIPT_DIR/common_env.sh"
activate_sds_env

CHR=""
RAW_VCF=""
OUT_ROOT="$SDS_SDS_INPUT_ROOT"
TEST_MODE=0
FORCE=0
EMIT_AUDIT_SIDECARS=1
USE_HEADER_INTERSECTION=0
INTERSECTION_DIR="$SDS_RAW_HEADER_INTERSECTION_ROOT/raw_joint"

usage() {
    cat >&2 <<'EOF'
Usage: run_raw_joint_chr_inputs.sh --chr N --vcf PATH [--out-root DIR] [--test] [--force] [--skip-audit-sidecars] [--use-header-intersection] [--intersection-dir DIR]
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --chr) CHR="$2"; shift 2 ;;
        --vcf) RAW_VCF="$2"; shift 2 ;;
        --out-root) OUT_ROOT="$2"; shift 2 ;;
        --test) TEST_MODE=1; shift ;;
        --force) FORCE=1; shift ;;
        --skip-audit-sidecars) EMIT_AUDIT_SIDECARS=0; shift ;;
        --use-header-intersection) USE_HEADER_INTERSECTION=1; shift ;;
        --intersection-dir) INTERSECTION_DIR="$2"; shift 2 ;;
        *) usage ;;
    esac
done

[[ -n "$CHR" && -n "$RAW_VCF" ]] || usage
[[ -f "$RAW_VCF" ]] || { echo "[Error] raw VCF not found: $RAW_VCF" >&2; exit 1; }

declare -a COMMON_ARGS=(--chr "$CHR" --vcf "$RAW_VCF" --out-root "$OUT_ROOT")
if [[ "$TEST_MODE" -eq 1 ]]; then
    COMMON_ARGS+=(--test)
fi
if [[ "$FORCE" -eq 1 ]]; then
    COMMON_ARGS+=(--force)
fi
if [[ "$EMIT_AUDIT_SIDECARS" -eq 0 ]]; then
    COMMON_ARGS+=(--skip-audit-sidecars)
fi

run_one() {
    local out_pop="$1"
    local sample_list="$2"
    local effective_sample_list="$sample_list"

    if [[ "$USE_HEADER_INTERSECTION" -eq 1 ]]; then
        effective_sample_list="$INTERSECTION_DIR/${out_pop}.raw_chr${CHR}_intersection.txt"
        bash "$SCRIPT_DIR/build_raw_header_intersection_sample_list.sh" \
            --vcf "$RAW_VCF" \
            --sample-list "$sample_list" \
            --output "$effective_sample_list" >/dev/null
    fi

    echo ">>> raw-joint input build: ${out_pop} chr${CHR}"
    bash "$SCRIPT_DIR/run_sds_input.sh" \
        "${COMMON_ARGS[@]}" \
        --sample-list "$effective_sample_list" \
        --out-pop "$out_pop"
}

run_one NCN "$(find_population_sample_list NCN)"
run_one SCN "$(find_population_sample_list SCN)"
run_one MERGED3971 "$(find_population_sample_list MERGED3971)"

echo "✓ DONE: raw-joint chr${CHR} inputs for NCN / SCN / MERGED3971"
