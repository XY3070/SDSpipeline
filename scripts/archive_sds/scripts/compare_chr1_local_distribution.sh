#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck source=/share/home/grp-wangyf/xuyuan/sds/scripts/common_env.sh
source "$SCRIPT_DIR/common_env.sh"
activate_sds_env
PYTHON_BIN="$SDS_ENV_PREFIX/bin/python"

BASELINE_RAW=""
COMPARE_RAW=""
OUT_DIR=""
BASELINE_LABEL="baseline"
COMPARE_LABEL="compare"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --baseline-raw) BASELINE_RAW="$2"; shift 2 ;;
        --compare-raw) COMPARE_RAW="$2"; shift 2 ;;
        --out-dir) OUT_DIR="$2"; shift 2 ;;
        --baseline-label) BASELINE_LABEL="$2"; shift 2 ;;
        --compare-label) COMPARE_LABEL="$2"; shift 2 ;;
        *)
            echo "Unknown parameter: $1" >&2
            exit 1
            ;;
    esac
done

if [[ -z "$BASELINE_RAW" || -z "$COMPARE_RAW" || -z "$OUT_DIR" ]]; then
    echo "Usage: $0 --baseline-raw PATH --compare-raw PATH --out-dir DIR [--baseline-label LABEL] [--compare-label LABEL]" >&2
    exit 1
fi

[[ -f "$BASELINE_RAW" ]] || { echo "[Error] Missing baseline raw file: $BASELINE_RAW" >&2; exit 1; }
[[ -f "$COMPARE_RAW" ]] || { echo "[Error] Missing compare raw file: $COMPARE_RAW" >&2; exit 1; }
mkdir -p "$OUT_DIR"

run_one() {
    local raw_path="$1"
    local label="$2"
    local prefix="$OUT_DIR/$label"

    "$PYTHON_BIN" "$SCRIPT_DIR/normalize_chr1_local_sds.py" \
        --input "$raw_path" \
        --output-prefix "$prefix"

    "$PYTHON_BIN" "$SCRIPT_DIR/report_tail_normality.py" \
        --input-normalized-tsv "${prefix}.normalized.tsv" \
        --output-prefix "${prefix}_tail"
}

run_one "$BASELINE_RAW" "$BASELINE_LABEL"
run_one "$COMPARE_RAW" "$COMPARE_LABEL"

"$PYTHON_BIN" "$SCRIPT_DIR/compare_sds_frequency_bins.py" \
    --label-a "$BASELINE_LABEL" \
    --label-b "$COMPARE_LABEL" \
    --bins-a "$OUT_DIR/${BASELINE_LABEL}.frequency_bins.tsv" \
    --bins-b "$OUT_DIR/${COMPARE_LABEL}.frequency_bins.tsv" \
    --output-tsv "$OUT_DIR/${BASELINE_LABEL}_vs_${COMPARE_LABEL}.frequency_bins_compare.tsv"

cat <<EOF
baseline_normalized	$OUT_DIR/${BASELINE_LABEL}.normalized.tsv
baseline_bins	$OUT_DIR/${BASELINE_LABEL}.frequency_bins.tsv
baseline_tail	$OUT_DIR/${BASELINE_LABEL}_tail.tsv
compare_normalized	$OUT_DIR/${COMPARE_LABEL}.normalized.tsv
compare_bins	$OUT_DIR/${COMPARE_LABEL}.frequency_bins.tsv
compare_tail	$OUT_DIR/${COMPARE_LABEL}_tail.tsv
bins_compare	$OUT_DIR/${BASELINE_LABEL}_vs_${COMPARE_LABEL}.frequency_bins_compare.tsv
EOF
