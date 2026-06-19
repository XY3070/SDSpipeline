#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

POP=""
CHR=""
ARRAY_QUEUE="normal"
SUBMITTER_QUEUE="smp"
SUBMITTER_SLOTS="4"
IN_ROOT="$BASE_DIR/data/processed/sds_input_rebuilt_main_contract_20260511"
OUT_ROOT="$BASE_DIR/data/processed/sds_output_gravel_chb_ne100k_newinput_20260511"
G_FILE="$BASE_DIR/tmp/gravel_chb_gamma_ne100k_20260502/gravel_chb_present100000.g_file.txt"
JOB_GROUP="via_smp"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pop) POP="$2"; shift 2 ;;
        --chr) CHR="$2"; shift 2 ;;
        --array-queue) ARRAY_QUEUE="$2"; shift 2 ;;
        --submitter-queue) SUBMITTER_QUEUE="$2"; shift 2 ;;
        --submitter-slots) SUBMITTER_SLOTS="$2"; shift 2 ;;
        --in-root) IN_ROOT="$2"; shift 2 ;;
        --out-root) OUT_ROOT="$2"; shift 2 ;;
        --g-file) G_FILE="$2"; shift 2 ;;
        --job-group) JOB_GROUP="$2"; shift 2 ;;
        *) echo "Unknown parameter: $1" >&2; exit 1 ;;
    esac
done

[[ -n "$POP" && -n "$CHR" ]] || {
    echo "Usage: $0 --pop POP --chr CHR [--array-queue normal] [--submitter-queue smp]" >&2
    exit 1
}

LOG_DIR="$OUT_ROOT/$POP/logs"
mkdir -p "$LOG_DIR"
WRAPPER_JOB="SDS_SUBMITTER_${POP}_chr${CHR}_${JOB_GROUP}"
INNER_CMD="bash \"$SCRIPT_DIR/submit_sds_compute_chunked_chr.sh\" --pop \"$POP\" --chr \"$CHR\" --queue \"$ARRAY_QUEUE\" --in-root \"$IN_ROOT\" --out-root \"$OUT_ROOT\" --g-file \"$G_FILE\" --job-group \"$JOB_GROUP\""

bsub -q "$SUBMITTER_QUEUE" -n "$SUBMITTER_SLOTS" -R "span[hosts=1]" -J "$WRAPPER_JOB" \
    -o "$LOG_DIR/${WRAPPER_JOB}.out" -e "$LOG_DIR/${WRAPPER_JOB}.err" \
    /bin/bash -lc "$INNER_CMD"
