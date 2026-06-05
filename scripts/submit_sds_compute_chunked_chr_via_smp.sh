#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$SCRIPT_DIR/common_env.sh"

POP=""
CHR=""
ARRAY_QUEUE="normal"
SUBMITTER_QUEUE="smp"
SUBMITTER_SLOTS="4"
IN_ROOT="$SDS_SDS_INPUT_ROOT"
OUT_ROOT="$SDS_SDS_OUTPUT_ROOT"
G_FILE=""
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

if [[ -z "$G_FILE" ]]; then
    G_FILE="$(find_default_g_file "$BASE_DIR" "$POP" || true)"
fi
[[ -n "$G_FILE" && -f "$G_FILE" ]] || { echo "[Error] g_file not found: $G_FILE" >&2; exit 1; }

LOG_DIR="$OUT_ROOT/$POP/logs"
mkdir -p "$LOG_DIR"
WRAPPER_JOB="SDS_SUBMITTER_${POP}_chr${CHR}_${JOB_GROUP}"
INNER_CMD="bash \"$SCRIPT_DIR/submit_sds_compute_chunked_chr.sh\" --pop \"$POP\" --chr \"$CHR\" --queue \"$ARRAY_QUEUE\" --in-root \"$IN_ROOT\" --out-root \"$OUT_ROOT\" --g-file \"$G_FILE\" --job-group \"$JOB_GROUP\""

bsub -q "$SUBMITTER_QUEUE" -n "$SUBMITTER_SLOTS" -R "span[hosts=1]" -J "$WRAPPER_JOB" \
    -o "$LOG_DIR/${WRAPPER_JOB}.out" -e "$LOG_DIR/${WRAPPER_JOB}.err" \
    /bin/bash -lc "$INNER_CMD"
