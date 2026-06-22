#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$SCRIPT_DIR/common_env.sh"
activate_sds_env

POP=""
CHR=""
OUT_ROOT="$SDS_SDS_OUTPUT_ROOT"
IN_ROOT="$SDS_SDS_INPUT_ROOT"
QUEUE="normal"
CHUNK_ROWS="5000"
ARRAY_PARALLEL="64"
G_FILE=""
INIT="0.00001"
S_FILE_NCOL="20000"
JOB_GROUP="track"
OUTPUT_PREFIX=""
declare -a WINDOWS=()

usage() {
    cat >&2 <<'EOF'
Usage: submit_chr_sds_track_with_review.sh --pop POP --chr N --job-group TAG --g-file PATH --output-prefix PREFIX [--out-root DIR] [--in-root DIR] [--queue QUEUE] [--chunk-rows N] [--array-parallel N] --window NAME=START-END [--window ...]
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pop) POP="$2"; shift 2 ;;
        --chr) CHR="$2"; shift 2 ;;
        --out-root) OUT_ROOT="$2"; shift 2 ;;
        --in-root) IN_ROOT="$2"; shift 2 ;;
        --queue) QUEUE="$2"; shift 2 ;;
        --chunk-rows) CHUNK_ROWS="$2"; shift 2 ;;
        --array-parallel) ARRAY_PARALLEL="$2"; shift 2 ;;
        --g-file) G_FILE="$2"; shift 2 ;;
        --init) INIT="$2"; shift 2 ;;
        --s-file-ncol) S_FILE_NCOL="$2"; shift 2 ;;
        --job-group) JOB_GROUP="$2"; shift 2 ;;
        --output-prefix) OUTPUT_PREFIX="$2"; shift 2 ;;
        --window) WINDOWS+=("$2"); shift 2 ;;
        *) usage ;;
    esac
done

[[ -n "$POP" && -n "$CHR" && -n "$G_FILE" && -n "$OUTPUT_PREFIX" ]] || usage
[[ "${#WINDOWS[@]}" -gt 0 ]] || usage
[[ -f "$G_FILE" ]] || { echo "[Error] g_file not found: $G_FILE" >&2; exit 1; }

SUBMIT_OUTPUT="$("$SCRIPT_DIR/submit_sds_compute_chunked_chr.sh" \
    --pop "$POP" \
    --chr "$CHR" \
    --out-root "$OUT_ROOT" \
    --in-root "$IN_ROOT" \
    --queue "$QUEUE" \
    --chunk-rows "$CHUNK_ROWS" \
    --array-parallel "$ARRAY_PARALLEL" \
    --g-file "$G_FILE" \
    --init "$INIT" \
    --s-file-ncol "$S_FILE_NCOL" \
    --job-group "$JOB_GROUP")"

printf '%s\n' "$SUBMIT_OUTPUT"

FINAL_JOB_ID="$(awk -F '\t' '$1=="FINAL_JOB"{print $2}' <<< "$SUBMIT_OUTPUT")"
[[ -n "$FINAL_JOB_ID" ]] || { echo "[Error] FINAL_JOB missing from submit output" >&2; exit 1; }

REVIEW_DIR="$OUT_ROOT/$POP/review"
LOGDIR="$OUT_ROOT/$POP/logs"
mkdir -p "$REVIEW_DIR" "$LOGDIR"

REVIEW_CMD=("$SDS_ENV_PREFIX/bin/python" "$SCRIPT_DIR/summarize_sds_windows.py" \
    --input "$OUT_ROOT/$POP/chr${CHR}.sds.tsv" \
    --output-prefix "$REVIEW_DIR/$OUTPUT_PREFIX")
for window_spec in "${WINDOWS[@]}"; do
    REVIEW_CMD+=(--window "$window_spec")
done

printf -v REVIEW_CMD_STR '%q ' "${REVIEW_CMD[@]}"
REVIEW_JOB_NAME="SDS_${POP}_chr${CHR}_review_${JOB_GROUP}"
REVIEW_SUBMIT_OUT="$(bsub -q "$QUEUE" -w "done(${FINAL_JOB_ID})" -n 1 -J "$REVIEW_JOB_NAME" -o "$LOGDIR/${REVIEW_JOB_NAME}.out" -e "$LOGDIR/${REVIEW_JOB_NAME}.err" /bin/bash -lc "$REVIEW_CMD_STR" < /dev/null)"
REVIEW_JOB_ID="$(sed -n 's/.*<\([0-9]\+\)>.*/\1/p' <<< "$REVIEW_SUBMIT_OUT")"
[[ -n "$REVIEW_JOB_ID" ]] || { echo "[Error] Failed to parse review bsub output: $REVIEW_SUBMIT_OUT" >&2; exit 1; }

printf 'REVIEW_JOB\t%s\n' "$REVIEW_JOB_ID"
