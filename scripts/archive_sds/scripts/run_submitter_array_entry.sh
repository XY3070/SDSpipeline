#!/bin/bash
set -euo pipefail

if [[ $# -lt 7 || $# -gt 11 ]]; then
    echo "Usage: $0 POP CHR_LIST_FILE IN_ROOT OUT_ROOT G_FILE JOB_GROUP_PREFIX RESULT_DIR [QUEUE] [CHUNK_ROWS] [ARRAY_PARALLEL] [CHUNK_JOB_SLOTS]" >&2
    exit 1
fi

POP="$1"
CHR_LIST_FILE="$2"
IN_ROOT="$3"
OUT_ROOT="$4"
G_FILE="$5"
JOB_GROUP_PREFIX="$6"
RESULT_DIR="$7"
QUEUE="${8:-normal}"
CHUNK_ROWS="${9:-}"
ARRAY_PARALLEL="${10:-}"
CHUNK_JOB_SLOTS="${11:-}"
SKIP_BOUNDARY_MISSING_FRACTION="${SDS_SKIP_BOUNDARY_MISSING_FRACTION:-}"
BOUNDARY_MISSING_MODE="${SDS_BOUNDARY_MISSING_MODE:-}"

CHR="$(sed -n "${LSB_JOBINDEX}p" "$CHR_LIST_FILE")"
[[ -n "$CHR" ]] || { echo "[Error] No chromosome for LSB_JOBINDEX=${LSB_JOBINDEX}" >&2; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$RESULT_DIR"

OUT_FILE="$RESULT_DIR/chr${CHR}.submitter.out"
FINAL_FILE="$RESULT_DIR/chr${CHR}.final_job.tsv"

cmd=(
    bash "$SCRIPT_DIR/submit_sds_compute_chunked_chr.sh"
    --pop "$POP"
    --chr "$CHR"
    --queue "$QUEUE"
    --in-root "$IN_ROOT"
    --out-root "$OUT_ROOT"
    --g-file "$G_FILE"
    --job-group "${JOB_GROUP_PREFIX}_chr${CHR}"
)

if [[ -n "$CHUNK_ROWS" ]]; then
    cmd+=(--chunk-rows "$CHUNK_ROWS")
fi
if [[ -n "$ARRAY_PARALLEL" ]]; then
    cmd+=(--array-parallel "$ARRAY_PARALLEL")
fi
if [[ -n "$CHUNK_JOB_SLOTS" ]]; then
    cmd+=(--chunk-job-slots "$CHUNK_JOB_SLOTS")
fi
if [[ -n "$SKIP_BOUNDARY_MISSING_FRACTION" ]]; then
    cmd+=(--skip-boundary-missing-fraction "$SKIP_BOUNDARY_MISSING_FRACTION")
fi
if [[ -n "$BOUNDARY_MISSING_MODE" ]]; then
    cmd+=(--boundary-missing-mode "$BOUNDARY_MISSING_MODE")
fi

"${cmd[@]}" | tee "$OUT_FILE"

FINAL_JOB="$(awk -F '\t' '$1=="FINAL_JOB"{print $2; exit}' "$OUT_FILE")"
[[ -n "$FINAL_JOB" ]] || { echo "[Error] Failed to parse FINAL_JOB for chr${CHR}" >&2; exit 1; }
printf 'chr\tfinal_job\n%s\t%s\n' "$CHR" "$FINAL_JOB" > "$FINAL_FILE"
