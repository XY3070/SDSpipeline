#!/bin/bash
set -euo pipefail

if [[ $# -lt 8 ]]; then
    echo "Usage: $0 MANIFEST_TSV S_FILE O_FILE B_FILE G_FILE INIT S_FILE_NCOL CACHE_DIR" >&2
    exit 1
fi

MANIFEST_TSV="$1"
INDEX="${LSB_JOBINDEX:-}"
S_FILE="$2"
O_FILE="$3"
B_FILE="$4"
G_FILE="$5"
INIT="$6"
S_FILE_NCOL="$7"
CACHE_DIR="$8"
SKIP_BOUNDARY_MISSING_FRACTION="${SDS_SKIP_BOUNDARY_MISSING_FRACTION:-}"
BOUNDARY_MISSING_MODE="${SDS_BOUNDARY_MISSING_MODE:-}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[[ -n "$INDEX" ]] || { echo "[Error] LSB_JOBINDEX is not set" >&2; exit 1; }
line="$(sed -n "${INDEX}p" "$MANIFEST_TSV")"
[[ -n "$line" ]] || { echo "[Error] No manifest entry at index $INDEX" >&2; exit 1; }
IFS=$'\t' read -r chunk_file out_tsv out_parquet out_summary <<< "$line"

bash "$SCRIPT_DIR/run_sds_compute_chunk.sh" \
  "$chunk_file" "$out_tsv" "$S_FILE" "$O_FILE" "$B_FILE" "$G_FILE" "$INIT" "$S_FILE_NCOL" "$out_parquet" "$out_summary" "$CACHE_DIR" \
  "${SKIP_BOUNDARY_MISSING_FRACTION}" "${BOUNDARY_MISSING_MODE}"
