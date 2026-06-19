#!/bin/bash
set -euo pipefail

if [[ $# -lt 4 ]]; then
    echo "Usage: $0 MANIFEST_TSV MERGED_TSV MERGED_PARQUET SUMMARY_CSV" >&2
    exit 1
fi

MANIFEST_TSV="$1"
MERGED_TSV="$2"
MERGED_PARQUET="$3"
SUMMARY_CSV="$4"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/common_env.sh"
activate_sds_env
PYTHON_BIN="$SDS_ENV_PREFIX/bin/python"

mkdir -p "$(dirname "$MERGED_TSV")" "$(dirname "$MERGED_PARQUET")" "$(dirname "$SUMMARY_CSV")"

header_written=0
: > "$MERGED_TSV"

while IFS=$'\t' read -r chunk_file out_tsv out_parquet out_summary; do
    [[ -n "$out_tsv" ]] || continue
    [[ -s "$out_tsv" ]] || continue
    if [[ "$header_written" -eq 0 ]]; then
        cat "$out_tsv" >> "$MERGED_TSV"
        header_written=1
    else
        tail -n +2 "$out_tsv" >> "$MERGED_TSV"
    fi
done < "$MANIFEST_TSV"

if [[ "$header_written" -eq 0 ]]; then
    echo "[Error] No chunk outputs were available to merge from $MANIFEST_TSV" >&2
    exit 1
fi

"$PYTHON_BIN" "$SCRIPT_DIR/archive_sds_output.py" \
    "$MERGED_TSV" \
    "$MERGED_PARQUET" \
    --summary-csv "$SUMMARY_CSV"
