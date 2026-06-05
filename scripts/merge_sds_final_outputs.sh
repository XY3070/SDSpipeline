#!/bin/bash
set -euo pipefail

if [[ $# -lt 4 ]]; then
    echo "Usage: $0 FINAL_TSV FINAL_PARQUET SUMMARY_CSV PART1 [PART2 ...]" >&2
    exit 1
fi

FINAL_TSV="$1"
FINAL_PARQUET="$2"
SUMMARY_CSV="$3"
shift 3

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/common_env.sh"
activate_sds_env
PYTHON_BIN="$SDS_ENV_PREFIX/bin/python"

mkdir -p "$(dirname "$FINAL_TSV")" "$(dirname "$FINAL_PARQUET")" "$(dirname "$SUMMARY_CSV")"

header_written=0
: > "$FINAL_TSV"
for part_file in "$@"; do
    [[ -s "$part_file" ]] || continue
    if [[ "$header_written" -eq 0 ]]; then
        cat "$part_file" >> "$FINAL_TSV"
        header_written=1
    else
        tail -n +2 "$part_file" >> "$FINAL_TSV"
    fi
done

if [[ "$header_written" -eq 0 ]]; then
    echo "[Error] No arm outputs were available to merge" >&2
    exit 1
fi

"$PYTHON_BIN" "$SCRIPT_DIR/archive_sds_output.py" \
    "$FINAL_TSV" \
    "$FINAL_PARQUET" \
    --summary-csv "$SUMMARY_CSV"
