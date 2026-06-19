#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 --daf DAF --chunk-root DIR --gamma-prefix PATH" >&2
    exit 1
}

DAF=""
CHUNK_ROOT=""
GAMMA_PREFIX=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --daf) DAF="$2"; shift 2 ;;
        --chunk-root) CHUNK_ROOT="$2"; shift 2 ;;
        --gamma-prefix) GAMMA_PREFIX="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; usage ;;
    esac
done

[[ -n "$DAF" && -n "$CHUNK_ROOT" && -n "$GAMMA_PREFIX" ]] || usage

OUTFILE="${GAMMA_PREFIX}.${DAF}"
TMP_VALUES="${CHUNK_ROOT}/all_stats_${DAF}.tmp"

find "$CHUNK_ROOT" -type f -name "res_${DAF}_*.tab" -print0 \
    | xargs -r -0 cat \
    | cut -f2 > "$TMP_VALUES"

if [[ ! -s "$TMP_VALUES" ]]; then
    echo "No non-empty branch-stat values found for DAF $DAF under $CHUNK_ROOT" >&2
    : > "$OUTFILE"
    exit 1
fi

awk -v daf="$DAF" '
    {sum += $1; n += 1}
    END {
        if (n < 1) {
            exit 1
        }
        printf("%s\t%.16g\n", daf, sum / n)
    }
' "$TMP_VALUES" > "$OUTFILE"

echo "gamma_file=$OUTFILE"
