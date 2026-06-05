#!/bin/bash
set -euo pipefail

SAMPLE_LIST="$1"
OUTDIR="$2"
PREFIX="${3:?Error: Must provide prefix}"
SFILE_WITH_NAMES="${4:-}"

mkdir -p "$OUTDIR" "$OUTDIR/tmp"
OFILE="$OUTDIR/${PREFIX}_o_file.txt"
ORDER_FILE="$OUTDIR/tmp/${PREFIX}_order_s.txt"

if [[ -n "$SFILE_WITH_NAMES" && -s "$SFILE_WITH_NAMES" ]]; then
    gawk '{print $1}' "$SFILE_WITH_NAMES" > "$ORDER_FILE"
else
    gawk 'NF > 0 && $1 !~ /^#/ { print $1 }' "$SAMPLE_LIST" > "$ORDER_FILE"
fi

N_SAMPLES=$(wc -l < "$ORDER_FILE")
if [[ "$N_SAMPLES" -le 0 ]]; then
    echo "[Error] No samples found to generate o_file" >&2
    exit 1
fi

gawk -v n="$N_SAMPLES" 'BEGIN {
    for (i = 1; i <= n; i++) {
        printf "1"
        if (i < n) printf "\t"
    }
    printf "\n"
}' > "$OFILE"
