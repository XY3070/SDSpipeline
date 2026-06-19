#!/bin/bash
set -euo pipefail

INPUT_VCF="$1"
SAMPLE_LIST="$2"
SFILE_WITH_NAMES="$3"
OUTDIR="$4"
PREFIX="${5:-chr1}"
# Arguments 6-9 are ignored in uniform mode

OFILE="$OUTDIR/${PREFIX}_o_file.txt"
ORDER_FILE="$OUTDIR/tmp/${PREFIX}_order_s.txt"

echo "Generate o_file: $PREFIX (Mode: Uniform/Fast)"

# 1. Get sample count from the s_file (which is already cleaned and sorted)
# We use awk to handle spaces/tabs and ensure we just count valid lines
if [[ -s "$SFILE_WITH_NAMES" ]]; then
    awk '{print $1}' "$SFILE_WITH_NAMES" > "$ORDER_FILE"
else
    # Fallback if s_file missing (should not happen in pipeline)
    grep -v '^#' "$SAMPLE_LIST" | awk '{print $1}' > "$ORDER_FILE"
fi

N_SAMPLES=$(wc -l < "$ORDER_FILE")

# 2. Generate a single line of "1"s separated by tabs
# 'yes 1' generates infinite 1s
# 'head' limits it to N_SAMPLES
# 'tr' converts newlines to tabs
if [[ $N_SAMPLES -gt 0 ]]; then
    yes 1 | head -n "$N_SAMPLES" | tr '\n' '\t' | sed 's/\t$//' > "$OFILE"
    # Add a newline at the end to be safe
    echo "" >> "$OFILE"
else
    echo "[Error] No samples found to generate o_file"
    exit 1
fi

echo "✓ o_file generated: $N_SAMPLES samples (All set to 1)"