#!/bin/bash
set -euo pipefail

INPUT_VCF="$1"
SAMPLE_LIST="$2"
SFILE_WITH_NAMES="$3"
OUTDIR="$4"
PREFIX="${5:?Error: Must provide prefix}"
SPATIAL_WINDOW="${6:-1000}"
REGION="${7:-}"

TMP_VCF="$OUTDIR/tmp/${PREFIX}_filtered.vcf.gz"
TFILE_FINAL="$OUTDIR/${PREFIX}_t_file.txt"
ORDER_FILE="$OUTDIR/tmp/${PREFIX}_order_s.txt"

mkdir -p "$OUTDIR/tmp"
if [[ -n "$REGION" ]]; then REG_OPT="-r $REGION"; else REG_OPT=""; fi

echo "[Step 1] Filter SNPs for $PREFIX..."

# --- FIX: Use awk to handle space-separated lists (cut fails on spaces) ---
if [[ -s "$SFILE_WITH_NAMES" ]]; then
    awk '{print $1}' "$SFILE_WITH_NAMES" > "$ORDER_FILE"
else
    grep -v '^#' "$SAMPLE_LIST" | awk '{print $1}' > "$ORDER_FILE"
fi
# ------------------------------------------------------------------------

bcftools view $REG_OPT \
    -m2 -M2 -v snps \
    -S "$ORDER_FILE" \
    --force-samples \
    --threads 1 \
    -Oz -o "$TMP_VCF" \
    "$INPUT_VCF"

bcftools query -f '%ID\t%REF\t%ALT\t%POS[\t%GT]\n' "$TMP_VCF" -S "$ORDER_FILE" --force-samples | \
awk -v OFS="\t" '
{
    n_samples = NF - 4
    n_missing = 0
    for (i = 5; i <= NF; i++) {
        if ($i ~ /^0[\/|]0$/) $i = 0
        else if ($i ~ /^0[\/|]1$/ || $i ~ /^1[\/|]0$/) $i = 1
        else if ($i ~ /^1[\/|]1$/) $i = 2
        else { $i = "NA"; n_missing++ }
    }
    if (n_missing < n_samples) print $0
}
' > "$TFILE_FINAL"

echo "✓ t_file generated for $PREFIX"
