#!/bin/bash
set -euo pipefail
VCF="$1"
SLIST="$2"
OUT_ARM="$3"
REGION="$4"

mkdir -p "$(dirname "$OUT_ARM")"
CLEAN_SLIST="$(dirname "$OUT_ARM")/samples.clean.txt"
gawk 'NF > 0 && $1 !~ /^#/ { print $1 }' "$SLIST" > "$CLEAN_SLIST"

TARGET_CHR="${REGION%%:*}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common_env.sh
source "$SCRIPT_DIR/common_env.sh"
ACTUAL_CHR="$(resolve_vcf_chr "$VCF" "$TARGET_CHR")"
if [[ -z "$ACTUAL_CHR" ]]; then
    if [[ "$TARGET_CHR" == chr* ]]; then
        ACTUAL_CHR="${TARGET_CHR#chr}"
    else
        ACTUAL_CHR="chr${TARGET_CHR}"
    fi
fi
RANGE_PART="${REGION#*:}"
CLEAN_REG="${ACTUAL_CHR}:${RANGE_PART}"

bcftools view -r "$CLEAN_REG" -S "$CLEAN_SLIST" -f PASS -m2 -M2 -v snps --force-samples -Ou "$VCF" | \
bcftools +fill-tags -Ou -- -t AF,F_MISSING | \
bcftools view -i 'AF>=0.05 && AF<=0.95 && F_MISSING<=0.05' -Ou | \
bcftools query -f '%ID\t%REF\t%ALT\t%POS[\t%GT]\n' | \
gawk -v OFS="\t" '
{
    n_samples = NF - 4
    n_missing = 0
    out = sprintf("%s\t%s\t%s\t%s", $1, $2, $3, $4)
    for (i = 5; i <= NF; i++) {
        if ($i ~ /^0[\/|]0$/) out = out "\t0"
        else if ($i ~ /^0[\/|]1$/ || $i ~ /^1[\/|]0$/) out = out "\t1"
        else if ($i ~ /^1[\/|]1$/) out = out "\t2"
        else {
            out = out "\tNA"
            n_missing++
        }
    }
    if (n_missing < n_samples) print out
}' > "$OUT_ARM"
