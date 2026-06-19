#!/bin/bash
set -euo pipefail

INPUT_VCF="$1"
SAMPLE_LIST="$2"
OUTDIR="$3"
PREFIX="$4"
REGION="${5:-}"

mkdir -p "$OUTDIR/tmp"
SFILE_WITH_NAMES="$OUTDIR/tmp/${PREFIX}_s_file.with_names.txt"
SFILE_FINAL="$OUTDIR/${PREFIX}_s_file.txt"
CLEAN_LIST="$OUTDIR/tmp/${PREFIX}_clean_samples.list"
MAX_SINGLETON_MISSING_SAMPLES="${MAX_SINGLETON_MISSING_SAMPLES:-10}"

# --- FIX: Extract 1st column to handle space-separated lists ---
gawk '{print $1}' "$SAMPLE_LIST" > "$CLEAN_LIST"
SAMPLE_COUNT="$(awk 'END{print NR}' "$CLEAN_LIST")"
if [[ -z "$SAMPLE_COUNT" || "$SAMPLE_COUNT" -le 0 ]]; then
    echo "[Error] No samples found in $SAMPLE_LIST" >&2
    exit 1
fi
MAX_SINGLETON_MISSING_FRACTION="$(awk -v n="$SAMPLE_COUNT" -v cap="$MAX_SINGLETON_MISSING_SAMPLES" '
BEGIN {
    frac = cap / n
    if (frac > 0.005) frac = 0.005
    printf "%.12f", frac
}')"
# ---------------------------------------------------------------

if [[ -n "$REGION" ]]; then REG_OPT="-r $REGION"; else REG_OPT=""; fi

echo "[Step 1] Extracting singletons for $PREFIX..."
echo "[s_file] prefix=${PREFIX} sample_count=${SAMPLE_COUNT} max_missing_samples=${MAX_SINGLETON_MISSING_SAMPLES} max_missing_fraction=${MAX_SINGLETON_MISSING_FRACTION}" >&2

bcftools view $REG_OPT -S "$CLEAN_LIST" --force-samples -Ou "$INPUT_VCF" | \
bcftools view -i "AC=1 && F_MISSING<=${MAX_SINGLETON_MISSING_FRACTION}" -Ou | \
bcftools query -f '%POS[\t%GT]\n' | \
gawk -v sample_list="$CLEAN_LIST" -v out_file="$SFILE_WITH_NAMES" '
BEGIN {
    while ((getline name < sample_list) > 0) {
        if (name !~ /^#/) {
            user_samples[++n] = name
            s_map[name] = ""
        }
    }
    close(sample_list)
}
{
    pos = $1
    for (i = 2; i <= NF; i++) {
        if ($i ~ /[1]/) {
            sample_idx = i - 1
            if (sample_idx <= n) {
                name = user_samples[sample_idx]
                if (s_map[name] == "") {
                    s_map[name] = pos
                } else {
                    s_map[name] = s_map[name] "\t" pos
                }
            }
            break 
        }
    }
}
END {
    for (i = 1; i <= n; i++) {
        name = user_samples[i]
        if (s_map[name] == "") {
            print name, "NA" > out_file
        } else {
            print name, s_map[name] > out_file
        }
    }
}
'

gawk '{$1=""; sub(/^[ \t]+/, ""); print}' OFS="\t" "$SFILE_WITH_NAMES" > "$SFILE_FINAL"
echo "✓ s_file generated for $PREFIX"
