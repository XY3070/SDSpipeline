#!/bin/bash
set -euo pipefail
VCF="$1"
SLIST="$2"
TMPDIR="$3"
ARM="$4"  # p 或 q
REGION="$5"

mkdir -p "$TMPDIR"
OUT_ARM="$TMPDIR/s_${ARM}.with_names.txt"
CLEAN_SLIST="$TMPDIR/samples.clean.txt"
MAX_SINGLETON_MISSING_SAMPLES="${MAX_SINGLETON_MISSING_SAMPLES:-10}"

gawk 'NF > 0 && $1 !~ /^#/ { print $1 }' "$SLIST" > "$CLEAN_SLIST"
SAMPLE_COUNT="$(awk 'END{print NR}' "$CLEAN_SLIST")"
if [[ -z "$SAMPLE_COUNT" || "$SAMPLE_COUNT" -le 0 ]]; then
    echo "[Error] No samples found in $SLIST" >&2
    exit 1
fi
MAX_SINGLETON_MISSING_FRACTION="$(awk -v n="$SAMPLE_COUNT" -v cap="$MAX_SINGLETON_MISSING_SAMPLES" '
BEGIN {
    frac = cap / n
    if (frac > 0.005) frac = 0.005
    printf "%.12f", frac
}')"

# 处理不同 VCF 的染色体命名习惯 (1 vs chr1)
TARGET_CHR="${REGION%%:*}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/share/home/grp-wangyf/xuyuan/sds/scripts/common_env.sh
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

# 先用 bcftools 在 subset 后补齐 AC/F_MISSING，再把候选收窄到真正的 singleton 位点。
# 缺失阈值同时满足原主契约 F_MISSING<=0.005 和新要求 missing sample count<=10。
echo "[s_file] arm=${ARM} region=${CLEAN_REG} sample_count=${SAMPLE_COUNT} max_missing_samples=${MAX_SINGLETON_MISSING_SAMPLES} max_missing_fraction=${MAX_SINGLETON_MISSING_FRACTION}" >&2
bcftools view -r "$CLEAN_REG" -S "$CLEAN_SLIST" -f PASS -m2 -M2 -v snps --force-samples -Ou "$VCF" | \
bcftools +fill-tags -Ou -- -t AC,F_MISSING | \
bcftools view -i "AC=1 && F_MISSING<=${MAX_SINGLETON_MISSING_FRACTION}" -Ou | \
bcftools query -f '%POS[\t%GT]\n' | \
gawk -v sample_list="$CLEAN_SLIST" -v out_file="$OUT_ARM" '
BEGIN {
    while ((getline name < sample_list) > 0) {
        samples[++n] = name
        s_map[name] = ""
    }
    close(sample_list)
}
{
    pos = $1
    singleton_sample = ""
    for (i = 2; i <= NF; i++) {
        gt = $i
        if (gt ~ /^0[\/|]1$/ || gt ~ /^1[\/|]0$/) {
            singleton_sample = samples[i-1]
            break
        }
    }

    if (singleton_sample != "") {
        s_map[singleton_sample] = (s_map[singleton_sample] == "" ? pos : s_map[singleton_sample] "\t" pos)
    }
}
END { for (i = 1; i <= n; i++) print samples[i] "\t" (s_map[samples[i]] == "" ? "NA" : s_map[samples[i]]) > out_file }'
