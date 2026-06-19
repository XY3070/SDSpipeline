#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$SCRIPT_DIR/common_env.sh"
activate_sds_env

usage() {
    cat >&2 <<'EOF'
Usage: merge_split_chr_vcfs.sh --chr N --pop-a POPA --pop-b POPB --pop-out OUTPOP
EOF
    exit 1
}

CHR=""
POP_A=""
POP_B=""
POP_OUT=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --chr) CHR="$2"; shift 2 ;;
        --pop-a) POP_A="$2"; shift 2 ;;
        --pop-b) POP_B="$2"; shift 2 ;;
        --pop-out) POP_OUT="$2"; shift 2 ;;
        *) usage ;;
    esac
done

[[ -n "$CHR" && -n "$POP_A" && -n "$POP_B" && -n "$POP_OUT" ]] || usage

VCF_A="$BASE_DIR/data/vcf/$POP_A/UKBQC_${POP_A}_chr${CHR}.vcf.gz"
VCF_B="$BASE_DIR/data/vcf/$POP_B/UKBQC_${POP_B}_chr${CHR}.vcf.gz"
OUT_DIR="$BASE_DIR/data/vcf/$POP_OUT"
OUT_VCF="$OUT_DIR/UKBQC_${POP_OUT}_chr${CHR}.vcf.gz"
OUT_TBI="${OUT_VCF}.tbi"
PART_VCF="${OUT_VCF}.part.vcf.gz"
PART_TBI="${PART_VCF}.tbi"
SAMPLE_OUT="$BASE_DIR/data/${POP_OUT}.txt"
TMP_SAMPLES="$SAMPLE_OUT.part"

[[ -f "$VCF_A" ]] || { echo "[Error] Missing VCF: $VCF_A" >&2; exit 1; }
[[ -f "$VCF_B" ]] || { echo "[Error] Missing VCF: $VCF_B" >&2; exit 1; }

mkdir -p "$OUT_DIR"

awk 'NF > 0 && $1 !~ /^#/' "$BASE_DIR/data/${POP_A}.txt" "$BASE_DIR/data/${POP_B}.txt" > "$TMP_SAMPLES"
mv "$TMP_SAMPLES" "$SAMPLE_OUT"

rm -f "$PART_VCF" "$PART_TBI"

echo ">>> Merging chr${CHR}: $POP_A + $POP_B -> $POP_OUT"
echo ">>> Output: $OUT_VCF"

bcftools merge \
    -m none \
    --threads 4 \
    -Oz \
    -o "$PART_VCF" \
    "$VCF_A" \
    "$VCF_B"

bcftools index -t "$PART_VCF"

mv "$PART_VCF" "$OUT_VCF"
mv "$PART_TBI" "$OUT_TBI"

echo ">>> Verifying merged VCF sample count"
bcftools query -l "$OUT_VCF" | wc -l

echo "✓ DONE: $OUT_VCF"
