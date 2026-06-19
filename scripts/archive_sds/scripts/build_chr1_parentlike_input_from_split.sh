#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$SCRIPT_DIR/common_env.sh"
activate_sds_env

NCN_S="$BASE_DIR/data/processed/sds_input/NCN/chr1_s_file.txt"
SCN_S="$BASE_DIR/data/processed/sds_input/SCN/chr1_s_file.txt"
NCN_VCF="$BASE_DIR/data/vcf/NCN/UKBQC_NCN_chr1.vcf.gz"
SCN_VCF="$BASE_DIR/data/vcf/SCN/UKBQC_SCN_chr1.vcf.gz"
SRC_MERGED_DIR="$BASE_DIR/data/processed/sds_input/MERGED3971"
OUT_POP="MERGED3971PL"
OUT_DIR="$BASE_DIR/data/processed/sds_input/${OUT_POP}"
TMP_DIR="$OUT_DIR/tmp_parentlike_chr1"

mkdir -p "$OUT_DIR" "$TMP_DIR"

POS_NCN="$TMP_DIR/ncn_singleton_pos.tsv"
POS_SCN="$TMP_DIR/scn_singleton_pos.tsv"
MASK_NCN="$TMP_DIR/ncn_masked_by_scn.pos"
MASK_SCN="$TMP_DIR/scn_masked_by_ncn.pos"

echo ">>> Building unique singleton position lists"
gawk 'BEGIN{OFS="\t"} $0!="NA"{n=split($0,a,/\t/); for(i=1;i<=n;i++) if(a[i]!="") print "chr1", a[i]}' "$NCN_S" | sort -u -k1,1 -k2,2n > "$POS_NCN"
gawk 'BEGIN{OFS="\t"} $0!="NA"{n=split($0,a,/\t/); for(i=1;i<=n;i++) if(a[i]!="") print "chr1", a[i]}' "$SCN_S" | sort -u -k1,1 -k2,2n > "$POS_SCN"

echo ">>> Querying opposite-cohort variant positions"
bcftools view -T "$POS_NCN" -f PASS -m2 -M2 -v snps -Ou "$SCN_VCF" | bcftools query -f '%POS\n' | sort -u > "$MASK_NCN"
bcftools view -T "$POS_SCN" -f PASS -m2 -M2 -v snps -Ou "$NCN_VCF" | bcftools query -f '%POS\n' | sort -u > "$MASK_SCN"

echo ">>> Building parent-like merged s_file"
"$SDS_ENV_PREFIX/bin/python" "$SCRIPT_DIR/build_parentlike_s_file_from_split_cohorts.py" \
    --ncn-s-file "$NCN_S" \
    --scn-s-file "$SCN_S" \
    --ncn-mask-pos "$MASK_NCN" \
    --scn-mask-pos "$MASK_SCN" \
    --output "$OUT_DIR/chr1_s_file.txt"

echo ">>> Copying t/o/b and sample list"
cp "$SRC_MERGED_DIR/chr1_t_file.txt" "$OUT_DIR/chr1_t_file.txt"
cp "$SRC_MERGED_DIR/chr1_o_file.txt" "$OUT_DIR/chr1_o_file.txt"
cp "$SRC_MERGED_DIR/chr1_b_file.txt" "$OUT_DIR/chr1_b_file.txt"
cp "$BASE_DIR/data/MERGED3971.txt" "$BASE_DIR/data/${OUT_POP}.txt"

non_na="$("$SDS_ENV_PREFIX/bin/python" - <<'PY'
from pathlib import Path
path = Path("/data/home/grp-wangyf/xuyuan/sds/data/processed/sds_input/MERGED3971PL/chr1_s_file.txt")
na = 0
non = 0
with path.open() as h:
    for line in h:
        if line.rstrip("\n") == "NA":
            na += 1
        else:
            non += 1
print(f"NA_rows={na} nonNA_rows={non}")
PY
)"
echo ">>> $non_na"
echo "✓ DONE: $OUT_DIR"
