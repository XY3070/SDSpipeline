#!/bin/bash
set -euo pipefail

echo "== Jobs =="
bjobs -w 484772 484773 484774 484775 484776 2>/dev/null || true

echo
echo "== SCN review outputs =="
ls -lh \
  /data/home/grp-wangyf/xuyuan/sds/data/processed/sds_output_inputfix_olddefault_chr1/SCN/review/SCN.chr1_inputfix_olddefault.summary.tsv \
  /data/home/grp-wangyf/xuyuan/sds/data/processed/sds_output_inputfix_olddefault_chr1/SCN/review/SCN.chr1_inputfix_olddefault.top_hits.tsv \
  2>/dev/null || true

echo
echo "== MERGED3971 review outputs =="
ls -lh \
  /data/home/grp-wangyf/xuyuan/sds/data/processed/sds_output_parent3971_olddefault_chr1/MERGED3971/review/MERGED3971.chr1_parent3971_olddefault.summary.tsv \
  /data/home/grp-wangyf/xuyuan/sds/data/processed/sds_output_parent3971_olddefault_chr1/MERGED3971/review/MERGED3971.chr1_parent3971_olddefault.top_hits.tsv \
  2>/dev/null || true

echo
echo "== MERGED3971 VCF =="
ls -lh \
  /data/home/grp-wangyf/xuyuan/sds/data/vcf/MERGED3971/UKBQC_MERGED3971_chr1.vcf.gz \
  /data/home/grp-wangyf/xuyuan/sds/data/vcf/MERGED3971/UKBQC_MERGED3971_chr1.vcf.gz.tbi \
  2>/dev/null || true
