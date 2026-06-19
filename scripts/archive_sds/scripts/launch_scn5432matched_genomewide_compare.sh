#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SDS_PYTHON="/data/home/grp-wangyf/intern/miniforge3/envs/sds/bin/python"

OLD_ROOT="$BASE_DIR/data/processed/sds_output_gravel_chb_ne100k_newinput_20260511/SCN"
NEW_ROOT="$BASE_DIR/data/processed/sds_output_gravel_chb_ne100k_scn5432matched_genomewide_run2_20260526/SCN"
OUT_DIR="$BASE_DIR/data/processed/scn5432matched_genomewide_compare_20260526"

mkdir -p "$OUT_DIR"

"$SDS_PYTHON" "$SCRIPT_DIR/compare_sds_diagnostics.py" \
  --track "old=$OLD_ROOT/diagnostics/SCN.gravel_chb_ne100k_manual" \
  --track "matched=$NEW_ROOT/diagnostics/SCN.gravel_chb_ne100k_submitter_array" \
  --output-prefix "$OUT_DIR/SCN_gravel_old_vs_scn5432matched"

"$SDS_PYTHON" "$SCRIPT_DIR/compare_sds_frequency_bins.py" \
  --label-a old \
  --label-b matched \
  --bins-a "$OLD_ROOT/SCN.frequency_bins.tsv" \
  --bins-b "$NEW_ROOT/SCN.frequency_bins.tsv" \
  --output-tsv "$OUT_DIR/SCN_gravel_old_vs_scn5432matched.frequency_bins_compare.tsv"

printf 'compare_dir\t%s\n' "$OUT_DIR"
