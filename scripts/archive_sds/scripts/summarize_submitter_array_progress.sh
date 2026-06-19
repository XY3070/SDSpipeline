#!/bin/bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: $0 STATE_DIR [OUT_ROOT]" >&2
    exit 1
fi

STATE_DIR="$1"
MANIFEST="$STATE_DIR/manifest.tsv"
RESULT_DIR="$STATE_DIR/results"

[[ -f "$MANIFEST" ]] || { echo "[Error] manifest not found: $MANIFEST" >&2; exit 1; }
[[ -d "$RESULT_DIR" ]] || { echo "[Error] result dir not found: $RESULT_DIR" >&2; exit 1; }

OUT_ROOT="${2:-$(awk -F '\t' '$1=="out_root"{print $2; exit}' "$MANIFEST")}"
[[ -n "$OUT_ROOT" ]] || { echo "[Error] failed to determine out_root" >&2; exit 1; }

POP_DIRS_FOUND="$(find "$OUT_ROOT" -mindepth 1 -maxdepth 1 -type d | wc -l | tr -d ' ')"
if [[ "$POP_DIRS_FOUND" -eq 1 ]]; then
    POP_OUTDIR="$(find "$OUT_ROOT" -mindepth 1 -maxdepth 1 -type d | head -n 1)"
else
    POP_OUTDIR="$OUT_ROOT"
fi

total_submitter_out="$(find "$RESULT_DIR" -maxdepth 1 -type f -name 'chr*.submitter.out' | wc -l | tr -d ' ')"
nonempty_submitter_out="$(find "$RESULT_DIR" -maxdepth 1 -type f -name 'chr*.submitter.out' -size +0c | wc -l | tr -d ' ')"
total_final_job="$(find "$RESULT_DIR" -maxdepth 1 -type f -name 'chr*.final_job.tsv' | wc -l | tr -d ' ')"

total_chr_sds="$(find "$POP_OUTDIR" -maxdepth 1 -type f -name 'chr*.sds.tsv' | wc -l | tr -d ' ')"
total_chr_archive="$(find "$POP_OUTDIR/logs" -maxdepth 1 -type f -name 'chr*.archive.csv' 2>/dev/null | wc -l | tr -d ' ')"
has_normalized=0
has_frequency_bins=0
has_manhattan=0
has_diag_tail=0

[[ -f "$POP_OUTDIR/$(basename "$POP_OUTDIR").normalized.tsv" ]] && has_normalized=1
[[ -f "$POP_OUTDIR/$(basename "$POP_OUTDIR").frequency_bins.tsv" ]] && has_frequency_bins=1
[[ -f "$POP_OUTDIR/$(basename "$POP_OUTDIR").manhattan.png" ]] && has_manhattan=1
find "$POP_OUTDIR/diagnostics" -maxdepth 1 -type f -name '*.chrom_tail_counts.tsv' >/dev/null 2>&1 && has_diag_tail=1 || true

printf 'state_dir\t%s\n' "$STATE_DIR"
printf 'out_root\t%s\n' "$OUT_ROOT"
printf 'pop_outdir\t%s\n' "$POP_OUTDIR"
printf 'submitter_out_total\t%s\n' "$total_submitter_out"
printf 'submitter_out_nonempty\t%s\n' "$nonempty_submitter_out"
printf 'final_job_tsv_count\t%s\n' "$total_final_job"
printf 'chr_sds_tsv_count\t%s\n' "$total_chr_sds"
printf 'chr_archive_csv_count\t%s\n' "$total_chr_archive"
printf 'has_normalized\t%s\n' "$has_normalized"
printf 'has_frequency_bins\t%s\n' "$has_frequency_bins"
printf 'has_manhattan_png\t%s\n' "$has_manhattan"
printf 'has_diagnostics_chrom_tail\t%s\n' "$has_diag_tail"
