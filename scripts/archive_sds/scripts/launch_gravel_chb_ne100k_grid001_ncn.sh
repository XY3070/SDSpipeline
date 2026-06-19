#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SDS_PYTHON="/data/home/grp-wangyf/intern/miniforge3/envs/sds/bin/python"

LOW_OUTDIR="${1:-$BASE_DIR/tmp/gravel_chb_gamma_ne100k_grid001_$(date +%Y%m%d_%H%M%S)}"
BASE_GFILE="$BASE_DIR/tmp/gravel_chb_gamma_ne100k_20260502/gravel_chb_present100000.g_file.txt"
COMBINED_GFILE="$LOW_OUTDIR/gravel_chb_present100000.grid001.g_file.txt"
OUTPUT_ROOT="${2:-$BASE_DIR/data/processed/sds_output_gravel_chb_ne100k_grid001_ncn_$(date +%Y%m%d_%H%M%S)}"
INPUT_ROOT="$BASE_DIR/data/processed/sds_input_rebuilt_main_contract_20260511"
QUEUE="${3:-normal}"
LOG_DIR="$OUTPUT_ROOT/logs"
MANIFEST="$LOW_OUTDIR/launch_manifest.tsv"

mkdir -p "$LOW_OUTDIR" "$OUTPUT_ROOT" "$LOG_DIR"
[[ -f "$BASE_GFILE" ]] || { echo "[Error] Missing baseline g_file: $BASE_GFILE" >&2; exit 1; }

printf 'key\tvalue\n' > "$MANIFEST"
printf 'baseline_gfile\t%s\n' "$BASE_GFILE" >> "$MANIFEST"
printf 'low_outdir\t%s\n' "$LOW_OUTDIR" >> "$MANIFEST"
printf 'combined_gfile\t%s\n' "$COMBINED_GFILE" >> "$MANIFEST"
printf 'output_root\t%s\n' "$OUTPUT_ROOT" >> "$MANIFEST"

"$SDS_PYTHON" "$SCRIPT_DIR/submit_gravel_chb_gamma_parallel.py" \
  --dry-run \
  --queue smp \
  --outdir "$LOW_OUTDIR" \
  --present-ne 100000 \
  --pop-model Gravel_CHB \
  --sim-reps 1000 \
  --chunk-size 100 \
  --daf-start 0.01 \
  --daf-end 0.04 \
  --daf-step 0.01 \
  > "$LOW_OUTDIR/gamma_dry_run.json"

CHUNK_TASK_FILE="$LOW_OUTDIR/parallel_submit/chunks_tasks.jsonl"
AGG_TASK_FILE="$LOW_OUTDIR/parallel_submit/aggregates_tasks.jsonl"
CHUNK_TASKS="$(wc -l < "$CHUNK_TASK_FILE")"
AGG_TASKS="$(wc -l < "$AGG_TASK_FILE")"
[[ "$CHUNK_TASKS" -gt 0 && "$AGG_TASKS" -gt 0 ]] || { echo "[Error] Dry-run did not create task files" >&2; exit 1; }

CHUNK_SUBMIT_OUT="$(bsub -q smp -J "gravel_chb_ne100k_grid001_chunks[1-${CHUNK_TASKS}]%96" -n 1 -R "span[hosts=1]" -cwd "$BASE_DIR/sds" -o "$BASE_DIR/sds/logs/gravel_chb_ne100k_grid001_chunks_%J_%I.out" -e "$BASE_DIR/sds/logs/gravel_chb_ne100k_grid001_chunks_%J_%I.err" "$SDS_PYTHON" "$SCRIPT_DIR/submit_gravel_chb_gamma_parallel.py" --run-task-file "$CHUNK_TASK_FILE" < /dev/null)"
CHUNK_JOB="$(printf '%s\n' "$CHUNK_SUBMIT_OUT" | sed -n 's/.*<\([0-9]\+\)>.*/\1/p')"
[[ -n "$CHUNK_JOB" ]] || { echo "[Error] Could not parse CHUNK_JOB" >&2; exit 1; }
printf 'chunk_job\t%s\n' "$CHUNK_JOB" >> "$MANIFEST"

AGG_SUBMIT_OUT="$(bsub -q smp -w "done(${CHUNK_JOB})" -J "gravel_chb_ne100k_grid001_aggs[1-${AGG_TASKS}]%24" -n 1 -R "span[hosts=1]" -cwd "$BASE_DIR/sds" -o "$BASE_DIR/sds/logs/gravel_chb_ne100k_grid001_aggs_%J_%I.out" -e "$BASE_DIR/sds/logs/gravel_chb_ne100k_grid001_aggs_%J_%I.err" "$SDS_PYTHON" "$SCRIPT_DIR/submit_gravel_chb_gamma_parallel.py" --run-task-file "$AGG_TASK_FILE" < /dev/null)"
AGG_JOB="$(printf '%s\n' "$AGG_SUBMIT_OUT" | sed -n 's/.*<\([0-9]\+\)>.*/\1/p')"
[[ -n "$AGG_JOB" ]] || { echo "[Error] Could not parse AGG_JOB" >&2; exit 1; }
printf 'aggregate_job\t%s\n' "$AGG_JOB" >> "$MANIFEST"

FINALIZE_CMD="cd \"$BASE_DIR\" && \"$SDS_PYTHON\" \"$SCRIPT_DIR/submit_gravel_chb_gamma_parallel.py\" --finalize-only --outdir \"$LOW_OUTDIR\" --present-ne 100000 --pop-model Gravel_CHB --dummy-npz \"$BASE_DIR/tmp/region_ne0_positive_20260426/scaled_ne0_100000/NCN_scaled_ne0_100000.npz\" --sim-reps 1000 --chunk-size 100 --daf-start 0.01 --daf-end 0.04 --daf-step 0.01 --ms-make-dir /data/home/grp-wangyf/xuyuan/ms/scripts --ms-binary /data/home/grp-wangyf/xuyuan/ms/msdir/ms --backward-script /data/home/grp-wangyf/xuyuan/ms/scripts/backward.py"
FINALIZE_SUBMIT_OUT="$(bsub -q smp -w "done(${AGG_JOB})" -n 1 -R "span[hosts=1]" -J gravel_chb_ne100k_grid001_finalize -o "$LOW_OUTDIR/finalize.out" -e "$LOW_OUTDIR/finalize.err" /bin/bash -lc "$FINALIZE_CMD" < /dev/null)"
FINAL_GAMMA_JOB="$(printf '%s\n' "$FINALIZE_SUBMIT_OUT" | sed -n 's/.*<\([0-9]\+\)>.*/\1/p')"
[[ -n "$FINAL_GAMMA_JOB" ]] || { echo "[Error] Could not parse FINAL_GAMMA_JOB" >&2; exit 1; }
printf 'final_gamma_job\t%s\n' "$FINAL_GAMMA_JOB" >> "$MANIFEST"

BUILD_CMD="cd \"$BASE_DIR\" && \"$SDS_PYTHON\" \"$SCRIPT_DIR/build_extended_gamma_file.py\" --baseline-gfile \"$BASE_GFILE\" --low-piece-dir \"$LOW_OUTDIR/final_gamma\" --output-gfile \"$COMBINED_GFILE\""
BUILD_SUBMIT_OUT="$(bsub -q smp -w "done(${FINAL_GAMMA_JOB})" -n 1 -R "span[hosts=1]" -J gravel_chb_grid001_build -o "$LOW_OUTDIR/build_extended_gamma.out" -e "$LOW_OUTDIR/build_extended_gamma.err" /bin/bash -lc "$BUILD_CMD" < /dev/null)"
BUILD_JOB="$(printf '%s\n' "$BUILD_SUBMIT_OUT" | sed -n 's/.*<\([0-9]\+\)>.*/\1/p')"
[[ -n "$BUILD_JOB" ]] || { echo "[Error] Could not parse BUILD_JOB" >&2; exit 1; }
printf 'build_job\t%s\n' "$BUILD_JOB" >> "$MANIFEST"

LAUNCH_CMD="bash \"$SCRIPT_DIR/launch_ncn_from_ready_gfile.sh\" NCN \"$INPUT_ROOT\" \"$OUTPUT_ROOT\" \"$COMBINED_GFILE\" \"$QUEUE\""
LAUNCH_SUBMIT_OUT="$(bsub -q "$QUEUE" -w "done(${BUILD_JOB})" -n 1 -R "span[hosts=1]" -J ncn_gravel_grid001_launcher -o "$LOG_DIR/ncn_launcher.out" -e "$LOG_DIR/ncn_launcher.err" /bin/bash -lc "$LAUNCH_CMD" < /dev/null)"
LAUNCH_JOB="$(printf '%s\n' "$LAUNCH_SUBMIT_OUT" | sed -n 's/.*<\([0-9]\+\)>.*/\1/p')"
[[ -n "$LAUNCH_JOB" ]] || { echo "[Error] Could not parse LAUNCH_JOB" >&2; exit 1; }
printf 'launch_job\t%s\n' "$LAUNCH_JOB" >> "$MANIFEST"

printf 'manifest\t%s\n' "$MANIFEST"
