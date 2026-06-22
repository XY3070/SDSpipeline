#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

OUT_ROOT="$SDS_DATA_PROCESSED_ROOT/sds_output_ncn_smcppquick_chr1_6_20260526"
IN_ROOT="$SDS_DATA_PROCESSED_ROOT/sds_input_rebuilt_main_contract_20260511"
G_FILE="$SDS_TMP_ROOT/ncn_smcppquick_smokegrid_20260526/gravel_chb_present87282.g_file.txt"
CHR_SPEC="1,6"
POP="NCN"
TS="$(date +%Y%m%d_%H%M%S)_$$"
STATE_DIR="$SDS_RUNS_ROOT/ncn_smcppquick_chr1_6_submitter_array_$TS"
CHR_LIST="$STATE_DIR/chrs.txt"
RESULT_DIR="$STATE_DIR/results"
MANIFEST="$STATE_DIR/manifest.tsv"
LOG_DIR="$OUT_ROOT/$POP/logs"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --out-root) OUT_ROOT="$2"; shift 2 ;;
        --in-root) IN_ROOT="$2"; shift 2 ;;
        --g-file) G_FILE="$2"; shift 2 ;;
        --chrs) CHR_SPEC="$2"; shift 2 ;;
        *) echo "Unknown parameter: $1" >&2; exit 1 ;;
    esac
done

[[ -f "$G_FILE" ]] || { echo "[Error] g_file not found: $G_FILE" >&2; exit 1; }

mkdir -p "$STATE_DIR" "$RESULT_DIR" "$LOG_DIR" "$OUT_ROOT/$POP"
tr ',' '\n' <<< "$CHR_SPEC" > "$CHR_LIST"

N="$(wc -l < "$CHR_LIST")"
[[ "$N" -gt 0 ]] || { echo "No chromosomes requested."; exit 1; }

ARRAY_OUT="$(bsub -q smp -n 4 -R "span[hosts=1]" -J "SDS_SUBMITTER_NCNSMCPPQ[1-${N}]%2" -o "$LOG_DIR/NCN_smcppquick_submitter_%I.out" -e "$LOG_DIR/NCN_smcppquick_submitter_%I.err" /bin/bash -lc "bash \"$SCRIPT_DIR/run_submitter_array_entry.sh\" \"$POP\" \"$CHR_LIST\" \"$IN_ROOT\" \"$OUT_ROOT\" \"$G_FILE\" ncn_smcppquick_submitter_array \"$RESULT_DIR\"" < /dev/null)"
ARRAY_JOB="$(sed -n 's/.*<\([0-9]\+\)>.*/\1/p' <<< "$ARRAY_OUT")"
[[ -n "$ARRAY_JOB" ]] || { echo "[Error] Failed to parse submitter array job id" >&2; exit 1; }

FINALIZER_OUT="$(bsub -q normal -w "done(${ARRAY_JOB})" -n 1 -R "span[hosts=1]" -J "SDS_SUBMITTER_NCNSMCPPQ_FINALIZE" -o "$LOG_DIR/NCN_smcppquick_submitter_finalize.out" -e "$LOG_DIR/NCN_smcppquick_submitter_finalize.err" /bin/bash -lc "bash \"$SCRIPT_DIR/submit_postprocess_from_final_jobs.sh\" \"$POP\" \"$RESULT_DIR\" \"$OUT_ROOT\" normal" < /dev/null)"
FINALIZER_JOB="$(sed -n 's/.*<\([0-9]\+\)>.*/\1/p' <<< "$FINALIZER_OUT")"
[[ -n "$FINALIZER_JOB" ]] || { echo "[Error] Failed to parse finalizer job id" >&2; exit 1; }

printf 'key\tvalue\n' > "$MANIFEST"
printf 'state_dir\t%s\nchr_list\t%s\nresult_dir\t%s\nout_root\t%s\ng_file\t%s\nsubmitter_array_job\t%s\nfinalizer_job\t%s\n' \
    "$STATE_DIR" "$CHR_LIST" "$RESULT_DIR" "$OUT_ROOT" "$G_FILE" "$ARRAY_JOB" "$FINALIZER_JOB" >> "$MANIFEST"

printf 'manifest\t%s\n' "$MANIFEST"
