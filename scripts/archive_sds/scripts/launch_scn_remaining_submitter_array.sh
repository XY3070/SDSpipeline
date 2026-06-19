#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
OUT_ROOT="$BASE_DIR/data/processed/sds_output_gravel_chb_ne100k_newinput_20260511"
IN_ROOT="$BASE_DIR/data/processed/sds_input_rebuilt_main_contract_20260511"
G_FILE="$BASE_DIR/tmp/gravel_chb_gamma_ne100k_20260502/gravel_chb_present100000.g_file.txt"
CHR_SPEC=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --chrs) CHR_SPEC="$2"; shift 2 ;;
        *) echo "Unknown parameter: $1" >&2; exit 1 ;;
    esac
done

TS="$(date +%Y%m%d_%H%M%S)_$$"
STATE_DIR="$BASE_DIR/data/processed/scn_submitter_array_$TS"
CHR_LIST="$STATE_DIR/chrs.txt"
RESULT_DIR="$STATE_DIR/results"
MANIFEST="$STATE_DIR/manifest.tsv"
LOG_DIR="$OUT_ROOT/SCN/logs"

mkdir -p "$STATE_DIR" "$RESULT_DIR" "$LOG_DIR"

if [[ -n "$CHR_SPEC" ]]; then
    tr ',' '\n' <<< "$CHR_SPEC" > "$CHR_LIST"
else
python3 - <<'PY' > "$CHR_LIST"
from pathlib import Path
out = Path('/data/home/grp-wangyf/xuyuan/sds/data/processed/sds_output_gravel_chb_ne100k_newinput_20260511/SCN')
for c in range(1, 23):
    if not (out / f'chr{c}.sds.tsv').exists():
        print(c)
PY
fi

N="$(wc -l < "$CHR_LIST")"
[[ "$N" -gt 0 ]] || { echo "No missing SCN chromosomes."; exit 0; }

ARRAY_OUT="$(bsub -q smp -n 4 -R "span[hosts=1]" -J "SDS_SUBMITTER_SCN_REMAINING[1-${N}]%2" -o "$LOG_DIR/SCN_submitter_remaining_%I.out" -e "$LOG_DIR/SCN_submitter_remaining_%I.err" /bin/bash -lc "bash \"$SCRIPT_DIR/run_submitter_array_entry.sh\" SCN \"$CHR_LIST\" \"$IN_ROOT\" \"$OUT_ROOT\" \"$G_FILE\" gravel_chb_ne100k_SCN_submitter_array \"$RESULT_DIR\"" < /dev/null)"
ARRAY_JOB="$(sed -n 's/.*<\([0-9]\+\)>.*/\1/p' <<< "$ARRAY_OUT")"
[[ -n "$ARRAY_JOB" ]] || { echo "[Error] Failed to parse submitter array job id" >&2; exit 1; }

FINALIZER_OUT="$(bsub -q normal -w "done(${ARRAY_JOB})" -n 1 -R "span[hosts=1]" -J "SDS_SUBMITTER_SCN_FINALIZE" -o "$LOG_DIR/SCN_submitter_finalize.out" -e "$LOG_DIR/SCN_submitter_finalize.err" /bin/bash -lc "bash \"$SCRIPT_DIR/submit_postprocess_from_final_jobs.sh\" SCN \"$RESULT_DIR\" \"$OUT_ROOT\" normal" < /dev/null)"
FINALIZER_JOB="$(sed -n 's/.*<\([0-9]\+\)>.*/\1/p' <<< "$FINALIZER_OUT")"
[[ -n "$FINALIZER_JOB" ]] || { echo "[Error] Failed to parse finalizer job id" >&2; exit 1; }

printf 'key\tvalue\n' > "$MANIFEST"
printf 'state_dir\t%s\nchr_list\t%s\nresult_dir\t%s\nsubmitter_array_job\t%s\nfinalizer_job\t%s\n' \
    "$STATE_DIR" "$CHR_LIST" "$RESULT_DIR" "$ARRAY_JOB" "$FINALIZER_JOB" >> "$MANIFEST"

printf 'manifest\t%s\n' "$MANIFEST"
