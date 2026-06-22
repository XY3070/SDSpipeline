#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$SCRIPT_DIR/common_env.sh"
activate_sds_env

SUBMIT_CHUNKED="$SCRIPT_DIR/submit_sds_compute_chunked_chr.sh"
POSTPROCESS="$SCRIPT_DIR/postprocess_sds_results.py"
DIAGNOSE="$SCRIPT_DIR/diagnose_sds_scan.py"

POP="${1:-NCN}"
INPUT_ROOT="${2:-$SDS_DATA_PROCESSED_ROOT/sds_input_rebuilt_main_contract_20260511}"
OUTPUT_ROOT="${3:-$SDS_DATA_PROCESSED_ROOT/sds_output_gravel_chb_ne100k_grid001_ncn_$(date +%Y%m%d_%H%M%S)}"
G_FILE="${4:-}"
QUEUE="${5:-normal}"

[[ -n "$G_FILE" && -f "$G_FILE" ]] || { echo "[Error] Missing g_file: $G_FILE" >&2; exit 1; }
mkdir -p "$OUTPUT_ROOT/$POP/logs"
LOGDIR="$OUTPUT_ROOT/$POP/logs"

FINAL_JOBS=()
for chr in $(seq 1 22); do
    submit_output="$(bash "$SUBMIT_CHUNKED" \
        --pop "$POP" \
        --chr "$chr" \
        --queue "$QUEUE" \
        --in-root "$INPUT_ROOT" \
        --out-root "$OUTPUT_ROOT" \
        --g-file "$G_FILE" \
        --job-group "gravel_chb_ne100k_grid001_chr${chr}")"
    printf '%s\n' "$submit_output" > "$LOGDIR/chr${chr}.submit.out"
    final_job="$(printf '%s\n' "$submit_output" | awk -F '\t' '$1=="FINAL_JOB"{print $2; exit}')"
    [[ -n "$final_job" ]] || { echo "[Error] Could not parse FINAL_JOB for chr${chr}" >&2; exit 1; }
    FINAL_JOBS+=("$final_job")
done

DEP_EXPR="$(printf 'done(%s) && ' "${FINAL_JOBS[@]}")"
DEP_EXPR="${DEP_EXPR% && }"

POST_CMD="cd \"$BASE_DIR\" && \"$SDS_PYTHON\" \"$POSTPROCESS\" --input-dir \"$OUTPUT_ROOT/$POP\" --pop $POP"
POST_SUBMIT_OUT="$(bsub -q "$QUEUE" -w "$DEP_EXPR" -n 1 -R "span[hosts=1]" -J "${POP}_gravel_grid001_post" -o "$LOGDIR/postprocess.out" -e "$LOGDIR/postprocess.err" /bin/bash -lc "$POST_CMD" < /dev/null)"
POST_JOB="$(printf '%s\n' "$POST_SUBMIT_OUT" | sed -n 's/.*<\([0-9]\+\)>.*/\1/p')"
[[ -n "$POST_JOB" ]] || { echo "[Error] Could not parse POST_JOB" >&2; exit 1; }

DIAG_CMD="mkdir -p \"$OUTPUT_ROOT/$POP/diagnostics\" && \"$SDS_PYTHON\" \"$DIAGNOSE\" --input-normalized-tsv \"$OUTPUT_ROOT/$POP/${POP}.normalized.tsv\" --output-prefix \"$OUTPUT_ROOT/$POP/diagnostics/${POP}.gravel_chb_ne100k_grid001\""
DIAG_SUBMIT_OUT="$(bsub -q "$QUEUE" -w "done(${POST_JOB})" -n 1 -R "span[hosts=1]" -J "${POP}_gravel_grid001_diag" -o "$LOGDIR/diag.out" -e "$LOGDIR/diag.err" /bin/bash -lc "$DIAG_CMD" < /dev/null)"
DIAG_JOB="$(printf '%s\n' "$DIAG_SUBMIT_OUT" | sed -n 's/.*<\([0-9]\+\)>.*/\1/p')"
[[ -n "$DIAG_JOB" ]] || { echo "[Error] Could not parse DIAG_JOB" >&2; exit 1; }

printf 'output_root\t%s\n' "$OUTPUT_ROOT"
printf 'post_job\t%s\n' "$POST_JOB"
printf 'diag_job\t%s\n' "$DIAG_JOB"
