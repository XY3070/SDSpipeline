#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

SDS_PYTHON="/data/home/grp-wangyf/intern/miniforge3/envs/sds/bin/python"
RUN_INPUT="$SCRIPT_DIR/run_sds_input.sh"
RUN_COMPUTE="$SCRIPT_DIR/run_sds_compute.sh"
POSTPROCESS="$SCRIPT_DIR/postprocess_sds_results.py"
DIAGNOSE="$SCRIPT_DIR/diagnose_sds_scan.py"
G_FILE="$BASE_DIR/tmp/gravel_chb_gamma_ne100k_20260502/gravel_chb_present100000.g_file.txt"

QUEUE="normal"
INPUT_ARRAY_CAP="4"
COMPUTE_ARRAY_CAP="6"
SMOKE_CHR="10"
INPUT_ROOT="$BASE_DIR/data/processed/sds_input_rebuilt_main_contract_20260511"
SMOKE_OUT_ROOT="$BASE_DIR/data/processed/sds_output_gravel_chb_ne100k_smoke_20260511"
OUT_ROOT="$BASE_DIR/data/processed/sds_output_gravel_chb_ne100k_newinput_20260511"
TS="$(date +%Y%m%d_%H%M%S)"
MANIFEST_DIR="$BASE_DIR/data/processed/gravel_chb_ne100k_newinput_launch_$TS"
MANIFEST="$MANIFEST_DIR/launch_manifest.tsv"
INPUT_LOG_DIR="$INPUT_ROOT/logs"
SMOKE_LOG_DIR="$SMOKE_OUT_ROOT/logs"
OUT_LOG_DIR="$OUT_ROOT/logs"

mkdir -p "$MANIFEST_DIR" "$INPUT_LOG_DIR" "$SMOKE_LOG_DIR" "$OUT_LOG_DIR"

[[ -f "$G_FILE" ]] || { echo "[Error] Missing Gravel_CHB g_file: $G_FILE" >&2; exit 1; }

submit_and_parse_job_id() {
    local submit_out="$1"
    local job_id
    job_id="$(sed -n 's/.*<\([0-9]\+\)>.*/\1/p' <<< "$submit_out")"
    [[ -n "$job_id" ]] || {
        echo "[Error] Failed to parse bsub output: $submit_out" >&2
        exit 1
    }
    printf '%s\n' "$job_id"
}

record_manifest() {
    local key="$1"
    local value="$2"
    printf '%s\t%s\n' "$key" "$value" >> "$MANIFEST"
}

printf 'key\tvalue\n' > "$MANIFEST"
record_manifest ts "$TS"
record_manifest queue "$QUEUE"
record_manifest smoke_chr "$SMOKE_CHR"
record_manifest g_file "$G_FILE"
record_manifest demographic_model "SDS-MODEL-003"
record_manifest present_day_ne0 "100000"
record_manifest input_root "$INPUT_ROOT"
record_manifest smoke_out_root "$SMOKE_OUT_ROOT"
record_manifest output_root "$OUT_ROOT"

submit_smoke_input() {
    local pop="$1"
    local cmd
    local submit_out
    cmd="bash \"$RUN_INPUT\" --pop $pop --chr $SMOKE_CHR --out-root \"$INPUT_ROOT\" --test --force"
    submit_out="$(bsub -q "$QUEUE" -n 1 -R "span[hosts=1]" -J "sds_${pop}_newinput_smoke_chr${SMOKE_CHR}" -o "$INPUT_LOG_DIR/${pop}_smoke_input_chr${SMOKE_CHR}.out" -e "$INPUT_LOG_DIR/${pop}_smoke_input_chr${SMOKE_CHR}.err" /bin/bash -lc "$cmd" < /dev/null)"
    submit_and_parse_job_id "$submit_out"
}

submit_smoke_compute() {
    local pop="$1"
    local dep_job="$2"
    local cmd
    local submit_out
    cmd="bash \"$RUN_COMPUTE\" --pop $pop --chr $SMOKE_CHR --test --in-root \"$INPUT_ROOT/test\" --out-root \"$SMOKE_OUT_ROOT\" --g-file \"$G_FILE\" --force"
    submit_out="$(bsub -q "$QUEUE" -w "done(${dep_job})" -n 2 -R "span[hosts=1]" -J "sds_${pop}_gravel_chb_ne100k_smoke_chr${SMOKE_CHR}" -o "$SMOKE_LOG_DIR/${pop}_smoke_compute_chr${SMOKE_CHR}.out" -e "$SMOKE_LOG_DIR/${pop}_smoke_compute_chr${SMOKE_CHR}.err" /bin/bash -lc "$cmd" < /dev/null)"
    submit_and_parse_job_id "$submit_out"
}

submit_full_input_array() {
    local pop="$1"
    local dep_expr="$2"
    local cmd
    local submit_out
    cmd="bash \"$RUN_INPUT\" --pop $pop --chr \${LSB_JOBINDEX} --out-root \"$INPUT_ROOT\" --force"
    submit_out="$(bsub -q "$QUEUE" -w "$dep_expr" -n 1 -R "span[hosts=1]" -J "sds_${pop}_newinput_full[1-22]%${INPUT_ARRAY_CAP}" -o "$INPUT_LOG_DIR/${pop}_full_input_%I.out" -e "$INPUT_LOG_DIR/${pop}_full_input_%I.err" /bin/bash -lc "$cmd" < /dev/null)"
    submit_and_parse_job_id "$submit_out"
}

submit_full_compute_array() {
    local pop="$1"
    local dep_job="$2"
    local cmd
    local submit_out
    cmd="bash \"$RUN_COMPUTE\" --pop $pop --chr \${LSB_JOBINDEX} --in-root \"$INPUT_ROOT\" --out-root \"$OUT_ROOT\" --g-file \"$G_FILE\" --force"
    submit_out="$(bsub -q "$QUEUE" -w "done(${dep_job})" -n 2 -R "span[hosts=1]" -J "sds_${pop}_gravel_chb_ne100k_newinput[1-22]%${COMPUTE_ARRAY_CAP}" -o "$OUT_LOG_DIR/${pop}_compute_%I.out" -e "$OUT_LOG_DIR/${pop}_compute_%I.err" /bin/bash -lc "$cmd" < /dev/null)"
    submit_and_parse_job_id "$submit_out"
}

submit_postprocess() {
    local pop="$1"
    local dep_job="$2"
    local cmd
    local submit_out
    cmd="cd \"$BASE_DIR\" && \"$SDS_PYTHON\" \"$POSTPROCESS\" --input-dir \"$OUT_ROOT/$pop\" --pop $pop"
    submit_out="$(bsub -q "$QUEUE" -w "done(${dep_job})" -n 1 -R "span[hosts=1]" -J "sds_${pop}_gravel_chb_ne100k_newinput_postprocess" -o "$OUT_LOG_DIR/${pop}_postprocess.out" -e "$OUT_LOG_DIR/${pop}_postprocess.err" /bin/bash -lc "$cmd" < /dev/null)"
    submit_and_parse_job_id "$submit_out"
}

submit_diagnostics() {
    local pop="$1"
    local dep_job="$2"
    local output_prefix="$OUT_ROOT/$pop/diagnostics/${pop}.gravel_chb_ne100k_newinput"
    local cmd
    local submit_out
    cmd="mkdir -p \"$OUT_ROOT/$pop/diagnostics\" && \"$SDS_PYTHON\" \"$DIAGNOSE\" --input-normalized-tsv \"$OUT_ROOT/$pop/${pop}.normalized.tsv\" --output-prefix \"$output_prefix\""
    submit_out="$(bsub -q "$QUEUE" -w "done(${dep_job})" -n 1 -R "span[hosts=1]" -J "sds_${pop}_gravel_chb_ne100k_newinput_diag" -o "$OUT_LOG_DIR/${pop}_diag.out" -e "$OUT_LOG_DIR/${pop}_diag.err" /bin/bash -lc "$cmd" < /dev/null)"
    submit_and_parse_job_id "$submit_out"
}

NCN_SMOKE_INPUT_JOB="$(submit_smoke_input NCN)"
SCN_SMOKE_INPUT_JOB="$(submit_smoke_input SCN)"
record_manifest ncn_smoke_input_job "$NCN_SMOKE_INPUT_JOB"
record_manifest scn_smoke_input_job "$SCN_SMOKE_INPUT_JOB"

NCN_SMOKE_COMPUTE_JOB="$(submit_smoke_compute NCN "$NCN_SMOKE_INPUT_JOB")"
SCN_SMOKE_COMPUTE_JOB="$(submit_smoke_compute SCN "$SCN_SMOKE_INPUT_JOB")"
record_manifest ncn_smoke_compute_job "$NCN_SMOKE_COMPUTE_JOB"
record_manifest scn_smoke_compute_job "$SCN_SMOKE_COMPUTE_JOB"

FULL_INPUT_DEP="done(${NCN_SMOKE_COMPUTE_JOB}) && done(${SCN_SMOKE_COMPUTE_JOB})"
NCN_FULL_INPUT_ARRAY_JOB="$(submit_full_input_array NCN "$FULL_INPUT_DEP")"
SCN_FULL_INPUT_ARRAY_JOB="$(submit_full_input_array SCN "$FULL_INPUT_DEP")"
record_manifest ncn_full_input_array_job "$NCN_FULL_INPUT_ARRAY_JOB"
record_manifest scn_full_input_array_job "$SCN_FULL_INPUT_ARRAY_JOB"

NCN_FULL_COMPUTE_ARRAY_JOB="$(submit_full_compute_array NCN "$NCN_FULL_INPUT_ARRAY_JOB")"
SCN_FULL_COMPUTE_ARRAY_JOB="$(submit_full_compute_array SCN "$SCN_FULL_INPUT_ARRAY_JOB")"
record_manifest ncn_full_compute_array_job "$NCN_FULL_COMPUTE_ARRAY_JOB"
record_manifest scn_full_compute_array_job "$SCN_FULL_COMPUTE_ARRAY_JOB"

NCN_POSTPROCESS_JOB="$(submit_postprocess NCN "$NCN_FULL_COMPUTE_ARRAY_JOB")"
SCN_POSTPROCESS_JOB="$(submit_postprocess SCN "$SCN_FULL_COMPUTE_ARRAY_JOB")"
record_manifest ncn_postprocess_job "$NCN_POSTPROCESS_JOB"
record_manifest scn_postprocess_job "$SCN_POSTPROCESS_JOB"

NCN_DIAG_JOB="$(submit_diagnostics NCN "$NCN_POSTPROCESS_JOB")"
SCN_DIAG_JOB="$(submit_diagnostics SCN "$SCN_POSTPROCESS_JOB")"
record_manifest ncn_diag_job "$NCN_DIAG_JOB"
record_manifest scn_diag_job "$SCN_DIAG_JOB"

printf 'manifest\t%s\n' "$MANIFEST"
printf 'status_script\t%s\n' "$SCRIPT_DIR/check_gravel_chb_ne100k_newinput_status.sh $MANIFEST\n"
