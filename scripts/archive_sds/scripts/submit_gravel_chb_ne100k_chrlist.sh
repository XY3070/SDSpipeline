#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

SDS_PYTHON="/data/home/grp-wangyf/intern/miniforge3/envs/sds/bin/python"
SUBMIT_CHUNKED="$SCRIPT_DIR/submit_sds_compute_chunked_chr.sh"
POSTPROCESS="$SCRIPT_DIR/postprocess_sds_results.py"
DIAGNOSE="$SCRIPT_DIR/diagnose_sds_scan.py"

POP=""
CHR_LIST=""
QUEUE="normal"
INPUT_ROOT="$BASE_DIR/data/processed/sds_input_rebuilt_main_contract_20260511"
OUT_ROOT="$BASE_DIR/data/processed/sds_output_gravel_chb_ne100k_newinput_20260511"
G_FILE="$BASE_DIR/tmp/gravel_chb_gamma_ne100k_20260502/gravel_chb_present100000.g_file.txt"
MODEL_ID="SDS-MODEL-003"
NE0_LABEL="100000"
RUN_LABEL="gravel_chb_ne100k"
TS="$(date +%Y%m%d_%H%M%S)_$$"
MANIFEST_DIR="$BASE_DIR/data/processed/${RUN_LABEL}_chrlist_${TS}"
MANIFEST="$MANIFEST_DIR/manifest.tsv"
LOG_DIR="$OUT_ROOT/logs"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pop) POP="$2"; shift 2 ;;
        --chrs) CHR_LIST="$2"; shift 2 ;;
        --queue) QUEUE="$2"; shift 2 ;;
        --in-root) INPUT_ROOT="$2"; shift 2 ;;
        --out-root) OUT_ROOT="$2"; shift 2 ;;
        --g-file) G_FILE="$2"; shift 2 ;;
        --model-id) MODEL_ID="$2"; shift 2 ;;
        --ne0) NE0_LABEL="$2"; shift 2 ;;
        --run-label) RUN_LABEL="$2"; shift 2 ;;
        *) echo "Unknown parameter: $1" >&2; exit 1 ;;
    esac
done

[[ -n "$POP" && -n "$CHR_LIST" ]] || {
    echo "Usage: $0 --pop POP --chrs 1,2,3 [--queue normal] [--in-root DIR] [--out-root DIR] [--g-file FILE] [--model-id ID] [--ne0 VALUE] [--run-label LABEL]" >&2
    exit 1
}

MANIFEST_DIR="$BASE_DIR/data/processed/${RUN_LABEL}_chrlist_${TS}"
MANIFEST="$MANIFEST_DIR/manifest.tsv"
LOG_DIR="$OUT_ROOT/logs"

mkdir -p "$MANIFEST_DIR" "$LOG_DIR"

record_manifest() {
    printf '%s\t%s\n' "$1" "$2" >> "$MANIFEST"
}

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

submit_one_chr() {
    local chr="$1"
    local output=""
    local attempt
    for attempt in 1 2 3 4 5; do
        if output="$(bash "$SUBMIT_CHUNKED" \
            --pop "$POP" \
            --chr "$chr" \
            --queue "$QUEUE" \
            --in-root "$INPUT_ROOT" \
            --out-root "$OUT_ROOT" \
            --g-file "$G_FILE" \
            --job-group "${RUN_LABEL}_${POP}_chr${chr}_manual" 2>&1)"; then
            printf '%s\n' "$output"
            return 0
        fi
        sleep 5
    done
    printf '%s\n' "$output" >&2
    return 1
}

submit_postprocess() {
    local dep_expr="$1"
    local cmd submit_out
    cmd="cd \"$BASE_DIR\" && \"$SDS_PYTHON\" \"$POSTPROCESS\" --input-dir \"$OUT_ROOT/$POP\" --pop $POP"
    submit_out="$(bsub -q "$QUEUE" -w "$dep_expr" -n 1 -R "span[hosts=1]" -J "sds_${POP}_${RUN_LABEL}_chrlist_postprocess" -o "$LOG_DIR/${POP}_${RUN_LABEL}_chrlist_postprocess.out" -e "$LOG_DIR/${POP}_${RUN_LABEL}_chrlist_postprocess.err" /bin/bash -lc "$cmd" < /dev/null)"
    submit_and_parse_job_id "$submit_out"
}

submit_diagnostics() {
    local dep_job="$1"
    local output_prefix="$OUT_ROOT/$POP/diagnostics/${POP}.${RUN_LABEL}.chrlist"
    local cmd submit_out
    cmd="mkdir -p \"$OUT_ROOT/$POP/diagnostics\" && \"$SDS_PYTHON\" \"$DIAGNOSE\" --input-normalized-tsv \"$OUT_ROOT/$POP/${POP}.normalized.tsv\" --output-prefix \"$output_prefix\""
    submit_out="$(bsub -q "$QUEUE" -w "done(${dep_job})" -n 1 -R "span[hosts=1]" -J "sds_${POP}_${RUN_LABEL}_chrlist_diag" -o "$LOG_DIR/${POP}_${RUN_LABEL}_chrlist_diag.out" -e "$LOG_DIR/${POP}_${RUN_LABEL}_chrlist_diag.err" /bin/bash -lc "$cmd" < /dev/null)"
    submit_and_parse_job_id "$submit_out"
}

printf 'key\tvalue\n' > "$MANIFEST"
record_manifest ts "$TS"
record_manifest pop "$POP"
record_manifest chr_list "$CHR_LIST"
record_manifest queue "$QUEUE"
record_manifest input_root "$INPUT_ROOT"
record_manifest output_root "$OUT_ROOT"
record_manifest g_file "$G_FILE"
record_manifest demographic_model "$MODEL_ID"
record_manifest present_day_ne0 "$NE0_LABEL"
record_manifest run_label "$RUN_LABEL"

IFS=',' read -r -a CHRS <<< "$CHR_LIST"
final_jobs=()
for chr in "${CHRS[@]}"; do
    submit_output="$(submit_one_chr "$chr")"
    final_job="$(printf '%s\n' "$submit_output" | awk -F '\t' '$1=="FINAL_JOB"{print $2; exit}')"
    [[ -n "$final_job" ]] || {
        echo "[Error] Could not parse FINAL_JOB for ${POP} chr${chr}" >&2
        printf '%s\n' "$submit_output" >&2
        exit 1
    }
    final_jobs+=("$final_job")
    record_manifest "chr${chr}_final_job" "$final_job"
    sleep 2
done

dep_expr="$(printf 'done(%s) && ' "${final_jobs[@]}")"
dep_expr="${dep_expr% && }"
post_job="$(submit_postprocess "$dep_expr")"
diag_job="$(submit_diagnostics "$post_job")"
record_manifest postprocess_job "$post_job"
record_manifest diag_job "$diag_job"

printf 'manifest\t%s\n' "$MANIFEST"
