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
OUT_ROOT="$BASE_DIR/data/processed/sds_output_olddefault_newinput_20260514"
G_FILE="$BASE_DIR/g_file.txt"
TS="$(date +%Y%m%d_%H%M%S)_$$"
MANIFEST_DIR="$BASE_DIR/data/processed/olddefault_newinput_chrlist_${TS}"
MANIFEST="$MANIFEST_DIR/manifest.tsv"
LOG_DIR="$OUT_ROOT/logs"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pop) POP="$2"; shift 2 ;;
        --chrs) CHR_LIST="$2"; shift 2 ;;
        --queue) QUEUE="$2"; shift 2 ;;
        *) echo "Unknown parameter: $1" >&2; exit 1 ;;
    esac
done

[[ -n "$POP" && -n "$CHR_LIST" ]] || {
    echo "Usage: $0 --pop POP --chrs 1,2,3 [--queue normal]" >&2
    exit 1
}

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
            --job-group "olddefault_newinput_${POP}_chr${chr}" 2>&1)"; then
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
    submit_out="$(bsub -q "$QUEUE" -w "$dep_expr" -n 1 -R "span[hosts=1]" -J "sds_${POP}_olddefault_newinput_postprocess" -o "$LOG_DIR/${POP}_olddefault_newinput_postprocess.out" -e "$LOG_DIR/${POP}_olddefault_newinput_postprocess.err" /bin/bash -lc "$cmd" < /dev/null)"
    submit_and_parse_job_id "$submit_out"
}

submit_diagnostics() {
    local dep_job="$1"
    local output_prefix="$OUT_ROOT/$POP/diagnostics/${POP}.olddefault_newinput"
    local cmd submit_out
    cmd="mkdir -p \"$OUT_ROOT/$POP/diagnostics\" && \"$SDS_PYTHON\" \"$DIAGNOSE\" --input-normalized-tsv \"$OUT_ROOT/$POP/${POP}.normalized.tsv\" --output-prefix \"$output_prefix\""
    submit_out="$(bsub -q "$QUEUE" -w "done(${dep_job})" -n 1 -R "span[hosts=1]" -J "sds_${POP}_olddefault_newinput_diag" -o "$LOG_DIR/${POP}_olddefault_newinput_diag.out" -e "$LOG_DIR/${POP}_olddefault_newinput_diag.err" /bin/bash -lc "$cmd" < /dev/null)"
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
record_manifest demographic_model "SDS-MODEL-001"
record_manifest present_day_ne0 "implicit / unknown from historical file provenance"

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
