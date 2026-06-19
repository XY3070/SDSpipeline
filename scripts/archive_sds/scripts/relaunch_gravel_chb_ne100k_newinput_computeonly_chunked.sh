#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

SDS_PYTHON="/data/home/grp-wangyf/intern/miniforge3/envs/sds/bin/python"
SUBMIT_CHUNKED="$SCRIPT_DIR/submit_sds_compute_chunked_chr.sh"
POSTPROCESS="$SCRIPT_DIR/postprocess_sds_results.py"
DIAGNOSE="$SCRIPT_DIR/diagnose_sds_scan.py"

QUEUE="normal"
INPUT_ROOT="$BASE_DIR/data/processed/sds_input_rebuilt_main_contract_20260511"
OUT_ROOT="$BASE_DIR/data/processed/sds_output_gravel_chb_ne100k_newinput_20260511"
G_FILE="$BASE_DIR/tmp/gravel_chb_gamma_ne100k_20260502/gravel_chb_present100000.g_file.txt"
POPS_RAW="${POPS:-NCN SCN}"
TS="$(date +%Y%m%d_%H%M%S)"
MANIFEST_DIR="$BASE_DIR/data/processed/gravel_chb_ne100k_newinput_relaunch_$TS"
MANIFEST="$MANIFEST_DIR/relaunch_manifest.tsv"
LOG_DIR="$OUT_ROOT/logs"

mkdir -p "$MANIFEST_DIR" "$LOG_DIR"
[[ -x "$SUBMIT_CHUNKED" ]] || { echo "[Error] Missing submitter: $SUBMIT_CHUNKED" >&2; exit 1; }
[[ -f "$G_FILE" ]] || { echo "[Error] Missing g_file: $G_FILE" >&2; exit 1; }

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

collect_missing_chrs() {
    local pop="$1"
    local outdir="$OUT_ROOT/$pop"
    local chr
    local missing=()
    for chr in $(seq 1 22); do
        if [[ ! -f "$outdir/chr${chr}.sds.tsv" ]]; then
            missing+=("$chr")
        fi
    done
    printf '%s\n' "${missing[@]}"
}

submit_compute_for_chr() {
    local pop="$1"
    local chr="$2"
    local output
    local attempt
    for attempt in 1 2 3; do
        if output="$(bash "$SUBMIT_CHUNKED" \
            --pop "$pop" \
            --chr "$chr" \
            --queue "$QUEUE" \
            --in-root "$INPUT_ROOT" \
            --out-root "$OUT_ROOT" \
            --g-file "$G_FILE" \
            --job-group "gravel_chb_ne100k_newinput_fast_chr${chr}" 2>&1)"; then
            printf '%s\n' "$output"
            return 0
        fi
        sleep 5
    done
    printf '%s\n' "$output" >&2
    return 1
}

final_job_from_submit_output() {
    awk -F '\t' '$1=="FINAL_JOB"{print $2; exit}'
}

submit_postprocess() {
    local pop="$1"
    local dep_expr="$2"
    local cmd
    local submit_out
    cmd="cd \"$BASE_DIR\" && \"$SDS_PYTHON\" \"$POSTPROCESS\" --input-dir \"$OUT_ROOT/$pop\" --pop $pop"
    submit_out="$(bsub -q "$QUEUE" -w "$dep_expr" -n 1 -R "span[hosts=1]" -J "sds_${pop}_gravel_chb_ne100k_newinput_postprocess_rerun" -o "$LOG_DIR/${pop}_postprocess_rerun.out" -e "$LOG_DIR/${pop}_postprocess_rerun.err" /bin/bash -lc "$cmd" < /dev/null)"
    submit_and_parse_job_id "$submit_out"
}

submit_diagnostics() {
    local pop="$1"
    local dep_job="$2"
    local output_prefix="$OUT_ROOT/$pop/diagnostics/${pop}.gravel_chb_ne100k_newinput_rerun"
    local cmd
    local submit_out
    cmd="mkdir -p \"$OUT_ROOT/$pop/diagnostics\" && \"$SDS_PYTHON\" \"$DIAGNOSE\" --input-normalized-tsv \"$OUT_ROOT/$pop/${pop}.normalized.tsv\" --output-prefix \"$output_prefix\""
    submit_out="$(bsub -q "$QUEUE" -w "done(${dep_job})" -n 1 -R "span[hosts=1]" -J "sds_${pop}_gravel_chb_ne100k_newinput_diag_rerun" -o "$LOG_DIR/${pop}_diag_rerun.out" -e "$LOG_DIR/${pop}_diag_rerun.err" /bin/bash -lc "$cmd" < /dev/null)"
    submit_and_parse_job_id "$submit_out"
}

printf 'key\tvalue\n' > "$MANIFEST"
record_manifest ts "$TS"
record_manifest queue "$QUEUE"
record_manifest input_root "$INPUT_ROOT"
record_manifest output_root "$OUT_ROOT"
record_manifest g_file "$G_FILE"
record_manifest demographic_model "SDS-MODEL-003"
record_manifest present_day_ne0 "100000"
record_manifest relaunch_mode "compute_only_missing_chromosomes"
record_manifest pops "$POPS_RAW"

for pop in $POPS_RAW; do
    mapfile -t missing_chrs < <(collect_missing_chrs "$pop")
    if [[ "${#missing_chrs[@]}" -eq 0 ]]; then
        record_manifest "${pop}_missing_chrs" "none"
        continue
    fi

    record_manifest "${pop}_missing_chrs" "$(printf '%s,' "${missing_chrs[@]}" | sed 's/,$//')"
    final_jobs=()
    for chr in "${missing_chrs[@]}"; do
        submit_output="$(submit_compute_for_chr "$pop" "$chr")"
        final_job="$(printf '%s\n' "$submit_output" | final_job_from_submit_output)"
        [[ -n "$final_job" ]] || {
            echo "[Error] Could not parse FINAL_JOB for ${pop} chr${chr}" >&2
            printf '%s\n' "$submit_output" >&2
            exit 1
        }
        final_jobs+=("$final_job")
        record_manifest "${pop}_chr${chr}_final_job" "$final_job"
    done

    dep_expr="$(printf 'done(%s) && ' "${final_jobs[@]}" | sed 's/ && $//')"
    post_job="$(submit_postprocess "$pop" "$dep_expr")"
    diag_job="$(submit_diagnostics "$pop" "$post_job")"
    record_manifest "${pop}_postprocess_job" "$post_job"
    record_manifest "${pop}_diag_job" "$diag_job"
done

printf 'manifest\t%s\n' "$MANIFEST"
