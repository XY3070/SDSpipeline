#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

OUT_ROOT="$BASE_DIR/data/processed/sds_output_scn_olddefault_genomewide_hybridpolicy_20260531"
IN_ROOT="$BASE_DIR/data/processed/sds_input_rebuilt_main_contract_20260511"
G_FILE="$BASE_DIR/g_file.txt"
RELAXED_CHRS="1,2,3,4,5,6,7,8,10,11,12,13,14,15,17,18,19,20,21,22"
CAP_CHRS="9,16"
POP="SCN"
QUEUE="normal"
SUBMITTER_PARALLEL="2"
CHUNK_ROWS="20000"
ARRAY_PARALLEL=""
CHUNK_JOB_SLOTS=""
SKIP_FINALIZER="0"
RELAXED_SKIP_BOUNDARY_MISSING_FRACTION="0.35"
CAP_SKIP_BOUNDARY_MISSING_FRACTION="0.35"
CAP_BOUNDARY_MISSING_MODE="cap_to_boundary"
TS="$(date +%Y%m%d_%H%M%S)_$$"
STATE_DIR="$BASE_DIR/data/processed/scn_olddefault_genomewide_hybridqarm_submitter_array_$TS"
RELAXED_CHR_LIST="$STATE_DIR/relaxed_chrs.txt"
CAP_CHR_LIST="$STATE_DIR/cap_chrs.txt"
RESULT_DIR="$STATE_DIR/results"
MANIFEST="$STATE_DIR/manifest.tsv"
LOG_DIR="$OUT_ROOT/$POP/logs"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --out-root) OUT_ROOT="$2"; shift 2 ;;
        --in-root) IN_ROOT="$2"; shift 2 ;;
        --g-file) G_FILE="$2"; shift 2 ;;
        --relaxed-chrs) RELAXED_CHRS="$2"; shift 2 ;;
        --cap-chrs) CAP_CHRS="$2"; shift 2 ;;
        --queue) QUEUE="$2"; shift 2 ;;
        --submitter-parallel) SUBMITTER_PARALLEL="$2"; shift 2 ;;
        --chunk-rows) CHUNK_ROWS="$2"; shift 2 ;;
        --array-parallel) ARRAY_PARALLEL="$2"; shift 2 ;;
        --chunk-job-slots) CHUNK_JOB_SLOTS="$2"; shift 2 ;;
        --relaxed-skip-boundary-missing-fraction) RELAXED_SKIP_BOUNDARY_MISSING_FRACTION="$2"; shift 2 ;;
        --cap-skip-boundary-missing-fraction) CAP_SKIP_BOUNDARY_MISSING_FRACTION="$2"; shift 2 ;;
        --cap-boundary-missing-mode) CAP_BOUNDARY_MISSING_MODE="$2"; shift 2 ;;
        --skip-finalizer) SKIP_FINALIZER="1"; shift ;;
        *) echo "Unknown parameter: $1" >&2; exit 1 ;;
    esac
done

[[ -f "$G_FILE" ]] || { echo "[Error] g_file not found: $G_FILE" >&2; exit 1; }

mkdir -p "$STATE_DIR" "$RESULT_DIR" "$LOG_DIR" "$OUT_ROOT/$POP"
tr ',' '\n' <<< "$RELAXED_CHRS" | sed '/^$/d' > "$RELAXED_CHR_LIST"
tr ',' '\n' <<< "$CAP_CHRS" | sed '/^$/d' > "$CAP_CHR_LIST"

submit_group() {
    local group_name="$1"
    local chr_list="$2"
    local env_skip="$3"
    local env_mode="$4"

    local n
    n="$(wc -l < "$chr_list")"
    if [[ "$n" -le 0 ]]; then
        printf '%s\n' "SKIP"
        return 0
    fi

    local env_prefix="env SDS_SKIP_BOUNDARY_MISSING_FRACTION=\"$env_skip\""
    if [[ -n "$env_mode" ]]; then
        env_prefix="$env_prefix SDS_BOUNDARY_MISSING_MODE=\"$env_mode\""
    fi

    local submit_out
    submit_out="$(bsub -q smp -n 4 -R "span[hosts=1]" \
        -J "SDS_SUBMITTER_SCNOLDHYB_${group_name}[1-${n}]%${SUBMITTER_PARALLEL}" \
        -o "$LOG_DIR/SCN_olddefault_hybrid_${group_name}_submitter_%I.out" \
        -e "$LOG_DIR/SCN_olddefault_hybrid_${group_name}_submitter_%I.err" \
        /bin/bash -lc "$env_prefix bash \"$SCRIPT_DIR/run_submitter_array_entry.sh\" \"$POP\" \"$chr_list\" \"$IN_ROOT\" \"$OUT_ROOT\" \"$G_FILE\" scn_olddefault_hybrid_${group_name}_submitter_array \"$RESULT_DIR\" \"$QUEUE\" \"${CHUNK_ROWS}\" \"${ARRAY_PARALLEL}\" \"${CHUNK_JOB_SLOTS}\"" < /dev/null)"
    local job_id
    job_id="$(sed -n 's/.*<\([0-9]\+\)>.*/\1/p' <<< "$submit_out")"
    [[ -n "$job_id" ]] || { echo "[Error] Failed to parse ${group_name} submitter array job id" >&2; exit 1; }
    printf '%s\n' "$job_id"
}

RELAXED_ARRAY_JOB="$(submit_group relaxed "$RELAXED_CHR_LIST" "$RELAXED_SKIP_BOUNDARY_MISSING_FRACTION" "")"
CAP_ARRAY_JOB="$(submit_group cap "$CAP_CHR_LIST" "$CAP_SKIP_BOUNDARY_MISSING_FRACTION" "$CAP_BOUNDARY_MISSING_MODE")"

FINALIZER_JOB="SKIPPED"
if [[ "$SKIP_FINALIZER" != "1" ]]; then
    FINALIZER_DEPS=()
    if [[ "$RELAXED_ARRAY_JOB" != "SKIP" ]]; then
        FINALIZER_DEPS+=("done(${RELAXED_ARRAY_JOB})")
    fi
    if [[ "$CAP_ARRAY_JOB" != "SKIP" ]]; then
        FINALIZER_DEPS+=("done(${CAP_ARRAY_JOB})")
    fi
    [[ "${#FINALIZER_DEPS[@]}" -gt 0 ]] || { echo "[Error] No chromosome groups requested." >&2; exit 1; }
    FINALIZER_DEP_EXPR="$(printf '%s && ' "${FINALIZER_DEPS[@]}")"
    FINALIZER_DEP_EXPR="${FINALIZER_DEP_EXPR% && }"
    FINALIZER_OUT="$(bsub -q "$QUEUE" -w "$FINALIZER_DEP_EXPR" -n 1 -R "span[hosts=1]" -J "SDS_SUBMITTER_SCNOLDHYB_FINALIZE" -o "$LOG_DIR/SCN_olddefault_hybrid_finalize.out" -e "$LOG_DIR/SCN_olddefault_hybrid_finalize.err" /bin/bash -lc "bash \"$SCRIPT_DIR/submit_postprocess_from_final_jobs.sh\" \"$POP\" \"$RESULT_DIR\" \"$OUT_ROOT\" \"$QUEUE\"" < /dev/null)"
    FINALIZER_JOB="$(sed -n 's/.*<\([0-9]\+\)>.*/\1/p' <<< "$FINALIZER_OUT")"
    [[ -n "$FINALIZER_JOB" ]] || { echo "[Error] Failed to parse finalizer job id" >&2; exit 1; }
fi

printf 'key\tvalue\n' > "$MANIFEST"
printf 'state_dir\t%s\nrelaxed_chr_list\t%s\ncap_chr_list\t%s\nresult_dir\t%s\nout_root\t%s\ng_file\t%s\nqueue\t%s\nsubmitter_parallel\t%s\nchunk_rows\t%s\narray_parallel\t%s\nchunk_job_slots\t%s\nrelaxed_skip_boundary_missing_fraction\t%s\ncap_skip_boundary_missing_fraction\t%s\ncap_boundary_missing_mode\t%s\nskip_finalizer\t%s\nrelaxed_submitter_array_job\t%s\ncap_submitter_array_job\t%s\nfinalizer_job\t%s\n' \
    "$STATE_DIR" "$RELAXED_CHR_LIST" "$CAP_CHR_LIST" "$RESULT_DIR" "$OUT_ROOT" "$G_FILE" "$QUEUE" "$SUBMITTER_PARALLEL" "$CHUNK_ROWS" "$ARRAY_PARALLEL" "$CHUNK_JOB_SLOTS" "$RELAXED_SKIP_BOUNDARY_MISSING_FRACTION" "$CAP_SKIP_BOUNDARY_MISSING_FRACTION" "$CAP_BOUNDARY_MISSING_MODE" "$SKIP_FINALIZER" "$RELAXED_ARRAY_JOB" "$CAP_ARRAY_JOB" "$FINALIZER_JOB" >> "$MANIFEST"

printf 'manifest\t%s\n' "$MANIFEST"
