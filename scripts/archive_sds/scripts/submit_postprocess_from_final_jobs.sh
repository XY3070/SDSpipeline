#!/bin/bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 POP RESULT_DIR OUT_ROOT QUEUE" >&2
    exit 1
fi

POP="$1"
RESULT_DIR="$2"
OUT_ROOT="$3"
QUEUE="$4"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SDS_PYTHON="/data/home/grp-wangyf/intern/miniforge3/envs/sds/bin/python"
POSTPROCESS="$SCRIPT_DIR/postprocess_sds_results.py"
DIAGNOSE="$SCRIPT_DIR/diagnose_sds_scan.py"
LOG_DIR="$OUT_ROOT/logs"
mkdir -p "$LOG_DIR"

mapfile -t CHROMS < <(awk -F '\t' 'FNR>1{print $1}' "$RESULT_DIR"/chr*.final_job.tsv 2>/dev/null | sed '/^$/d' | sort -u)
[[ "${#CHROMS[@]}" -gt 0 ]] || { echo "[Error] No final jobs found under $RESULT_DIR" >&2; exit 1; }

WAIT_TIMEOUT_SEC="${POSTPROCESS_WAIT_TIMEOUT_SEC:-14400}"
WAIT_INTERVAL_SEC="${POSTPROCESS_WAIT_INTERVAL_SEC:-30}"
deadline=$(( $(date +%s) + WAIT_TIMEOUT_SEC ))

declare -A REQUIRE_P_ARM=()
declare -A REQUIRE_Q_ARM=()
for chr in "${CHROMS[@]}"; do
    submitter_out="$RESULT_DIR/chr${chr}.submitter.out"
    p_job=""
    q_job=""
    if [[ -f "$submitter_out" ]]; then
        p_job="$(awk -F '\t' '$1=="P_ARRAY_JOB"{print $2; exit}' "$submitter_out")"
        q_job="$(awk -F '\t' '$1=="Q_ARRAY_JOB"{print $2; exit}' "$submitter_out")"
    fi
    if [[ "$p_job" != "SKIP" ]]; then
        REQUIRE_P_ARM["$chr"]=1
    fi
    if [[ "$q_job" != "SKIP" ]]; then
        REQUIRE_Q_ARM["$chr"]=1
    fi
done

missing_files=()
while true; do
    missing_files=()
    for chr in "${CHROMS[@]}"; do
        if [[ -n "${REQUIRE_P_ARM[$chr]:-}" ]]; then
            path="$OUT_ROOT/$POP/chr${chr}_p.sds.tsv"
            [[ -s "$path" ]] || missing_files+=("$path")
        fi
        if [[ -n "${REQUIRE_Q_ARM[$chr]:-}" ]]; then
            path="$OUT_ROOT/$POP/chr${chr}_q.sds.tsv"
            [[ -s "$path" ]] || missing_files+=("$path")
        fi
    done

    if [[ "${#missing_files[@]}" -eq 0 ]]; then
        break
    fi

    now="$(date +%s)"
    if (( now >= deadline )); then
        printf '[Error] Timed out waiting for postprocess inputs under %s\n' "$OUT_ROOT/$POP" >&2
        printf '[Error] Missing %d chr*_p/q.sds.tsv files:\n' "${#missing_files[@]}" >&2
        printf '  %s\n' "${missing_files[@]}" >&2
        exit 1
    fi

    sleep "$WAIT_INTERVAL_SEC"
done

POST_CMD="cd \"$BASE_DIR\" && \"$SDS_PYTHON\" \"$POSTPROCESS\" --input-dir \"$OUT_ROOT/$POP\" --pop $POP"
POST_OUT="$(bsub -q "$QUEUE" -n 1 -R "span[hosts=1]" -J "sds_${POP}_gravel_chb_ne100k_submitter_array_postprocess" -o "$LOG_DIR/${POP}_submitter_array_postprocess.out" -e "$LOG_DIR/${POP}_submitter_array_postprocess.err" /bin/bash -lc "$POST_CMD" < /dev/null)"
POST_JOB="$(sed -n 's/.*<\([0-9]\+\)>.*/\1/p' <<< "$POST_OUT")"
[[ -n "$POST_JOB" ]] || { echo "[Error] Failed to parse postprocess job id" >&2; exit 1; }

DIAG_CMD="mkdir -p \"$OUT_ROOT/$POP/diagnostics\" && \"$SDS_PYTHON\" \"$DIAGNOSE\" --input-normalized-tsv \"$OUT_ROOT/$POP/${POP}.normalized.tsv\" --output-prefix \"$OUT_ROOT/$POP/diagnostics/${POP}.gravel_chb_ne100k_submitter_array\""
DIAG_OUT="$(bsub -q "$QUEUE" -w "done(${POST_JOB})" -n 1 -R "span[hosts=1]" -J "sds_${POP}_gravel_chb_ne100k_submitter_array_diag" -o "$LOG_DIR/${POP}_submitter_array_diag.out" -e "$LOG_DIR/${POP}_submitter_array_diag.err" /bin/bash -lc "$DIAG_CMD" < /dev/null)"
DIAG_JOB="$(sed -n 's/.*<\([0-9]\+\)>.*/\1/p' <<< "$DIAG_OUT")"
[[ -n "$DIAG_JOB" ]] || { echo "[Error] Failed to parse diagnostics job id" >&2; exit 1; }

printf 'postprocess_job\t%s\ndiag_job\t%s\n' "$POST_JOB" "$DIAG_JOB"
