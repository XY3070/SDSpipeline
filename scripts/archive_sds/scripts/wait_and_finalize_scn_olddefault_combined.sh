#!/bin/bash
set -euo pipefail

if [[ $# -lt 4 ]]; then
    echo "Usage: $0 OUT_ROOT EXPECTED_FINAL_JOB_COUNT COMPARE_OUT_DIR STATE_DIR..." >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

OUT_ROOT="$1"
EXPECTED_COUNT="$2"
COMPARE_OUT_DIR="$3"
shift 3
STATE_DIRS=("$@")

TIMEOUT_SEC="${TIMEOUT_SEC:-172800}"
SLEEP_SEC="${SLEEP_SEC:-60}"

count_final_jobs() {
    local total=0
    local state_dir
    for state_dir in "${STATE_DIRS[@]}"; do
        if [[ -d "$state_dir/results" ]]; then
            local n
            n="$(find "$state_dir/results" -maxdepth 1 -type f -name 'chr*.final_job.tsv' | wc -l | tr -d ' ')"
            total=$((total + n))
        fi
    done
    printf '%s\n' "$total"
}

start_ts="$(date +%s)"
while true; do
    total_final_jobs="$(count_final_jobs)"
    if [[ "$total_final_jobs" -ge "$EXPECTED_COUNT" ]]; then
        break
    fi
    now_ts="$(date +%s)"
    if (( now_ts - start_ts > TIMEOUT_SEC )); then
        echo "[Error] timed out waiting for final_job tsv count; have=$total_final_jobs expected=$EXPECTED_COUNT" >&2
        exit 1
    fi
    sleep "$SLEEP_SEC"
done

TS="$(date +%Y%m%d_%H%M%S)_$$"
COMBINED_DIR="$BASE_DIR/data/processed/scn_olddefault_combined_final_jobs_$TS"
mkdir -p "$COMBINED_DIR"

declare -A seen
for state_dir in "${STATE_DIRS[@]}"; do
    result_dir="$state_dir/results"
    [[ -d "$result_dir" ]] || continue
    while IFS= read -r f; do
        base="$(basename "$f")"
        if [[ -n "${seen[$base]:-}" ]]; then
            echo "[Error] duplicate final_job file name detected: $base" >&2
            echo "  first: ${seen[$base]}" >&2
            echo "  second: $f" >&2
            exit 1
        fi
        seen["$base"]="$f"
        cp "$f" "$COMBINED_DIR/$base"
    done < <(find "$result_dir" -maxdepth 1 -type f -name 'chr*.final_job.tsv' | sort)
done

actual_combined="$(find "$COMBINED_DIR" -maxdepth 1 -type f -name 'chr*.final_job.tsv' | wc -l | tr -d ' ')"
if [[ "$actual_combined" -lt "$EXPECTED_COUNT" ]]; then
    echo "[Error] combined final_job dir incomplete: have=$actual_combined expected=$EXPECTED_COUNT" >&2
    exit 1
fi

printf 'combined_dir\t%s\n' "$COMBINED_DIR"
bash "$SCRIPT_DIR/submit_postprocess_from_final_jobs.sh" SCN "$COMBINED_DIR" "$OUT_ROOT" normal
bash "$SCRIPT_DIR/wait_and_compare_scn_olddefault_genomewide.sh" \
    --root "$OUT_ROOT/SCN" \
    --out-dir "$COMPARE_OUT_DIR"
