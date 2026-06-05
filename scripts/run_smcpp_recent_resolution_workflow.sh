#!/bin/bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 POP" >&2
    exit 1
fi

POP="$(printf '%s' "$1" | tr '[:lower:]' '[:upper:]')"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
SDS_PYTHON="${SDS_PYTHON:-/data/home/grp-wangyf/intern/miniforge3/envs/sds/bin/python}"

REFERENCE_BASE="${REFERENCE_BASE:-${POP}_fine}"
COARSE_BASE="${COARSE_BASE:-$POP}"
WINDOW_SIZE="${WINDOW_SIZE:-10}"
KNOTS="${KNOTS:-16}"
TIMEPOINT_END="${TIMEPOINT_END:-100000}"
VARIANTS="${VARIANTS:-tp20_w10_k16:20,tp50_w10_k16:50}"
WAIT_FOR_COMPLETION="${WAIT_FOR_COMPLETION:-0}"
AUTO_SUBMIT_SCN="${AUTO_SUBMIT_SCN:-0}"
FOLLOWUP_POP="${FOLLOWUP_POP:-}"

if [[ -z "$FOLLOWUP_POP" && "$AUTO_SUBMIT_SCN" == "1" && "$POP" == "NCN" ]]; then
    FOLLOWUP_POP="SCN"
fi

IFS=',' read -r -a VARIANT_ITEMS <<<"$VARIANTS"
JOB_IDS=()
CANDIDATE_BASES=()
FIRST_TAG=""
FIRST_START=""

wait_for_job() {
    local job_id="$1"
    while bjobs "$job_id" >/dev/null 2>&1; do
        sleep 60
    done
    if ! bhist -l "$job_id" | grep -q "Done successfully"; then
        echo "[smc++] job $job_id did not complete successfully" >&2
        return 1
    fi
}

for item in "${VARIANT_ITEMS[@]}"; do
    tag="${item%%:*}"
    start="${item##*:}"
    if [[ -z "$FIRST_TAG" ]]; then
        FIRST_TAG="$tag"
        FIRST_START="$start"
    fi
    CANDIDATE_BASES+=("${POP}_${tag}")
    SUBMIT_OUTPUT="$("$SCRIPT_DIR/submit_smcpp_recent_job.sh" "$POP" "$tag" "$WINDOW_SIZE" "$KNOTS" "$start" "$TIMEPOINT_END")"
    printf '%s\n' "$SUBMIT_OUTPUT"
    JOB_ID="$(printf '%s\n' "$SUBMIT_OUTPUT" | sed -n 's/.*job_id=\([0-9][0-9]*\).*/\1/p' | tail -n 1)"
    if [[ -n "$JOB_ID" ]]; then
        JOB_IDS+=("$JOB_ID")
    fi
done

if [[ "$WAIT_FOR_COMPLETION" != "1" ]]; then
    printf '[smc++] submitted candidate bases: %s\n' "${CANDIDATE_BASES[*]}"
    printf '[smc++] job ids: %s\n' "${JOB_IDS[*]:-none}"
    exit 0
fi

for job_id in "${JOB_IDS[@]}"; do
    wait_for_job "$job_id"
done

REPORT_OUTPUT="$(
    "$SDS_PYTHON" \
        "$PROJECT_ROOT/benchmark/demography/evaluate_smcpp_recent_resolution.py" \
        --pop "$POP" \
        --coarse-base "$COARSE_BASE" \
        --reference-base "$REFERENCE_BASE" \
        --candidate-bases "${CANDIDATE_BASES[@]}"
)"
printf '%s\n' "$REPORT_OUTPUT"

REPORT_JSON="$PROJECT_ROOT/benchmark/demography/$POP/smcpp/recent_sensitivity/${POP}_recent_resolution_report.json"
"$SDS_PYTHON" "$PROJECT_ROOT/benchmark/demography/plot_smcpp_recent_sensitivity.py" --report "$REPORT_JSON"

if [[ -n "$FOLLOWUP_POP" ]]; then
    SHOULD_SUBMIT="$(
        "$SDS_PYTHON" - <<PY
import json
from pathlib import Path
payload = json.loads(Path("$REPORT_JSON").read_text())
print("1" if payload["decision"]["passes"] else "0")
PY
    )"
    if [[ "$SHOULD_SUBMIT" == "1" ]]; then
        "$SCRIPT_DIR/submit_smcpp_recent_job.sh" \
            "$FOLLOWUP_POP" \
            "$FIRST_TAG" \
            "$WINDOW_SIZE" \
            "$KNOTS" \
            "$FIRST_START" \
            "$TIMEPOINT_END"
    else
        printf '[smc++] %s decision did not pass; skipping follow-up submit for %s\n' "$POP" "$FOLLOWUP_POP"
    fi
fi
