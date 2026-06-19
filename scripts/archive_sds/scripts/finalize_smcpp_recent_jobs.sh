#!/bin/bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 POP JOB_ID [JOB_ID ...]" >&2
    exit 1
fi

POP="$(printf '%s' "$1" | tr '[:lower:]' '[:upper:]')"
shift
JOB_IDS=("$@")

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
SDS_PYTHON="${SDS_PYTHON:-/data/home/grp-wangyf/intern/miniforge3/envs/sds/bin/python}"
COARSE_BASE="${COARSE_BASE:-$POP}"
REFERENCE_BASE="${REFERENCE_BASE:-${POP}_fine}"
CANDIDATE_BASES_CSV="${CANDIDATE_BASES_CSV:-${POP}_tp20_w10_k16,${POP}_tp50_w10_k16}"
FOLLOWUP_POP="${FOLLOWUP_POP:-}"
FOLLOWUP_TAG="${FOLLOWUP_TAG:-tp20_w10_k16}"
FOLLOWUP_WINDOW_SIZE="${FOLLOWUP_WINDOW_SIZE:-10}"
FOLLOWUP_KNOTS="${FOLLOWUP_KNOTS:-16}"
FOLLOWUP_TIMEPOINT_START="${FOLLOWUP_TIMEPOINT_START:-20}"
FOLLOWUP_TIMEPOINT_END="${FOLLOWUP_TIMEPOINT_END:-100000}"

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

for job_id in "${JOB_IDS[@]}"; do
    printf '[smc++] waiting for job %s\n' "$job_id"
    wait_for_job "$job_id"
done

IFS=',' read -r -a CANDIDATE_BASES <<<"$CANDIDATE_BASES_CSV"
if [[ "${#CANDIDATE_BASES[@]}" -eq 0 ]]; then
    echo "[smc++] no candidate bases provided" >&2
    exit 1
fi

"$SDS_PYTHON" \
    "$PROJECT_ROOT/benchmark/demography/evaluate_smcpp_recent_resolution.py" \
    --pop "$POP" \
    --coarse-base "$COARSE_BASE" \
    --reference-base "$REFERENCE_BASE" \
    --candidate-bases "${CANDIDATE_BASES[@]}"

REPORT_JSON="$PROJECT_ROOT/benchmark/demography/$POP/smcpp/recent_sensitivity/${POP}_recent_resolution_report.json"
"$SDS_PYTHON" \
    "$PROJECT_ROOT/benchmark/demography/plot_smcpp_recent_sensitivity.py" \
    --report "$REPORT_JSON"

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
            "$FOLLOWUP_TAG" \
            "$FOLLOWUP_WINDOW_SIZE" \
            "$FOLLOWUP_KNOTS" \
            "$FOLLOWUP_TIMEPOINT_START" \
            "$FOLLOWUP_TIMEPOINT_END"
    else
        printf '[smc++] %s report did not pass; skipping follow-up submit for %s\n' "$POP" "$FOLLOWUP_POP"
    fi
fi
