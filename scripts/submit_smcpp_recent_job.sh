#!/bin/bash
set -euo pipefail

if [[ $# -lt 6 || $# -gt 7 ]]; then
    cat >&2 <<'EOF'
Usage: submit_smcpp_recent_job.sh POP TAG WINDOW_SIZE KNOTS TIMEPOINT_START TIMEPOINT_END [EM_ITERATIONS]

Example:
  submit_smcpp_recent_job.sh NCN tp20_w10_k16 10 16 20 100000
EOF
    exit 1
fi

POP="$(printf '%s' "$1" | tr '[:lower:]' '[:upper:]')"
TAG="$2"
WINDOW_SIZE="$3"
KNOTS="$4"
TIMEPOINT_START="$5"
TIMEPOINT_END="$6"
EM_ITERATIONS="${7:-}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common_env.sh
source "$SCRIPT_DIR/common_env.sh"
LOG_DIR="${SMCPP_LOG_DIR:-$SDS_RUNS_ROOT/smcpp/logs}"
BENCH_ROOT="${SMCPP_BENCH_ROOT:-$SDS_DEMOGRAPHY_ROOT}"
RUNNER_SCRIPT="$SCRIPT_DIR/run_smcpp_benchmark_lsf.sh"
QUEUE="${SMCPP_LSF_QUEUE:-normal}"
SLOTS="${SMCPP_LSF_SLOTS:-1}"
JOB_NAME="${SMCPP_JOB_NAME:-smcpp_bench_${POP}_${TAG}}"
OUT_LOG="$LOG_DIR/smcpp_bench_${POP}_${TAG}.%J.out"
ERR_LOG="$LOG_DIR/smcpp_bench_${POP}_${TAG}.%J.err"
BASE_NAME="${POP}_${TAG}"

mkdir -p "$LOG_DIR"
[[ -f "$RUNNER_SCRIPT" ]] || { echo "[Error] Missing dependency: $RUNNER_SCRIPT" >&2; exit 1; }

EXTRA_ARGS=(--knots "$KNOTS" --timepoints "$TIMEPOINT_START" "$TIMEPOINT_END")
CMD=(
    env
    "BENCH_ROOT=$BENCH_ROOT"
    "SMCPP_BASE=$BASE_NAME"
    "SMCPP_WINDOW_SIZE=$WINDOW_SIZE"
    "SMCPP_EXTRA_ARGS=${EXTRA_ARGS[*]}"
)

if [[ -n "$EM_ITERATIONS" ]]; then
    CMD+=("SMCPP_EM_ITERATIONS=$EM_ITERATIONS")
fi

CMD+=(
    /bin/bash
    "$RUNNER_SCRIPT"
    "$POP"
)

printf '[smc++] submit job=%s queue=%s slots=%s\n' "$JOB_NAME" "$QUEUE" "$SLOTS"
printf '[smc++] command:'
printf ' %q' "${CMD[@]}"
printf '\n'

if [[ "${DRY_RUN:-0}" == "1" ]]; then
    exit 0
fi

BSUB_OUTPUT="$(
    bsub \
        -q "$QUEUE" \
        -n "$SLOTS" \
        -J "$JOB_NAME" \
        -o "$OUT_LOG" \
        -e "$ERR_LOG" \
        "${CMD[@]}"
)"

printf '%s\n' "$BSUB_OUTPUT"
JOB_ID="$(printf '%s\n' "$BSUB_OUTPUT" | sed -n 's/.*<\([0-9][0-9]*\)>.*/\1/p')"
if [[ -n "$JOB_ID" ]]; then
    printf '[smc++] submitted job_id=%s log_prefix=%s\n' "$JOB_ID" "$LOG_DIR/smcpp_bench_${POP}_${TAG}.$JOB_ID"
fi
