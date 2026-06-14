#!/bin/bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 POP" >&2
    exit 1
fi

POP="$1"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# shellcheck source=common_env.sh
source "$SCRIPT_DIR/common_env.sh"
activate_sds_env

export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/matplotlib}"

BENCH_ROOT="${BENCH_ROOT:-$SDS_DEMOGRAPHY_ROOT}"
SMCPP_IMAGE="${SMCPP_IMAGE:-$SDS_SMCPP_IMAGE}"
SMCPP_CORES="${SMCPP_CORES:-8}"
SMCPP_WINDOW_SIZE="${SMCPP_WINDOW_SIZE:-20}"
SMCPP_KNOTS="${SMCPP_KNOTS:-8}"
SMCPP_SPLINE="${SMCPP_SPLINE:-piecewise}"
SMCPP_TIMEPOINTS="${SMCPP_TIMEPOINTS:-10.0 100000.0}"
SMCPP_MU="${SMCPP_MU:-1.25e-8}"
SMCPP_MISSING_CUTOFF="${SMCPP_MISSING_CUTOFF:-1000000}"
SMCPP_NONSEG_CUTOFF="${SMCPP_NONSEG_CUTOFF:-1000000}"
SMCPP_BASE="${SMCPP_BASE:-$POP}"
SMCPP_EXTRA_ARGS="${SMCPP_EXTRA_ARGS:-}"
SMCPP_SINGULARITY="${SMCPP_SINGULARITY:-}"
SMCPP_SINGULARITY_CONFDIR="${SMCPP_SINGULARITY_CONFDIR:-$SDS_SMCPP_SINGULARITY_CONFDIR}"
SMCPP_ROOTFS="${SMCPP_ROOTFS:-$SDS_SMCPP_ROOTFS}"
SMCPP_MPLCONFIGDIR="${SMCPP_MPLCONFIGDIR:-$SDS_SMCPP_MPLCONFIGDIR}"
SMCPP_PROBE_CATEGORY=""

resolve_runtime_candidate() {
    local candidate="$1"
    local resolved=""
    if [[ "$candidate" == */* ]]; then
        [[ -x "$candidate" ]] || return 1
        printf '%s\n' "$candidate"
        return 0
    fi
    resolved="$(command -v "$candidate" 2>/dev/null || true)"
    [[ -n "$resolved" ]] || return 1
    printf '%s\n' "$resolved"
}

detect_smcpp_runtime_for_log() {
    if [[ -d "$SMCPP_ROOTFS" ]]; then
        SMCPP_RUNTIME_PATH="$SMCPP_ROOTFS"
        SMCPP_RUNTIME_SOURCE="shared smc++ rootfs"
        return 0
    fi
    if [[ -n "$SMCPP_SINGULARITY" ]]; then
        SMCPP_RUNTIME_PATH="$(resolve_runtime_candidate "$SMCPP_SINGULARITY" || true)"
        if [[ -n "$SMCPP_RUNTIME_PATH" ]]; then
            SMCPP_RUNTIME_SOURCE="SMCPP_SINGULARITY"
        else
            SMCPP_RUNTIME_PATH="$SMCPP_SINGULARITY"
            SMCPP_RUNTIME_SOURCE="SMCPP_SINGULARITY (unresolved override)"
        fi
        return 0
    fi
    if [[ -n "${SINGULARITY:-}" ]]; then
        SMCPP_RUNTIME_PATH="$(resolve_runtime_candidate "$SINGULARITY" || true)"
        if [[ -n "$SMCPP_RUNTIME_PATH" ]]; then
            SMCPP_RUNTIME_SOURCE="SINGULARITY environment variable"
            return 0
        fi
    fi
    SMCPP_RUNTIME_PATH="$(command -v singularity 2>/dev/null || true)"
    if [[ -n "$SMCPP_RUNTIME_PATH" ]]; then
        SMCPP_RUNTIME_SOURCE="PATH lookup for singularity"
        return 0
    fi
    SMCPP_RUNTIME_PATH="$(command -v apptainer 2>/dev/null || true)"
    if [[ -n "$SMCPP_RUNTIME_PATH" ]]; then
        SMCPP_RUNTIME_SOURCE="PATH lookup for apptainer"
        return 0
    fi
    if [[ -x "$SDS_SMCPP_SINGULARITY_BIN" ]]; then
        SMCPP_RUNTIME_PATH="$SDS_SMCPP_SINGULARITY_BIN"
        SMCPP_RUNTIME_SOURCE="project fallback"
        return 0
    fi
    return 1
}

classify_smcpp_probe_failure() {
    local probe_log="$1"
    if grep -q "Couldn't not parse configuration file" "$probe_log"; then
        SMCPP_PROBE_CATEGORY="config-missing"
    elif grep -q "No setuid installation found" "$probe_log"; then
        SMCPP_PROBE_CATEGORY="no-setuid-runtime"
    elif grep -q "socket communication error" "$probe_log"; then
        SMCPP_PROBE_CATEGORY="userns-socket-permission"
    elif grep -q "newuidmap was not found" "$probe_log"; then
        SMCPP_PROBE_CATEGORY="fakeroot-missing-newuidmap"
    elif grep -q "Could not locate an smc++ matplotlib config directory" "$probe_log"; then
        SMCPP_PROBE_CATEGORY="matplotlib-config-missing"
    else
        SMCPP_PROBE_CATEGORY="runner-probe-failed"
    fi
}

CMD=(
    "$SDS_PYTHON"
    "$SCRIPT_DIR/run_smcpp_benchmark.py"
    --pop "$POP"
    --subset-samples "$BENCH_ROOT/$POP/subset_100.samples.txt"
    --distinguished-pairs "$BENCH_ROOT/$POP/smcpp_distinguished_pairs.tsv"
    --subset-vcf-dir "$BENCH_ROOT/$POP/subset_vcf"
    --output-dir "$BENCH_ROOT/$POP/smcpp"
    --smcpp-image "$SMCPP_IMAGE"
    --cores "$SMCPP_CORES"
    --window-size "$SMCPP_WINDOW_SIZE"
    --knots "$SMCPP_KNOTS"
    --spline "$SMCPP_SPLINE"
    --timepoints $SMCPP_TIMEPOINTS
    --mu "$SMCPP_MU"
    --missing-cutoff "$SMCPP_MISSING_CUTOFF"
    --nonseg-cutoff "$SMCPP_NONSEG_CUTOFF"
    --base "$SMCPP_BASE"
)

if [[ -d "$SMCPP_ROOTFS" ]]; then
    CMD+=(--smcpp-rootfs "$SMCPP_ROOTFS")
fi

if [[ -d "$SMCPP_MPLCONFIGDIR" ]]; then
    CMD+=(--smcpp-mplconfigdir "$SMCPP_MPLCONFIGDIR")
fi

if [[ -n "$SMCPP_SINGULARITY" ]]; then
    CMD+=(--singularity "$SMCPP_SINGULARITY")
fi

if [[ -d "$SMCPP_SINGULARITY_CONFDIR" ]]; then
    CMD+=(--singularity-confdir "$SMCPP_SINGULARITY_CONFDIR")
fi

if [[ -n "${SMCPP_EM_ITERATIONS:-}" ]]; then
    CMD+=(--em-iterations "$SMCPP_EM_ITERATIONS")
fi

if [[ -n "$SMCPP_EXTRA_ARGS" ]]; then
    # Intentionally split extra args on shell whitespace for simple overrides.
    # This wrapper is only used by the submitting user on the same machine.
    # shellcheck disable=SC2206
    EXTRA_ARR=($SMCPP_EXTRA_ARGS)
    CMD+=("${EXTRA_ARR[@]}")
fi

printf '[smc++] host=%s pop=%s cores=%s\n' "$(hostname)" "$POP" "$SMCPP_CORES"
printf '[smc++] benchmark root=%s\n' "$BENCH_ROOT"
printf '[smc++] image=%s\n' "$SMCPP_IMAGE"
if detect_smcpp_runtime_for_log; then
    printf '[smc++] runner target=%s source=%s\n' "$SMCPP_RUNTIME_PATH" "$SMCPP_RUNTIME_SOURCE"
else
    printf '[smc++] runner target=<unresolved> source=auto-detect\n'
fi
if [[ -d "$SMCPP_SINGULARITY_CONFDIR" ]]; then
    printf '[smc++] singularity confdir=%s\n' "$SMCPP_SINGULARITY_CONFDIR"
fi
if [[ -d "$SMCPP_MPLCONFIGDIR" ]]; then
    printf '[smc++] matplotlib config=%s\n' "$SMCPP_MPLCONFIGDIR"
fi
printf '[smc++] command:'
printf ' %q' "${CMD[@]}"
printf '\n'

PROBE_LOG="$(mktemp /tmp/smcpp_probe.XXXXXX.log)"
trap 'rm -f "$PROBE_LOG"' EXIT
PROBE_CMD=("${CMD[@]}" --probe-only)
printf '[smc++] probe command:'
printf ' %q' "${PROBE_CMD[@]}"
printf '\n'
if ! "${PROBE_CMD[@]}" >"$PROBE_LOG" 2>&1; then
    cat "$PROBE_LOG" >&2
    classify_smcpp_probe_failure "$PROBE_LOG"
    printf '[smc++] probe failed category=%s\n' "$SMCPP_PROBE_CATEGORY" >&2
    exit 1
fi
cat "$PROBE_LOG"
printf '[smc++] probe succeeded\n'

"${CMD[@]}"
