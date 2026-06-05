#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck source=common_env.sh
source "$SCRIPT_DIR/common_env.sh"
activate_sds_env

T_FILE="$1"
OUT_FILE="$2"
S_FILE="$3"
O_FILE="$4"
B_FILE="$5"
G_FILE="$6"
INIT="${7:-0.00001}"
S_FILE_NCOL="${8:-20000}"
OUT_PARQUET="${9:-}"
SUMMARY_CSV="${10:-}"
CACHE_DIR="${11:-}"
SKIP_BOUNDARY_MISSING_FRACTION="${12:-}"
BOUNDARY_MISSING_MODE="${13:-}"
PYTHON_BIN="$SDS_ENV_PREFIX/bin/python"

[[ -x "$PYTHON_BIN" ]] || { echo "[Error] Python not found in SDS env: $PYTHON_BIN" >&2; exit 1; }

EXTRA_ARGS=()
if [[ -n "$OUT_PARQUET" ]]; then
    EXTRA_ARGS+=(--output-parquet "$OUT_PARQUET")
fi
if [[ -n "$SUMMARY_CSV" ]]; then
    EXTRA_ARGS+=(--summary-csv "$SUMMARY_CSV")
fi
if [[ -n "$CACHE_DIR" ]]; then
    EXTRA_ARGS+=(--pickle-cache-dir "$CACHE_DIR")
fi
if [[ -n "$SKIP_BOUNDARY_MISSING_FRACTION" ]]; then
    EXTRA_ARGS+=(--skip-boundary-missing-fraction "$SKIP_BOUNDARY_MISSING_FRACTION")
fi
if [[ -n "$BOUNDARY_MISSING_MODE" ]]; then
    EXTRA_ARGS+=(--boundary-missing-mode "$BOUNDARY_MISSING_MODE")
fi

NUMBA_CACHE_DIR="$REPO_ROOT/.cache/numba"
mkdir -p "$NUMBA_CACHE_DIR"

THREAD_BUDGET="${SDS_NUMBA_THREADS:-${LSB_DJOB_NUMPROC:-1}}"
if ! [[ "$THREAD_BUDGET" =~ ^[0-9]+$ ]] || [[ "$THREAD_BUDGET" -lt 1 ]]; then
    THREAD_BUDGET=1
fi

export PYTHONPATH="$REPO_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export NUMBA_CACHE_DIR
export SDS_NUMBA_THREADS="$THREAD_BUDGET"
export NUMBA_NUM_THREADS="$THREAD_BUDGET"
export OMP_NUM_THREADS="$THREAD_BUDGET"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

"$PYTHON_BIN" "$SCRIPT_DIR/compute_SDS.py" \
    "$S_FILE" \
    "$T_FILE" \
    "$O_FILE" \
    "$B_FILE" \
    "$G_FILE" \
    "$INIT" \
    "$S_FILE_NCOL" \
    --output "$OUT_FILE" \
    "${EXTRA_ARGS[@]}"
