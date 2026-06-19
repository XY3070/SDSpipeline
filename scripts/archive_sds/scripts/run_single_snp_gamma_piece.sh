#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 --pop POP --daf DAF --scenario-npz PATH --present-ne INT --gamma-prefix PATH --workdir PATH [--sim-reps N]" >&2
    exit 1
}

POP=""
DAF=""
SCENARIO_NPZ=""
PRESENT_NE=""
GAMMA_PREFIX=""
WORKDIR=""
SIM_REPS="1000"
TIP_BRANCH_DUMP="0"
MS_MAKE_DIR="/data/home/grp-wangyf/xuyuan/ms/scripts"
MS_BINARY="/data/home/grp-wangyf/xuyuan/ms/msdir/ms"
BACKWARD_SCRIPT="/data/home/grp-wangyf/xuyuan/ms/scripts/backward.py"
SAMPLE_SIZE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pop) POP="$2"; shift 2 ;;
        --daf) DAF="$2"; shift 2 ;;
        --scenario-npz) SCENARIO_NPZ="$2"; shift 2 ;;
        --present-ne) PRESENT_NE="$2"; shift 2 ;;
        --gamma-prefix) GAMMA_PREFIX="$2"; shift 2 ;;
        --workdir) WORKDIR="$2"; shift 2 ;;
        --sim-reps) SIM_REPS="$2"; shift 2 ;;
        --dump-tip-branches) TIP_BRANCH_DUMP="1"; shift 1 ;;
        --ms-make-dir) MS_MAKE_DIR="$2"; shift 2 ;;
        --ms-binary) MS_BINARY="$2"; shift 2 ;;
        --backward-script) BACKWARD_SCRIPT="$2"; shift 2 ;;
        --sample-size) SAMPLE_SIZE="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; usage ;;
    esac
done

[[ -n "$POP" && -n "$DAF" && -n "$SCENARIO_NPZ" && -n "$PRESENT_NE" && -n "$GAMMA_PREFIX" && -n "$WORKDIR" ]] || usage

mkdir -p "$(dirname "$GAMMA_PREFIX")" "$(dirname "$WORKDIR")"
/bin/rm -rf "$WORKDIR"

make -C "$MS_MAKE_DIR" sim_single_daf \
    DAF="$DAF" \
    POP_MODEL="$POP" \
    NE_STAT="median" \
    NPZ_PATH="$SCENARIO_NPZ" \
    SIM_NUM_REPLICATIONS="$SIM_REPS" \
    PRESENT_DIPLOID_POPULATION_SIZE="$PRESENT_NE" \
    GAMMA_PREFIX="$GAMMA_PREFIX" \
    WORKDIR="$WORKDIR" \
    TIP_BRANCH_DUMP="$TIP_BRANCH_DUMP" \
    MS="$MS_BINARY" \
    ${SAMPLE_SIZE:+SAMPLE_SIZE="$SAMPLE_SIZE"} \
    "SIMUPOP_BACKWARD=$BACKWARD_SCRIPT -m $POP --npz_path $SCENARIO_NPZ --ne_stat median"

OUTPUT_PATH="${GAMMA_PREFIX}.${DAF}"
if [[ ! -s "$OUTPUT_PATH" ]]; then
    echo "ERROR: expected non-empty gamma output at $OUTPUT_PATH" >&2
    exit 1
fi
