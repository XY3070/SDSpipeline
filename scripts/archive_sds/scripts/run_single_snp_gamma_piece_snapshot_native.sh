#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 --pop POP --daf DAF --scenario-npz PATH --present-ne INT --out-path PATH --workdir PATH [--sim-reps N] [--snapshot-makefile PATH]" >&2
    exit 1
}

POP=""
DAF=""
SCENARIO_NPZ=""
PRESENT_NE=""
OUT_PATH=""
WORKDIR=""
SIM_REPS="1000"
SNAPSHOT_MAKEFILE="/data/home/grp-wangyf/xuyuan/ms/scripts/snapshots/20260405_1528/Makefile.bak"
MS_BINARY="/data/home/grp-wangyf/xuyuan/ms/msdir/ms"
BACKWARD_SCRIPT="/data/home/grp-wangyf/xuyuan/ms/scripts/backward.py"
SAMPLE_SIZE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pop) POP="$2"; shift 2 ;;
        --daf) DAF="$2"; shift 2 ;;
        --scenario-npz) SCENARIO_NPZ="$2"; shift 2 ;;
        --present-ne) PRESENT_NE="$2"; shift 2 ;;
        --out-path) OUT_PATH="$2"; shift 2 ;;
        --workdir) WORKDIR="$2"; shift 2 ;;
        --sim-reps) SIM_REPS="$2"; shift 2 ;;
        --snapshot-makefile) SNAPSHOT_MAKEFILE="$2"; shift 2 ;;
        --ms-binary) MS_BINARY="$2"; shift 2 ;;
        --backward-script) BACKWARD_SCRIPT="$2"; shift 2 ;;
        --sample-size) SAMPLE_SIZE="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; usage ;;
    esac
done

[[ -n "$POP" && -n "$DAF" && -n "$SCENARIO_NPZ" && -n "$PRESENT_NE" && -n "$OUT_PATH" && -n "$WORKDIR" ]] || usage

RUNROOT="$WORKDIR/runroot"
/bin/rm -rf "$WORKDIR"
mkdir -p "$RUNROOT" "$(dirname "$OUT_PATH")"
cp "$SNAPSHOT_MAKEFILE" "$RUNROOT/Makefile"

make -C "$RUNROOT" sim_single_daf \
    DAF="$DAF" \
    POP_MODEL="$POP" \
    NPZ_PATH="$SCENARIO_NPZ" \
    NE_STAT="median" \
    SIM_NUM_REPLICATIONS="$SIM_REPS" \
    PRESENT_DIPLOID_POPULATION_SIZE="$PRESENT_NE" \
    ${SAMPLE_SIZE:+SAMPLE_SIZE="$SAMPLE_SIZE"} \
    MS="$MS_BINARY" \
    "SIMUPOP_BACKWARD=$BACKWARD_SCRIPT -m $POP --npz_path $SCENARIO_NPZ --ne_stat median"

SOURCE_PATH="$RUNROOT/sds_input.gamma_shapes.$DAF"
if [[ ! -s "$SOURCE_PATH" ]]; then
    echo "ERROR: expected non-empty gamma output at $SOURCE_PATH" >&2
    exit 1
fi
cp "$SOURCE_PATH" "$OUT_PATH"
