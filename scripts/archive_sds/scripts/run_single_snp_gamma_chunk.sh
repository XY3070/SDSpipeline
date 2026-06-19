#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 --pop POP --daf DAF --start-rep N --end-rep N --scenario-npz PATH --present-ne INT --workdir PATH [--ms-binary PATH]" >&2
    exit 1
}

POP=""
DAF=""
START_REP=""
END_REP=""
SCENARIO_NPZ=""
PRESENT_NE=""
WORKDIR=""
TIP_BRANCH_DUMP="0"
MS_MAKE_DIR="/data/home/grp-wangyf/xuyuan/ms/scripts"
MS_BINARY="/data/home/grp-wangyf/xuyuan/ms/msdir/ms"
BACKWARD_SCRIPT="/data/home/grp-wangyf/xuyuan/ms/scripts/backward.py"
SAMPLE_SIZE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pop) POP="$2"; shift 2 ;;
        --daf) DAF="$2"; shift 2 ;;
        --start-rep) START_REP="$2"; shift 2 ;;
        --end-rep) END_REP="$2"; shift 2 ;;
        --scenario-npz) SCENARIO_NPZ="$2"; shift 2 ;;
        --present-ne) PRESENT_NE="$2"; shift 2 ;;
        --workdir) WORKDIR="$2"; shift 2 ;;
        --dump-tip-branches) TIP_BRANCH_DUMP="1"; shift 1 ;;
        --ms-make-dir) MS_MAKE_DIR="$2"; shift 2 ;;
        --ms-binary) MS_BINARY="$2"; shift 2 ;;
        --backward-script) BACKWARD_SCRIPT="$2"; shift 2 ;;
        --sample-size) SAMPLE_SIZE="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; usage ;;
    esac
done

[[ -n "$POP" && -n "$DAF" && -n "$START_REP" && -n "$END_REP" && -n "$SCENARIO_NPZ" && -n "$PRESENT_NE" && -n "$WORKDIR" ]] || usage

/bin/rm -rf "$WORKDIR"
mkdir -p "$WORKDIR"

fail_count=0
for i in $(seq "$START_REP" "$END_REP"); do
    make -C "$MS_MAKE_DIR" my_simulation1 \
        DAF="$DAF" \
        I="$i" \
        WORKDIR="$WORKDIR" \
        TIP_BRANCH_DUMP="$TIP_BRANCH_DUMP" \
        POP_MODEL="$POP" \
        NE_STAT="median" \
        NPZ_PATH="$SCENARIO_NPZ" \
        PRESENT_DIPLOID_POPULATION_SIZE="$PRESENT_NE" \
        MS="$MS_BINARY" \
        ${SAMPLE_SIZE:+SAMPLE_SIZE="$SAMPLE_SIZE"} \
        "SIMUPOP_BACKWARD=$BACKWARD_SCRIPT -m $POP --npz_path $SCENARIO_NPZ --ne_stat median" \
        > /dev/null 2>&1 || { echo "Warning: Sim $i failed" >&2; fail_count=$((fail_count + 1)); }
done

if [[ "$fail_count" -gt 0 ]]; then
    echo "ERROR: $fail_count simulation replicates failed in $WORKDIR" >&2
    exit 1
fi

echo "chunk_complete workdir=$WORKDIR reps=${START_REP}-${END_REP}"
