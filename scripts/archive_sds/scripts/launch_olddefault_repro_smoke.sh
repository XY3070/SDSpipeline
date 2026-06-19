#!/bin/bash
set -euo pipefail

ROOT="/data/home/grp-wangyf/xuyuan"
SDS_ROOT="$ROOT/sds"
OUTDIR="$SDS_ROOT/tmp/olddefault_repro_smoke_20260525"
LOGDIR="$SDS_ROOT/logs"
SCENARIO_NPZ="$SDS_ROOT/tmp/region_ne0_positive_20260426/scaled_ne0_100000/NCN_scaled_ne0_100000.npz"
CURRENT_MAKE_DIR="$ROOT/ms/scripts"
SNAPSHOT_MAKE_DIR="$ROOT/ms/scripts/snapshots/20260405_1528"
PRESENT_NE="219253"
SIM_REPS="100"
DAFS=("0.05" "0.07" "0.10")

mkdir -p "$OUTDIR" "$LOGDIR" "$OUTDIR/results"

MANIFEST="$OUTDIR/manifest.tsv"
{
    printf "label\tmakefile_source\tdaf\tpresent_ne\tsim_reps\tgamma_path\n"
} > "$MANIFEST"

job_ids=()

submit_one() {
    local label="$1"
    local make_dir="$2"
    local daf="$3"
    local gamma_prefix="$OUTDIR/results/$label/gamma"
    local workdir="$OUTDIR/results/$label/workdir"
    mkdir -p "$(dirname "$gamma_prefix")"
    printf "%s\t%s\t%s\t%s\t%s\t%s.%s\n" \
        "$label" "$(basename "$make_dir")" "$daf" "$PRESENT_NE" "$SIM_REPS" "$gamma_prefix" "$daf" >> "$MANIFEST"
    local outlog="$LOGDIR/${label}.%J.out"
    local errlog="$LOGDIR/${label}.%J.err"
    local cmd=(
        bsub -q smp -J "$label" -n 1 -R "span[hosts=1]"
        -o "$outlog" -e "$errlog"
        /bin/bash "$SDS_ROOT/scripts/run_single_snp_gamma_piece.sh"
        --pop Gravel_CHB
        --daf "$daf"
        --scenario-npz "$SCENARIO_NPZ"
        --present-ne "$PRESENT_NE"
        --gamma-prefix "$gamma_prefix"
        --workdir "$workdir"
        --sim-reps "$SIM_REPS"
        --ms-make-dir "$make_dir"
    )
    local submit_output
    submit_output="$("${cmd[@]}")"
    echo "$submit_output"
    local job_id
    job_id="$(sed -n 's/.*Job <\([0-9]\+\)>.*/\1/p' <<<"$submit_output")"
    job_ids+=("$job_id")
}

for daf in "${DAFS[@]}"; do
    submit_one "olddefault_repro_current_daf${daf/./p}" "$CURRENT_MAKE_DIR" "$daf"
    submit_one "olddefault_repro_snapshot_daf${daf/./p}" "$SNAPSHOT_MAKE_DIR" "$daf"
done

dep=""
for job_id in "${job_ids[@]}"; do
    if [[ -n "$dep" ]]; then
        dep="${dep}&&"
    fi
    dep="${dep}done(${job_id})"
done

bsub -q smp -w "$dep" -J olddefault_repro_summary -n 1 -R "span[hosts=1]" \
    -o "$LOGDIR/olddefault_repro_summary.%J.out" \
    -e "$LOGDIR/olddefault_repro_summary.%J.err" \
    /data/home/grp-wangyf/intern/miniforge3/envs/sds/bin/python \
    "$SDS_ROOT/scripts/summarize_olddefault_repro_smoke.py" \
    --root "$OUTDIR"

echo "manifest: $MANIFEST"
