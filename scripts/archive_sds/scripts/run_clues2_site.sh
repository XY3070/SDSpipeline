#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
RELATE_DIR="$(cd "$BASE_DIR/../relate" && pwd)"
# shellcheck source=/data/home/grp-wangyf/xuyuan/sds/scripts/common_env.sh
source "$SCRIPT_DIR/common_env.sh"
activate_relate_runtime

POP=""
CHR=""
POS=""
SITE_ID=""
OUT_ROOT="$SDS_RELATE_OUTPUT_ROOT"
NUM_SAMPLES="${NUM_SAMPLES:-200}"
MU="${MU:-1.25e-8}"
T_CUTOFF="${T_CUTOFF:-600}"
DF="${DF:-600}"
TIME_BINS="${TIME_BINS:-50 200}"
COAL_PATH=""
WITH_TRAJ=0

normalize_annot_file() {
    local annot_path="$1"
    local tmp_path

    [[ -f "$annot_path" ]] || return 0

    tmp_path="$(mktemp "${annot_path}.tmp.XXXXXX")"
    awk '{ sub(/;+$/, "", $0); print }' "$annot_path" > "$tmp_path"
    if ! cmp -s "$annot_path" "$tmp_path"; then
        mv "$tmp_path" "$annot_path"
    else
        rm -f "$tmp_path"
    fi
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pop)
            POP="$2"
            shift 2
            ;;
        --chr)
            CHR="$2"
            shift 2
            ;;
        --pos)
            POS="$2"
            shift 2
            ;;
        --site-id)
            SITE_ID="$2"
            shift 2
            ;;
        --out-root)
            OUT_ROOT="$2"
            shift 2
            ;;
        --num-samples)
            NUM_SAMPLES="$2"
            shift 2
            ;;
        --mu)
            MU="$2"
            shift 2
            ;;
        --coal)
            COAL_PATH="$2"
            shift 2
            ;;
        --t-cutoff)
            T_CUTOFF="$2"
            shift 2
            ;;
        --df)
            DF="$2"
            shift 2
            ;;
        --time-bins)
            shift
            TIME_BINS=""
            while [[ $# -gt 0 && "$1" != --* ]]; do
                TIME_BINS+="${TIME_BINS:+ }$1"
                shift
            done
            ;;
        --with-traj)
            WITH_TRAJ=1
            shift
            ;;
        *)
            echo "Unknown parameter: $1" >&2
            exit 1
            ;;
    esac
done

if [[ -z "$POP" || -z "$CHR" || -z "$POS" ]]; then
    echo "Usage: $0 --pop POP --chr CHR --pos POS [--site-id NAME] [--coal FILE] [--with-traj]" >&2
    exit 1
fi

if [[ -z "$SITE_ID" ]]; then
    SITE_ID="chr${CHR}_${POS}"
fi

if [[ -z "$COAL_PATH" ]]; then
    if [[ -f "$OUT_ROOT/$POP/popsize/${POP}.popsize.coal" ]]; then
        COAL_PATH="$OUT_ROOT/$POP/popsize/${POP}.popsize.coal"
    elif [[ -f "$OUT_ROOT/refs/hg38_generic_1000G.coal" ]]; then
        COAL_PATH="$OUT_ROOT/refs/hg38_generic_1000G.coal"
    fi
fi

[[ -n "$COAL_PATH" ]] || { echo "[Error] No .coal file found. Provide --coal explicitly or prepare refs/popsize first." >&2; exit 1; }
[[ -f "$COAL_PATH" ]] || { echo "[Error] coal file not found: $COAL_PATH" >&2; exit 1; }

PYTHON_BIN="$SDS_ENV_PREFIX/bin/python"
export MPLCONFIGDIR="$BASE_DIR/.cache/matplotlib"
mkdir -p "$MPLCONFIGDIR"

RELATE_PREFIX="$OUT_ROOT/$POP/relate/UKBQC_${POP}_chr${CHR}"
PREP_PREFIX="$OUT_ROOT/$POP/prepared/UKBQC_${POP}_chr${CHR}"
CLUES_DIR="$OUT_ROOT/$POP/clues/$SITE_ID"
mkdir -p "$CLUES_DIR"

[[ -f "${RELATE_PREFIX}.anc" ]] || { echo "[Error] Relate anc missing: ${RELATE_PREFIX}.anc" >&2; exit 1; }
[[ -f "${RELATE_PREFIX}.mut" ]] || { echo "[Error] Relate mut missing: ${RELATE_PREFIX}.mut" >&2; exit 1; }
[[ -f "${PREP_PREFIX}.haps.gz" ]] || { echo "[Error] Prepared haps missing: ${PREP_PREFIX}.haps.gz" >&2; exit 1; }

normalize_annot_file "${PREP_PREFIX}.annot"

BRANCH_PREFIX="$CLUES_DIR/$SITE_ID"
BRANCH_BASENAME="$SITE_ID"

if [[ -f "${BRANCH_PREFIX}.anc" && -f "${BRANCH_PREFIX}.mut" && -f "${BRANCH_PREFIX}.newick" && -f "${BRANCH_PREFIX}.sites" ]]; then
    printf '[Skip] SampleBranchLengths output already exists: %s\n' "$BRANCH_PREFIX"
else
    (
        cd "$CLUES_DIR"
        "$RELATE_DIR/scripts/SampleBranchLengths/SampleBranchLengths.sh" \
            -i "$RELATE_PREFIX" \
            -o "$BRANCH_BASENAME" \
            -m "$MU" \
            --coal "$COAL_PATH" \
            --format n \
            --num_samples "$NUM_SAMPLES" \
            --first_bp "$POS" \
            --last_bp "$POS"
    )
fi

printf '[Run] Building CLUES inputs: %s\n' "$BRANCH_PREFIX"
"$PYTHON_BIN" -u "$SCRIPT_DIR/prepare_clues_inputs.py" \
    --sites "${BRANCH_PREFIX}.sites" \
    --mut "${BRANCH_PREFIX}.mut" \
    --haps "${PREP_PREFIX}.haps.gz" \
    --position "$POS" \
    --out-prefix "$BRANCH_PREFIX"

printf '[Run] Converting Relate samples for CLUES2: %s\n' "$BRANCH_PREFIX"
"$PYTHON_BIN" -u "$BASE_DIR/../CLUES2/RelateToCLUES.py" \
    --RelateSamples "${BRANCH_PREFIX}.newick" \
    --DerivedFile "${BRANCH_PREFIX}_derived.txt" \
    --out "$BRANCH_PREFIX"

POP_FREQ="$(tr -d '\n' < "${BRANCH_PREFIX}_popfreq.txt")"

CLUES_ARGS=(
    "$BASE_DIR/../CLUES2/inference.py"
    --coal "$COAL_PATH"
    --times "${BRANCH_PREFIX}_times.txt"
    --popFreq "$POP_FREQ"
    --out "${BRANCH_PREFIX}.clues"
    --tCutoff "$T_CUTOFF"
    --df "$DF"
    --CI 0.95
)

if [[ -n "$TIME_BINS" ]]; then
    read -r -a TIME_BINS_ARR <<< "$TIME_BINS"
    CLUES_ARGS+=(--timeBins "${TIME_BINS_ARR[@]}")
fi

if [[ "$WITH_TRAJ" -eq 0 ]]; then
    CLUES_ARGS+=(--noAlleleTraj)
fi

printf '[Run] Running CLUES2 inference: %s\n' "${BRANCH_PREFIX}.clues"
"$PYTHON_BIN" -u "${CLUES_ARGS[@]}"

if [[ "$WITH_TRAJ" -eq 1 ]]; then
    printf '[Run] Plotting allele trajectory: %s\n' "${BRANCH_PREFIX}.clues"
    "$PYTHON_BIN" -u "$BASE_DIR/../CLUES2/plot_traj.py" \
        --freqs "${BRANCH_PREFIX}.clues_freqs.txt" \
        --post "${BRANCH_PREFIX}.clues_post.txt" \
        --figure "${BRANCH_PREFIX}.clues" \
        --generation_time 28
fi

printf '[Done] CLUES2 site finished: %s\n' "$BRANCH_PREFIX"
