#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RELATE_DIR="$PROJECT_ROOT/relate"
BENCH_ROOT="$PROJECT_ROOT/benchmark/demography"
# shellcheck source=/data/home/grp-wangyf/xuyuan/sds/scripts/common_env.sh
source "$SCRIPT_DIR/common_env.sh"
activate_relate_runtime

POP=""
CHROMOSOMES="${CHROMOSOMES:-1-22}"
THREADS="${THREADS:-16}"
RELATE_INIT_NE="${RELATE_INIT_NE:-30000}"
MU="${MU:-1.25e-8}"
YEARS_PER_GEN="${YEARS_PER_GEN:-28}"
BINS="${BINS:-1.3,5,0.1}"
NUM_ITER="${NUM_ITER:-5}"
THRESHOLD="${THRESHOLD:-0.0}"
SEED="${SEED:-1}"
TAG="${TAG:-relate_recent}"
SMCPP_SUFFIX="${SMCPP_SUFFIX:-_fine_smcpp.csv}"
FORCE=0
SKIP_COMPARE=0

usage() {
    cat <<'EOF' >&2
Usage: run_relate_recent_ne.sh --pop POP [options]

Options:
  --pop POP               Population label, e.g. NCN or SCN.
  --chromosomes SPEC      Chromosome list/ranges. Default: 1-22.
  --threads N             Relate thread count. Default: 16.
  --init-ne N             Initial Relate -N value. Default: 30000.
  --mu MU                 Mutation rate. Default: 1.25e-8.
  --years-per-gen N       Years per generation. Default: 28.
  --bins SPEC             Relate epoch bins. Default: 1.3,5,0.1.
  --num-iter N            EstimatePopulationSize iterations. Default: 5.
  --threshold X           EstimatePopulationSize threshold. Default: 0.0.
  --seed N                EstimatePopulationSize seed. Default: 1.
  --tag TAG               Output tag. Default: relate_recent.
  --smcpp-suffix SUFFIX   SMC++ CSV suffix for comparison. Default: _fine_smcpp.csv.
  --skip-compare          Do not regenerate comparison figures.
  --force                 Rebuild chromosome trees and popsize outputs.
EOF
    exit "${1:-1}"
}

expand_chromosomes() {
    python3 - "$1" <<'PY'
import sys

values = []
for token in sys.argv[1].split(","):
    item = token.strip()
    if not item:
        continue
    if "-" in item:
        start_text, end_text = item.split("-", 1)
        start = int(start_text)
        end = int(end_text)
        if end < start:
            raise SystemExit(f"Invalid chromosome range: {item}")
        values.extend(range(start, end + 1))
    else:
        values.append(int(item))
print("\n".join(str(v) for v in sorted(dict.fromkeys(values))))
PY
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage 0
            ;;
        --pop)
            POP="$2"
            shift 2
            ;;
        --chromosomes)
            CHROMOSOMES="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --init-ne)
            RELATE_INIT_NE="$2"
            shift 2
            ;;
        --mu)
            MU="$2"
            shift 2
            ;;
        --years-per-gen)
            YEARS_PER_GEN="$2"
            shift 2
            ;;
        --bins)
            BINS="$2"
            shift 2
            ;;
        --num-iter)
            NUM_ITER="$2"
            shift 2
            ;;
        --threshold)
            THRESHOLD="$2"
            shift 2
            ;;
        --seed)
            SEED="$2"
            shift 2
            ;;
        --tag)
            TAG="$2"
            shift 2
            ;;
        --smcpp-suffix)
            SMCPP_SUFFIX="$2"
            shift 2
            ;;
        --skip-compare)
            SKIP_COMPARE=1
            shift
            ;;
        --force)
            FORCE=1
            shift
            ;;
        *)
            echo "Unknown parameter: $1" >&2
            usage
            ;;
    esac
done

[[ -n "$POP" ]] || usage 1
POP="${POP^^}"

mapfile -t CHR_LIST < <(expand_chromosomes "$CHROMOSOMES")
[[ "${#CHR_LIST[@]}" -gt 0 ]] || { echo "[Error] No chromosomes expanded from $CHROMOSOMES" >&2; exit 1; }

SUBSET_FILE="$BENCH_ROOT/$POP/subset_100.samples.txt"
[[ -f "$SUBSET_FILE" ]] || { echo "[Error] Missing subset file: $SUBSET_FILE" >&2; exit 1; }
SUBSET_SIZE="$(awk 'NF {count++} END {print count+0}' "$SUBSET_FILE")"
BASE_PREFIX="UKBQC_${POP}_subset${SUBSET_SIZE}"

PREP_ARGS=(--pop "$POP" --chromosomes "$CHROMOSOMES")
if [[ "$FORCE" -eq 1 ]]; then
    PREP_ARGS+=(--force)
fi
"$SCRIPT_DIR/prepare_recent_ne_inputs.sh" "${PREP_ARGS[@]}"

RELATE_ROOT="$BENCH_ROOT/$POP/relate_recent"
PREP_DIR="$RELATE_ROOT/prepared"
TREE_DIR="$RELATE_ROOT/trees"
POPSIZE_DIR="$RELATE_ROOT/popsize"
SHARED_REF_ROOT="$BENCH_ROOT/relate_shared"
POPLABELS="$RELATE_ROOT/${POP}.poplabels"
mkdir -p "$TREE_DIR" "$POPSIZE_DIR"

CHR_NAMES_FILE="$RELATE_ROOT/${TAG}_chr_names.txt"
printf '%s\n' "${CHR_LIST[@]}" > "$CHR_NAMES_FILE"

for CHR in "${CHR_LIST[@]}"; do
    PREP_PREFIX="$PREP_DIR/${BASE_PREFIX}_chr${CHR}"
    OUT_BASENAME="${BASE_PREFIX}_chr${CHR}"
    OUT_PREFIX="$TREE_DIR/$OUT_BASENAME"
    MAP_PATH="$SHARED_REF_ROOT/refs/relate_input_files_zenodo_15801307/Relate_input_files/GRCh38/recomb_map/genetic_map_chr${CHR}.txt"

    [[ -f "${PREP_PREFIX}.haps.gz" ]] || { echo "[Error] Missing prepared haps: ${PREP_PREFIX}.haps.gz" >&2; exit 1; }
    [[ -f "${PREP_PREFIX}.sample.gz" ]] || { echo "[Error] Missing prepared sample: ${PREP_PREFIX}.sample.gz" >&2; exit 1; }
    [[ -f "${PREP_PREFIX}.dist.gz" ]] || { echo "[Error] Missing prepared dist: ${PREP_PREFIX}.dist.gz" >&2; exit 1; }
    [[ -f "${PREP_PREFIX}.annot" ]] || { echo "[Error] Missing prepared annot: ${PREP_PREFIX}.annot" >&2; exit 1; }
    if [[ "$FORCE" -eq 1 ]]; then
        rm -f "${OUT_PREFIX}.anc" "${OUT_PREFIX}.mut"
    fi

    if [[ ! -f "${OUT_PREFIX}.anc" || ! -f "${OUT_PREFIX}.mut" ]]; then
        [[ -f "$MAP_PATH" ]] || { echo "[Error] Missing recombination map: $MAP_PATH" >&2; exit 1; }
        (
            cd "$TREE_DIR"
            env -u SHELLOPTS "$RELATE_DIR/scripts/RelateParallel/RelateParallel.sh" \
                --haps "${PREP_PREFIX}.haps.gz" \
                --sample "${PREP_PREFIX}.sample.gz" \
                --map "$MAP_PATH" \
                --dist "${PREP_PREFIX}.dist.gz" \
                --annot "${PREP_PREFIX}.annot" \
                -m "$MU" \
                -N "$RELATE_INIT_NE" \
                -o "$OUT_BASENAME" \
                --threads "$THREADS"
        )
    else
        printf '[Skip] Relate trees already exist: %s\n' "$OUT_PREFIX"
    fi
done

POPSIZE_PREFIX="$POPSIZE_DIR/${POP}_${TAG}"
if [[ "$FORCE" -eq 1 ]]; then
    rm -f "${POPSIZE_PREFIX}.coal" "${POPSIZE_PREFIX}.pairwise.coal" "${POPSIZE_PREFIX}.pairwise.bin"
    rm -f "${POPSIZE_PREFIX}_ne.csv" "${POPSIZE_PREFIX}.png" "${POPSIZE_PREFIX}_manifest.json"
fi

if [[ ! -f "${POPSIZE_PREFIX}.coal" ]]; then
    env -u SHELLOPTS "$RELATE_DIR/scripts/EstimatePopulationSize/EstimatePopulationSize.sh" \
        -i "$TREE_DIR/${BASE_PREFIX}" \
        --chr "$CHR_NAMES_FILE" \
        -m "$MU" \
        --poplabels "$POPLABELS" \
        --pop_of_interest "$POP" \
        -o "$POPSIZE_PREFIX" \
        --threads "$THREADS" \
        --bins "$BINS" \
        --years_per_gen "$YEARS_PER_GEN" \
        --num_iter "$NUM_ITER" \
        --threshold "$THRESHOLD" \
        --seed "$SEED" \
        --noplot
else
    printf '[Skip] Relate population size output already exists: %s.coal\n' "$POPSIZE_PREFIX"
fi

"$SDS_ENV_PREFIX/bin/python" \
    "$BENCH_ROOT/relate_coal_to_ne.py" \
    --coal "${POPSIZE_PREFIX}.coal" \
    --population "$POP" \
    --output-csv "${POPSIZE_PREFIX}_ne.csv" \
    --output-json "${POPSIZE_PREFIX}_manifest.json" \
    --output-png "${POPSIZE_PREFIX}.png" \
    --years-per-gen "$YEARS_PER_GEN"

if [[ "$SKIP_COMPARE" -eq 0 ]]; then
    AVAILABLE_POPS=()
    for candidate in NCN SCN; do
        if [[ -f "$BENCH_ROOT/$candidate/relate_recent/popsize/${candidate}_${TAG}_ne.csv" ]]; then
            AVAILABLE_POPS+=("$candidate")
        fi
    done
    if [[ "${#AVAILABLE_POPS[@]}" -gt 0 ]]; then
        POP_CSV="$(IFS=,; echo "${AVAILABLE_POPS[*]}")"
        OUTPUT_PREFIX="Population_History_Recent_Method_Comparison"
        if [[ "${#AVAILABLE_POPS[@]}" -eq 1 ]]; then
            OUTPUT_PREFIX="${AVAILABLE_POPS[0]}_Recent_Method_Comparison"
        fi
        "$SDS_ENV_PREFIX/bin/python" \
            "$BENCH_ROOT/evaluate_recent_ne_methods.py" \
            --root "$BENCH_ROOT" \
            --pops "$POP_CSV" \
            --smcpp-suffix "$SMCPP_SUFFIX" \
            --relate-tag "$TAG" \
            --output-prefix "$OUTPUT_PREFIX"
    fi
fi

printf '[Done] Relate recent Ne finished for %s with tag %s\n' "$POP" "$TAG"
