#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
BENCH_ROOT="$PROJECT_ROOT/benchmark/demography"
SINGER_BIN_DIR="$PROJECT_ROOT/SINGER-0.1.8-beta/releases/singer-0.1.8-beta-linux-x86_64"
PHLASH_PYTHON="${PHLASH_PYTHON:-/data/home/grp-wangyf/intern/miniforge3/envs/phlash/bin/python}"

POP=""
CHROMOSOMES="${CHROMOSOMES:-20-22}"
MU="${MU:-1.25e-8}"
NE="${NE:-}"
RATIO="${RATIO:-1}"
BLOCK_LENGTH="${BLOCK_LENGTH:-1000000}"
POSTERIOR_SAMPLES="${POSTERIOR_SAMPLES:-100}"
THIN="${THIN:-20}"
POLAR="${POLAR:-0.99}"
FREQ="${FREQ:-1}"
NUM_CORES="${NUM_CORES:-16}"
TMRCA_WINDOW_SIZE="${TMRCA_WINDOW_SIZE:-1000000}"
FORCE=0
SKIP_EVAL=0

usage() {
    cat <<'EOF' >&2
Usage: run_singer_recent_pilot.sh --pop POP [options]

Options:
  --pop POP                 Population label. First intended use is NCN.
  --chromosomes SPEC        Chromosome list/ranges. Default: 20-22.
  --mu MU                   Mutation rate. Default: 1.25e-8.
  --ne N                    Explicit diploid Ne prior. Default: auto-estimate from clean VCFs.
  --ratio X                 Recombination/mutation ratio. Default: 1.
  --block-length N          parallel_singer block length. Default: 1000000.
  --posterior-samples N     Number of posterior samples. Default: 100.
  --thin N                  Thinning interval. Default: 20.
  --polar X                 Polarization probability. Default: 0.99.
  --freq N                  Convert to tskit every N samples. Default: 1.
  --num-cores N             parallel_singer core count. Default: 16.
  --tmrca-window-size N     Window size for diagnostic summaries. Default: 1000000.
  --skip-eval               Skip tskit-based diagnostic summaries.
  --force                   Re-run SINGER and overwrite diagnostics.
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
        --mu)
            MU="$2"
            shift 2
            ;;
        --ne)
            NE="$2"
            shift 2
            ;;
        --ratio)
            RATIO="$2"
            shift 2
            ;;
        --block-length)
            BLOCK_LENGTH="$2"
            shift 2
            ;;
        --posterior-samples)
            POSTERIOR_SAMPLES="$2"
            shift 2
            ;;
        --thin)
            THIN="$2"
            shift 2
            ;;
        --polar)
            POLAR="$2"
            shift 2
            ;;
        --freq)
            FREQ="$2"
            shift 2
            ;;
        --num-cores)
            NUM_CORES="$2"
            shift 2
            ;;
        --tmrca-window-size)
            TMRCA_WINDOW_SIZE="$2"
            shift 2
            ;;
        --skip-eval)
            SKIP_EVAL=1
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
[[ -x "$SINGER_BIN_DIR/parallel_singer" ]] || { echo "[Error] parallel_singer not found: $SINGER_BIN_DIR/parallel_singer" >&2; exit 1; }
[[ -x "$PHLASH_PYTHON" ]] || { echo "[Error] phlash Python not found: $PHLASH_PYTHON" >&2; exit 1; }

mapfile -t CHR_LIST < <(expand_chromosomes "$CHROMOSOMES")
[[ "${#CHR_LIST[@]}" -gt 0 ]] || { echo "[Error] No chromosomes expanded from $CHROMOSOMES" >&2; exit 1; }

PREP_ARGS=(--pop "$POP" --chromosomes "$CHROMOSOMES" --skip-relate)
if [[ "$FORCE" -eq 1 ]]; then
    PREP_ARGS+=(--force)
fi
"$SCRIPT_DIR/prepare_recent_ne_inputs.sh" "${PREP_ARGS[@]}"

SUBSET_FILE="$BENCH_ROOT/$POP/subset_100.samples.txt"
[[ -f "$SUBSET_FILE" ]] || { echo "[Error] Missing subset file: $SUBSET_FILE" >&2; exit 1; }
SUBSET_SIZE="$(awk 'NF {count++} END {print count+0}' "$SUBSET_FILE")"
BASE_PREFIX="UKBQC_${POP}_subset${SUBSET_SIZE}"
SINGER_ROOT="$BENCH_ROOT/$POP/singer_recent"
INPUT_DIR="$SINGER_ROOT/input"
OUTPUT_DIR="$SINGER_ROOT/output"
DIAG_DIR="$SINGER_ROOT/diagnostics"
mkdir -p "$OUTPUT_DIR" "$DIAG_DIR"

if [[ -z "$NE" ]]; then
    PRIOR_JSON="$INPUT_DIR/${POP}_singer_prior_ne.json"
    VCFS=()
    for CHR in "${CHR_LIST[@]}"; do
        VCFS+=("$INPUT_DIR/${BASE_PREFIX}_chr${CHR}.phlash.vcf.gz")
    done
    "$PHLASH_PYTHON" \
        "$BENCH_ROOT/estimate_singer_prior_ne.py" \
        --population "$POP" \
        --mu "$MU" \
        --output "$PRIOR_JSON" \
        "${VCFS[@]}"
    NE="$(python3 - "$PRIOR_JSON" <<'PY'
import json, sys
data = json.load(open(sys.argv[1]))
print(int(round(data["ne_diploid_estimate"])))
PY
)"
    printf '[SINGER] auto-estimated diploid Ne prior: %s\n' "$NE"
fi

for CHR in "${CHR_LIST[@]}"; do
    INPUT_VCF="$INPUT_DIR/${BASE_PREFIX}_chr${CHR}.phlash.vcf.gz"
    [[ -f "$INPUT_VCF" ]] || { echo "[Error] Missing SINGER input VCF: $INPUT_VCF" >&2; exit 1; }
    INPUT_PREFIX="${INPUT_VCF%.vcf.gz}"
    OUT_PREFIX="$OUTPUT_DIR/${BASE_PREFIX}_chr${CHR}"

    if [[ "$FORCE" -eq 1 ]]; then
        rm -f "${OUT_PREFIX}"_*.trees "${OUT_PREFIX}"_nodes_*.txt "${OUT_PREFIX}"_branches_*.txt \
            "${OUT_PREFIX}"_muts_*.txt "${OUT_PREFIX}"_recombs_*.txt
    fi

    if ! compgen -G "${OUT_PREFIX}_*.trees" > /dev/null; then
        "$SINGER_BIN_DIR/parallel_singer" \
            -Ne "$NE" \
            -m "$MU" \
            -ratio "$RATIO" \
            -L "$BLOCK_LENGTH" \
            -vcf "$INPUT_PREFIX" \
            -output "$OUT_PREFIX" \
            -n "$POSTERIOR_SAMPLES" \
            -thin "$THIN" \
            -polar "$POLAR" \
            -freq "$FREQ" \
            -num_cores "$NUM_CORES"
    else
        printf '[Skip] SINGER trees already exist for chr%s: %s\n' "$CHR" "$OUT_PREFIX"
    fi

    if [[ "$SKIP_EVAL" -eq 0 ]]; then
        "$PHLASH_PYTHON" \
            "$BENCH_ROOT/evaluate_singer_recent_pilot.py" \
            --glob "${OUT_PREFIX}_*.trees" \
            --mutation-rate "$MU" \
            --output-prefix "$DIAG_DIR/${BASE_PREFIX}_chr${CHR}" \
            --window-size "$TMRCA_WINDOW_SIZE"
    fi
done

printf '[Done] SINGER recent pilot finished for %s across chromosomes %s\n' "$POP" "$CHROMOSOMES"
