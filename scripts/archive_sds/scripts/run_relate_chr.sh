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
THREADS="${THREADS:-16}"
NE="${NE:-30000}"
MU="${MU:-1.25e-8}"
OUT_ROOT="$BASE_DIR/data/processed/relate_clues2"
FORCE=0
MAP_PATH=""

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
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --Ne)
            NE="$2"
            shift 2
            ;;
        --mu)
            MU="$2"
            shift 2
            ;;
        --out-root)
            OUT_ROOT="$2"
            shift 2
            ;;
        --map)
            MAP_PATH="$2"
            shift 2
            ;;
        --force)
            FORCE=1
            shift
            ;;
        *)
            echo "Unknown parameter: $1" >&2
            exit 1
            ;;
    esac
done

if [[ -z "$POP" || -z "$CHR" ]]; then
    echo "Usage: $0 --pop POP --chr CHR [--threads N] [--Ne HAPLOID_NE] [--mu MU] [--map FILE] [--out-root DIR] [--force]" >&2
    exit 1
fi

PREP_PREFIX="$OUT_ROOT/$POP/prepared/UKBQC_${POP}_chr${CHR}"
OUT_DIR="$OUT_ROOT/$POP/relate"
OUT_BASENAME="UKBQC_${POP}_chr${CHR}"
OUT_PREFIX="$OUT_DIR/$OUT_BASENAME"
if [[ -z "$MAP_PATH" ]]; then
    MAP_PATH="$OUT_ROOT/refs/relate_input_files_zenodo_15801307/Relate_input_files/GRCh38/recomb_map/genetic_map_chr${CHR}.txt"
fi
mkdir -p "$OUT_DIR"

[[ -f "${PREP_PREFIX}.haps.gz" ]] || { echo "[Error] Missing prepared haps: ${PREP_PREFIX}.haps.gz" >&2; exit 1; }
[[ -f "${PREP_PREFIX}.sample.gz" ]] || { echo "[Error] Missing prepared sample: ${PREP_PREFIX}.sample.gz" >&2; exit 1; }
[[ -f "${PREP_PREFIX}.annot" ]] || { echo "[Error] Missing prepared annot: ${PREP_PREFIX}.annot" >&2; exit 1; }
[[ -f "$MAP_PATH" ]] || { echo "[Error] Missing map file: $MAP_PATH" >&2; exit 1; }

normalize_annot_file "${PREP_PREFIX}.annot"

if [[ "$FORCE" -eq 1 ]]; then
    rm -f "${OUT_PREFIX}.anc" "${OUT_PREFIX}.mut"
fi

if [[ -f "${OUT_PREFIX}.anc" && -f "${OUT_PREFIX}.mut" && "$FORCE" -ne 1 ]]; then
    printf '[Skip] Relate output already exists: %s\n' "$OUT_PREFIX"
    exit 0
fi

(
    cd "$OUT_DIR"
    "$RELATE_DIR/scripts/RelateParallel/RelateParallel.sh" \
        --haps "${PREP_PREFIX}.haps.gz" \
        --sample "${PREP_PREFIX}.sample.gz" \
        --map "$MAP_PATH" \
        --dist "${PREP_PREFIX}.dist.gz" \
        --annot "${PREP_PREFIX}.annot" \
        -m "$MU" \
        -N "$NE" \
        -o "$OUT_BASENAME" \
        --threads "$THREADS"
)

printf '[Done] Relate chromosome run finished: %s\n' "$OUT_PREFIX"
