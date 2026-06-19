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
OUT_ROOT="$BASE_DIR/data/processed/relate_clues2"
FORCE=0
ANCESTOR_PATH=""
MASK_PATH=""

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
        --out-root)
            OUT_ROOT="$2"
            shift 2
            ;;
        --ancestor)
            ANCESTOR_PATH="$2"
            shift 2
            ;;
        --mask)
            MASK_PATH="$2"
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
    echo "Usage: $0 --pop POP --chr CHR [--out-root DIR] [--ancestor FILE] [--mask FILE] [--force]" >&2
    exit 1
fi

PHASED_PREFIX="$BASE_DIR/data/vcf/$POP/shapeit5/UKBQC_${POP}_chr${CHR}.phased"
PHASED_VCF="${PHASED_PREFIX}.vcf.gz"
SAMPLE_LIST="$BASE_DIR/data/${POP}.txt"
REF_DIR="$OUT_ROOT/refs"
POP_DIR="$OUT_ROOT/$POP"
RAW_DIR="$POP_DIR/raw_input"
PREP_DIR="$POP_DIR/prepared"
LOG_DIR="$BASE_DIR/logs/relate_clues2"
TMP_DIR="$POP_DIR/tmp/chr${CHR}_prepare"
POPLABELS_PATH="$POP_DIR/${POP}.poplabels"
RAW_PREFIX="$RAW_DIR/UKBQC_${POP}_chr${CHR}"
PREP_PREFIX="$PREP_DIR/UKBQC_${POP}_chr${CHR}"

mkdir -p "$RAW_DIR" "$PREP_DIR" "$LOG_DIR" "$TMP_DIR" "$POP_DIR"
trap '/bin/rm -rf "$TMP_DIR"' EXIT

[[ -f "$PHASED_VCF" ]] || { echo "[Error] phased VCF not found: $PHASED_VCF" >&2; exit 1; }
[[ -f "$SAMPLE_LIST" ]] || { echo "[Error] sample list not found: $SAMPLE_LIST" >&2; exit 1; }

if [[ -z "$ANCESTOR_PATH" ]]; then
    for candidate in "$REF_DIR"/hg38_ancestor_chr"${CHR}".fa.gz "$REF_DIR"/hg38_ancestor_chr"${CHR}".fa; do
        if [[ -f "$candidate" || -L "$candidate" ]]; then
            ANCESTOR_PATH="$candidate"
            break
        fi
    done
fi

if [[ -z "$MASK_PATH" ]]; then
    for candidate in "$REF_DIR"/hg38_mask_chr"${CHR}".fa.gz "$REF_DIR"/hg38_mask_chr"${CHR}".fa "$REF_DIR"/hg38_mask_chr"${CHR}".fasta.gz "$REF_DIR"/hg38_mask_chr"${CHR}".fasta; do
        if [[ -f "$candidate" || -L "$candidate" ]]; then
            MASK_PATH="$candidate"
            break
        fi
    done
fi

[[ -n "$ANCESTOR_PATH" ]] || { echo "[Error] ancestor reference not found under $REF_DIR" >&2; exit 1; }
[[ -n "$MASK_PATH" ]] || { echo "[Error] mask reference not found under $REF_DIR" >&2; exit 1; }

materialize_reference() {
    local input_path="$1"
    if [[ "$input_path" == *.gz ]]; then
        local out_path="${input_path%.gz}"
        if [[ ! -f "$out_path" || "$FORCE" -eq 1 ]]; then
            gzip -cd "$input_path" > "$out_path"
        fi
        printf '%s\n' "$out_path"
        return 0
    fi
    printf '%s\n' "$input_path"
}

ANCESTOR_PATH="$(materialize_reference "$ANCESTOR_PATH")"
MASK_PATH="$(materialize_reference "$MASK_PATH")"

build_poplabels() {
    local vcf_path="$1"
    local sample_list_path="$2"
    local pop="$3"
    local out_path="$4"

    local vcf_ids="$TMP_DIR/vcf_ids.txt"
    local vcf_sorted="$TMP_DIR/vcf_ids.sorted.txt"
    local sample_sorted="$TMP_DIR/sample_ids.sorted.txt"

    bcftools query -l "$vcf_path" > "$vcf_ids"
    sort "$vcf_ids" > "$vcf_sorted"
    sort "$sample_list_path" > "$sample_sorted"

    if ! diff -u "$sample_sorted" "$vcf_sorted" > "$TMP_DIR/${pop}_ids.diff"; then
        echo "[Error] Sample list and VCF IDs differ for $pop. See $TMP_DIR/${pop}_ids.diff" >&2
        exit 1
    fi

    {
        printf 'sample\tpopulation\tgroup\tsex\n'
        awk -v pop="$pop" '{print $1 "\t" pop "\t" pop "\t0"}' "$vcf_ids"
    } > "$out_path"
}

if [[ ! -f "$POPLABELS_PATH" || "$FORCE" -eq 1 ]]; then
    build_poplabels "$PHASED_VCF" "$SAMPLE_LIST" "$POP" "$POPLABELS_PATH"
fi

if [[ "$FORCE" -eq 1 ]]; then
    rm -f "${RAW_PREFIX}.haps" "${RAW_PREFIX}.sample"
    rm -f "${PREP_PREFIX}.haps.gz" "${PREP_PREFIX}.sample.gz" "${PREP_PREFIX}.dist.gz" "${PREP_PREFIX}.annot"
fi

if [[ ! -f "${RAW_PREFIX}.haps" || ! -f "${RAW_PREFIX}.sample" ]]; then
    "$RELATE_DIR/bin/RelateFileFormats" \
        --mode ConvertFromVcf \
        --haps "${RAW_PREFIX}.haps" \
        --sample "${RAW_PREFIX}.sample" \
        -i "$PHASED_PREFIX"
fi

if [[ ! -f "${PREP_PREFIX}.haps.gz" || ! -f "${PREP_PREFIX}.sample.gz" || ! -f "${PREP_PREFIX}.annot" ]]; then
    "$RELATE_DIR/scripts/PrepareInputFiles/PrepareInputFiles.sh" \
        --haps "${RAW_PREFIX}.haps" \
        --sample "${RAW_PREFIX}.sample" \
        --ancestor "$ANCESTOR_PATH" \
        --mask "$MASK_PATH" \
        --poplabels "$POPLABELS_PATH" \
        -o "$PREP_PREFIX"
fi

normalize_annot_file "${PREP_PREFIX}.annot"

printf '[Done] Relate input prepared: %s\n' "$PREP_PREFIX"
