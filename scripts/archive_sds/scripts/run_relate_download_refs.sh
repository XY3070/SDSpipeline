#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck source=/data/home/grp-wangyf/xuyuan/sds/scripts/common_env.sh
source "$SCRIPT_DIR/common_env.sh"
activate_sds_env

OUT_ROOT="${SDS_RELATE_OUTPUT_ROOT:-$SDS_RESULTS_ROOT/legacy/relate_clues2}"
FORCE=0
DOWNLOAD_URL="https://zenodo.org/api/records/15801307/files/Relate_input_files.tgz/content"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --out-root)
            OUT_ROOT="$2"
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

REF_DIR="$OUT_ROOT/refs"
ARCHIVE_PATH="$REF_DIR/Relate_input_files.tgz"
EXTRACT_DIR="$REF_DIR/relate_input_files_zenodo_15801307"
MANIFEST_PATH="$REF_DIR/relate_input_files_manifest.txt"
LINKS_PATH="$REF_DIR/resolved_hg38_paths.tsv"

mkdir -p "$REF_DIR"

if [[ "$FORCE" -eq 1 ]]; then
    rm -rf "$EXTRACT_DIR"
    rm -f "$ARCHIVE_PATH" "$MANIFEST_PATH" "$LINKS_PATH"
fi

if [[ ! -f "$ARCHIVE_PATH" ]]; then
    tmp_archive="${ARCHIVE_PATH}.tmp"
    rm -f "$tmp_archive"
    curl -L --fail --output "$tmp_archive" "$DOWNLOAD_URL"
    mv "$tmp_archive" "$ARCHIVE_PATH"
fi

if [[ ! -d "$EXTRACT_DIR" ]]; then
    mkdir -p "$EXTRACT_DIR"
    tar -xzf "$ARCHIVE_PATH" -C "$EXTRACT_DIR"
fi

tar -tzf "$ARCHIVE_PATH" > "$MANIFEST_PATH"

resolve_chr_name() {
    local chr="$1"
    if [[ "$chr" == "23" ]]; then
        printf 'X\n'
    elif [[ "$chr" == "24" ]]; then
        printf 'Y\n'
    else
        printf '%s\n' "$chr"
    fi
}

ancestor_dir="$EXTRACT_DIR/Relate_input_files/GRCh38/human_ancestor_GRCh38"
mask_dir="$EXTRACT_DIR/Relate_input_files/GRCh38/20160622_genome_mask_GRCh38/StrictMask"
generic_coal_path="$EXTRACT_DIR/Relate_input_files/coal_rates/1000G_auto.coal"

[[ -d "$ancestor_dir" ]] || { echo "[Error] Missing ancestor dir: $ancestor_dir" >&2; exit 1; }
[[ -d "$mask_dir" ]] || { echo "[Error] Missing mask dir: $mask_dir" >&2; exit 1; }
[[ -f "$generic_coal_path" ]] || { echo "[Error] Missing generic coal file: $generic_coal_path" >&2; exit 1; }

{
    printf 'key\tpath\n'
    printf 'generic_coal\t%s\n' "$generic_coal_path"
} > "$LINKS_PATH"

ln -sfn "$generic_coal_path" "$REF_DIR/hg38_generic_1000G.coal"

for chr in $(seq 1 22) X Y; do
    chr_name="$(resolve_chr_name "$chr")"
    ancestor_path="$ancestor_dir/homo_sapiens_ancestor_${chr_name}.fa.gz"
    mask_path="$mask_dir/20160622.chr${chr_name}.mask.fasta.gz"

    [[ -f "$ancestor_path" ]] || { echo "[Error] Missing ancestor file: $ancestor_path" >&2; exit 1; }
    [[ -f "$mask_path" ]] || { echo "[Error] Missing mask file: $mask_path" >&2; exit 1; }

    ln -sfn "$ancestor_path" "$REF_DIR/hg38_ancestor_chr${chr_name}.fa.gz"
    ln -sfn "$mask_path" "$REF_DIR/hg38_mask_chr${chr_name}.fa.gz"

    {
        printf 'ancestor_chr%s\t%s\n' "$chr_name" "$ancestor_path"
        printf 'mask_chr%s\t%s\n' "$chr_name" "$mask_path"
    } >> "$LINKS_PATH"
done

printf '[Done] refs ready under %s\n' "$REF_DIR"
