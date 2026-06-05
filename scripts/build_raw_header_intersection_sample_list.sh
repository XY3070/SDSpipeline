#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$SCRIPT_DIR/common_env.sh"
activate_sds_env

RAW_VCF=""
TARGET_SAMPLE_LIST=""
OUT_FILE=""

usage() {
    cat >&2 <<'EOF'
Usage: build_raw_header_intersection_sample_list.sh --vcf PATH --sample-list PATH --output PATH
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf) RAW_VCF="$2"; shift 2 ;;
        --sample-list) TARGET_SAMPLE_LIST="$2"; shift 2 ;;
        --output) OUT_FILE="$2"; shift 2 ;;
        *) usage ;;
    esac
done

[[ -n "$RAW_VCF" && -n "$TARGET_SAMPLE_LIST" && -n "$OUT_FILE" ]] || usage
[[ -f "$RAW_VCF" ]] || { echo "[Error] raw VCF not found: $RAW_VCF" >&2; exit 1; }
[[ -f "$TARGET_SAMPLE_LIST" ]] || { echo "[Error] target sample list not found: $TARGET_SAMPLE_LIST" >&2; exit 1; }

TMPDIR="$(mktemp -d)"
trap '/bin/rm -rf "$TMPDIR"' EXIT

RAW_SORTED="$TMPDIR/raw.samples.sorted.txt"
TARGET_SORTED="$TMPDIR/target.samples.sorted.txt"

bcftools query -l "$RAW_VCF" | sort > "$RAW_SORTED"
gawk 'NF > 0 && $1 !~ /^#/ { print $1 }' "$TARGET_SAMPLE_LIST" | sort > "$TARGET_SORTED"

mkdir -p "$(dirname "$OUT_FILE")"
comm -12 "$RAW_SORTED" "$TARGET_SORTED" > "$OUT_FILE"

echo "output=$OUT_FILE"
echo "count=$(wc -l < "$OUT_FILE")"
