#!/bin/bash
set -euo pipefail

# ========= Customer configuration ========
VCF_DIR="/share/home/grp-wangyf/xuyuan/sds/plink/smc_vcfs"
VCF_PREFIX="UKBQC.for_smc"
SAMPLE_LIST="/share/home/grp-wangyf/xuyuan/sds/plink/final_samples_to_keep.list"
OUTDIR="/share/home/grp-wangyf/xuyuan/sds/data/processed/sds_input"
LOGDIR="$OUTDIR/logs"
SCRIPT_DIR="/share/home/grp-wangyf/xuyuan/sds/scripts"

# SDS parameters
SPATIAL_WINDOW=1000
MIN_DP=10
MIN_GQ=20
N_PARALLEL=8
FORCE=0
CHROMOSOMES=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--chr) CHROMOSOMES="$2"; shift 2 ;;
        -j|--parallel) N_PARALLEL="$2"; shift 2 ;;
        -t|--threads) shift 2 ;;
        -f|--force) FORCE=1; shift ;;
        --skip-env) shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

[[ -z "$CHROMOSOMES" ]] && CHROMOSOMES=$(seq 1 22)
mkdir -p "$OUTDIR/tmp" "$LOGDIR"

# ========= Process Function =========
process_chr() {
    local CHR=$1
    local PREFIX="chr${CHR}"
    local LOGFILE="$LOGDIR/${PREFIX}.log"
    
    # --- FIX: Strict file extension check ---
    local VCF_FILE=""
    if [[ -f "$VCF_DIR/${VCF_PREFIX}-chr${CHR}.vcf.gz" ]]; then
        VCF_FILE="$VCF_DIR/${VCF_PREFIX}-chr${CHR}.vcf.gz"
    elif [[ -f "$VCF_DIR/${VCF_PREFIX}-chr${CHR}.vcf.bgz" ]]; then
        VCF_FILE="$VCF_DIR/${VCF_PREFIX}-chr${CHR}.vcf.bgz"
    else
        echo "Error: VCF file not found for $CHR" >&2
        return 1
    fi
    # ----------------------------------------

    exec > >(tee "$LOGFILE") 2>&1
    echo "[$(date)] Processing $PREFIX using $(basename "$VCF_FILE")..." >&2

    # 1. b_file
    if [[ $FORCE -eq 1 ]] || [[ ! -f "$OUTDIR/${PREFIX}_b_file.txt" ]]; then
        bash "$SCRIPT_DIR/generate_b_file.sh" "$CHR" "$OUTDIR" "$PREFIX"
    fi

    # 2. s_file
    if [[ $FORCE -eq 1 ]] || [[ ! -f "$OUTDIR/${PREFIX}_s_file.txt" ]]; then
        echo "  > Generating s_file..." >&2
        # Pass empty region "" to force sequential read (bypassing broken/missing index)
        bash "$SCRIPT_DIR/generate_s_file.sh" \
            "$VCF_FILE" "$SAMPLE_LIST" "$OUTDIR" "$PREFIX" ""
    fi

    # 3. t_file
    if [[ $FORCE -eq 1 ]] || [[ ! -f "$OUTDIR/${PREFIX}_t_file.txt" ]]; then
        echo "  > Generating t_file..." >&2
        local SFILE_WITH_NAMES="$OUTDIR/tmp/${PREFIX}_s_file.with_names.txt"
        bash "$SCRIPT_DIR/generate_t_file.sh" \
            "$VCF_FILE" "$SAMPLE_LIST" "$SFILE_WITH_NAMES" \
            "$OUTDIR" "$PREFIX" "$SPATIAL_WINDOW" ""
    fi

    # 4. o_file
    if [[ $FORCE -eq 1 ]] || [[ ! -f "$OUTDIR/${PREFIX}_o_file.txt" ]]; then
         echo "  > Generating o_file..." >&2
         local SFILE_WITH_NAMES="$OUTDIR/tmp/${PREFIX}_s_file.with_names.txt"
         bash "$SCRIPT_DIR/generate_o_file.sh" \
            "$VCF_FILE" "$SAMPLE_LIST" "$SFILE_WITH_NAMES" \
            "$OUTDIR" "$PREFIX" "$MIN_DP" "$MIN_GQ" "" ""
    fi

    echo "[$(date)] $PREFIX Done" >&2
}
export -f process_chr
export VCF_DIR VCF_PREFIX SAMPLE_LIST OUTDIR LOGDIR SPATIAL_WINDOW MIN_DP MIN_GQ FORCE SCRIPT_DIR

# ========= Run =========
echo "Starting parallel processing on: $CHROMOSOMES"
echo $CHROMOSOMES | tr ' ' '\n' | \
    parallel -j $N_PARALLEL --line-buffer --tagstring "{}" \
    'process_chr {}'