#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../scripts/common_env.sh
source "$SCRIPT_DIR/../scripts/common_env.sh"

# ========= Customer config =========
VCF_DIR="${SDS_VCF_ROOT:-$SDS_INPUT_ROOT/raw/vcf}"
VCF_PREFIX="UKBQC.for_smc"
SAMPLE_LIST="${SDS_SAMPLE_LIST_ROOT:-$SDS_INPUT_ROOT/freeze/sample_lists}/final_samples_to_keep.list"
OUTDIR="${SDS_SDS_INPUT_ROOT:-$SDS_RESULTS_ROOT/production/sds_input}"
LOGDIR="$OUTDIR/logs"

SPATIAL_WINDOW=1000
MIN_DP=10
MIN_GQ=20
O_FILE_MODE="fraction"
N_PARALLEL=8
NICE_VALUE=15

# ========= Customer parameters =========
FORCE=0
CHROMOSOMES=""
SKIP_ENV=0

while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--chr) CHROMOSOMES="$2"; shift 2 ;;
        -j|--parallel) N_PARALLEL="$2"; shift 2 ;;
        -w|--window) SPATIAL_WINDOW="$2"; shift 2 ;;
        -t|--threads) N_THREADS_PER_CHR="$2"; shift 2 ;;
        -f|--force) FORCE=1; shift ;;
        -n|--nice) NICE_VALUE="$2"; shift 2 ;;
        --vcf_dir) VCF_DIR="$2"; shift 2 ;;
        --vcf_prefix) VCF_PREFIX="$2"; shift 2 ;;
        --o-mode) O_FILE_MODE="$2"; shift 2 ;;
        --skip-env) SKIP_ENV=1; shift ;;
        -h|--help)
            echo "Usage: $0 [options]"
            echo "  -c, --chr CHRS      Chromosomes (default: 1-22)"
            echo "  -j, --parallel N    Number of parallel jobs (default: 8)"
            echo "  -f, --force         Force regeneration"
            echo "      --vcf_dir PATH  Override VCF directory (default: $VCF_DIR)"
            echo "      --vcf_prefix PX Override VCF prefix (default: $VCF_PREFIX)"
            echo "      --o-mode MODE   o_file mode: fraction|count|uniform (default: fraction)"
            echo "      --skip-env      Skip environment check"
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

[[ -z "$CHROMOSOMES" ]] && CHROMOSOMES=$(seq 1 22 | tr '\n' ' ')

# ========= Init environment =========
mkdir -p "$OUTDIR/tmp" "$LOGDIR"

if [[ -f ~/.bashrc ]]; then
    source ~/.bashrc
    mamba activate SDS 2>/dev/null || true
fi

cd "$SCRIPT_DIR"

# ========= Simplified environment check =========
echo "=========================================="
echo "Environment check"
echo "=========================================="

# Check tools
echo "Check required tools..."
command -v bcftools >/dev/null || { echo "[Error] bcftools is not installed"; exit 1; }
command -v gawk >/dev/null || { echo "[Error] gawk is not installed"; exit 1; }
command -v parallel >/dev/null || { echo "[Error] GNU parallel is not installed"; exit 1; }
echo "  ✓ All tools are installed"

# Check inputs
echo "Check input files..."
[[ -d "$VCF_DIR" ]] || { echo "[Error] VCF directory does not exist: $VCF_DIR"; exit 1; }
[[ -f "$SAMPLE_LIST" ]] || { echo "[Error] Sample list does not exist: $SAMPLE_LIST"; exit 1; }
echo "  ✓ Input files check passed"

# Check VCF files
if [[ $SKIP_ENV -eq 1 ]]; then
    echo "[Warning] Skipping VCF check"
else
echo "Check VCF files..."
echo "Directory: $VCF_DIR"
echo "Prefix: $VCF_PREFIX"
echo "Chromosomes: $CHROMOSOMES"
n_found=0
for chr in $CHROMOSOMES; do
    vcf_gz="$VCF_DIR/${VCF_PREFIX}-chr${chr}.vcf.gz"
    vcf_bgz="$VCF_DIR/${VCF_PREFIX}-chr${chr}.vcf.bgz"
    if [[ -f "$vcf_gz" ]] || [[ -f "$vcf_bgz" ]]; then
        n_found=$((n_found + 1))
    fi
done

if [[ $n_found -eq 0 ]]; then
    for f in "$VCF_DIR"/*-chr*.vcf.gz "$VCF_DIR"/*-chr*.vcf.bgz; do
        [[ -f "$f" ]] || continue
        base=$(basename "$f")
        auto_prefix=${base%%-chr*}
        echo "  ✓ Found VCF file: $base (prefix: $auto_prefix)"
        VCF_PREFIX="$auto_prefix"
        n_found=0
        for chr in $CHROMOSOMES; do
            vcf_gz="$VCF_DIR/${VCF_PREFIX}-chr${chr}.vcf.gz"
            vcf_bgz="$VCF_DIR/${VCF_PREFIX}-chr${chr}.vcf.bgz"
            if [[ -f "$vcf_gz" ]] || [[ -f "$vcf_bgz" ]]; then
                n_found=$((n_found + 1))
            fi
        done
        [[ $n_found -gt 0 ]] && { echo "  ✓ Auto-detected prefix: $VCF_PREFIX (matched $n_found chromosomes)"; break; }
    done
    if [[ $n_found -eq 0 ]]; then
        echo "[Error] No VCF files found"
        echo "Expected format: $VCF_DIR/${VCF_PREFIX}-chr*.vcf.gz or .vcf.bgz"
        echo "Actual files (first 5 .vcf.gz):"
        ls -lh "$VCF_DIR"/*.vcf.gz 2>/dev/null | head -5 || echo "  (none)"
        echo "Actual files (first 5 .vcf.bgz):"
        ls -lh "$VCF_DIR"/*.vcf.bgz 2>/dev/null | head -5 || echo "  (none)"
        echo "Directory listing (first 10 items):"
        ls -lh "$VCF_DIR" 2>/dev/null | head -10 || echo "  (directory not readable or empty)"
        exit 1
    fi
fi
echo "  ✓ Found $n_found VCF files"
fi

# 样本数
n_samples=$(grep -c -v '^#' "$SAMPLE_LIST" || true)
echo "  ✓ Sample count: $n_samples"

echo ""
echo "Configuration:"
            echo "  - Chromosomes: $CHROMOSOMES"
            echo "  - Parallel jobs: $N_PARALLEL"
            echo "  - Threads per chromosome: ${N_THREADS_PER_CHR:-1}"
            echo "  - Nice value: $NICE_VALUE"
            echo "  - o_file mode: $O_FILE_MODE"
            echo ""

PROGRESS_LOG="$LOGDIR/progress.log"
echo "[$(date)] === SDS input generation started ===" > "$PROGRESS_LOG"

# ========= process_chr function =========
process_chr() {
    local CHR=$1
    local VCF_FILE="$VCF_DIR/${VCF_PREFIX}-chr${CHR}.vcf.gz"
    local PREFIX="chr${CHR}"
    local LOGFILE="$LOGDIR/${PREFIX}.log"
    
    {
        echo "=========================================="
        echo "[$(date)] Start processing $PREFIX"
        echo "=========================================="
        
        if [[ ! -f "$VCF_FILE" ]]; then
            echo "[Error] VCF does not exist: $VCF_FILE"
            echo "[$(date)] $PREFIX FAILED (no VCF)" >> "$PROGRESS_LOG"
            return 1
        fi
        
        # Check if already processed
        if [[ $FORCE -eq 0 ]]; then
            if [[ -f "$OUTDIR/${PREFIX}_s_file.txt" ]] && \
               [[ -f "$OUTDIR/${PREFIX}_t_file.txt" ]] && \
               [[ -f "$OUTDIR/${PREFIX}_o_file.txt" ]] && \
               [[ -f "$OUTDIR/${PREFIX}_b_file.txt" ]]; then
                echo "[Skip] All files already exist"
                echo "[$(date)] $PREFIX SKIPPED" >> "$PROGRESS_LOG"
                return 0
            fi
        fi
        
        local START_TIME=$(date +%s)

        echo ""
        echo "[1/5] b_file"
        if [[ $FORCE -eq 0 ]] && [[ -f "$OUTDIR/${PREFIX}_b_file.txt" ]]; then
            echo "[Skip] b_file already exists"
        else
            nice -n $NICE_VALUE bash "$SCRIPT_DIR/generate_b_file.sh" \
                "$CHR" "$OUTDIR" "$PREFIX" || {
                echo "[Error] b_file generation failed"
                echo "[$(date)] $PREFIX FAILED (b_file)" >> "$PROGRESS_LOG"
                return 1
            }
        fi

        # read p/q arm boundaries
        if [[ -f "$OUTDIR/${PREFIX}_b_file.txt" ]]; then
            read P_START P_END < <(head -1 "$OUTDIR/${PREFIX}_b_file.txt")
            read Q_START Q_END < <(sed -n '2p' "$OUTDIR/${PREFIX}_b_file.txt")
        else
            echo "[Error] b_file not found after generation"
            echo "[$(date)] $PREFIX FAILED (b_file missing)" >> "$PROGRESS_LOG"
            return 1
        fi

        # detect contig id
        CONTIG_ID=$(bcftools view -h "$VCF_FILE" | gawk -v chr="$CHR" 'match($0, /##contig=<ID=([^,>]+)/, a){id=a[1]; if(id==chr){print id; exit} else if(id=="chr" chr){print id; exit}} END{ }')
        [[ -z "$CONTIG_ID" ]] && CONTIG_ID="$CHR"
        REG_P="$CONTIG_ID:${P_START}-${P_END}"
        REG_Q="$CONTIG_ID:${Q_START}-${Q_END}"

        echo ""
        echo "[2/5] s_file (p/q arms in parallel)"
        if [[ $FORCE -eq 0 ]] && [[ -f "$OUTDIR/${PREFIX}_s_file.txt" ]] && [[ -f "$OUTDIR/tmp/${PREFIX}_s_file.with_names.txt" ]]; then
            echo "[Skip] s_file already exists"
        else
            nice -n $NICE_VALUE bash "$SCRIPT_DIR/generate_s_file.sh" \
                "$VCF_FILE" "$SAMPLE_LIST" "$OUTDIR" "${PREFIX}_p" "$REG_P" &
            pid_p=$!
            nice -n $NICE_VALUE bash "$SCRIPT_DIR/generate_s_file.sh" \
                "$VCF_FILE" "$SAMPLE_LIST" "$OUTDIR" "${PREFIX}_q" "$REG_Q" &
            pid_q=$!
            wait $pid_p || { echo "[Error] s_file p-arm failed"; echo "[$(date)] $PREFIX FAILED (s_file p)" >> "$PROGRESS_LOG"; return 1; }
            wait $pid_q || { echo "[Error] s_file q-arm failed"; echo "[$(date)] $PREFIX FAILED (s_file q)" >> "$PROGRESS_LOG"; return 1; }

            gawk -v OFS="\t" 'FNR==NR{names[NR]=$1; next} NR==FNR+0{}' "$SAMPLE_LIST" >/dev/null 2>&1 || true

            gawk -v OFS="\t" -v slist="$SAMPLE_LIST" '
            FNR==NR{p[$1]=($2=="NA"?"NA":substr($0, index($0,$2))); next}
            {
                q[$1]=($2=="NA"?"NA":substr($0, index($0,$2)))
            }
            END{
                while((getline n<slist)>0){ if(n~ /^\s*$/ || n~ /^#/){continue} name=n; 
                    ps=(name in p)?p[name]:"NA"; qs=(name in q)?q[name]:"NA";
                    if(ps=="NA" && qs=="NA"){ print name, "NA" }
                    else if(ps=="NA"){ print name, qs }
                    else if(qs=="NA"){ print name, ps }
                    else { print name, ps "\t" qs }
                }
            }
            ' "$OUTDIR/tmp/${PREFIX}_p_s_file.with_names.txt" "$OUTDIR/tmp/${PREFIX}_q_s_file.with_names.txt" > "$OUTDIR/tmp/${PREFIX}_s_file.with_names.txt"

            gawk '{$1=""; sub(/^[ \t]+/, ""); print}' OFS="\t" "$OUTDIR/tmp/${PREFIX}_s_file.with_names.txt" > "$OUTDIR/${PREFIX}_s_file.txt"
        fi
        
        echo ""
        echo "[3/5] t_file (p/q arms in parallel)"
        if [[ $FORCE -eq 0 ]] && [[ -f "$OUTDIR/${PREFIX}_t_file.txt" ]]; then
            echo "[Skip] t_file already exists"
        else
            nice -n $NICE_VALUE bash "$SCRIPT_DIR/generate_t_file.sh" \
                "$VCF_FILE" "$SAMPLE_LIST" "$OUTDIR/tmp/${PREFIX}_s_file.with_names.txt" \
                "$OUTDIR" "${PREFIX}_p" "$SPATIAL_WINDOW" "$REG_P" &
            pid_tp=$!
            nice -n $NICE_VALUE bash "$SCRIPT_DIR/generate_t_file.sh" \
                "$VCF_FILE" "$SAMPLE_LIST" "$OUTDIR/tmp/${PREFIX}_s_file.with_names.txt" \
                "$OUTDIR" "${PREFIX}_q" "$SPATIAL_WINDOW" "$REG_Q" &
            pid_tq=$!
            wait $pid_tp || { echo "[Error] t_file p-arm failed"; echo "[$(date)] $PREFIX FAILED (t_file p)" >> "$PROGRESS_LOG"; return 1; }
            wait $pid_tq || { echo "[Error] t_file q-arm failed"; echo "[$(date)] $PREFIX FAILED (t_file q)" >> "$PROGRESS_LOG"; return 1; }
            cat "$OUTDIR/${PREFIX}_p_t_file.txt" "$OUTDIR/${PREFIX}_q_t_file.txt" > "$OUTDIR/${PREFIX}_t_file.txt"
        fi
        
        echo ""
        echo "[4/5] o_file"
        if [[ $FORCE -eq 0 ]] && [[ -f "$OUTDIR/${PREFIX}_o_file.txt" ]]; then
            echo "[Skip] o_file already exists"
        else
            nice -n $NICE_VALUE bash "$SCRIPT_DIR/generate_o_file.sh" \
                "$VCF_FILE" "$SAMPLE_LIST" "$OUTDIR/tmp/${PREFIX}_s_file.with_names.txt" \
                "$OUTDIR" "$PREFIX" "$MIN_DP" "$MIN_GQ" "$O_FILE_MODE" || {
                echo "[Error] o_file generation failed"
                echo "[$(date)] $PREFIX FAILED (o_file)" >> "$PROGRESS_LOG"
                return 1
            }
        fi
        
        echo ""
        echo "[5/5] b_file (already generated)"
        
        local END_TIME=$(date +%s)
        local DURATION=$((END_TIME - START_TIME))
        
        echo ""
        echo "✓ $PREFIX processed successfully (${DURATION} s)"
        echo "[$(date)] $PREFIX SUCCESS (${DURATION}s)" >> "$PROGRESS_LOG"
        
    } >"$LOGFILE" 2>&1
}

export -f process_chr
export VCF_DIR VCF_PREFIX SAMPLE_LIST OUTDIR LOGDIR SPATIAL_WINDOW MIN_DP MIN_GQ O_FILE_MODE FORCE NICE_VALUE SCRIPT_DIR N_THREADS_PER_CHR

# ========= Parallel processing =========
echo "=========================================="
echo "Parallel processing started"
echo "=========================================="

echo $CHROMOSOMES | tr ' ' '\n' | \
    parallel -j $N_PARALLEL \
        --progress \
        --joblog "$LOGDIR/parallel.log" \
        --halt-on-error 0 \
        'process_chr {}'

# ========= Result summary =========
echo ""
echo "=========================================="
echo "Result summary"
echo "=========================================="

n_success=0
n_failed=0
n_skipped=0

for chr in $CHROMOSOMES; do
    if grep -q "chr${chr} SUCCESS" "$PROGRESS_LOG" 2>/dev/null; then
        n_success=$((n_success + 1))
    elif grep -q "chr${chr} SKIPPED" "$PROGRESS_LOG" 2>/dev/null; then
        n_skipped=$((n_skipped + 1))
    else
        n_failed=$((n_failed + 1))
    fi
done

echo "Success: $n_success"
echo "Skipped: $n_skipped"
echo "Failed: $n_failed"

if [[ $n_failed -gt 0 ]]; then
    echo ""
    echo "see $LOGDIR/parallel.log for details"
    exit 1
fi

echo ""
echo "✓ All chromosomes processed successfully! Output directory: $OUTDIR"
