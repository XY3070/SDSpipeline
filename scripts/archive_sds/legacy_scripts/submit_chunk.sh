#BSUB -J SDS_Final[1-22]
#BSUB -q q8358
#BSUB -n 14                  # 申请 14 个核 (14 * 22 = 308，填满队列)
#BSUB -R "span[hosts=1]"     # 确保这 14 个核在同一节点
#BSUB -o /share/home/grp-wangyf/xuyuan/sdSPY/scripts/logs/sds_%J_%I.out
#BSUB -e /share/home/grp-wangyf/xuyuan/sdSPY/scripts/logs/sds_%J_%I.err

# === 环境加载 ===
set +e 
echo "Job started on $(hostname) at $(date)"
[ -f ~/.bash_profile ] && source ~/.bash_profile || true
[ -f ~/.profile ] && source ~/.profile || true
export PS1="dummy"
[ -f ~/.bashrc ] && source ~/.bashrc || true
export PATH="$HOME/.cargo/bin:$HOME/.local/bin:$PATH"

if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate SDS || true
fi
set -euo pipefail

# ================= 配置 =================
N_THREADS=$LSB_DJOB_NUMPROC
[ -z "$N_THREADS" ] && N_THREADS=1

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

SCRIPT_DIR="/share/home/grp-wangyf/xuyuan/sdSPY/scripts"
IN_DIR="/share/home/grp-wangyf/xuyuan/sds/data/processed/sds_input"
OUT_DIR="/share/home/grp-wangyf/xuyuan/sds/data/processed/sds_output"
VCF_DIR="/share/home/grp-wangyf/xuyuan/sds/plink/smc_vcfs"
VCF_PREFIX="UKBQC.for_smc"
SAMPLE_LIST="/share/home/grp-wangyf/xuyuan/sds/plink/clean_sample_list.txt"

CHR=$LSB_JOBINDEX
PREFIX="chr${CHR}"

mkdir -p "$OUT_DIR/chunks/${PREFIX}"
cd "$SCRIPT_DIR"

echo "=========================================================="
echo "Processing ${PREFIX} with ${N_THREADS} cores"
echo "=========================================================="

# 1. 检查输入 VCF
if [[ -f "$VCF_DIR/${VCF_PREFIX}-${PREFIX}.vcf.gz" ]]; then
    VCF_FILE="$VCF_DIR/${VCF_PREFIX}-${PREFIX}.vcf.gz"
elif [[ -f "$VCF_DIR/${VCF_PREFIX}-${PREFIX}.vcf.bgz" ]]; then
    VCF_FILE="$VCF_DIR/${VCF_PREFIX}-${PREFIX}.vcf.bgz"
else
    echo "Error: VCF file not found"
    exit 1
fi

# 2. 检查 s_file 中间文件
SFILE_WITH_NAMES="$IN_DIR/tmp/${PREFIX}_s_file.with_names.txt"
if [[ ! -f "$SFILE_WITH_NAMES" ]]; then
    echo "Warning: '$SFILE_WITH_NAMES' missing. Regenerating..."
    bash generate_s_file.sh "$VCF_FILE" "$SAMPLE_LIST" "$IN_DIR" "$PREFIX" ""
fi

# 3. 重新生成 FULL t_file
# 注意：如果确定 t_file 已经存在且正确，可以注释掉下面这行以节省时间
echo "[Step 1] Regenerating FULL t_file..."
bash generate_t_file.sh "$VCF_FILE" "$SAMPLE_LIST" "$SFILE_WITH_NAMES" "$IN_DIR" "$PREFIX" "0" ""

# ================= 切分与并行 =================
T_FILE="$IN_DIR/${PREFIX}_t_file.txt"
S_FILE="$IN_DIR/${PREFIX}_s_file.txt"
O_FILE="$IN_DIR/${PREFIX}_o_file.txt"
B_FILE="$IN_DIR/${PREFIX}_b_file.txt"
G_FILE="$IN_DIR/g_file.txt"
RUN_SCRIPT="$SCRIPT_DIR/run_chunk_sds.sh"

CHUNK_DIR="$IN_DIR/chunks/${PREFIX}"
rm -f "$CHUNK_DIR"/* 
mkdir -p "$CHUNK_DIR"

echo "[Step 2] Splitting t_file into chunks (5000 lines)..."
split -l 5000 -d -a 4 "$T_FILE" "$CHUNK_DIR/chunk_"
CHUNKS_COUNT=$(ls "$CHUNK_DIR" | wc -l)
echo "Created $CHUNKS_COUNT chunks."

echo "[Step 3] Running SDS computation..."

# 使用 parallel 调用外部脚本
# {.} 表示去除扩展名的文件名 (这里 chunk_xxxx)
ls "$CHUNK_DIR"/chunk_* | parallel -j "$N_THREADS" \
    "$RUN_SCRIPT" \
    {} \
    "$OUT_DIR/chunks/${PREFIX}/{/.}.res" \
    "$S_FILE" \
    "$O_FILE" \
    "$B_FILE" \
    "$G_FILE"

# ================= 合并结果 =================
echo "[Step 4] Merging results..."
FINAL_OUT="$OUT_DIR/${PREFIX}_sds_res_FULL.txt"

# 写入 Header
echo -e "ID\tAA\tDA\tPOS\tDAF\tnG0\tnG1\tnG2\trSDS\tSuggestedInitPoint" > "$FINAL_OUT"

# 合并内容 (跳过第一行)
for res_file in $(ls "$OUT_DIR/chunks/${PREFIX}"/*.res | sort); do
    if [[ -s "$res_file" ]]; then
        tail -n +2 "$res_file" >> "$FINAL_OUT"
    fi
done

echo "Done. Final output: $FINAL_OUT"
# 清理
rm -rf "$CHUNK_DIR" "$OUT_DIR/chunks/${PREFIX}"