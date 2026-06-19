#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck source=/share/home/grp-wangyf/xuyuan/sds/scripts/common_env.sh
source "$SCRIPT_DIR/common_env.sh"
activate_sds_env

require_nonempty_file() {
    local path="$1"
    local label="$2"
    if [[ ! -s "$path" ]]; then
        echo "[Error] ${label} missing or empty: $path" >&2
        exit 1
    fi
}

count_non_na_rows() {
    local path="$1"
    gawk 'BEGIN{n=0} $0 != "NA" {n++} END{print n}' "$path"
}

write_input_audit_sidecars() {
    local sample_order="$1"
    local s_file="$2"
    local t_file="$3"
    local o_file="$4"
    local outdir="$5"
    local chr="$6"

    local sample_order_out="$outdir/chr${chr}_sample_order.txt"
    local s_audit="$outdir/chr${chr}_s_file.audit.tsv"
    local t_audit="$outdir/chr${chr}_t_file.audit.tsv"
    local o_audit="$outdir/chr${chr}_o_file.audit.tsv"
    local t_header="$outdir/tmp/chr${chr}_t_header.tmp"

    mkdir -p "$outdir/tmp"
    cp "$sample_order" "$sample_order_out"

    {
        printf 'sample_id\tsingleton_positions...\n'
        paste "$sample_order_out" "$s_file"
    } > "$s_audit"

    {
        printf 'ID\tAA\tDA\tPOS'
        while IFS= read -r sample; do
            printf '\t%s' "$sample"
        done < "$sample_order_out"
        printf '\n'
    } > "$t_header"
    cat "$t_header" "$t_file" > "$t_audit"

    {
        paste -sd'\t' "$sample_order_out"
        cat "$o_file"
    } > "$o_audit"
}

# ========= 默认配置 =========
POP=""
CHR=""
OUT_ROOT="$BASE_DIR/data/processed/sds_input"
VCF_OVERRIDE=""
SAMPLE_LIST_OVERRIDE=""
OUT_POP=""
DO_S=true; DO_T=true; DO_O=true; DO_B=true; TEST_MODE=false
FORCE=0
EMIT_AUDIT_SIDECARS=1

# ========= 参数解析 =========
while [[ $# -gt 0 ]]; do
    case $1 in
        --pop) POP="$2"; shift 2 ;;
        --chr) CHR="$2"; shift 2 ;;
        --vcf) VCF_OVERRIDE="$2"; shift 2 ;;
        --sample-list) SAMPLE_LIST_OVERRIDE="$2"; shift 2 ;;
        --out-pop) OUT_POP="$2"; shift 2 ;;
        --out-root) OUT_ROOT="$2"; shift 2 ;;
        --skip-s) DO_S=false; shift ;;
        --skip-t) DO_T=false; shift ;;
        --skip-o) DO_O=false; shift ;;
        --skip-b) DO_B=false; shift ;;
        --test) TEST_MODE=true; shift ;;
        --force) FORCE=1; shift ;;
        --skip-audit-sidecars) EMIT_AUDIT_SIDECARS=0; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
done

# ========= 路径推断 =========
if [[ -z "$CHR" ]]; then
    echo "Usage: $0 --chr CHR [--pop POP] [--vcf PATH --sample-list PATH --out-pop POP] [--out-root DIR] [--skip-s] [--skip-t] [--skip-o] [--skip-b] [--test] [--force] [--skip-audit-sidecars]" >&2
    exit 1
fi

if [[ -z "$OUT_POP" ]]; then
    OUT_POP="$POP"
fi

if [[ -z "$OUT_POP" ]]; then
    echo "[Error] Need --pop POP or --out-pop POP" >&2
    exit 1
fi

if [[ -n "$VCF_OVERRIDE" ]]; then
    VCF_FILE="$VCF_OVERRIDE"
else
    [[ -n "$POP" ]] || { echo "[Error] --pop POP is required when --vcf is not provided" >&2; exit 1; }
    VCF_FILE="$BASE_DIR/data/vcf/$POP/UKBQC_${POP}_chr${CHR}.vcf.gz"
    if [[ ! -f "$VCF_FILE" ]]; then
        ALT_VCF="$BASE_DIR/plink/vcf_output/$POP/UKBQC_${POP}_chr${CHR}.vcf.gz"
        [[ -f "$ALT_VCF" ]] && VCF_FILE="$ALT_VCF"
    fi
fi

if [[ -n "$SAMPLE_LIST_OVERRIDE" ]]; then
    SAMPLE_LIST="$SAMPLE_LIST_OVERRIDE"
else
    [[ -n "$POP" ]] || { echo "[Error] --sample-list PATH is required when --pop is not provided" >&2; exit 1; }
    SAMPLE_LIST="$BASE_DIR/data/${POP}.txt"
    if [[ ! -f "$SAMPLE_LIST" ]]; then
        ALT_SLIST="$BASE_DIR/data/metadata/${POP}.txt"
        [[ -f "$ALT_SLIST" ]] && SAMPLE_LIST="$ALT_SLIST"
    fi
fi

if [ "$TEST_MODE" = true ]; then
    OUTDIR="$OUT_ROOT/test/$OUT_POP"
else
    OUTDIR="$OUT_ROOT/$OUT_POP"
fi
LOGDIR="$OUTDIR/logs"
MY_TMP="$OUTDIR/tmp/chr${CHR}_$(date +%s)_$RANDOM"
mkdir -p "$MY_TMP" "$OUTDIR" "$LOGDIR"

trap '/bin/rm -rf "$MY_TMP"' EXIT

command -v gawk >/dev/null || { echo "[Error] gawk is not installed" >&2; exit 1; }
[[ -f "$VCF_FILE" ]] || { echo "[Error] VCF not found: $VCF_FILE" >&2; exit 1; }
[[ -f "$SAMPLE_LIST" ]] || { echo "[Error] Sample list not found: $SAMPLE_LIST" >&2; exit 1; }
if [[ "$DO_S" = true || "$DO_T" = true ]]; then
    command -v bcftools >/dev/null || { echo "[Error] bcftools is not installed" >&2; exit 1; }
fi

validate_sample_list_against_vcf() {
    local sample_list="$1"
    local vcf_file="$2"
    local tmp_vcf_samples="$MY_TMP/vcf_samples.txt"
    local tmp_missing="$MY_TMP/missing_samples.txt"

    bcftools query -l "$vcf_file" > "$tmp_vcf_samples"
    require_nonempty_file "$tmp_vcf_samples" "VCF sample list"

    gawk '
    NR == FNR { have[$1] = 1; next }
    { if (!($1 in have)) print $1 }
    ' "$tmp_vcf_samples" "$sample_list" > "$tmp_missing"

    if [[ -s "$tmp_missing" ]]; then
        echo "[Error] Sample list contains names absent from VCF header: $sample_list" >&2
        sed -n '1,10p' "$tmp_missing" >&2
        exit 1
    fi
}

echo ">>> SDS env: /data/home/grp-wangyf/intern/miniforge3/envs/sds"
echo ">>> Population label: $OUT_POP | Chromosome: $CHR"
echo ">>> VCF: $VCF_FILE"
echo ">>> Sample list: $SAMPLE_LIST"

CLEAN_SAMPLE_ORDER="$MY_TMP/sample_order.txt"
gawk 'NF > 0 && $1 !~ /^#/ { print $1 }' "$SAMPLE_LIST" > "$CLEAN_SAMPLE_ORDER"
require_nonempty_file "$CLEAN_SAMPLE_ORDER" "sample order"
validate_sample_list_against_vcf "$CLEAN_SAMPLE_ORDER" "$VCF_FILE"

S_OUT="$OUTDIR/chr${CHR}_s_file.txt"
T_OUT="$OUTDIR/chr${CHR}_t_file.txt"
O_OUT="$OUTDIR/chr${CHR}_o_file.txt"

# 1. 生成边界文件 (b_file)
if [ "$DO_B" = true ] && { [ "$FORCE" -eq 1 ] || [ ! -f "$OUTDIR/chr${CHR}_b_file.txt" ]; }; then
    bash "$SCRIPT_DIR/generate_b_file.sh" "$CHR" "$OUTDIR" "chr${CHR}"
elif [ "$DO_B" = true ]; then
    echo "[Skip] b_file exists. Use --force to regenerate."
fi

# 获取分臂区间 (从 b_file 读取)
B_FILE="$OUTDIR/chr${CHR}_b_file.txt"
[[ -f "$B_FILE" ]] || { echo "[Error] b_file not found: $B_FILE" >&2; exit 1; }
if [[ "$(wc -l < "$B_FILE")" -lt 2 ]]; then
    echo "[Error] b_file malformed: $B_FILE" >&2
    exit 1
fi
read -r P_START P_END < <(sed -n '1p' "$B_FILE")
read -r Q_START Q_END < <(sed -n '2p' "$B_FILE")
P_REG="chr${CHR}:${P_START}-${P_END}"
Q_REG="chr${CHR}:${Q_START}-${Q_END}"

# 如果是测试模式，覆盖区间为前 2Mb
if [ "$TEST_MODE" = true ]; then
    TEST_START="$(find_first_variant_pos "$VCF_FILE" "chr${CHR}")"
    if [[ -z "$TEST_START" ]]; then
        echo "[Error] No variants found on chr${CHR} in $VCF_FILE" >&2
        exit 1
    fi
    TEST_END=$((TEST_START + 2000000 - 1))
    if (( TEST_END > Q_END )); then
        TEST_END="$Q_END"
    fi
    echo "[TEST] Mode enabled: processing 2Mb from first variant"
    echo "[TEST] Region: chr${CHR}:${TEST_START}-${TEST_END}"
    P_REG="chr${CHR}:${TEST_START}-${TEST_END}"; Q_REG="NONE"
fi

# 2. 并行处理 p 臂和 q 臂 (Trunking)
# 这里启动 & 后台运行，实现单个 Job 内的 2 CPU 并行
echo ">>> Extracting s_file and t_file in parallel arms..."

if [ "$DO_S" = true ] && { [ "$FORCE" -eq 1 ] || [ ! -f "$S_OUT" ]; }; then
    s_p_pid=""
    s_q_pid=""
    bash "$SCRIPT_DIR/generate_s_file.sh" "$VCF_FILE" "$SAMPLE_LIST" "$MY_TMP" "p" "$P_REG" &
    s_p_pid=$!
    if [ "$Q_REG" != "NONE" ]; then
        bash "$SCRIPT_DIR/generate_s_file.sh" "$VCF_FILE" "$SAMPLE_LIST" "$MY_TMP" "q" "$Q_REG" &
        s_q_pid=$!
    fi
    wait "$s_p_pid"
    if [ -n "$s_q_pid" ]; then
        wait "$s_q_pid"
    fi

    if [ "$Q_REG" = "NONE" ]; then
        cp "$MY_TMP/s_p.with_names.txt" "$MY_TMP/chr${CHR}_s_file.with_names.txt"
    else
        gawk -v OFS="\t" -v slist="$SAMPLE_LIST" '
        FNR == NR { p[$1] = ($2 == "NA" ? "NA" : substr($0, index($0, $2))); next }
        { q[$1] = ($2 == "NA" ? "NA" : substr($0, index($0, $2))) }
        END {
            while ((getline line < slist) > 0) {
                if (line ~ /^\s*$/ || line ~ /^#/) continue
                split(line, fields, /[[:space:]]+/)
                name = fields[1]
                ps = (name in p) ? p[name] : "NA"
                qs = (name in q) ? q[name] : "NA"
                if (ps == "NA" && qs == "NA") print name, "NA"
                else if (ps == "NA") print name, qs
                else if (qs == "NA") print name, ps
                else print name, ps "\t" qs
            }
            close(slist)
        }' "$MY_TMP/s_p.with_names.txt" "$MY_TMP/s_q.with_names.txt" > "$MY_TMP/chr${CHR}_s_file.with_names.txt"
    fi

    gawk '{$1=""; sub(/^[ \t]+/, ""); print}' OFS="\t" \
        "$MY_TMP/chr${CHR}_s_file.with_names.txt" > "$S_OUT"
    require_nonempty_file "$MY_TMP/chr${CHR}_s_file.with_names.txt" "s_file intermediate"
    require_nonempty_file "$S_OUT" "s_file"
elif [ "$DO_S" = true ]; then
    echo "[Skip] s_file exists. Use --force to regenerate."
    require_nonempty_file "$S_OUT" "s_file"
fi

if [ -f "$S_OUT" ]; then
    non_na_rows="$(count_non_na_rows "$S_OUT")"
    if [ "${non_na_rows:-0}" -eq 0 ]; then
        echo "[Error] s_file contains no singleton positions: $S_OUT" >&2
        echo "[Error] Current VCF/site universe is unsuitable for SDS singleton input at chr${CHR} for $POP." >&2
        exit 1
    fi
fi

if [ "$DO_T" = true ] && { [ "$FORCE" -eq 1 ] || [ ! -f "$T_OUT" ]; }; then
    t_p_pid=""
    t_q_pid=""
    bash "$SCRIPT_DIR/generate_t_file.sh" "$VCF_FILE" "$SAMPLE_LIST" "$MY_TMP/t_p.txt" "$P_REG" &
    t_p_pid=$!
    if [ "$Q_REG" != "NONE" ]; then
        bash "$SCRIPT_DIR/generate_t_file.sh" "$VCF_FILE" "$SAMPLE_LIST" "$MY_TMP/t_q.txt" "$Q_REG" &
        t_q_pid=$!
    fi
    wait "$t_p_pid"
    if [ -n "$t_q_pid" ]; then
        wait "$t_q_pid"
    fi

    if [ "$Q_REG" = "NONE" ]; then
        cp "$MY_TMP/t_p.txt" "$T_OUT"
    else
        cat "$MY_TMP/t_p.txt" "$MY_TMP/t_q.txt" > "$T_OUT"
    fi
    require_nonempty_file "$T_OUT" "t_file"
elif [ "$DO_T" = true ]; then
    echo "[Skip] t_file exists. Use --force to regenerate."
    require_nonempty_file "$T_OUT" "t_file"
fi

# 3. 生成观测性文件 (o_file)
if [ "$DO_O" = true ] && { [ "$FORCE" -eq 1 ] || [ ! -f "$O_OUT" ]; }; then
    bash "$SCRIPT_DIR/generate_o_file.sh" "$SAMPLE_LIST" "$OUTDIR" "chr${CHR}" "$MY_TMP/chr${CHR}_s_file.with_names.txt"
    require_nonempty_file "$O_OUT" "o_file"
elif [ "$DO_O" = true ]; then
    echo "[Skip] o_file exists. Use --force to regenerate."
    require_nonempty_file "$O_OUT" "o_file"
fi

if [ "$DO_B" = true ]; then
    require_nonempty_file "$B_FILE" "b_file"
fi

if [ "$EMIT_AUDIT_SIDECARS" -eq 1 ] && [ -f "$S_OUT" ] && [ -f "$T_OUT" ] && [ -f "$O_OUT" ]; then
    write_input_audit_sidecars "$CLEAN_SAMPLE_ORDER" "$S_OUT" "$T_OUT" "$O_OUT" "$OUTDIR" "$CHR"
fi

echo "✓ DONE: chr${CHR} for $POP"
