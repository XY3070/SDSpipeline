#!/bin/bash
set -euo pipefail

# ========= Parameter reading =========
CHR="$1"
OUTDIR="$2"
PREFIX="${3:-chr${CHR}}"

# 输出文件
mkdir -p "$OUTDIR"
BFILE="$OUTDIR/${PREFIX}_b_file.txt"

echo "=========================================="
echo "Generate b_file (skip centromere)"
echo "=========================================="

# ========= GRCh38 centromere coordinates (UCSC) =========
# Format: CENTROMERE_START CENTROMERE_END
declare -A CENTROMERE_COORDS=(
    ["1"]="122026459 124932724"
    ["2"]="92188145 94090557"
    ["3"]="90772458 93655574"
    ["4"]="49712061 51743951"
    ["5"]="46485900 50059807"
    ["6"]="58553888 59829934"
    ["7"]="58169653 61528020"
    ["8"]="44033744 45877265"
    ["9"]="43389635 45518558"
    ["10"]="39686682 41593521"
    ["11"]="51078348 54425074"
    ["12"]="34769407 37185252"
    ["13"]="16000000 18051248"
    ["14"]="16000000 18173523"
    ["15"]="17083673 19725254"
    ["16"]="36311158 38265669"
    ["17"]="22813679 26616164"
    ["18"]="15460899 20861206"
    ["19"]="24498980 27190874"
    ["20"]="26436232 30038348"
    ["21"]="10864560 12915808"
    ["22"]="12954788 15054318"
)

echo "[Step 1/2] Get chr$CHR coordinates (reference: GRCh38/UCSC)..."

if [[ -z "${CENTROMERE_COORDS[$CHR]:-}" ]]; then
    echo "[Error] Chromosome $CHR does not have centromere coordinates defined"
    exit 1
fi

read CENT_START CENT_END <<< "${CENTROMERE_COORDS[$CHR]}"

# Define chromosome total length (approximate, from previous version)
# This is only for display purposes, the actual processing uses regions
declare -A CHR_LENGTHS=(
    ["1"]=248956422 ["2"]=242193529 ["3"]=198295559 ["4"]=190214555
    ["5"]=181538259 ["6"]=170805979 ["7"]=159345973 ["8"]=145138636
    ["9"]=138394717 ["10"]=133797422 ["11"]=135086622 ["12"]=133275309
    ["13"]=114364328 ["14"]=107043718 ["15"]=101991189 ["16"]=90338345
    ["17"]=83257441 ["18"]=80373285 ["19"]=58617616 ["20"]=64444167
    ["21"]=46709983 ["22"]=50818468
)
CHR_LEN=${CHR_LENGTHS[$CHR]}

# Define two continuous regions for parallel processing
# Region 1: 1 - centromere_end
# Region 2: centromere_end - chromosome_end
P_START=1
P_END=$CENT_END
Q_START=$((CENT_END + 1))
Q_END=$CHR_LEN

echo "  chr$CHR structure (for parallel trunking):"
echo "    Trunk 1: $P_START - $P_END"
echo "    Trunk 2: $Q_START - $Q_END"

echo "[Step 2/2] Generate b_file (two-line format for trunking)..."

# b_file format:
# Line 1: Trunk 1 region
# Line 2: Trunk 2 region
cat > "$BFILE" <<EOF
$P_START	$P_END
$Q_START	$Q_END
EOF

# Calculate statistics
CENTROMERE_SIZE=$((CENT_END - CENT_START))
P_ARM_SIZE=$((P_END - P_START + 1))
Q_ARM_SIZE=$((Q_END - Q_START + 1))
TOTAL_ANALYZED=$((P_ARM_SIZE + Q_ARM_SIZE))
CHROMOSOME_SIZE=$CHR_LEN

echo ""
echo "=========================================="
echo "Statistics"
echo "=========================================="
echo "Chromosome total length: $(echo "scale=1; $CHROMOSOME_SIZE / 1000000" | bc) Mb"
echo "  - P-arm: $(echo "scale=1; $P_ARM_SIZE / 1000000" | bc) Mb"
echo "  - centromere: $(echo "scale=1; $CENTROMERE_SIZE / 1000000" | bc) Mb (skip)"
echo "  - q-arm: $(echo "scale=1; $Q_ARM_SIZE / 1000000" | bc) Mb"
echo ""
echo "Analyzed region: $(echo "scale=1; $TOTAL_ANALYZED / 1000000" | bc) Mb"
echo "Skipped proportion: $(echo "scale=1; ($CHROMOSOME_SIZE - $TOTAL_ANALYZED) * 100 / $CHROMOSOME_SIZE" | bc)%"
echo ""
echo "✓ b_file generated: $BFILE"
echo ""
echo "File content:"
cat "$BFILE"