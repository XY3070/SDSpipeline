#!/bin/bash
# SDS input file validation script

set -euo pipefail

OUTDIR="${1:?Please provide output directory}"
CHROMOSOMES="${2:-$(seq 1 22)}"

echo "=========================================="
echo "SDS input file validation"
echo "=========================================="
echo "Output directory: $OUTDIR"
echo "Chromosomes: $CHROMOSOMES"
echo ""

TOTAL_ERRORS=0

for CHR in $CHROMOSOMES; do
    PREFIX="chr${CHR}"
    echo "----------------------------------------"
    echo "Validate $PREFIX"
    echo "----------------------------------------"
    
    ERRORS=0
    
    # Validate file existence
    for suffix in s_file.txt t_file.txt o_file.txt b_file.txt; do
        FILE="$OUTDIR/${PREFIX}_${suffix}"
        if [[ ! -f "$FILE" ]]; then
            echo "  ✗ Missing: $suffix"
            ((ERRORS++))
        else
            echo "  ✓ Exists: $suffix"
        fi
    done
    
    # If all files exist, perform detailed validation
    if [[ $ERRORS -eq 0 ]]; then
        echo ""
        echo "Detailed validation:"
        
        # Validate s_file
        S_FILE="$OUTDIR/${PREFIX}_s_file.txt"
        N_SAMPLES_S=$(wc -l < "$S_FILE")
        echo "  s_file: $N_SAMPLES_S 个样本"
        
        # Check for empty lines or format errors
        if grep -q '^[[:space:]]*$' "$S_FILE"; then
            echo "    ✗ Warning: s_file contains empty lines"
            ((ERRORS++))
        fi
        
        # Validate t_file
        T_FILE="$OUTDIR/${PREFIX}_t_file.txt"
        N_SNPS=$(wc -l < "$T_FILE")
        N_COLS=$(head -1 "$T_FILE" | awk '{print NF}')
        N_SAMPLES_T=$((N_COLS - 4))
        
        echo "  t_file: $N_SNPS SNPs, $N_SAMPLES_T samples"
        
        if [[ $N_SAMPLES_S -ne $N_SAMPLES_T ]]; then
            echo "    ✗ Error: s_file and t_file have different number of samples!"
            ((ERRORS++))
        else
            echo "    ✓ Samples match"
        fi
        
        # Check genotype encoding
        INVALID_GT=$(awk '{for(i=5;i<=NF;i++) if($i!="0"&&$i!="1"&&$i!="2"&&$i!="NA") print $i}' "$T_FILE" | head -5)
        if [[ -n "$INVALID_GT" ]]; then
            echo "    ✗ Warning: Invalid genotype encoding found: $INVALID_GT"
        else
            echo "    ✓ Genotype encoding is correct"
        fi
        
        # Calculate missing rate
        TOTAL_GT=$((N_SNPS * N_SAMPLES_T))
        MISSING_GT=$(awk '{for(i=5;i<=NF;i++) if($i=="NA") n++} END{print n}' "$T_FILE")
        MISSING_RATE=$(echo "scale=2; $MISSING_GT * 100 / $TOTAL_GT" | bc)
        echo "    Missing rate: ${MISSING_RATE}%"
        
        # Validate o_file
        O_FILE="$OUTDIR/${PREFIX}_o_file.txt"
        N_VALUES=$(awk '{print NF}' "$O_FILE")
        
        echo "  o_file: $N_VALUES values"
        
        if [[ $N_VALUES -ne $N_SAMPLES_S ]]; then
            echo "    ✗ Error: o_file values count does not match sample count!"
            ((ERRORS++))
        else
            echo "    ✓ Values count matches samples"
        fi
        
        # Calculate o_file range
        O_STATS=$(awk '{
            for(i=1;i<=NF;i++) {
                sum+=$i
                if(NR==1&&i==1) min=max=$i
                if($i<min) min=$i
                if($i>max) max=$i
            }
        }
        END {
            printf "min=%d max=%d avg=%.0f", min, max, sum/NF
        }' "$O_FILE")
        echo "    Callable bases: $O_STATS"
        
        # Validate b_file
        B_FILE="$OUTDIR/${PREFIX}_b_file.txt"
        N_REGIONS=$(wc -l < "$B_FILE")
        
        echo "  b_file: $N_REGIONS regions"
        
        if [[ $N_REGIONS -ne 2 ]]; then
            echo "    ✗ Error: b_file should have 2 lines (p and q arms)"
            ((ERRORS++))
        else
            echo "    ✓ Region count matches arms"
            
            # Show region bounds
            echo "    Region bounds:"
            awk '{printf "      %d - %d (%.1f Mb)\n", $1, $2, ($2-$1)/1000000}' "$B_FILE"
        fi
        
        # Validate SNP positions
        OUTSIDE_SNP=$(awk '
        NR==FNR {
            start[NR]=$1; end[NR]=$2; n_regions=NR
            next
        }
        {
            pos=$4
            in_region=0
            for(i=1; i<=n_regions; i++) {
                if(pos>=start[i] && pos<=end[i]) {
                    in_region=1
                    break
                }
            }
            if(!in_region) outside++
        }
        END {
            print outside
        }
        ' "$B_FILE" "$T_FILE")
        
        if [[ $OUTSIDE_SNP -gt 0 ]]; then
            echo "    ✗ Error: $OUTSIDE_SNP SNPs outside defined regions"
            ((ERRORS++))
        else
            echo "    ✓ All SNPs inside analysis regions"
        fi
        
    fi
    
    if [[ $ERRORS -gt 0 ]]; then
        echo ""
        echo "  Status: ✗ Found $ERRORS errors"
        ((TOTAL_ERRORS += ERRORS))
    else
        echo ""
        echo "  Status: ✓ Validation passed"
    fi
    
    echo ""
done

echo "=========================================="
echo "Validation Summary"
echo "=========================================="

if [[ $TOTAL_ERRORS -gt 0 ]]; then
    echo "✗ Total errors found: $TOTAL_ERRORS"
    echo "Please check the above error messages and regenerate the files"
    exit 1
else
    echo "✓ All files validated successfully!"
    echo ""
    echo "Files are ready for SDS calculation"
    echo ""
    echo "Next steps:"
    echo "  1. Confirm all files are in correct format"
    echo "  2. Run the SDS main program"
    echo "  3. Suggested command:"
    echo "     sds --s_file chr1_s_file.txt \\"
    echo "         --t_file chr1_t_file.txt \\"
    echo "         --o_file chr1_o_file.txt \\"
    echo "         --b_file chr1_b_file.txt \\"
    echo "         --output chr1_sds_results.txt"
fi
