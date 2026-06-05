#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck source=common_env.sh
source "$SCRIPT_DIR/common_env.sh"
activate_sds_env
PYTHON_BIN="$SDS_PYTHON"

POP=""
CHR=""
OUT_ROOT="$SDS_SDS_OUTPUT_ROOT"
IN_ROOT="$SDS_SDS_INPUT_ROOT"
VCF_OVERRIDE=""
TEST_MODE=false
FORCE=0
KEEP_TEMP=0
G_FILE=""
INIT="0.00001"
S_FILE_NCOL="20000"
CHUNK_ROWS="0"
MAX_PARALLEL_CHUNKS="0"
REQUESTED_SLOTS="${LSB_DJOB_NUMPROC:-0}"

while [[ $# -gt 0 ]]; do
    case $1 in
        --pop) POP="$2"; shift 2 ;;
        --chr) CHR="$2"; shift 2 ;;
        --in-root) IN_ROOT="$2"; shift 2 ;;
        --test) TEST_MODE=true; shift ;;
        --force) FORCE=1; shift ;;
        --keep-temp) KEEP_TEMP=1; shift ;;
        --out-root) OUT_ROOT="$2"; shift 2 ;;
        --vcf) VCF_OVERRIDE="$2"; shift 2 ;;
        --g-file) G_FILE="$2"; shift 2 ;;
        --init) INIT="$2"; shift 2 ;;
        --s-file-ncol) S_FILE_NCOL="$2"; shift 2 ;;
        --chunk-rows) CHUNK_ROWS="$2"; shift 2 ;;
        --max-parallel-chunks) MAX_PARALLEL_CHUNKS="$2"; shift 2 ;;
        *)
            echo "Unknown parameter: $1" >&2
            exit 1
            ;;
    esac
done

if [[ -z "$POP" || -z "$CHR" ]]; then
    echo "Usage: $0 --pop POP --chr CHR [--test] [--force] [--keep-temp] [--out-root DIR] [--g-file PATH] [--init VALUE] [--s-file-ncol N] [--chunk-rows N] [--max-parallel-chunks N]" >&2
    exit 1
fi

if [[ -n "$VCF_OVERRIDE" ]]; then
    VCF_FILE="$VCF_OVERRIDE"
else
    VCF_FILE="$(find_population_vcf "$POP" "$CHR" || true)"
fi

if [[ -z "$G_FILE" ]]; then
    G_FILE="$(find_default_g_file "$BASE_DIR" "$POP" || true)"
fi

[[ -n "$G_FILE" ]] || { echo "[Error] g_file.txt not found. Provide --g-file PATH." >&2; exit 1; }
[[ -f "$G_FILE" ]] || { echo "[Error] g_file not found: $G_FILE" >&2; exit 1; }

INDIR="$IN_ROOT/$POP"
OUTDIR="$OUT_ROOT/$POP"
LOGDIR="$OUTDIR/logs"
TMPDIR="$OUTDIR/tmp/chr${CHR}_$(date +%s)_$RANDOM"
CACHE_DIR="$OUTDIR/cache"
CHUNK_LOG_DIR="$LOGDIR/chr${CHR}_chunks"

mkdir -p "$OUTDIR" "$LOGDIR" "$TMPDIR" "$CACHE_DIR" "$CHUNK_LOG_DIR"

if [[ "$KEEP_TEMP" -eq 0 ]]; then
    trap '/bin/rm -rf "$TMPDIR"' EXIT
fi

S_FILE="$INDIR/chr${CHR}_s_file.txt"
T_FILE="$INDIR/chr${CHR}_t_file.txt"
O_FILE="$INDIR/chr${CHR}_o_file.txt"
B_FILE="$INDIR/chr${CHR}_b_file.txt"

[[ -f "$S_FILE" ]] || { echo "[Error] Missing s_file: $S_FILE" >&2; exit 1; }
[[ -f "$T_FILE" ]] || { echo "[Error] Missing t_file: $T_FILE" >&2; exit 1; }
[[ -f "$O_FILE" ]] || { echo "[Error] Missing o_file: $O_FILE" >&2; exit 1; }
[[ -f "$B_FILE" ]] || { echo "[Error] Missing b_file: $B_FILE" >&2; exit 1; }
[[ -f "$VCF_FILE" ]] || { echo "[Error] VCF not found: $VCF_FILE" >&2; exit 1; }

read -r P_START P_END < <(sed -n '1p' "$B_FILE")
read -r Q_START Q_END < <(sed -n '2p' "$B_FILE")

P_T="$TMPDIR/chr${CHR}_p.t.feather"
Q_T="$TMPDIR/chr${CHR}_q.t.feather"
P_OUT="$OUTDIR/chr${CHR}_p.sds.tsv"
Q_OUT="$OUTDIR/chr${CHR}_q.sds.tsv"
FINAL_OUT="$OUTDIR/chr${CHR}.sds.tsv"
P_PARQUET="$OUTDIR/chr${CHR}_p.sds.parquet"
Q_PARQUET="$OUTDIR/chr${CHR}_q.sds.parquet"
FINAL_PARQUET="$OUTDIR/chr${CHR}.sds.parquet"
P_SUMMARY="$LOGDIR/chr${CHR}_p.compute.csv"
Q_SUMMARY="$LOGDIR/chr${CHR}_q.compute.csv"
ARCHIVE_SUMMARY="$LOGDIR/chr${CHR}.archive.csv"

launch_with_limit() {
    local limit="$1"
    shift
    local oldest_pid

    "$@" &
    if [[ "$limit" -gt 0 ]]; then
        while [[ "$(jobs -pr | wc -l)" -ge "$limit" ]]; do
            oldest_pid="$(jobs -pr | head -n 1)"
            [[ -n "$oldest_pid" ]] || break
            wait "$oldest_pid"
        done
    fi
}

merge_chunk_outputs() {
    local out_file="$1"
    shift
    local header_written=0

    : > "$out_file"
    for part_file in "$@"; do
        [[ -s "$part_file" ]] || continue
        if [[ "$header_written" -eq 0 ]]; then
            cat "$part_file" >> "$out_file"
            header_written=1
        else
            tail -n +2 "$part_file" >> "$out_file"
        fi
    done

    [[ "$header_written" -eq 1 ]] || {
        echo "[Error] No chunk outputs were merged into $out_file" >&2
        exit 1
    }
}

compute_arm_chunks() {
    local arm_label="$1"
    local input_feather="$2"
    local arm_out="$3"
    local arm_parquet="$4"
    local arm_summary="$5"

    local chunk_dir="$TMPDIR/chunks_${arm_label}"
    local split_summary="$CHUNK_LOG_DIR/chr${CHR}_${arm_label}.split.csv"
    mkdir -p "$chunk_dir"

    if [[ "$CHUNK_ROWS" -gt 0 ]]; then
        "$PYTHON_BIN" "$SCRIPT_DIR/split_feather_by_rows.py" \
            "$input_feather" "$chunk_dir" \
            --chunk-rows "$CHUNK_ROWS" \
            --prefix "chr${CHR}_${arm_label}" \
            --summary-csv "$split_summary"
    else
        cp "$input_feather" "$chunk_dir/chr${CHR}_${arm_label}_0000.feather"
    fi

    mapfile -t chunk_files < <(find "$chunk_dir" -maxdepth 1 -type f -name '*.feather' | sort)
    [[ "${#chunk_files[@]}" -gt 0 ]] || {
        echo "[Error] No compute chunks found for arm ${arm_label}" >&2
        exit 1
    }

    local chunk_outputs=()
    local idx=0
    for chunk_file in "${chunk_files[@]}"; do
        local part_out="$chunk_dir/${arm_label}_chunk_${idx}.sds.tsv"
        local part_parquet="$chunk_dir/${arm_label}_chunk_${idx}.parquet"
        local part_summary="$CHUNK_LOG_DIR/chr${CHR}_${arm_label}_chunk_${idx}.csv"
        chunk_outputs+=("$part_out")
        launch_with_limit "$MAX_PARALLEL_CHUNKS" \
            bash "$SCRIPT_DIR/run_sds_compute_chunk.sh" \
            "$chunk_file" "$part_out" "$S_FILE" "$O_FILE" "$B_FILE" "$G_FILE" "$INIT" "$S_FILE_NCOL" "$part_parquet" "$part_summary" "$CACHE_DIR"
        idx=$((idx + 1))
    done

    wait
    merge_chunk_outputs "$arm_out" "${chunk_outputs[@]}"
    "$PYTHON_BIN" "$SCRIPT_DIR/archive_sds_output.py" \
        "$arm_out" \
        "$arm_parquet" \
        --summary-csv "$arm_summary"
}

if [[ "$FORCE" -eq 0 && -f "$FINAL_OUT" ]]; then
    echo "[Skip] Final SDS output already exists: $FINAL_OUT"
    exit 0
fi

echo ">>> SDS runtime: $(sds_env_label)"
echo ">>> Compute population: $POP | Chromosome: $CHR"
echo ">>> Input dir: $INDIR"
echo ">>> g_file: $G_FILE"

if [[ "$TEST_MODE" = true ]]; then
    TEST_START="$(find_first_variant_pos "$VCF_FILE" "chr${CHR}")"
    if [[ -z "$TEST_START" ]]; then
        echo "[Error] No variants found on chr${CHR} in VCF for test mode" >&2
        exit 1
    fi
    TEST_END=$((TEST_START + 2000000 - 1))
    if (( TEST_END > Q_END )); then
        TEST_END="$Q_END"
    fi
    echo "[TEST] Mode enabled: computing 2Mb from first variant"
    echo "[TEST] Region: chr${CHR}:${TEST_START}-${TEST_END}"
    "$PYTHON_BIN" "$SCRIPT_DIR/filter_t_to_feather.py" \
        "$T_FILE" "$P_T" \
        --start "$TEST_START" \
        --end "$TEST_END" \
        --summary-csv "$P_SUMMARY"
    /bin/rm -f "$Q_T" "$Q_SUMMARY"
else
    "$PYTHON_BIN" "$SCRIPT_DIR/filter_t_to_feather.py" \
        "$T_FILE" "$P_T" \
        --start "$P_START" \
        --end "$P_END" \
        --summary-csv "$P_SUMMARY"
    "$PYTHON_BIN" "$SCRIPT_DIR/filter_t_to_feather.py" \
        "$T_FILE" "$Q_T" \
        --start "$Q_START" \
        --end "$Q_END" \
        --summary-csv "$Q_SUMMARY"
fi

if [[ ! -f "$P_T" && ! -f "$Q_T" ]]; then
    echo "[Error] No compute chunks were generated from $T_FILE" >&2
    exit 1
fi

ARM_CONCURRENCY=0
if [[ -f "$P_T" ]]; then
    ARM_CONCURRENCY=$((ARM_CONCURRENCY + 1))
fi
if [[ "$TEST_MODE" = false && -f "$Q_T" ]]; then
    ARM_CONCURRENCY=$((ARM_CONCURRENCY + 1))
fi
if [[ "$ARM_CONCURRENCY" -lt 1 ]]; then
    ARM_CONCURRENCY=1
fi

if [[ "$MAX_PARALLEL_CHUNKS" -eq 0 && "$CHUNK_ROWS" -gt 0 && "$REQUESTED_SLOTS" -gt 0 ]]; then
    MAX_PARALLEL_CHUNKS=$((REQUESTED_SLOTS / ARM_CONCURRENCY))
    if [[ "$MAX_PARALLEL_CHUNKS" -lt 1 ]]; then
        MAX_PARALLEL_CHUNKS=1
    fi
fi

THREAD_BUDGET="${SDS_NUMBA_THREADS:-}"
if [[ -z "$THREAD_BUDGET" ]]; then
    if [[ "$REQUESTED_SLOTS" -gt 0 ]]; then
        if [[ "$MAX_PARALLEL_CHUNKS" -gt 1 ]]; then
            THREAD_BUDGET=1
        else
            THREAD_BUDGET=$((REQUESTED_SLOTS / ARM_CONCURRENCY))
            if [[ "$THREAD_BUDGET" -lt 1 ]]; then
                THREAD_BUDGET=1
            fi
        fi
    elif [[ "$MAX_PARALLEL_CHUNKS" -gt 1 ]]; then
        THREAD_BUDGET=1
    fi
fi

if [[ -n "$THREAD_BUDGET" ]]; then
    export SDS_NUMBA_THREADS="$THREAD_BUDGET"
fi

echo ">>> Requested LSF slots: ${REQUESTED_SLOTS}"
echo ">>> Arm concurrency: ${ARM_CONCURRENCY}"
echo ">>> Chunk rows: ${CHUNK_ROWS}"
echo ">>> Max parallel chunks per arm: ${MAX_PARALLEL_CHUNKS}"
echo ">>> Numba threads per compute process: ${SDS_NUMBA_THREADS:-auto}"

if [[ -f "$P_T" ]]; then
    compute_arm_chunks "p" "$P_T" "$P_OUT" "$P_PARQUET" "$P_SUMMARY" &
fi

if [[ "$TEST_MODE" = false && -f "$Q_T" ]]; then
    compute_arm_chunks "q" "$Q_T" "$Q_OUT" "$Q_PARQUET" "$Q_SUMMARY" &
fi

wait

{
    header_written=0
    for part_file in "$P_OUT" "$Q_OUT"; do
        [[ -s "$part_file" ]] || continue
        if [[ "$header_written" -eq 0 ]]; then
            cat "$part_file"
            header_written=1
        else
            tail -n +2 "$part_file"
        fi
    done
} > "$FINAL_OUT"

[[ -s "$FINAL_OUT" ]] || { echo "[Error] Final SDS output was not created: $FINAL_OUT" >&2; exit 1; }

"$PYTHON_BIN" "$SCRIPT_DIR/archive_sds_output.py" \
    "$FINAL_OUT" \
    "$FINAL_PARQUET" \
    --summary-csv "$ARCHIVE_SUMMARY"

echo "✓ DONE: chr${CHR} compute for $POP"
