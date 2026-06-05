#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$SCRIPT_DIR/common_env.sh"
activate_sds_env
PYTHON_BIN="$SDS_PYTHON"

POP=""
CHR=""
OUT_ROOT="$SDS_SDS_OUTPUT_ROOT"
IN_ROOT="$SDS_SDS_INPUT_ROOT"
QUEUE="normal"
CHUNK_ROWS=""
ARRAY_PARALLEL=""
CHUNK_JOB_SLOTS=""
G_FILE=""
INIT="0.00001"
S_FILE_NCOL="20000"
JOB_GROUP="chunked"
WAIT_FOR_FINAL=0
SKIP_BOUNDARY_MISSING_FRACTION=""
BOUNDARY_MISSING_MODE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pop) POP="$2"; shift 2 ;;
        --chr) CHR="$2"; shift 2 ;;
        --out-root) OUT_ROOT="$2"; shift 2 ;;
        --in-root) IN_ROOT="$2"; shift 2 ;;
        --queue) QUEUE="$2"; shift 2 ;;
        --chunk-rows) CHUNK_ROWS="$2"; shift 2 ;;
        --array-parallel) ARRAY_PARALLEL="$2"; shift 2 ;;
        --chunk-job-slots) CHUNK_JOB_SLOTS="$2"; shift 2 ;;
        --g-file) G_FILE="$2"; shift 2 ;;
        --init) INIT="$2"; shift 2 ;;
        --s-file-ncol) S_FILE_NCOL="$2"; shift 2 ;;
        --skip-boundary-missing-fraction) SKIP_BOUNDARY_MISSING_FRACTION="$2"; shift 2 ;;
        --boundary-missing-mode) BOUNDARY_MISSING_MODE="$2"; shift 2 ;;
        --job-group) JOB_GROUP="$2"; shift 2 ;;
        --wait) WAIT_FOR_FINAL=1; shift ;;
        *) echo "Unknown parameter: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "$POP" || -z "$CHR" ]]; then
    echo "Usage: $0 --pop POP --chr CHR [--out-root DIR] [--queue QUEUE] [--chunk-rows N] [--array-parallel N] [--chunk-job-slots N] [--g-file PATH] [--wait]" >&2
    exit 1
fi

QUEUE="$(sds_effective_queue "$QUEUE")"
if [[ -z "$CHUNK_ROWS" ]]; then
    CHUNK_ROWS="$(sds_queue_default_chunk_rows "$QUEUE")"
fi
if [[ -z "$ARRAY_PARALLEL" ]]; then
    ARRAY_PARALLEL="$(sds_queue_default_array_parallel "$QUEUE")"
fi
if [[ -z "$CHUNK_JOB_SLOTS" ]]; then
    CHUNK_JOB_SLOTS="$(sds_queue_default_chunk_job_slots "$QUEUE")"
fi
CHUNK_NUMBA_THREADS="$(sds_queue_default_numba_threads "$QUEUE" "$CHUNK_JOB_SLOTS")"

if [[ -z "$G_FILE" ]]; then
    G_FILE="$(find_default_g_file "$BASE_DIR" "$POP" || true)"
fi
[[ -n "$G_FILE" && -f "$G_FILE" ]] || { echo "[Error] g_file not found: $G_FILE" >&2; exit 1; }

INDIR="$IN_ROOT/$POP"
OUTDIR="$OUT_ROOT/$POP"
LOGDIR="$OUTDIR/logs"
WORKDIR="$OUTDIR/chunk_submit/chr${CHR}_${JOB_GROUP}_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$LOGDIR" "$WORKDIR"

S_FILE="$INDIR/chr${CHR}_s_file.txt"
T_FILE="$INDIR/chr${CHR}_t_file.txt"
O_FILE="$INDIR/chr${CHR}_o_file.txt"
B_FILE="$INDIR/chr${CHR}_b_file.txt"
[[ -f "$S_FILE" && -f "$T_FILE" && -f "$O_FILE" && -f "$B_FILE" ]] || {
    echo "[Error] Missing input files under $INDIR for chr${CHR}" >&2
    exit 1
}

S_NON_NA_ROWS="$(gawk 'BEGIN{n=0} $0 != "NA" {n++} END{print n}' "$S_FILE")"
if [[ "${S_NON_NA_ROWS:-0}" -eq 0 ]]; then
    echo "[Error] s_file contains no singleton positions: $S_FILE" >&2
    echo "[Error] Refusing to submit chunked SDS compute for chr${CHR} ${POP}." >&2
    exit 1
fi

read -r P_START P_END < <(sed -n '1p' "$B_FILE")
read -r Q_START Q_END < <(sed -n '2p' "$B_FILE")

P_T="$WORKDIR/chr${CHR}_p.t.feather"
Q_T="$WORKDIR/chr${CHR}_q.t.feather"
"$PYTHON_BIN" "$SCRIPT_DIR/filter_t_to_feather.py" "$T_FILE" "$P_T" --start "$P_START" --end "$P_END" --summary-csv "$WORKDIR/chr${CHR}_p.filter.csv"
"$PYTHON_BIN" "$SCRIPT_DIR/filter_t_to_feather.py" "$T_FILE" "$Q_T" --start "$Q_START" --end "$Q_END" --summary-csv "$WORKDIR/chr${CHR}_q.filter.csv"

create_arm_manifest() {
    local arm="$1"
    local arm_feather="$2"
    local chunk_dir="$WORKDIR/${arm}_chunks"
    local manifest="$WORKDIR/chr${CHR}_${arm}.manifest.tsv"

    mkdir -p "$chunk_dir"
    "$PYTHON_BIN" "$SCRIPT_DIR/split_feather_by_rows.py" \
        "$arm_feather" "$chunk_dir" \
        --chunk-rows "$CHUNK_ROWS" \
        --prefix "chr${CHR}_${arm}" \
        --summary-csv "$WORKDIR/chr${CHR}_${arm}.split.csv"

    : > "$manifest"
    local idx=0
    local chunk_file
    while IFS= read -r chunk_file; do
        local out_tsv="$chunk_dir/${arm}_chunk_${idx}.sds.tsv"
        local out_parquet="$chunk_dir/${arm}_chunk_${idx}.parquet"
        local out_summary="$chunk_dir/${arm}_chunk_${idx}.csv"
        printf '%s\t%s\t%s\t%s\n' "$chunk_file" "$out_tsv" "$out_parquet" "$out_summary" >> "$manifest"
        idx=$((idx + 1))
    done < <(find "$chunk_dir" -maxdepth 1 -type f -name '*.feather' | sort)

    printf '%s\n' "$manifest"
}

submit_arm_array() {
    local arm="$1"
    local manifest="$2"
    local chunk_count
    local job_name="SDS_${POP}_chr${CHR}_${arm}_${JOB_GROUP}"
    chunk_count="$(wc -l < "$manifest")"
    if [[ "$chunk_count" -le 0 ]]; then
        printf '%s\n' "SKIP"
        return 0
    fi

    local submit_out
    submit_out="$(bsub -q "$QUEUE" -n "$CHUNK_JOB_SLOTS" -R "span[hosts=1]" -J "${job_name}[1-${chunk_count}]%${ARRAY_PARALLEL}" -o "$LOGDIR/${job_name}_%I.out" -e "$LOGDIR/${job_name}_%I.err" env SDS_NUMBA_THREADS="$CHUNK_NUMBA_THREADS" SDS_SKIP_BOUNDARY_MISSING_FRACTION="$SKIP_BOUNDARY_MISSING_FRACTION" SDS_BOUNDARY_MISSING_MODE="$BOUNDARY_MISSING_MODE" bash "$SCRIPT_DIR/run_sds_compute_manifest_entry.sh" "$manifest" "$S_FILE" "$O_FILE" "$B_FILE" "$G_FILE" "$INIT" "$S_FILE_NCOL" "$OUTDIR/cache" < /dev/null)"
    local job_id
    job_id="$(sed -n 's/.*<\([0-9]\+\)>.*/\1/p' <<< "$submit_out")"
    [[ -n "$job_id" ]] || { echo "[Error] Failed to parse bsub output: $submit_out" >&2; exit 1; }
    printf '%s\n' "$job_id"
}

submit_arm_merge() {
    local arm="$1"
    local array_job_id="$2"
    local manifest="$3"
    if [[ -z "$array_job_id" || "$array_job_id" == "SKIP" ]]; then
        printf '%s\n' "SKIP"
        return 0
    fi
    local arm_tsv="$OUTDIR/chr${CHR}_${arm}.sds.tsv"
    local arm_parquet="$OUTDIR/chr${CHR}_${arm}.sds.parquet"
    local arm_summary="$LOGDIR/chr${CHR}_${arm}.compute.csv"
    local job_name="SDS_${POP}_chr${CHR}_${arm}_merge_${JOB_GROUP}"
    local submit_out
    submit_out="$(bsub -q "$QUEUE" -w "done(${array_job_id})" -n 1 -J "$job_name" -o "$LOGDIR/${job_name}.out" -e "$LOGDIR/${job_name}.err" /bin/bash -lc "bash \"$SCRIPT_DIR/merge_sds_chunk_outputs.sh\" \"$manifest\" \"$arm_tsv\" \"$arm_parquet\" \"$arm_summary\"" < /dev/null)"
    local job_id
    job_id="$(sed -n 's/.*<\([0-9]\+\)>.*/\1/p' <<< "$submit_out")"
    [[ -n "$job_id" ]] || { echo "[Error] Failed to parse merge bsub output: $submit_out" >&2; exit 1; }
    printf '%s\n' "$job_id"
}

P_MANIFEST="$(create_arm_manifest p "$P_T")"
Q_MANIFEST="$(create_arm_manifest q "$Q_T")"
P_ARRAY_JOB="$(submit_arm_array p "$P_MANIFEST")"
Q_ARRAY_JOB="$(submit_arm_array q "$Q_MANIFEST")"
P_MERGE_JOB="$(submit_arm_merge p "$P_ARRAY_JOB" "$P_MANIFEST")"
Q_MERGE_JOB="$(submit_arm_merge q "$Q_ARRAY_JOB" "$Q_MANIFEST")"

FINAL_DEPS=()
if [[ "$P_MERGE_JOB" != "SKIP" ]]; then
    FINAL_DEPS+=("done(${P_MERGE_JOB})")
fi
if [[ "$Q_MERGE_JOB" != "SKIP" ]]; then
    FINAL_DEPS+=("done(${Q_MERGE_JOB})")
fi
if [[ "${#FINAL_DEPS[@]}" -eq 0 ]]; then
    echo "[Error] No non-empty arms were available for chr${CHR}" >&2
    exit 1
fi
FINAL_DEP_EXPR="$(printf '%s && ' "${FINAL_DEPS[@]}")"
FINAL_DEP_EXPR="${FINAL_DEP_EXPR% && }"

FINAL_TSV="$OUTDIR/chr${CHR}.sds.tsv"
FINAL_PARQUET="$OUTDIR/chr${CHR}.sds.parquet"
FINAL_SUMMARY="$LOGDIR/chr${CHR}.archive.csv"
FINAL_JOB_NAME="SDS_${POP}_chr${CHR}_final_${JOB_GROUP}"
FINAL_CMD="bash \"$SCRIPT_DIR/merge_sds_final_outputs.sh\" \"$FINAL_TSV\" \"$FINAL_PARQUET\" \"$FINAL_SUMMARY\" \"$OUTDIR/chr${CHR}_p.sds.tsv\" \"$OUTDIR/chr${CHR}_q.sds.tsv\""
FINAL_SUBMIT_OUT="$(bsub -q "$QUEUE" -w "$FINAL_DEP_EXPR" -n 1 -J "$FINAL_JOB_NAME" -o "$LOGDIR/${FINAL_JOB_NAME}.out" -e "$LOGDIR/${FINAL_JOB_NAME}.err" /bin/bash -lc "$FINAL_CMD" < /dev/null)"
FINAL_JOB_ID="$(sed -n 's/.*<\([0-9]\+\)>.*/\1/p' <<< "$FINAL_SUBMIT_OUT")"
[[ -n "$FINAL_JOB_ID" ]] || { echo "[Error] Failed to parse final bsub output: $FINAL_SUBMIT_OUT" >&2; exit 1; }

if [[ "$WAIT_FOR_FINAL" -eq 1 ]]; then
    bwait -w "ended(${FINAL_JOB_ID})"
fi

cat <<EOF
WORKDIR	$WORKDIR
QUEUE	$QUEUE
CHUNK_ROWS	$CHUNK_ROWS
ARRAY_PARALLEL	$ARRAY_PARALLEL
CHUNK_JOB_SLOTS	$CHUNK_JOB_SLOTS
NUMBA_THREADS	$CHUNK_NUMBA_THREADS
P_ARRAY_JOB	$P_ARRAY_JOB
Q_ARRAY_JOB	$Q_ARRAY_JOB
P_MERGE_JOB	$P_MERGE_JOB
Q_MERGE_JOB	$Q_MERGE_JOB
FINAL_JOB	$FINAL_JOB_ID
FINAL_TSV	$FINAL_TSV
FINAL_PARQUET	$FINAL_PARQUET
SKIP_BOUNDARY_MISSING_FRACTION	$SKIP_BOUNDARY_MISSING_FRACTION
BOUNDARY_MISSING_MODE	$BOUNDARY_MISSING_MODE
EOF
