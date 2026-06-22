#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$SCRIPT_DIR/common_env.sh"
activate_sds_env

POP="NCN"
CHR="1"
G_FILE="$BASE_DIR/g_file.txt"
SPARSE_T="$SDS_TMP_ROOT/chr1_tfile_sparse_like_legacy_20260505/chr1_t_file.txt"
INPUT_DIR="$SDS_SDS_INPUT_ROOT/$POP"
OUT_DIR="$SDS_TMP_ROOT/chr1_sparse_shadow_q_20260505"
LOG_DIR="$OUT_DIR/logs"
CACHE_DIR="$OUT_DIR/cache"
INIT="0.00001"
S_FILE_NCOL="20000"
PYTHON_BIN="$SDS_ENV_PREFIX/bin/python"

mkdir -p "$OUT_DIR" "$LOG_DIR" "$CACHE_DIR"

S_FILE="$INPUT_DIR/chr${CHR}_s_file.txt"
O_FILE="$INPUT_DIR/chr${CHR}_o_file.txt"
B_FILE="$INPUT_DIR/chr${CHR}_b_file.txt"

[[ -s "$SPARSE_T" ]] || { echo "[Error] Missing sparse t_file: $SPARSE_T" >&2; exit 1; }
[[ -s "$S_FILE" ]] || { echo "[Error] Missing s_file: $S_FILE" >&2; exit 1; }
[[ -s "$O_FILE" ]] || { echo "[Error] Missing o_file: $O_FILE" >&2; exit 1; }
[[ -s "$B_FILE" ]] || { echo "[Error] Missing b_file: $B_FILE" >&2; exit 1; }
[[ -s "$G_FILE" ]] || { echo "[Error] Missing g_file: $G_FILE" >&2; exit 1; }

read -r _P_START _P_END < <(sed -n '1p' "$B_FILE")
read -r Q_START Q_END < <(sed -n '2p' "$B_FILE")

Q_T="$OUT_DIR/chr1_q.t.feather"
Q_OUT="$OUT_DIR/chr1_q.sds.tsv"
Q_PARQUET="$OUT_DIR/chr1_q.sds.parquet"
Q_SUMMARY="$LOG_DIR/chr1_q.compute.csv"
ARCHIVE_SUMMARY="$LOG_DIR/chr1_q.archive.csv"

"$PYTHON_BIN" "$SCRIPT_DIR/filter_t_to_feather.py" \
    "$SPARSE_T" "$Q_T" \
    --start "$Q_START" \
    --end "$Q_END" \
    --summary-csv "$Q_SUMMARY"

bash "$SCRIPT_DIR/run_sds_compute_chunk.sh" \
    "$Q_T" "$Q_OUT" "$S_FILE" "$O_FILE" "$B_FILE" "$G_FILE" "$INIT" "$S_FILE_NCOL" "$Q_PARQUET" "$Q_SUMMARY" "$CACHE_DIR"

"$PYTHON_BIN" "$SCRIPT_DIR/archive_sds_output.py" \
    "$Q_OUT" \
    "$Q_PARQUET" \
    --summary-csv "$ARCHIVE_SUMMARY"

echo "shadow_q_output=$Q_OUT"
