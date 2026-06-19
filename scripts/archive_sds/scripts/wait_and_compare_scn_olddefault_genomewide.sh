#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

ROOT="$BASE_DIR/data/processed/sds_output_scn_olddefault_genomewide_20260527/SCN"
OUT_DIR="$BASE_DIR/data/processed/scn_olddefault_genomewide_compare_20260528"
TIMEOUT_SEC="${TIMEOUT_SEC:-43200}"
SLEEP_SEC="${SLEEP_SEC:-60}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --root) ROOT="$2"; shift 2 ;;
        --out-dir) OUT_DIR="$2"; shift 2 ;;
        *) echo "Unknown parameter: $1" >&2; exit 1 ;;
    esac
done

FREQ="$ROOT/SCN.frequency_bins.tsv"
DIAG="$ROOT/diagnostics/SCN.gravel_chb_ne100k_submitter_array.chrom_tail_counts.tsv"

start_ts="$(date +%s)"
while true; do
    if [[ -f "$FREQ" && -f "$DIAG" ]]; then
        break
    fi
    now_ts="$(date +%s)"
    if (( now_ts - start_ts > TIMEOUT_SEC )); then
        echo "[Error] timed out waiting for SCN olddefault genomewide outputs" >&2
        echo "missing_freq=$([[ -f "$FREQ" ]] && echo 0 || echo 1)" >&2
        echo "missing_diag=$([[ -f "$DIAG" ]] && echo 0 || echo 1)" >&2
        exit 1
    fi
    sleep "$SLEEP_SEC"
done

bash "$SCRIPT_DIR/launch_scn_olddefault_genomewide_compare.sh" \
    --olddefault-root "$ROOT" \
    --out-dir "$OUT_DIR"
