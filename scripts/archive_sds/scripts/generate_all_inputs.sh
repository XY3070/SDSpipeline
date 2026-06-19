#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

NPROC_TOTAL="$(nproc)"
PARALLEL_JOBS=$((NPROC_TOTAL / 2))
if (( PARALLEL_JOBS < 1 )); then
    PARALLEL_JOBS=1
fi

cd "$BASE_DIR"

for pop in NCN SCN MERGED3971; do
    for chr in $(seq 1 22); do
        printf '%s %s\n' "$pop" "$chr"
    done
done | xargs -n 2 -P "$PARALLEL_JOBS" bash -lc '
    set -euo pipefail
    cd "'"$BASE_DIR"'" || exit 1
    printf "[%s] START %s chr%s\n" "$(date +%F\ %T)" "$1" "$2"
    bash scripts/run_sds_input.sh --pop "$1" --chr "$2" --force
    printf "[%s] DONE %s chr%s\n" "$(date +%F\ %T)" "$1" "$2"
' _
