#!/bin/bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 /path/to/relaunch_manifest.tsv" >&2
    exit 1
fi

MANIFEST="$1"
[[ -f "$MANIFEST" ]] || { echo "[Error] Manifest not found: $MANIFEST" >&2; exit 1; }

value_of() {
    local key="$1"
    gawk -F '\t' -v k="$key" '$1 == k {print $2; exit}' "$MANIFEST"
}

show_job() {
    local label="$1"
    local job_id="$2"
    [[ -n "$job_id" ]] || return 0
    local status
    status="$(bjobs "$job_id" 2>/dev/null | awk 'NR==2{print $3}')"
    if [[ -z "$status" ]]; then
        status="$(bhist -n 1 "$job_id" 2>/dev/null | awk 'NR==2{print $4}')"
    fi
    printf '%s\t%s\t%s\n' "$label" "$job_id" "${status:-unknown}"
}

printf 'manifest\t%s\n' "$MANIFEST"
printf 'output_root\t%s\n' "$(value_of output_root)"
for pop in NCN SCN; do
    printf '%s_missing\t%s\n' "$pop" "$(value_of "${pop}_missing_chrs")"
    for chr in $(seq 1 22); do
        show_job "${pop}_chr${chr}_final_job" "$(value_of "${pop}_chr${chr}_final_job")"
    done
    show_job "${pop}_postprocess_job" "$(value_of "${pop}_postprocess_job")"
    show_job "${pop}_diag_job" "$(value_of "${pop}_diag_job")"
done
