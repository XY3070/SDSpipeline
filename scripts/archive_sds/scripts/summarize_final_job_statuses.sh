#!/bin/bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 STATE_DIR..." >&2
    exit 1
fi

declare -A seen
ids=()
for state_dir in "$@"; do
    result_dir="$state_dir/results"
    [[ -d "$result_dir" ]] || continue
    while IFS=$'\t' read -r chr final_job; do
        [[ "$chr" == "chr" ]] && continue
        [[ -n "$final_job" ]] || continue
        if [[ -z "${seen[$final_job]:-}" ]]; then
            seen["$final_job"]=1
            ids+=("$final_job")
        fi
    done < <(cat "$result_dir"/chr*.final_job.tsv 2>/dev/null || true)
done

printf 'final_job_count\t%s\n' "${#ids[@]}"
[[ "${#ids[@]}" -gt 0 ]] || exit 0

declare -A counts
for job_id in "${ids[@]}"; do
    stat="$(rtk bjobs "$job_id" 2>/dev/null | awk 'NR==2{print $3}')" || true
    if [[ -z "$stat" ]]; then
        stat="UNKNOWN"
    fi
    counts["$stat"]=$(( ${counts["$stat"]:-0} + 1 ))
    printf 'job\t%s\t%s\n' "$job_id" "$stat"
done

for key in PEND RUN DONE EXIT UNKNOWN; do
    printf 'count_%s\t%s\n' "$key" "${counts[$key]:-0}"
done
