#!/bin/bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 /path/to/launch_manifest.tsv" >&2
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
    if [[ -z "$job_id" ]]; then
        printf '%s\tmissing\n' "$label"
        return
    fi
    local status
    status="$(bjobs "$job_id" 2>/dev/null | awk 'NR==2{print $3}')"
    if [[ -z "$status" ]]; then
        status="$(bhist -n 1 "$job_id" 2>/dev/null | awk 'NR==2{print $4}')"
    fi
    printf '%s\t%s\t%s\n' "$label" "$job_id" "${status:-unknown}"
}

printf 'manifest\t%s\n' "$MANIFEST"
printf 'g_file\t%s\n' "$(value_of g_file)"
printf 'input_root\t%s\n' "$(value_of input_root)"
printf 'output_root\t%s\n' "$(value_of output_root)"

for key in \
    ncn_smoke_input_job \
    scn_smoke_input_job \
    ncn_smoke_compute_job \
    scn_smoke_compute_job \
    ncn_full_input_array_job \
    scn_full_input_array_job \
    ncn_full_compute_array_job \
    scn_full_compute_array_job \
    ncn_postprocess_job \
    scn_postprocess_job \
    ncn_diag_job \
    scn_diag_job; do
    show_job "$key" "$(value_of "$key")"
done

for path in \
    "$(value_of output_root)/NCN/NCN.normalized.tsv" \
    "$(value_of output_root)/SCN/SCN.normalized.tsv" \
    "$(value_of output_root)/NCN/NCN.manhattan.png" \
    "$(value_of output_root)/SCN/SCN.manhattan.png"; do
    if [[ -f "$path" ]]; then
        printf 'file\tpresent\t%s\n' "$path"
    else
        printf 'file\tmissing\t%s\n' "$path"
    fi
done
