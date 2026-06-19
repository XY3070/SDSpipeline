#!/bin/bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SDS_PYTHON="/data/home/grp-wangyf/intern/miniforge3/envs/sds/bin/python"
LOGDIR="$REPO_ROOT/logs"
mkdir -p "$LOGDIR"

SIM_REPS="1000"
OUTDIR=""
POP="NCN"
CHROM="4"
SNP_ID="chr4:99317841:T:C"
POS="99317841"
CONSTANT_NE="152000"
ARTICLE_LABEL="Nyuwa"
ARTICLE_AF="0.7669"
ARTICLE_SDS="10.3947110809814"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sim-reps) SIM_REPS="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --pop) POP="$2"; shift 2 ;;
        --chrom) CHROM="$2"; shift 2 ;;
        --snp-id) SNP_ID="$2"; shift 2 ;;
        --pos) POS="$2"; shift 2 ;;
        --constant-ne) CONSTANT_NE="$2"; shift 2 ;;
        --article-label) ARTICLE_LABEL="$2"; shift 2 ;;
        --article-af) ARTICLE_AF="$2"; shift 2 ;;
        --article-sds) ARTICLE_SDS="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "$OUTDIR" ]]; then
    OUTDIR="$REPO_ROOT/tmp/gamma_single_snp_$(date +%Y%m%d_%H%M%S)"
fi
mkdir -p "$OUTDIR"

SCAN_ARGS=(
    --pop "$POP"
    --chrom "$CHROM"
    --snp-id "$SNP_ID"
    --pos "$POS"
)
if [[ -n "$CONSTANT_NE" ]]; then
    SCAN_ARGS+=(--constant-ne "$CONSTANT_NE")
fi
if [[ -n "$ARTICLE_LABEL" ]]; then
    SCAN_ARGS+=(--article-label "$ARTICLE_LABEL")
fi
if [[ -n "$ARTICLE_AF" ]]; then
    SCAN_ARGS+=(--article-af "$ARTICLE_AF")
fi
if [[ -n "$ARTICLE_SDS" ]]; then
    SCAN_ARGS+=(--article-sds "$ARTICLE_SDS")
fi

echo "[submit] preparing scenario inputs under $OUTDIR"
"$SDS_PYTHON" "$REPO_ROOT/scripts/scan_single_snp_gamma_sensitivity.py" --prepare-only --outdir "$OUTDIR" "${SCAN_ARGS[@]}" >/dev/null

read -r DAF_COMPLEMENT DAF_MAIN <<EOF2
$(
python3 - <<'PY' "$OUTDIR/target_metadata.json"
import json, sys
with open(sys.argv[1]) as handle:
    data = json.load(handle)
print(data["daf_complement_exact"], data["daf_exact"])
PY
)
EOF2

JOB_IDS=()
MANIFEST="$OUTDIR/submit_jobs.tsv"
printf "kind\tscenario\tfreq_label\tdaf\tjob_id\n" >"$MANIFEST"
SCENARIO_MANIFEST="$OUTDIR/scenario_manifest.tsv"

if [[ ! -f "$SCENARIO_MANIFEST" ]]; then
    echo "[submit] missing scenario manifest: $SCENARIO_MANIFEST" >&2
    exit 1
fi

parse_job_id() {
    sed -n 's/Job <\([0-9]\+\)>.*/\1/p'
}

while IFS=$'\t' read -r scenario_label scenario_file_label scenario_source scenario_percentile scenario_dir scenario_npz present_ne; do
    present_ne="${present_ne//$'\r'/}"
    for freq_label in complement daf; do
        if [[ "$freq_label" == "complement" ]]; then
            daf="$DAF_COMPLEMENT"
        else
            daf="$DAF_MAIN"
        fi
        gamma_prefix="$scenario_dir/${scenario_file_label}_${freq_label}.gamma"
        workdir="$scenario_dir/${scenario_file_label}_${freq_label}.work"
        log_out="$LOGDIR/${scenario_file_label}_${freq_label}_%J.out"
        log_err="$LOGDIR/${scenario_file_label}_${freq_label}_%J.err"
        submit_out="$(
            bsub -q normal \
                -J "gamma_${scenario_file_label}_${freq_label}" \
                -n 1 \
                -R "span[hosts=1]" \
                -cwd "$REPO_ROOT" \
                -o "$log_out" \
                -e "$log_err" \
                bash "$REPO_ROOT/scripts/run_single_snp_gamma_piece.sh" \
                    --pop "$POP" \
                    --daf "$daf" \
                    --scenario-npz "$scenario_npz" \
                    --present-ne "$(printf '%.0f' "$present_ne")" \
                    --gamma-prefix "$gamma_prefix" \
                    --workdir "$workdir" \
                    --sim-reps "$SIM_REPS"
        )"
        job_id="$(printf "%s\n" "$submit_out" | parse_job_id)"
        if [[ -z "$job_id" ]]; then
            echo "[submit] failed to parse job id from: $submit_out" >&2
            exit 1
        fi
        JOB_IDS+=("$job_id")
        printf "gamma\t%s\t%s\t%s\t%s\n" "$scenario_label" "$freq_label" "$daf" "$job_id" >>"$MANIFEST"
        echo "[submit] gamma job $scenario_label/$freq_label -> $job_id"
    done
done < <(tail -n +2 "$SCENARIO_MANIFEST")

dep_expr=""
for job_id in "${JOB_IDS[@]}"; do
    if [[ -z "$dep_expr" ]]; then
        dep_expr="done($job_id)"
    else
        dep_expr="${dep_expr}&&done($job_id)"
    fi
done

final_log_out="$LOGDIR/gamma_single_snp_finalize_%J.out"
final_log_err="$LOGDIR/gamma_single_snp_finalize_%J.err"
final_submit_out="$(
    bsub -q normal \
        -J "gamma_single_snp_finalize" \
        -n 1 \
        -cwd "$REPO_ROOT" \
        -w "$dep_expr" \
        -o "$final_log_out" \
        -e "$final_log_err" \
        "$SDS_PYTHON" "$REPO_ROOT/scripts/scan_single_snp_gamma_sensitivity.py" --outdir "$OUTDIR" "${SCAN_ARGS[@]}"
)"
final_job_id="$(printf "%s\n" "$final_submit_out" | parse_job_id)"
if [[ -z "$final_job_id" ]]; then
    echo "[submit] failed to parse finalize job id from: $final_submit_out" >&2
    exit 1
fi
printf "finalize\tall\tall\tNA\t%s\n" "$final_job_id" >>"$MANIFEST"

report_log_out="$LOGDIR/gamma_single_snp_report_%J.out"
report_log_err="$LOGDIR/gamma_single_snp_report_%J.err"
report_submit_out="$(
    bsub -q normal \
        -J "gamma_single_snp_report" \
        -n 1 \
        -cwd "$REPO_ROOT" \
        -w "done($final_job_id)" \
        -o "$report_log_out" \
        -e "$report_log_err" \
        "$SDS_PYTHON" "$REPO_ROOT/scripts/report_single_snp_gamma_sensitivity.py" --outdir "$OUTDIR"
)"
report_job_id="$(printf "%s\n" "$report_submit_out" | parse_job_id)"
if [[ -z "$report_job_id" ]]; then
    echo "[submit] failed to parse report job id from: $report_submit_out" >&2
    exit 1
fi
printf "report\tall\tall\tNA\t%s\n" "$report_job_id" >>"$MANIFEST"

echo "[submit] finalize job -> $final_job_id"
echo "[submit] report job -> $report_job_id"
echo "[submit] outdir: $OUTDIR"
echo "[submit] manifest: $MANIFEST"
