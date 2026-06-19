#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$SCRIPT_DIR/common_env.sh"
activate_sds_env

CHR="1"
IN_ROOT="$BASE_DIR/data/processed/sds_input_raw_joint_chr1_full"
OUT_ROOT="$BASE_DIR/data/processed/sds_output_raw_joint_olddefault_chr1"
COMPARE_DIR="$BASE_DIR/data/processed/chr1_adjudication_compare_raw_joint"
QUEUE="normal"
CHUNK_ROWS="20000"
ARRAY_PARALLEL="64"
G_FILE="$BASE_DIR/g_file.txt"
COMPARE_LABEL="NCN_SCN_MERGED3971_raw_joint"
TRACKS="NCN,SCN,MERGED3971"

usage() {
    cat >&2 <<'EOF'
Usage: launch_raw_joint_chr1_olddefault_compare.sh [--in-root DIR] [--out-root DIR] [--compare-dir DIR] [--queue QUEUE] [--chunk-rows N] [--array-parallel N] [--g-file PATH] [--compare-label LABEL] [--tracks CSV]
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --in-root) IN_ROOT="$2"; shift 2 ;;
        --out-root) OUT_ROOT="$2"; shift 2 ;;
        --compare-dir) COMPARE_DIR="$2"; shift 2 ;;
        --queue) QUEUE="$2"; shift 2 ;;
        --chunk-rows) CHUNK_ROWS="$2"; shift 2 ;;
        --array-parallel) ARRAY_PARALLEL="$2"; shift 2 ;;
        --g-file) G_FILE="$2"; shift 2 ;;
        --compare-label) COMPARE_LABEL="$2"; shift 2 ;;
        --tracks) TRACKS="$2"; shift 2 ;;
        *) usage ;;
    esac
done

[[ -f "$G_FILE" ]] || { echo "[Error] g_file not found: $G_FILE" >&2; exit 1; }

WINDOW_125="peak125=125100001-125200000"
WINDOW_143="peak143=143200001-143300000"
WINDOW_35="control35=35735220-36063259"

submit_track() {
    local pop="$1"
    local prefix="$2"

    bash "$SCRIPT_DIR/submit_chr_sds_track_with_review.sh" \
        --pop "$pop" \
        --chr "$CHR" \
        --in-root "$IN_ROOT" \
        --out-root "$OUT_ROOT" \
        --queue "$QUEUE" \
        --chunk-rows "$CHUNK_ROWS" \
        --array-parallel "$ARRAY_PARALLEL" \
        --g-file "$G_FILE" \
        --job-group "$prefix" \
        --output-prefix "$prefix" \
        --window "$WINDOW_125" \
        --window "$WINDOW_143" \
        --window "$WINDOW_35"
}

contains_track() {
    local needle="$1"
    local item
    IFS=',' read -r -a _track_list <<< "$TRACKS"
    for item in "${_track_list[@]}"; do
        [[ "$item" == "$needle" ]] && return 0
    done
    return 1
}

NCN_SUBMIT=""
SCN_SUBMIT=""
MERGED_SUBMIT=""
if contains_track NCN; then
    NCN_SUBMIT="$(submit_track NCN ncn_rawjoint_olddefault)"
    printf '%s\n' "$NCN_SUBMIT"
fi
if contains_track SCN; then
    SCN_SUBMIT="$(submit_track SCN scn_rawjoint_olddefault)"
    printf '%s\n' "$SCN_SUBMIT"
fi
if contains_track MERGED3971; then
    MERGED_SUBMIT="$(submit_track MERGED3971 merged3971_rawjoint_olddefault)"
    printf '%s\n' "$MERGED_SUBMIT"
fi

NCN_REVIEW=""
SCN_REVIEW=""
MERGED_REVIEW=""
if [[ -n "$NCN_SUBMIT" ]]; then
    NCN_REVIEW="$(awk -F '\t' '$1=="REVIEW_JOB"{print $2}' <<< "$NCN_SUBMIT")"
fi
if [[ -n "$SCN_SUBMIT" ]]; then
    SCN_REVIEW="$(awk -F '\t' '$1=="REVIEW_JOB"{print $2}' <<< "$SCN_SUBMIT")"
fi
if [[ -n "$MERGED_SUBMIT" ]]; then
    MERGED_REVIEW="$(awk -F '\t' '$1=="REVIEW_JOB"{print $2}' <<< "$MERGED_SUBMIT")"
fi

declare -a WAIT_IDS=()
[[ -n "$NCN_REVIEW" ]] && WAIT_IDS+=("$NCN_REVIEW")
[[ -n "$SCN_REVIEW" ]] && WAIT_IDS+=("$SCN_REVIEW")
[[ -n "$MERGED_REVIEW" ]] && WAIT_IDS+=("$MERGED_REVIEW")
[[ "${#WAIT_IDS[@]}" -gt 0 ]] || { echo "[Error] no review jobs to wait for" >&2; exit 1; }

mkdir -p "$COMPARE_DIR" "$OUT_ROOT/logs"
COMPARE_OUTPUT="$COMPARE_DIR/${COMPARE_LABEL}.summary.tsv"
COMPARE_JOB_NAME="SDS_chr1_compare_${COMPARE_LABEL}"
COMPARE_CMD=(
    "$SDS_ENV_PREFIX/bin/python" "$SCRIPT_DIR/compare_chr1_review_tracks.py"
    --track "NCN=$OUT_ROOT/NCN/review/ncn_rawjoint_olddefault.summary.tsv"
    --track "SCN=$OUT_ROOT/SCN/review/scn_rawjoint_olddefault.summary.tsv"
    --track "MERGED3971=$OUT_ROOT/MERGED3971/review/merged3971_rawjoint_olddefault.summary.tsv"
    --output "$COMPARE_OUTPUT"
)
printf -v COMPARE_CMD_STR '%q ' "${COMPARE_CMD[@]}"
WAIT_EXPR=""
for jid in "${WAIT_IDS[@]}"; do
    if [[ -n "$WAIT_EXPR" ]]; then
        WAIT_EXPR+=" && "
    fi
    WAIT_EXPR+="done(${jid})"
done
COMPARE_SUBMIT="$(bsub -q "$QUEUE" -w "$WAIT_EXPR" -n 1 -J "$COMPARE_JOB_NAME" -o "$OUT_ROOT/logs/${COMPARE_JOB_NAME}.out" -e "$OUT_ROOT/logs/${COMPARE_JOB_NAME}.err" /bin/bash -lc "$COMPARE_CMD_STR" < /dev/null)"
COMPARE_JOB_ID="$(sed -n 's/.*<\([0-9]\+\)>.*/\1/p' <<< "$COMPARE_SUBMIT")"
[[ -n "$COMPARE_JOB_ID" ]] || { echo "[Error] failed to parse compare job id: $COMPARE_SUBMIT" >&2; exit 1; }

cat <<EOF
NCN_REVIEW_JOB	$NCN_REVIEW
SCN_REVIEW_JOB	$SCN_REVIEW
MERGED_REVIEW_JOB	$MERGED_REVIEW
COMPARE_JOB	$COMPARE_JOB_ID
COMPARE_OUTPUT	$COMPARE_OUTPUT
EOF
