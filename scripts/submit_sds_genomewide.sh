#!/bin/bash
# submit_sds_genomewide.sh — Genome-wide SDS scan orchestrator with dynamic
# cluster-aware concurrency.
#
# Submits all (population × chromosome) pairs as managed batches, re-sensing
# cluster state between waves to maximize throughput without saturating queues.
#
# Usage:
#   submit_sds_genomewide.sh \
#       --pops "NCN SCN" \
#       [--chromosomes "1-22"] \
#       [--queue QUEUE] \
#       [--max-concurrent-chrs 3] \
#       [--re-sense-interval 120] \
#       [--poll-interval 30] \
#       [--dry-run] \
#       [--g-file PATH] \
#       [--job-group LABEL] \
#       [--in-root DIR] \
#       [--out-root DIR] \
#       [--max-retries 2]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$SCRIPT_DIR/common_env.sh"
# shellcheck source=sense_cluster.sh
source "$SCRIPT_DIR/sense_cluster.sh"

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
POPS=""
CHR_RANGE="1-22"
QUEUE=""
MAX_CONCURRENT_CHRS=3
RE_SENSE_INTERVAL=120
POLL_INTERVAL=30
DRY_RUN=0
G_FILE=""
JOB_GROUP="genomewide"
IN_ROOT="$SDS_SDS_INPUT_ROOT"
OUT_ROOT="$SDS_SDS_OUTPUT_ROOT"
MAX_RETRIES=2

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --pops) POPS="$2"; shift 2 ;;
        --chromosomes) CHR_RANGE="$2"; shift 2 ;;
        --queue) QUEUE="$2"; shift 2 ;;
        --max-concurrent-chrs) MAX_CONCURRENT_CHRS="$2"; shift 2 ;;
        --re-sense-interval) RE_SENSE_INTERVAL="$2"; shift 2 ;;
        --poll-interval) POLL_INTERVAL="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        --g-file) G_FILE="$2"; shift 2 ;;
        --job-group) JOB_GROUP="$2"; shift 2 ;;
        --in-root) IN_ROOT="$2"; shift 2 ;;
        --out-root) OUT_ROOT="$2"; shift 2 ;;
        --max-retries) MAX_RETRIES="$2"; shift 2 ;;
        -h|--help)
            sed -n '2,20p' "$0" | sed 's/^# \?//' >&2
            exit 0
            ;;
        *) echo "[Error] Unknown parameter: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "$POPS" ]]; then
    echo "Usage: $0 --pops \"NCN SCN\" [--chromosomes 1-22] [--dry-run]" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Build work list
# ---------------------------------------------------------------------------
declare -a CHROMOSOMES=()
if [[ "$CHR_RANGE" =~ ^([0-9]+)-([0-9]+)$ ]]; then
    for (( c=${BASH_REMATCH[1]}; c<=${BASH_REMATCH[2]}; c++ )); do
        CHROMOSOMES+=("$c")
    done
elif [[ "$CHR_RANGE" =~ ^([0-9]+(,[0-9]+)*)$ ]]; then
    IFS=',' read -ra CHROMOSOMES <<< "$CHR_RANGE"
else
    echo "[Error] Invalid --chromosomes format: $CHR_RANGE (expected N-M or N,M,...)" >&2
    exit 1
fi

# Work list: "POP CHR" pairs
declare -a WORK_LIST=()
for pop in $POPS; do
    for chr in "${CHROMOSOMES[@]}"; do
        WORK_LIST+=("$pop $chr")
    done
done
TOTAL_TASKS=${#WORK_LIST[@]}

# ---------------------------------------------------------------------------
# State tracking
# ---------------------------------------------------------------------------
declare -A TASK_STATUS=()      # "POP:CHR" -> pending|submitted|done|failed
declare -A TASK_JOB_ID=()      # "POP:CHR" -> FINAL_JOB_ID
declare -A TASK_SUBMIT_TIME=() # "POP:CHR" -> epoch
declare -A TASK_RETRIES=()     # "POP:CHR" -> count
declare -A TASK_RETRY_AFTER=() # "POP:CHR" -> epoch (earliest retry time)

for entry in "${WORK_LIST[@]}"; do
    read -r pop chr <<< "$entry"
    key="${pop}:${chr}"
    TASK_STATUS["$key"]="pending"
    TASK_RETRIES["$key"]=0
    TASK_RETRY_AFTER["$key"]=0
done

# Sort pending queue: larger chromosomes first (chr1 = most chunks = longest)
# (already in natural order since CHROMOSOMES is 1..22)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
log() {
    echo "[GENOMEWIDE] $*" >&2
}

count_by_status() {
    local target="$1" count=0
    for key in "${!TASK_STATUS[@]}"; do
        [[ "${TASK_STATUS[$key]}" == "$target" ]] && (( count++ ))
    done
    echo "$count"
}

retry_backoff_seconds() {
    local attempt="$1"
    case "$attempt" in
        0) echo 0 ;;
        1) echo 300 ;;
        2) echo 900 ;;
        *) echo -1 ;;
    esac
}

poll_job_status() {
    local job_id="$1"
    local output
    output="$(bjobs -noheader -o 'stat' "$job_id" 2>/dev/null)" || {
        if bhist -l "$job_id" 2>/dev/null | grep -q "Done successfully" 2>/dev/null; then
            echo "DONE"
        elif bhist -l "$job_id" 2>/dev/null | grep -q "Exited with exit code" 2>/dev/null; then
            echo "FAILED"
        else
            echo "GONE"
        fi
        return
    }
    case "$output" in
        DONE)  echo "DONE" ;;
        EXIT)  echo "FAILED" ;;
        RUN|SSUSP|USUSP) echo "RUNNING" ;;
        PEND)  echo "PENDING" ;;
        *)     echo "UNKNOWN" ;;
    esac
}

submit_chr_task() {
    local pop="$1"
    local chr="$2"
    local queue_arg="${3:-}"

    local -a cmd=(
        bash "$SCRIPT_DIR/submit_sds_compute_chunked_chr.sh"
        --pop "$pop"
        --chr "$chr"
        --auto
        --in-root "$IN_ROOT"
        --out-root "$OUT_ROOT"
        --job-group "$JOB_GROUP"
    )
    [[ -n "$queue_arg" ]] && cmd+=(--queue "$queue_arg")
    [[ -n "$G_FILE" ]] && cmd+=(--g-file "$G_FILE")

    local submit_output
    submit_output="$("${cmd[@]}" 2>&1)" || {
        log "Submit failed for $pop:chr$chr: $submit_output"
        echo ""
        return 1
    }

    local final_job
    final_job="$(echo "$submit_output" | awk -F'\t' '$1 == "FINAL_JOB" { print $2 }')"
    echo "$final_job"
}

print_progress() {
    local done_count active_count pending_count failed_count now
    done_count="$(count_by_status done)"
    active_count="$(count_by_status submitted)"
    pending_count="$(count_by_status pending)"
    failed_count="$(count_by_status failed)"
    now="$(date '+%H:%M:%S')"

    log "$now | Done: ${done_count}/${TOTAL_TASKS} | Active: ${active_count} | Pending: ${pending_count} | Failed: ${failed_count}"
    if [[ -n "${SENSE_QUEUE:-}" ]]; then
        log "Cluster: $SENSE_QUEUE (score=${SENSE_SCORE:-?}, free=${SENSE_FREE_SLOTS:-?}, fairshare=${SENSE_FAIRSHARE:-?})"
    fi

    local active_keys=""
    for key in "${!TASK_STATUS[@]}"; do
        [[ "${TASK_STATUS[$key]}" == "submitted" ]] && active_keys+=" ${key}(${TASK_JOB_ID[$key]:-?})"
    done
    [[ -n "$active_keys" ]] && log "Active:${active_keys}"
}

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------
activate_sds_env

if [[ "$DRY_RUN" -eq 1 ]]; then
    log "DRY RUN: would submit $TOTAL_TASKS tasks:"
    for entry in "${WORK_LIST[@]}"; do
        read -r pop chr <<< "$entry"
        log "  $pop chr$chr"
    done
    eval "$(sense_cluster_emit "$QUEUE")"
    log "Cluster state: queue=$SENSE_QUEUE score=$SENSE_SCORE free=$SENSE_FREE_SLOTS fairshare=$SENSE_FAIRSHARE"
    log "Recommended: chunk_rows=$SENSE_CHUNK_ROWS array_parallel=$SENSE_ARRAY_PARALLEL slots=$SENSE_CHUNK_JOB_SLOTS threads=$SENSE_NUMBA_THREADS"
    exit 0
fi

for pop in $POPS; do
    for chr in "${CHROMOSOMES[@]}"; do
        local_s="$IN_ROOT/$pop/chr${chr}_s_file.txt"
        if [[ ! -f "$local_s" ]]; then
            echo "[Error] Missing s_file: $local_s (run SDS input pipeline first)" >&2
            exit 1
        fi
    done
done

# ---------------------------------------------------------------------------
# Main orchestration loop
# ---------------------------------------------------------------------------
log "Starting genome-wide SDS orchestrator: $TOTAL_TASKS tasks across $(echo "$POPS" | wc -w) population(s)"
LAST_SENSE_TIME=0
WAVE=0

while true; do
    done_count="$(count_by_status done)"
    active_count="$(count_by_status submitted)"
    pending_count="$(count_by_status pending)"
    failed_count="$(count_by_status failed)"

    # Check termination
    if (( done_count + failed_count >= TOTAL_TASKS )); then
        break
    fi

    # Re-sense cluster if interval elapsed
    now="$(date +%s)"
    if (( now - LAST_SENSE_TIME >= RE_SENSE_INTERVAL )); then
        eval "$(sense_cluster_emit "$QUEUE")"
        LAST_SENSE_TIME="$now"
    fi

    # --- Poll submitted jobs ---
    for key in "${!TASK_STATUS[@]}"; do
        [[ "${TASK_STATUS[$key]}" != "submitted" ]] && continue
        job_id="${TASK_JOB_ID[$key]}"
        [[ -z "$job_id" ]] && continue

        status="$(poll_job_status "$job_id")"
        case "$status" in
            DONE)
                TASK_STATUS["$key"]="done"
                (( active_count-- )) || true
                ;;
            FAILED|GONE)
                retries="${TASK_RETRIES[$key]}"
                if (( retries < MAX_RETRIES )); then
                    backoff="$(retry_backoff_seconds "$retries")"
                    if (( backoff >= 0 )); then
                        TASK_STATUS["$key"]="pending"
                        TASK_RETRIES["$key"]=$(( retries + 1 ))
                        TASK_RETRY_AFTER["$key"]=$(( now + backoff ))
                        TASK_JOB_ID["$key"]=""
                        log "Retry scheduled for $key (attempt $(( retries + 1 )), backoff ${backoff}s)"
                        (( active_count-- )) || true
                    else
                        TASK_STATUS["$key"]="failed"
                        (( active_count-- )) || true
                        log "FAILED permanently: $key (max retries exceeded)"
                    fi
                else
                    TASK_STATUS["$key"]="failed"
                    (( active_count-- )) || true
                    log "FAILED permanently: $key"
                fi
                ;;
        esac
    done

    # --- Compute wave budget ---
    free_slots="${SENSE_FREE_SLOTS:-0}"
    slot_budget=$(( free_slots / 10 ))
    (( slot_budget < 1 )) && slot_budget=0
    wave_budget=$(( MAX_CONCURRENT_CHRS - active_count ))
    (( wave_budget < 0 )) && wave_budget=0
    (( slot_budget < wave_budget )) && wave_budget=$slot_budget
    (( free_slots == 0 )) && wave_budget=0

    # --- Submit pending tasks up to budget ---
    submitted_this_wave=0
    for entry in "${WORK_LIST[@]}"; do
        (( submitted_this_wave >= wave_budget )) && break

        read -r pop chr <<< "$entry"
        key="${pop}:${chr}"
        [[ "${TASK_STATUS[$key]}" != "pending" ]] && continue

        # Check retry backoff
        retry_after="${TASK_RETRY_AFTER[$key]:-0}"
        (( now < retry_after )) && continue

        (( WAVE++ )) || true
        log "Submitting $key (wave $WAVE)"

        final_job="$(submit_chr_task "$pop" "$chr" "$QUEUE")"
        if [[ -n "$final_job" ]]; then
            TASK_STATUS["$key"]="submitted"
            TASK_JOB_ID["$key"]="$final_job"
            TASK_SUBMIT_TIME["$key"]="$(date +%s)"
            (( active_count++ )) || true
            (( submitted_this_wave++ )) || true
            log "Submitted $key -> FINAL_JOB=$final_job"
        else
            log "Submit returned empty job ID for $key, will retry next wave"
        fi
    done

    # --- Progress ---
    print_progress

    # --- Check if everything is done or no progress possible ---
    active_count="$(count_by_status submitted)"
    pending_count="$(count_by_status pending)"
    if (( active_count == 0 && pending_count == 0 )); then
        break
    fi
    if (( active_count == 0 && pending_count > 0 && wave_budget == 0 )); then
        log "Cluster saturated or no free slots. Waiting ${POLL_INTERVAL}s..."
    fi

    sleep "$POLL_INTERVAL"
done

# ---------------------------------------------------------------------------
# Final report
# ---------------------------------------------------------------------------
done_count="$(count_by_status done)"
failed_count="$(count_by_status failed)"

echo ""
log "=== GENOME-WIDE SDS COMPLETE ==="
log "Done: $done_count / $TOTAL_TASKS"
log "Failed: $failed_count / $TOTAL_TASKS"

if (( failed_count > 0 )); then
    log "Failed tasks:"
    for key in "${!TASK_STATUS[@]}"; do
        [[ "${TASK_STATUS[$key]}" == "failed" ]] && log "  $key"
    done
    exit 1
fi

# Emit job IDs for downstream tracking
echo "POP	CHR	FINAL_JOB"
for entry in "${WORK_LIST[@]}"; do
    read -r pop chr <<< "$entry"
    key="${pop}:${chr}"
    echo "${pop}	${chr}	${TASK_JOB_ID[$key]:-}"
done
