#!/bin/bash
# sense_cluster.sh — LSF cluster resource sensing library for SDS pipeline.
#
# Usage:
#   source sense_cluster.sh        # import functions
#   bash sense_cluster.sh --eval   # emit sourceable SENSE_* variables
#   bash sense_cluster.sh --report # human-readable diagnostic to stderr
#
# All functions are safe to call when LSF is unavailable; they return
# conservative defaults instead of failing.

set -euo pipefail

SENSE_QUEUE_LIST="${SENSE_QUEUE_LIST:-normal smp}"
SENSE_USER="${SENSE_USER:-$(whoami 2>/dev/null || echo "${USER:-unknown}")}"

SENSE_MAX_ARRAY_PARALLEL="${SENSE_MAX_ARRAY_PARALLEL:-64}"
SENSE_MIN_ARRAY_PARALLEL="${SENSE_MIN_ARRAY_PARALLEL:-4}"
SENSE_MAX_CHUNK_JOB_SLOTS="${SENSE_MAX_CHUNK_JOB_SLOTS:-8}"
SENSE_MIN_CHUNK_JOB_SLOTS="${SENSE_MIN_CHUNK_JOB_SLOTS:-1}"
SENSE_MIN_CHUNK_ROWS="${SENSE_MIN_CHUNK_ROWS:-1000}"
SENSE_MAX_CHUNK_ROWS="${SENSE_MAX_CHUNK_ROWS:-20000}"

# ---------------------------------------------------------------------------
# LSF availability check
# ---------------------------------------------------------------------------

_sense_lsf_available() {
    command -v bqueues >/dev/null 2>&1 \
        && command -v bhosts >/dev/null 2>&1 \
        && command -v bjobs >/dev/null 2>&1
}

# ---------------------------------------------------------------------------
# Queue introspection
# ---------------------------------------------------------------------------

sense_queue_hostgroup() {
    local queue="$1"
    _sense_lsf_available || { echo ""; return 1; }
    bqueues -l "$queue" 2>/dev/null \
        | awk '/^HOSTS:/{gsub(/\//, "", $2); print $2; exit}'
}

sense_queue_free_slots() {
    local queue="$1"
    local hostgroup
    hostgroup="$(sense_queue_hostgroup "$queue" 2>/dev/null)" || { echo 0; return; }
    [[ -z "$hostgroup" ]] && { echo 0; return; }
    _sense_lsf_available || { echo 0; return; }
    bhosts -w "$hostgroup" 2>/dev/null \
        | awk 'NR>1 && $2 == "ok" { free = $4 - $6; if (free > 0) total += free } END { print total+0 }'
}

sense_queue_total_slots() {
    local queue="$1"
    local hostgroup
    hostgroup="$(sense_queue_hostgroup "$queue" 2>/dev/null)" || { echo 0; return; }
    [[ -z "$hostgroup" ]] && { echo 0; return; }
    _sense_lsf_available || { echo 0; return; }
    bhosts -w "$hostgroup" 2>/dev/null \
        | awk 'NR>1 { total += $4 } END { print total+0 }'
}

sense_queue_load() {
    local queue="$1"
    _sense_lsf_available || { echo "0 0 0"; return; }
    bqueues -l "$queue" 2>/dev/null \
        | awk '
            /^PRIO NICE/ { header_seen=1; next }
            header_seen && /^ *[0-9]/ {
                print $8, $9, $10
                exit
            }
          '
}

sense_queue_status() {
    local queue="$1"
    _sense_lsf_available || { echo "unknown"; return 1; }
    bqueues -w "$queue" 2>/dev/null \
        | awk 'NR>1 { print $3 }'
}

sense_user_fairshare_priority() {
    local queue="$1"
    local user="${2:-$SENSE_USER}"
    _sense_lsf_available || { echo "0"; return; }
    bqueues -l "$queue" 2>/dev/null \
        | awk -v user="$user" '
            /^SHARE_INFO_FOR:/ { in_section=1; next }
            in_section && /^$/ { in_section=0 }
            in_section && $1 == user { print $3; exit }
          '
}

sense_queue_pending_jobs() {
    local queue="$1"
    local user="${2:-$SENSE_USER}"
    _sense_lsf_available || { echo 0; return; }
    local count
    count="$(bjobs -p -u "$user" -q "$queue" 2>/dev/null | tail -n +2 | wc -l)" || true
    echo "${count:-0}" | tr -d ' '
}

sense_queue_running_jobs() {
    local queue="$1"
    local user="${2:-$SENSE_USER}"
    _sense_lsf_available || { echo 0; return; }
    local count
    count="$(bjobs -r -u "$user" -q "$queue" 2>/dev/null | tail -n +2 | wc -l)" || true
    echo "${count:-0}" | tr -d ' '
}

# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

sense_score_queue() {
    local queue="$1"
    local user="${2:-$SENSE_USER}"

    local status
    status="$(sense_queue_status "$queue" 2>/dev/null || echo "unknown")"
    if [[ "$status" != Open:* ]]; then
        echo "0"
        return
    fi

    local fairshare free_slots pending_user
    fairshare="$(sense_user_fairshare_priority "$queue" "$user" 2>/dev/null || echo 0)"
    fairshare="${fairshare:-0}"
    free_slots="$(sense_queue_free_slots "$queue" 2>/dev/null || echo 0)"
    free_slots="${free_slots:-0}"
    pending_user="$(sense_queue_pending_jobs "$queue" "$user" 2>/dev/null || echo 0)"
    pending_user="${pending_user:-0}"

    awk -v fs="$fairshare" -v free="$free_slots" -v pend="$pending_user" \
        'BEGIN { printf "%.6f\n", fs * free / (1 + pend) }'
}

sense_pick_best_queue() {
    local user="${1:-$SENSE_USER}"
    local best_queue=""
    local best_score="0"

    for queue in $SENSE_QUEUE_LIST; do
        local score
        score="$(sense_score_queue "$queue" "$user" 2>/dev/null || echo 0)"
        if awk -v s="$score" -v b="$best_score" 'BEGIN { exit !(s > b) }' 2>/dev/null; then
            best_score="$score"
            best_queue="$queue"
        fi
    done

    if [[ -z "$best_queue" ]]; then
        best_queue="${SENSE_QUEUE_LIST%% *}"
        best_score="0"
    fi

    echo "$best_queue $best_score"
}

# ---------------------------------------------------------------------------
# Parameter derivation
# ---------------------------------------------------------------------------

sense_derive_sds_params() {
    local queue="$1"
    local free_slots
    free_slots="$(sense_queue_free_slots "$queue" 2>/dev/null || echo 0)"
    free_slots="${free_slots:-0}"

    local chunk_job_slots array_parallel chunk_rows numba_threads

    if [[ "$queue" == "smp" ]]; then
        chunk_job_slots=4
        array_parallel=$(( free_slots / (chunk_job_slots * 2) ))
        chunk_rows=10000
        numba_threads=$chunk_job_slots
    else
        chunk_job_slots=1
        array_parallel=$(( free_slots * 3 / 4 ))
        chunk_rows=10000
        numba_threads=1
    fi

    (( array_parallel < SENSE_MIN_ARRAY_PARALLEL )) && array_parallel=$SENSE_MIN_ARRAY_PARALLEL
    (( array_parallel > SENSE_MAX_ARRAY_PARALLEL )) && array_parallel=$SENSE_MAX_ARRAY_PARALLEL
    (( chunk_job_slots < SENSE_MIN_CHUNK_JOB_SLOTS )) && chunk_job_slots=$SENSE_MIN_CHUNK_JOB_SLOTS
    (( chunk_job_slots > SENSE_MAX_CHUNK_JOB_SLOTS )) && chunk_job_slots=$SENSE_MAX_CHUNK_JOB_SLOTS

    echo "$chunk_job_slots $array_parallel $chunk_rows $numba_threads"
}

# ---------------------------------------------------------------------------
# Top-level emitter
# ---------------------------------------------------------------------------

sense_cluster_emit() {
    local requested_queue="${1:-}"
    local user="${2:-$SENSE_USER}"

    if ! _sense_lsf_available; then
        local fallback_queue="${requested_queue:-normal}"
        cat <<EOF
SENSE_QUEUE="$fallback_queue"
SENSE_SCORE="0"
SENSE_FREE_SLOTS="0"
SENSE_QUEUE_NJOBS="0"
SENSE_QUEUE_PEND="0"
SENSE_QUEUE_RUN="0"
SENSE_FAIRSHARE="0"
SENSE_CHUNK_JOB_SLOTS="1"
SENSE_ARRAY_PARALLEL="32"
SENSE_CHUNK_ROWS="5000"
SENSE_NUMBA_THREADS="1"
SENSE_USER="$user"
SENSE_TIMESTAMP="$(date -Iseconds 2>/dev/null || date '+%Y-%m-%dT%H:%M:%S')"
SENSE_LSF_AVAILABLE="0"
EOF
        return 0
    fi

    local queue score
    if [[ -n "$requested_queue" ]]; then
        queue="$requested_queue"
        score="$(sense_score_queue "$queue" "$user" 2>/dev/null || echo 0)"
    else
        read -r queue score <<< "$(sense_pick_best_queue "$user")"
    fi

    local chunk_job_slots array_parallel chunk_rows numba_threads
    read -r chunk_job_slots array_parallel chunk_rows numba_threads \
        <<< "$(sense_derive_sds_params "$queue")"

    local free_slots queue_njobs queue_pend queue_run fairshare
    free_slots="$(sense_queue_free_slots "$queue" 2>/dev/null || echo 0)"
    read -r queue_njobs queue_pend queue_run <<< "$(sense_queue_load "$queue" 2>/dev/null || echo "0 0 0")"
    fairshare="$(sense_user_fairshare_priority "$queue" "$user" 2>/dev/null || echo 0)"

    cat <<EOF
SENSE_QUEUE="$queue"
SENSE_SCORE="${score:-0}"
SENSE_FREE_SLOTS="${free_slots:-0}"
SENSE_QUEUE_NJOBS="${queue_njobs:-0}"
SENSE_QUEUE_PEND="${queue_pend:-0}"
SENSE_QUEUE_RUN="${queue_run:-0}"
SENSE_FAIRSHARE="${fairshare:-0}"
SENSE_CHUNK_JOB_SLOTS="$chunk_job_slots"
SENSE_ARRAY_PARALLEL="$array_parallel"
SENSE_CHUNK_ROWS="$chunk_rows"
SENSE_NUMBA_THREADS="$numba_threads"
SENSE_USER="$user"
SENSE_TIMESTAMP="$(date -Iseconds 2>/dev/null || date '+%Y-%m-%dT%H:%M:%S')"
SENSE_LSF_AVAILABLE="1"
EOF
}

# ---------------------------------------------------------------------------
# Human-readable report
# ---------------------------------------------------------------------------

_sense_report() {
    local user="${1:-$SENSE_USER}"

    if ! _sense_lsf_available; then
        echo "[sense_cluster] LSF commands not found. Cannot generate report." >&2
        return 1
    fi

    echo "[sense_cluster] Cluster resource report for user: $user" >&2
    echo "[sense_cluster] Timestamp: $(date -Iseconds 2>/dev/null || date)" >&2
    echo "" >&2
    printf '%-12s %-12s %-10s %-10s %-10s %-10s %-10s\n' \
        "QUEUE" "STATUS" "FREE" "TOTAL" "FAIRSHARE" "PENDING" "SCORE" >&2
    printf '%-12s %-12s %-10s %-10s %-10s %-10s %-10s\n' \
        "-----" "------" "----" "-----" "---------" "-------" "-----" >&2

    local best_queue="" best_score="0"
    for queue in $SENSE_QUEUE_LIST; do
        local status free total fairshare pending score
        status="$(sense_queue_status "$queue" 2>/dev/null || echo "unknown")"
        free="$(sense_queue_free_slots "$queue" 2>/dev/null || echo 0)"
        total="$(sense_queue_total_slots "$queue" 2>/dev/null || echo 0)"
        fairshare="$(sense_user_fairshare_priority "$queue" "$user" 2>/dev/null || echo 0)"
        fairshare="${fairshare:-0}"
        pending="$(sense_queue_pending_jobs "$queue" "$user" 2>/dev/null || echo 0)"
        score="$(sense_score_queue "$queue" "$user" 2>/dev/null || echo 0)"

        printf '%-12s %-12s %-10s %-10s %-10s %-10s %-10s\n' \
            "$queue" "$status" "$free" "$total" "$fairshare" "$pending" "$score" >&2

        if awk -v s="$score" -v b="$best_score" 'BEGIN { exit !(s > b) }' 2>/dev/null; then
            best_score="$score"
            best_queue="$queue"
        fi
    done

    echo "" >&2
    echo "[sense_cluster] Best queue: ${best_queue:-none} (score: $best_score)" >&2

    if [[ -n "$best_queue" ]]; then
        local params
        params="$(sense_derive_sds_params "$best_queue")"
        read -r cjs ap cr nt <<< "$params"
        echo "[sense_cluster] Recommended params: chunk_job_slots=$cjs array_parallel=$ap chunk_rows=$cr numba_threads=$nt" >&2
    fi
}

# ---------------------------------------------------------------------------
# CLI entry point (only runs when executed, not sourced)
# ---------------------------------------------------------------------------

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    MODE="eval"
    REQUESTED_QUEUE=""

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --eval)    MODE="eval"; shift ;;
            --report)  MODE="report"; shift ;;
            --queue)   REQUESTED_QUEUE="$2"; shift 2 ;;
            --user)    SENSE_USER="$2"; shift 2 ;;
            -h|--help)
                echo "Usage: sense_cluster.sh [--eval|--report] [--queue QUEUE] [--user USER]" >&2
                exit 0
                ;;
            *) echo "Unknown option: $1" >&2; exit 1 ;;
        esac
    done

    case "$MODE" in
        eval)   sense_cluster_emit "$REQUESTED_QUEUE" "$SENSE_USER" ;;
        report) _sense_report "$SENSE_USER" ;;
    esac
fi
