#!/bin/bash
#BSUB -J SDS_INPUT[1-22]          # Job Array 处理 22 条染色体
#BSUB -q normal                   # 普通队列
#BSUB -n 2                        # 默认仅覆盖输入阶段；normal compute 默认转发到 chunked LSF array
#BSUB -o logs/sds_%J_%I.out
#BSUB -e logs/sds_%J_%I.err

set -euo pipefail

# When submitted via `bsub < script`, LSF copies the script to ~/.lsbatch,
# so BASH_SOURCE no longer points at the repo. Fall back to the job cwd.
if [[ -n "${BASH_SOURCE[0]:-}" ]]; then
    CANDIDATE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
else
    CANDIDATE_DIR=""
fi

if [[ -n "$CANDIDATE_DIR" && -f "$CANDIDATE_DIR/common_env.sh" ]]; then
    SCRIPT_DIR="$CANDIDATE_DIR"
elif [[ -f "$(pwd)/scripts/common_env.sh" ]]; then
    SCRIPT_DIR="$(pwd)/scripts"
else
    echo "[Error] Could not locate scripts/common_env.sh" >&2
    exit 1
fi
# shellcheck source=/share/home/grp-wangyf/xuyuan/sds/scripts/common_env.sh
source "$SCRIPT_DIR/common_env.sh"
activate_sds_env

POP="${POP:-NCN}"
STAGE="${STAGE:-input}"
TEST_MODE="${TEST_MODE:-0}"
CHR="${LSB_JOBINDEX:-}"
OUT_ROOT="${OUT_ROOT:-}"
G_FILE="${G_FILE:-}"
INIT="${INIT:-}"
S_FILE_NCOL="${S_FILE_NCOL:-}"
CHUNK_ROWS="${CHUNK_ROWS:-}"
MAX_PARALLEL_CHUNKS="${MAX_PARALLEL_CHUNKS:-}"
ARRAY_PARALLEL="${ARRAY_PARALLEL:-}"
CHUNK_JOB_SLOTS="${CHUNK_JOB_SLOTS:-}"
QUEUE="$(sds_effective_queue "${QUEUE:-${LSB_QUEUE:-normal}}")"
SDS_COMPUTE_MODE="${SDS_COMPUTE_MODE:-auto}"

if [[ -z "$CHR" ]]; then
    echo "[Error] LSB_JOBINDEX is not set" >&2
    exit 1
fi

TEST_ARG=()
if [[ "$TEST_MODE" = "1" ]]; then
    TEST_ARG+=(--test)
fi

COMPUTE_ARG=()
if [[ -n "$OUT_ROOT" ]]; then
    COMPUTE_ARG+=(--out-root "$OUT_ROOT")
fi
if [[ -n "$G_FILE" ]]; then
    COMPUTE_ARG+=(--g-file "$G_FILE")
fi
if [[ -n "$INIT" ]]; then
    COMPUTE_ARG+=(--init "$INIT")
fi
if [[ -n "$S_FILE_NCOL" ]]; then
    COMPUTE_ARG+=(--s-file-ncol "$S_FILE_NCOL")
fi
if [[ -n "$CHUNK_ROWS" ]]; then
    COMPUTE_ARG+=(--chunk-rows "$CHUNK_ROWS")
fi
if [[ -n "$MAX_PARALLEL_CHUNKS" ]]; then
    COMPUTE_ARG+=(--max-parallel-chunks "$MAX_PARALLEL_CHUNKS")
fi

CHUNKED_ARG=(
    --pop "$POP"
    --chr "$CHR"
    --queue "$QUEUE"
    --wait
)
if [[ -n "$OUT_ROOT" ]]; then
    CHUNKED_ARG+=(--out-root "$OUT_ROOT")
fi
if [[ -n "$G_FILE" ]]; then
    CHUNKED_ARG+=(--g-file "$G_FILE")
fi
if [[ -n "$INIT" ]]; then
    CHUNKED_ARG+=(--init "$INIT")
fi
if [[ -n "$S_FILE_NCOL" ]]; then
    CHUNKED_ARG+=(--s-file-ncol "$S_FILE_NCOL")
fi
if [[ -n "$CHUNK_ROWS" ]]; then
    CHUNKED_ARG+=(--chunk-rows "$CHUNK_ROWS")
fi
if [[ -n "$ARRAY_PARALLEL" ]]; then
    CHUNKED_ARG+=(--array-parallel "$ARRAY_PARALLEL")
fi
if [[ -n "$CHUNK_JOB_SLOTS" ]]; then
    CHUNKED_ARG+=(--chunk-job-slots "$CHUNK_JOB_SLOTS")
fi

should_use_chunked_submit() {
    if [[ "$TEST_MODE" = "1" ]]; then
        return 1
    fi
    case "$SDS_COMPUTE_MODE" in
        chunked) return 0 ;;
        local) return 1 ;;
        auto)
            [[ "$QUEUE" == "normal" ]]
            return
            ;;
        *)
            echo "[Error] Unsupported SDS_COMPUTE_MODE: $SDS_COMPUTE_MODE" >&2
            exit 1
            ;;
    esac
}

case "$STAGE" in
    input)
        if [[ ${#TEST_ARG[@]} -gt 0 ]]; then
            bash "$SCRIPT_DIR/run_sds_input.sh" --pop "$POP" --chr "$CHR" "${TEST_ARG[@]}"
        else
            bash "$SCRIPT_DIR/run_sds_input.sh" --pop "$POP" --chr "$CHR"
        fi
        ;;
    compute)
        if should_use_chunked_submit; then
            bash "$SCRIPT_DIR/submit_sds_compute_chunked_chr.sh" "${CHUNKED_ARG[@]}"
        elif [[ ${#TEST_ARG[@]} -gt 0 ]]; then
            bash "$SCRIPT_DIR/run_sds_compute.sh" --pop "$POP" --chr "$CHR" "${TEST_ARG[@]}" "${COMPUTE_ARG[@]}"
        else
            bash "$SCRIPT_DIR/run_sds_compute.sh" --pop "$POP" --chr "$CHR" "${COMPUTE_ARG[@]}"
        fi
        ;;
    full)
        if [[ ${#TEST_ARG[@]} -gt 0 ]]; then
            bash "$SCRIPT_DIR/run_sds_input.sh" --pop "$POP" --chr "$CHR" "${TEST_ARG[@]}"
            bash "$SCRIPT_DIR/run_sds_compute.sh" --pop "$POP" --chr "$CHR" "${TEST_ARG[@]}" "${COMPUTE_ARG[@]}"
        else
            bash "$SCRIPT_DIR/run_sds_input.sh" --pop "$POP" --chr "$CHR"
            if should_use_chunked_submit; then
                bash "$SCRIPT_DIR/submit_sds_compute_chunked_chr.sh" "${CHUNKED_ARG[@]}"
            else
                bash "$SCRIPT_DIR/run_sds_compute.sh" --pop "$POP" --chr "$CHR" "${COMPUTE_ARG[@]}"
            fi
        fi
        ;;
    *)
        echo "[Error] Unsupported STAGE: $STAGE" >&2
        echo "Valid values: input, compute, full" >&2
        exit 1
        ;;
esac
