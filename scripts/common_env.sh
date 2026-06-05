#!/bin/bash
set -euo pipefail

SDS_ENV_PREFIX="/data/home/grp-wangyf/intern/miniforge3/envs/sds"
MINIFORGE_ROOT="/data/home/grp-wangyf/intern/miniforge3"
MAMBA_BIN="$MINIFORGE_ROOT/bin/mamba"

ensure_utf8_locale() {
    if [[ "${LC_ALL:-}" == "C.UTF-8" && ! -d /usr/lib/locale/C.UTF-8 ]]; then
        unset LC_ALL
    fi

    export LANG="${LANG:-C.UTF-8}"
    export LC_CTYPE="${LC_CTYPE:-${LANG}}"
}

activate_sds_env() {
    ensure_utf8_locale

    if [[ ! -d "$SDS_ENV_PREFIX" ]]; then
        echo "[Error] SDS env not found: $SDS_ENV_PREFIX" >&2
        return 1
    fi

    export PATH="$SDS_ENV_PREFIX/bin:$MINIFORGE_ROOT/bin:$PATH"
    hash -r

    if [[ ! -x "$SDS_ENV_PREFIX/bin/python" ]]; then
        echo "[Error] Python not found in SDS env: $SDS_ENV_PREFIX/bin/python" >&2
        return 1
    fi
}

activate_relate_runtime() {
    activate_sds_env

    # Relate binaries require a newer libstdc++ than the system default.
    export LD_LIBRARY_PATH="$SDS_ENV_PREFIX/lib:$MINIFORGE_ROOT/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
}

mamba_run_in_sds_env() {
    activate_sds_env

    if [[ -x "$MAMBA_BIN" ]]; then
        "$MAMBA_BIN" run -p "$SDS_ENV_PREFIX" "$@"
        return $?
    fi

    "$@"
}

sds_effective_queue() {
    local queue="${1:-${LSB_QUEUE:-normal}}"
    if [[ -z "$queue" ]]; then
        queue="normal"
    fi
    printf '%s\n' "$queue"
}

sds_queue_default_chunk_rows() {
    local queue
    queue="$(sds_effective_queue "${1:-}")"
    case "$queue" in
        smp) printf '%s\n' 10000 ;;
        *) printf '%s\n' 5000 ;;
    esac
}

sds_queue_default_array_parallel() {
    local queue
    queue="$(sds_effective_queue "${1:-}")"
    case "$queue" in
        smp) printf '%s\n' 24 ;;
        *) printf '%s\n' 32 ;;
    esac
}

sds_queue_default_chunk_job_slots() {
    local queue
    queue="$(sds_effective_queue "${1:-}")"
    case "$queue" in
        smp) printf '%s\n' 4 ;;
        *) printf '%s\n' 1 ;;
    esac
}

sds_queue_default_numba_threads() {
    local queue="${1:-}"
    local slots="${2:-}"

    if [[ -z "$slots" ]]; then
        slots="$(sds_queue_default_chunk_job_slots "$queue")"
    fi
    if [[ -z "$slots" || "$slots" -lt 1 ]]; then
        slots=1
    fi
    printf '%s\n' "$slots"
}

find_default_g_file() {
    local base_dir="$1"
    local pop="${2:-}"
    local candidate
    local -a candidates=()

    if [[ -n "$pop" ]]; then
        candidates+=(
            "$base_dir/data/g_file.${pop}.txt"
            "$base_dir/data/processed/sds_input/g_file.${pop}.txt"
            "$base_dir/data/processed/g_file.${pop}.txt"
            "/data/home/grp-wangyf/xuyuan/ms/scripts/sds_input.gamma_shapes.${pop}.final"
        )
    fi

    candidates+=(
        "$base_dir/data/g_file.txt"
        "$base_dir/data/processed/sds_input/g_file.txt"
        "$base_dir/data/processed/g_file.txt"
        "/data/home/grp-wangyf/xuyuan/ms/scripts/sds_input.gamma_shapes.final"
    )

    for candidate in "${candidates[@]}"; do
        if [[ -f "$candidate" ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done

    return 1
}

resolve_vcf_chr() {
    local vcf_file="$1"
    local target_chr="$2"

    bcftools view -h "$vcf_file" | \
        gawk -v target="$target_chr" '
            match($0, /^##contig=<ID=([^,>]+)/, a) {
                id = a[1]
                if (id == target) {
                    print id
                    exit
                }
            }
        '
}

find_first_variant_pos() {
    local vcf_file="$1"
    local target_chr="$2"
    local actual_chr

    actual_chr="$(resolve_vcf_chr "$vcf_file" "$target_chr")"
    if [[ -z "$actual_chr" ]]; then
        if [[ "$target_chr" == chr* ]]; then
            actual_chr="${target_chr#chr}"
        else
            actual_chr="chr${target_chr}"
        fi
    fi

    bcftools query -r "${actual_chr}" -f '%POS\n' "$vcf_file" | sed -n '1p'
}
