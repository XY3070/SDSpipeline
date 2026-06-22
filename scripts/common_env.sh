#!/bin/bash
set -euo pipefail

COMMON_ENV_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SDS_PIPELINE_BOOTSTRAP_ROOT="$(cd "$COMMON_ENV_DIR/.." && pwd)"

SDS_PIPELINE_ROOT="${SDS_PIPELINE_ROOT:-$SDS_PIPELINE_BOOTSTRAP_ROOT}"
SDS_LOCAL_CONFIG="${SDS_LOCAL_CONFIG:-$SDS_PIPELINE_ROOT/config/paths.env}"
if [[ -f "$SDS_LOCAL_CONFIG" ]]; then
    # shellcheck source=/dev/null
    source "$SDS_LOCAL_CONFIG"
fi

SDS_PIPELINE_ROOT="${SDS_PIPELINE_ROOT:-$SDS_PIPELINE_BOOTSTRAP_ROOT}"
SDSLOG_ROOT="${SDSLOG_ROOT:-$SDS_PIPELINE_ROOT/../SDSlog}"
SDS_WORKSPACE_ROOT="${SDS_WORKSPACE_ROOT:-$SDS_PIPELINE_ROOT/../SDSworkspace}"

SDS_INPUT_ROOT="${SDS_INPUT_ROOT:-$SDS_WORKSPACE_ROOT/input}"
SDS_RESULTS_ROOT="${SDS_RESULTS_ROOT:-$SDS_WORKSPACE_ROOT/results}"
SDS_RUNS_ROOT="${SDS_RUNS_ROOT:-$SDS_WORKSPACE_ROOT/runs}"
SDS_ONEOFF_ROOT="${SDS_ONEOFF_ROOT:-$SDS_WORKSPACE_ROOT/oneoff}"
SDS_PROVENANCE_ROOT="${SDS_PROVENANCE_ROOT:-$SDS_WORKSPACE_ROOT/provenance}"
SDS_CACHE_ROOT="${SDS_CACHE_ROOT:-$SDS_WORKSPACE_ROOT/cache}"
SDS_TMP_ROOT="${SDS_TMP_ROOT:-$SDS_WORKSPACE_ROOT/tmp}"
SDS_EXTERNAL_ROOT="${SDS_EXTERNAL_ROOT:-$SDS_WORKSPACE_ROOT/external}"
SDS_DATA_PROCESSED_ROOT="${SDS_DATA_PROCESSED_ROOT:-$SDS_RESULTS_ROOT/legacy}"

SDS_VCF_ROOT="${SDS_VCF_ROOT:-$SDS_INPUT_ROOT/raw/vcf}"
SDS_SAMPLE_LIST_ROOT="${SDS_SAMPLE_LIST_ROOT:-$SDS_INPUT_ROOT/freeze/sample_lists}"
SDS_FREEZE_ROOT="${SDS_FREEZE_ROOT:-$SDS_INPUT_ROOT/freeze}"
SDS_REFERENCE_ROOT="${SDS_REFERENCE_ROOT:-$SDS_INPUT_ROOT/reference}"

SDS_PRODUCTION_ROOT="${SDS_PRODUCTION_ROOT:-$SDS_RESULTS_ROOT/production}"
SDS_AUDIT_ROOT="${SDS_AUDIT_ROOT:-$SDS_RESULTS_ROOT/audit}"
SDS_SDS_INPUT_ROOT="${SDS_SDS_INPUT_ROOT:-$SDS_PRODUCTION_ROOT/sds_input}"
SDS_SDS_OUTPUT_ROOT="${SDS_SDS_OUTPUT_ROOT:-$SDS_PRODUCTION_ROOT/sds_output}"
SDS_DEMOGRAPHY_ROOT="${SDS_DEMOGRAPHY_ROOT:-$SDS_PRODUCTION_ROOT/demography}"
SDS_GAMMA_ROOT="${SDS_GAMMA_ROOT:-$SDS_PRODUCTION_ROOT/gamma}"
SDS_RELATE_OUTPUT_ROOT="${SDS_RELATE_OUTPUT_ROOT:-$SDS_PRODUCTION_ROOT/relate_clues2}"
SDS_RAW_HEADER_INTERSECTION_ROOT="${SDS_RAW_HEADER_INTERSECTION_ROOT:-$SDS_AUDIT_ROOT/raw_header_intersections}"

SDS_MS_ROOT="${SDS_MS_ROOT:-$SDS_EXTERNAL_ROOT/ms}"
SDS_MS_SCRIPTS_DIR="${SDS_MS_SCRIPTS_DIR:-$SDS_MS_ROOT/scripts}"
SDS_MS_BINARY="${SDS_MS_BINARY:-$SDS_MS_ROOT/msdir/ms}"
SDS_BACKWARD_SCRIPT="${SDS_BACKWARD_SCRIPT:-$SDS_MS_SCRIPTS_DIR/backward.py}"
SDS_RELATE_DIR="${SDS_RELATE_DIR:-$SDS_EXTERNAL_ROOT/relate}"
SDS_SMCPP_ROOT="${SDS_SMCPP_ROOT:-$SDS_EXTERNAL_ROOT/smcpp}"
SDS_SMCPP_IMAGE="${SDS_SMCPP_IMAGE:-$SDS_SMCPP_ROOT/smcpp_latest.sif}"
SDS_SMCPP_ROOTFS="${SDS_SMCPP_ROOTFS:-$SDS_SMCPP_ROOT/smcpp_latest.sandbox}"
SDS_SMCPP_MPLCONFIGDIR="${SDS_SMCPP_MPLCONFIGDIR:-$SDS_SMCPP_ROOT/smcpp_matplotlib}"
SDS_SMCPP_RUNTIME_DIR="${SDS_SMCPP_RUNTIME_DIR:-$SDS_CACHE_ROOT/smcpp_runtime}"
SDS_SMCPP_SINGULARITY_CONFDIR="${SDS_SMCPP_SINGULARITY_CONFDIR:-$SDS_SMCPP_ROOT/singularity_conf}"
SDS_SMCPP_SINGULARITY_BIN="${SDS_SMCPP_SINGULARITY_BIN:-$SDS_EXTERNAL_ROOT/singularity/bin/singularity}"

SDS_LEGACY_CODE_ROOT="${SDS_LEGACY_CODE_ROOT:-$SDS_PIPELINE_ROOT/../sds}"
SDS_LEGACY_BENCHMARK_ROOT="${SDS_LEGACY_BENCHMARK_ROOT:-$SDS_PIPELINE_ROOT/../benchmark}"

SDS_ENV_PREFIX="${SDS_ENV_PREFIX:-${CONDA_PREFIX:-}}"
if [[ -z "${MINIFORGE_ROOT:-}" && -n "$SDS_ENV_PREFIX" ]]; then
    MINIFORGE_ROOT="$(cd "$(dirname "$(dirname "$SDS_ENV_PREFIX")")" && pwd 2>/dev/null || true)"
fi
MINIFORGE_ROOT="${MINIFORGE_ROOT:-}"

SDS_PYTHON="${SDS_PYTHON:-}"
if [[ -z "$SDS_PYTHON" && -n "$SDS_ENV_PREFIX" ]]; then
    SDS_PYTHON="$SDS_ENV_PREFIX/bin/python"
fi
if [[ -z "$SDS_PYTHON" ]]; then
    SDS_PYTHON="$(command -v python3 2>/dev/null || true)"
fi

MAMBA_BIN="${MAMBA_BIN:-}"
if [[ -z "$MAMBA_BIN" && -n "$MINIFORGE_ROOT" ]]; then
    MAMBA_BIN="$MINIFORGE_ROOT/bin/mamba"
fi

export SDS_PIPELINE_ROOT SDS_LOCAL_CONFIG SDSLOG_ROOT SDS_WORKSPACE_ROOT
export SDS_INPUT_ROOT SDS_RESULTS_ROOT SDS_RUNS_ROOT SDS_ONEOFF_ROOT SDS_PROVENANCE_ROOT
export SDS_CACHE_ROOT SDS_TMP_ROOT SDS_EXTERNAL_ROOT SDS_VCF_ROOT SDS_SAMPLE_LIST_ROOT
export SDS_FREEZE_ROOT SDS_REFERENCE_ROOT SDS_PRODUCTION_ROOT SDS_AUDIT_ROOT
export SDS_SDS_INPUT_ROOT SDS_SDS_OUTPUT_ROOT SDS_DEMOGRAPHY_ROOT SDS_GAMMA_ROOT
export SDS_RELATE_OUTPUT_ROOT SDS_RAW_HEADER_INTERSECTION_ROOT SDS_MS_ROOT SDS_MS_SCRIPTS_DIR
export SDS_MS_BINARY SDS_BACKWARD_SCRIPT SDS_RELATE_DIR SDS_SMCPP_ROOT SDS_SMCPP_IMAGE
export SDS_SMCPP_ROOTFS SDS_SMCPP_MPLCONFIGDIR SDS_SMCPP_RUNTIME_DIR
export SDS_SMCPP_SINGULARITY_CONFDIR SDS_SMCPP_SINGULARITY_BIN
export SDS_LEGACY_CODE_ROOT SDS_LEGACY_BENCHMARK_ROOT
export SDS_DATA_PROCESSED_ROOT
export SDS_ENV_PREFIX MINIFORGE_ROOT SDS_PYTHON MAMBA_BIN

ensure_utf8_locale() {
    if [[ "${LC_ALL:-}" == "C.UTF-8" && ! -d /usr/lib/locale/C.UTF-8 ]]; then
        unset LC_ALL
    fi

    export LANG="${LANG:-C.UTF-8}"
    export LC_CTYPE="${LC_CTYPE:-${LANG}}"
}

activate_sds_env() {
    ensure_utf8_locale

    if [[ -n "$SDS_ENV_PREFIX" && -d "$SDS_ENV_PREFIX" ]]; then
        export PATH="$SDS_ENV_PREFIX/bin${MINIFORGE_ROOT:+:$MINIFORGE_ROOT/bin}:$PATH"
        hash -r
        if [[ -z "$SDS_PYTHON" || ! -x "$SDS_PYTHON" ]]; then
            SDS_PYTHON="$SDS_ENV_PREFIX/bin/python"
        fi
    fi

    if [[ -z "$SDS_PYTHON" || ! -x "$SDS_PYTHON" ]]; then
        echo "[Error] SDS python runtime not found. Set SDS_PYTHON or SDS_ENV_PREFIX in $SDS_LOCAL_CONFIG" >&2
        return 1
    fi

    export SDS_PYTHON
}

activate_relate_runtime() {
    activate_sds_env

    # Relate binaries require a newer libstdc++ than the system default.
    if [[ -n "$SDS_ENV_PREFIX" && -d "$SDS_ENV_PREFIX/lib" ]]; then
        export LD_LIBRARY_PATH="$SDS_ENV_PREFIX/lib${MINIFORGE_ROOT:+:$MINIFORGE_ROOT/lib}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
    fi
}

mamba_run_in_sds_env() {
    activate_sds_env

    if [[ -n "$MAMBA_BIN" && -x "$MAMBA_BIN" && -n "$SDS_ENV_PREFIX" && -d "$SDS_ENV_PREFIX" ]]; then
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
        *) printf '%s\n' 10000 ;;
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

first_existing_file() {
    local candidate
    for candidate in "$@"; do
        if [[ -n "$candidate" && -f "$candidate" ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    return 1
}

sds_env_label() {
    if [[ -n "$SDS_ENV_PREFIX" && -d "$SDS_ENV_PREFIX" ]]; then
        printf '%s\n' "$SDS_ENV_PREFIX"
        return 0
    fi
    printf '%s\n' "$SDS_PYTHON"
}

find_default_g_file() {
    local base_dir="$1"
    local pop="${2:-}"
    local -a candidates=()

    candidates+=(
        "$SDS_GAMMA_ROOT/g_file.Gravel_EAS.txt"
    )

    if [[ -n "$pop" ]]; then
        candidates+=(
            "$SDS_GAMMA_ROOT/g_file.${pop}.txt"
            "$SDS_GAMMA_ROOT/${pop}/g_file.txt"
            "$SDS_GAMMA_ROOT/sds_input.gamma_shapes.${pop}.final"
            "$SDS_MS_SCRIPTS_DIR/sds_input.gamma_shapes.${pop}.final"
            "$SDS_LEGACY_CODE_ROOT/data/g_file.${pop}.txt"
            "$SDS_LEGACY_CODE_ROOT/data/processed/sds_input/g_file.${pop}.txt"
            "$SDS_LEGACY_CODE_ROOT/data/processed/g_file.${pop}.txt"
            "$base_dir/data/g_file.${pop}.txt"
            "$base_dir/data/processed/sds_input/g_file.${pop}.txt"
            "$base_dir/data/processed/g_file.${pop}.txt"
        )
    fi

    candidates+=(
        "$SDS_GAMMA_ROOT/g_file.txt"
        "$SDS_GAMMA_ROOT/default/g_file.txt"
        "$SDS_GAMMA_ROOT/sds_input.gamma_shapes.final"
        "$SDS_MS_SCRIPTS_DIR/sds_input.gamma_shapes.final"
        "$SDS_LEGACY_CODE_ROOT/data/g_file.txt"
        "$SDS_LEGACY_CODE_ROOT/data/processed/sds_input/g_file.txt"
        "$SDS_LEGACY_CODE_ROOT/data/processed/g_file.txt"
        "$base_dir/data/g_file.txt"
        "$base_dir/data/processed/sds_input/g_file.txt"
        "$base_dir/data/processed/g_file.txt"
    )

    first_existing_file "${candidates[@]}"
}

find_population_vcf() {
    local pop="$1"
    local chr="$2"
    local -a candidates=(
        "$SDS_VCF_ROOT/$pop/UKBQC_${pop}_chr${chr}.vcf.gz"
        "$SDS_VCF_ROOT/$pop/UKBQC_${pop}_chr${chr}.phased.vcf.gz"
        "$SDS_VCF_ROOT/$pop/shapeit5/UKBQC_${pop}_chr${chr}.phased.vcf.gz"
        "$SDS_LEGACY_CODE_ROOT/data/vcf/$pop/UKBQC_${pop}_chr${chr}.vcf.gz"
        "$SDS_LEGACY_CODE_ROOT/plink/vcf_output/$pop/UKBQC_${pop}_chr${chr}.vcf.gz"
        "$pop/UKBQC_${pop}_chr${chr}.vcf.gz"
    )
    first_existing_file "${candidates[@]}"
}

find_population_sample_list() {
    local pop="$1"
    local -a candidates=(
        "$SDS_SAMPLE_LIST_ROOT/${pop}.txt"
        "$SDS_FREEZE_ROOT/${pop}.txt"
        "$SDS_FREEZE_ROOT/sample_lists/${pop}.txt"
        "$SDS_FREEZE_ROOT/cohorts/${pop}.txt"
        "$SDS_LEGACY_CODE_ROOT/data/${pop}.txt"
        "$SDS_LEGACY_CODE_ROOT/data/metadata/${pop}.txt"
    )
    first_existing_file "${candidates[@]}"
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
