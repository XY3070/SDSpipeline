#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# shellcheck source=common_env.sh
source "$SCRIPT_DIR/common_env.sh"

SMCPP_IMAGE="${SMCPP_IMAGE:-$SDS_SMCPP_IMAGE}"
SMCPP_ROOTFS="${SMCPP_ROOTFS:-$SDS_SMCPP_ROOTFS}"
SMCPP_SINGULARITY_BIN="${SMCPP_SINGULARITY_BIN:-$SDS_SMCPP_SINGULARITY_BIN}"
SMCPP_SINGULARITY_CONFDIR="${SMCPP_SINGULARITY_CONFDIR:-$SDS_SMCPP_SINGULARITY_CONFDIR}"
SMCPP_MPLCONFIGDIR="${SMCPP_MPLCONFIGDIR:-$SDS_SMCPP_MPLCONFIGDIR}"
SMCPP_SIF_PARTITION_ID="${SMCPP_SIF_PARTITION_ID:-4}"
SMCPP_TMP_ROOT="${SDS_TMP_ROOT:-/tmp}"

resolve_runtime_candidate() {
    local candidate="$1"
    local resolved=""
    if [[ "$candidate" == */* ]]; then
        [[ -x "$candidate" ]] || return 1
        printf '%s\n' "$candidate"
        return 0
    fi
    resolved="$(command -v "$candidate" 2>/dev/null || true)"
    [[ -n "$resolved" ]] || return 1
    printf '%s\n' "$resolved"
}

cleanup() {
    if [[ -n "${TMP_SQFS:-}" ]]; then
        rm -f "$TMP_SQFS"
    fi
    if [[ -n "${TMP_ROOTFS:-}" && -d "${TMP_ROOTFS:-}" ]]; then
        rm -rf "$TMP_ROOTFS"
    fi
    if [[ -n "${TMP_CONFDIR:-}" && -d "${TMP_CONFDIR:-}" ]]; then
        rm -rf "$TMP_CONFDIR"
    fi
}
trap cleanup EXIT

if [[ ! -f "$SMCPP_IMAGE" ]]; then
    echo "[smc++] image not found: $SMCPP_IMAGE" >&2
    exit 1
fi

SMCPP_SINGULARITY_BIN="$(resolve_runtime_candidate "$SMCPP_SINGULARITY_BIN" || true)"
if [[ -z "$SMCPP_SINGULARITY_BIN" ]]; then
    echo "[smc++] singularity sif helper not found: $SMCPP_SINGULARITY_BIN" >&2
    exit 1
fi

if ! command -v unsquashfs >/dev/null 2>&1; then
    echo "[smc++] unsquashfs is required but was not found on PATH" >&2
    exit 1
fi

if [[ ! -d "$SMCPP_SINGULARITY_CONFDIR" ]]; then
    mkdir -p "$SMCPP_TMP_ROOT"
    TMP_CONFDIR="$(mktemp -d "$SMCPP_TMP_ROOT/singularity_conf_tmp.XXXXXX")"
    cp -a /etc/singularity/. "$TMP_CONFDIR/"
    mkdir -p "$(dirname "$SMCPP_SINGULARITY_CONFDIR")"
    mv "$TMP_CONFDIR" "$SMCPP_SINGULARITY_CONFDIR"
    unset TMP_CONFDIR
fi

if [[ ! -d "$SMCPP_ROOTFS" ]]; then
    mkdir -p "$SMCPP_TMP_ROOT"
    TMP_SQFS="$(mktemp "$SMCPP_TMP_ROOT/smcpp_rootfs.XXXXXX.squashfs")"
    TMP_ROOTFS="$(mktemp -d "$SMCPP_TMP_ROOT/smcpp_rootfs_tmp.XXXXXX")"
    "$SMCPP_SINGULARITY_BIN" sif dump "$SMCPP_SIF_PARTITION_ID" "$SMCPP_IMAGE" > "$TMP_SQFS"
    unsquashfs -f -d "$TMP_ROOTFS" "$TMP_SQFS"
    mv "$TMP_ROOTFS" "$SMCPP_ROOTFS"
    unset TMP_ROOTFS
fi

TEMPLATE_RC="$SMCPP_ROOTFS/usr/share/matplotlib/matplotlib.conf/matplotlibrc.template"
if [[ ! -f "$TEMPLATE_RC" ]]; then
    echo "[smc++] matplotlib template not found in rootfs: $TEMPLATE_RC" >&2
    exit 1
fi

mkdir -p "$SMCPP_MPLCONFIGDIR"
if [[ ! -f "$SMCPP_MPLCONFIGDIR/matplotlibrc" ]]; then
    cp "$TEMPLATE_RC" "$SMCPP_MPLCONFIGDIR/matplotlibrc"
fi

printf '[smc++] prepared singularity confdir=%s\n' "$SMCPP_SINGULARITY_CONFDIR"
printf '[smc++] prepared rootfs=%s\n' "$SMCPP_ROOTFS"
printf '[smc++] prepared matplotlib config=%s\n' "$SMCPP_MPLCONFIGDIR"
