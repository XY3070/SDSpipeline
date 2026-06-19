#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

SMCPP_IMAGE="${SMCPP_IMAGE:-$PROJECT_ROOT/smcpp_latest.sif}"
SMCPP_ROOTFS="${SMCPP_ROOTFS:-$PROJECT_ROOT/smcpp_latest.sandbox}"
SMCPP_SINGULARITY_BIN="${SMCPP_SINGULARITY_BIN:-$PROJECT_ROOT/bin/singularity}"
SMCPP_SINGULARITY_CONFDIR="${SMCPP_SINGULARITY_CONFDIR:-$PROJECT_ROOT/third_party/singularity_conf}"
SMCPP_MPLCONFIGDIR="${SMCPP_MPLCONFIGDIR:-$PROJECT_ROOT/third_party/smcpp_matplotlib}"
SMCPP_SIF_PARTITION_ID="${SMCPP_SIF_PARTITION_ID:-4}"

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

if [[ ! -x "$SMCPP_SINGULARITY_BIN" ]]; then
    echo "[smc++] singularity sif helper not found: $SMCPP_SINGULARITY_BIN" >&2
    exit 1
fi

if ! command -v unsquashfs >/dev/null 2>&1; then
    echo "[smc++] unsquashfs is required but was not found on PATH" >&2
    exit 1
fi

if [[ ! -d "$SMCPP_SINGULARITY_CONFDIR" ]]; then
    TMP_CONFDIR="$(mktemp -d "$PROJECT_ROOT/.singularity_conf_tmp.XXXXXX")"
    cp -a /etc/singularity/. "$TMP_CONFDIR/"
    mkdir -p "$(dirname "$SMCPP_SINGULARITY_CONFDIR")"
    mv "$TMP_CONFDIR" "$SMCPP_SINGULARITY_CONFDIR"
    unset TMP_CONFDIR
fi

if [[ ! -d "$SMCPP_ROOTFS" ]]; then
    TMP_SQFS="$(mktemp /tmp/smcpp_rootfs.XXXXXX.squashfs)"
    TMP_ROOTFS="$(mktemp -d "$PROJECT_ROOT/.smcpp_rootfs_tmp.XXXXXX")"
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
