#!/usr/bin/env python3
"""Single-figure ancestry plot: all samples on ONE stacked bar, with Region
labels on the x-axis marking where each Region's block starts/ends.

Sort order (left -> right):
    South (Guangzhou -> Shenzhen -> Hong Kong)
    East  (Shanghai)
    North (Jinan -> JiNing)
    Unknown
Within each Region block: sort by component-0 proportion (monotonic gradient).

X-axis: one tick per Region, label = "Region (n=...)", placed at the block
center. Vertical separators between Region blocks.

Usage:
    python plot_ancestry_one_figure.py [--k 2] [--out path.png]
"""
import argparse, os, sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = "/data/home/grp-wangyf/xuyuan/SDSworkspace/runs/admixture_9k_v1_20260619"
Q_DIR = os.path.join(ROOT, "unsupervised")
PSAM = os.path.join(ROOT, "input/9k_ldpruned.psam")
METADATA_TSV = "/data/home/grp-wangyf/xuyuan/raw/9k/sample_qc/cohort_freeze.tsv"

PALETTE = ["#0072B2", "#E69F00", "#009E73", "#CC79A7",
           "#56B4E9", "#D55E00", "#F0E442", "#949494"]

REGION_ORDER = ["Guangzhou", "Shenzhen", "HongKong",
                "Shanghai", "Jinan", "JiNing", "-"]
REGION_LABEL = {
    "Guangzhou": "South·Guangzhou",
    "Shenzhen":  "South·Shenzhen",
    "HongKong":  "South·Hong Kong",
    "Shanghai":  "East·Shanghai",
    "Jinan":     "North·Jinan",
    "JiNing":    "North·JiNing",
    "-":         "Unknown",
}


def load_q(path):
    rows = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            try:
                float(parts[0]); start = 0
            except ValueError:
                start = 1
            rows.append([float(x) for x in parts[start:]])
    return np.asarray(rows, dtype=float)


def load_labels_in_q_order(tsv_path, psam_path):
    LOW_INFO = {"", "0", "China", "EastAsia", "East Asian Ancestry"}
    sample_to_region = {}
    with open(tsv_path) as f:
        next(f, None)
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            sid = parts[0].strip() if parts else ""
            if not sid:
                continue
            reg = parts[4].strip() if len(parts) >= 5 else "-"
            sample_to_region[sid] = reg if reg not in LOW_INFO else "-"
    out, sids = [], []
    with open(psam_path) as f:
        next(f, None)
        for line in f:
            if not line.strip():
                continue
            sid = line.split("\t", 1)[0].strip()
            out.append(sample_to_region.get(sid, "-"))
            sids.append(sid)
    return out, sids


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=2)
    ap.add_argument("--out",
                    default=os.path.join(ROOT, "figures/ancestry_one_figure_K2.png"))
    args = ap.parse_args()

    q = load_q(os.path.join(Q_DIR, f"9k.{args.k}.Q"))
    labels, sids = load_labels_in_q_order(METADATA_TSV, PSAM)
    n = min(q.shape[0], len(labels))
    q = q[:n]
    labels = labels[:n]
    sids = sids[:n]

    # Build the sorted sample order: Region blocks, intra-Region sort.
    by_reg = defaultdict(list)
    for i in range(n):
        by_reg[labels[i]].append(i)

    ordered_idx = []
    region_blocks = []  # list of (region_label, n_in_block, start_in_ordered)
    for reg in REGION_ORDER:
        if reg not in by_reg:
            continue
        idxs = by_reg[reg]
        # Sort by component 0 proportion → monotonic gradient across components
        idxs.sort(key=lambda i: q[i, 0])
        start = len(ordered_idx)
        ordered_idx.extend(idxs)
        region_blocks.append((reg, len(idxs), start))

    Q_sorted = q[ordered_idx]
    n_total = Q_sorted.shape[0]
    k = Q_sorted.shape[1]

    # Figure: one long stacked bar
    fig_w = max(16.0, n_total / 200)
    fig, ax = plt.subplots(figsize=(fig_w, 3.8))

    x = np.arange(n_total)
    bottoms = np.zeros(n_total)
    for c in range(k):
        ax.bar(x, Q_sorted[:, c], bottom=bottoms,
               color=PALETTE[c % len(PALETTE)],
               width=1.0, edgecolor="none", linewidth=0,
               align="center")
        bottoms += Q_sorted[:, c]

    # Region separators + x-axis ticks/labels
    tick_pos, tick_lab = [], []
    for reg, n_in_block, start in region_blocks:
        # separator line at the start of each block (except the first)
        if start > 0:
            ax.axvline(start - 0.5, color="black", lw=1.2, alpha=0.9,
                       ymin=0, ymax=1)
        center = start + (n_in_block - 1) / 2
        tick_pos.append(center)
        tick_lab.append(f"{REGION_LABEL.get(reg, reg)}\n(n={n_in_block})")

    ax.set_xlim(-0.5, n_total - 0.5)
    ax.set_ylim(0, 1.02)
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_lab, fontsize=9, fontweight="bold",
                       rotation=45, ha="right", rotation_mode="anchor")
    ax.tick_params(axis="x", length=6, width=1.2, pad=8)
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0", "0.25", "0.5", "0.75", "1"], fontsize=9)
    ax.set_ylabel("Ancestry proportion", fontsize=11)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.2)
    ax.spines["left"].set_linewidth(1.0)

    # Legend
    handles = [Patch(color=PALETTE[c], label=f"k{c+1}") for c in range(k)]
    ax.legend(handles=handles, loc="upper right", frameon=True,
              facecolor="white", edgecolor="#cccccc",
              ncol=k, handlelength=1.0, fontsize=10)

    ax.set_title(
        f"ADAMIXTURE K={args.k}, 9k Chinese cohort — {n_total} samples, "
        f"sorted by Region",
        fontsize=13, pad=10,
    )

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {args.out}  ({n_total} samples, K={args.k}, "
          f"{len(region_blocks)} regions)")


if __name__ == "__main__":
    main()
