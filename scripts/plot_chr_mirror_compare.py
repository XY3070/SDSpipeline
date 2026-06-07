#!/usr/bin/env python3
"""Plot a mirrored single-chromosome Manhattan comparison from plot-points TSVs."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--top-points-tsv", required=True)
    parser.add_argument("--bottom-points-tsv", required=True)
    parser.add_argument("--top-label", required=True)
    parser.add_argument("--bottom-label", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--title", default="chr1 normalized SDS Manhattan comparison")
    parser.add_argument("--chrom-label", default="chr1")
    parser.add_argument("--top-threshold", type=float, default=None)
    parser.add_argument("--bottom-threshold", type=float, default=None)
    parser.add_argument("--point-size", type=float, default=3.0)
    return parser.parse_args()


def load_points(path: Path) -> list[tuple[float, float, int]]:
    rows = []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            x = float(row["plot_x"])
            y = float(row["neg_log10_p"])
            sig = int(row["is_significant"])
            if not math.isfinite(x) or not math.isfinite(y):
                continue
            rows.append((x, y, sig))
    if not rows:
        raise SystemExit(f"No points loaded from {path}")
    return rows


def main() -> int:
    args = parse_args()
    top_rows = load_points(Path(args.top_points_tsv))
    bottom_rows = load_points(Path(args.bottom_points_tsv))

    xmax = max(max(x for x, _y, _sig in top_rows), max(x for x, _y, _sig in bottom_rows))
    ymax = max(max(y for _x, y, _sig in top_rows), max(y for _x, y, _sig in bottom_rows))
    ymax = max(5.0, math.ceil(ymax * 1.05))

    out_prefix = Path(args.output_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(14, 6), dpi=220)
    top_x = [x / 1e6 for x, _y, _sig in top_rows]
    top_y = [y for _x, y, _sig in top_rows]
    top_sig = [sig for _x, _y, sig in top_rows]
    bottom_x = [x / 1e6 for x, _y, _sig in bottom_rows]
    bottom_y = [-y for _x, y, _sig in bottom_rows]
    bottom_sig = [sig for _x, _y, sig in bottom_rows]

    top_colors = ["#b2182b" if sig else "#4e79a7" for sig in top_sig]
    bottom_colors = ["#2166ac" if sig else "#9c755f" for sig in bottom_sig]

    ax.scatter(top_x, top_y, s=args.point_size, c=top_colors, alpha=0.65, linewidths=0, rasterized=True)
    ax.scatter(bottom_x, bottom_y, s=args.point_size, c=bottom_colors, alpha=0.65, linewidths=0, rasterized=True)

    ax.axhline(0.0, color="black", linewidth=1.0)
    if args.top_threshold is not None and args.top_threshold > 0:
        y = -math.log10(args.top_threshold)
        ax.axhline(y, color="#b2182b", linestyle="--", linewidth=1.0)
    if args.bottom_threshold is not None and args.bottom_threshold > 0:
        y = -math.log10(args.bottom_threshold)
        ax.axhline(-y, color="#2166ac", linestyle="--", linewidth=1.0)

    ax.set_xlim(0, xmax / 1e6)
    ax.set_ylim(-ymax, ymax)
    ax.set_xlabel(f"{args.chrom_label} position (Mb)")
    ax.set_ylabel(r"$-\log_{10}(p)$")
    ax.set_title(args.title)

    yticks = ax.get_yticks()
    ax.yaxis.set_major_locator(FixedLocator(yticks))
    ax.set_yticklabels([f"{abs(v):.0f}" if abs(v) >= 1 else "0" for v in yticks])

    ax.text(0.01, 0.96, args.top_label, transform=ax.transAxes, ha="left", va="top", fontsize=11, weight="bold")
    ax.text(0.01, 0.04, args.bottom_label, transform=ax.transAxes, ha="left", va="bottom", fontsize=11, weight="bold")

    # Highlight the two artifact windows from the QC test.
    for start, end, label in [
        (125_099_999, 125_200_000, "peak125"),
        (143_199_999, 143_300_000, "peak143"),
    ]:
        left = start / 1e6
        right = end / 1e6
        ax.axvspan(left, right, color="#d9d9d9", alpha=0.25, linewidth=0)
        ax.text((left + right) / 2, ymax * 0.92, label, ha="center", va="top", fontsize=9)

    fig.tight_layout()
    fig.savefig(out_prefix.with_suffix(".png"), facecolor="white")
    fig.savefig(out_prefix.with_suffix(".pdf"), facecolor="white")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
