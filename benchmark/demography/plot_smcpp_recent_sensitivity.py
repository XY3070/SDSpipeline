#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot recent-time SMC++ sensitivity overlays from a resolution report."
    )
    parser.add_argument("--report", required=True, help="Path to *_recent_resolution_report.json.")
    parser.add_argument(
        "--zoom-max",
        type=float,
        default=5000.0,
        help="Right bound for the zoomed recent-time plot. Default: 5000 generations.",
    )
    return parser.parse_args()


def load_smcpp_curve(path: Path) -> tuple[np.ndarray, np.ndarray]:
    xs: list[float] = []
    ys: list[float] = []
    with path.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            xs.append(float(row["x"]))
            ys.append(float(row["y"]))
    return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)


def configure_axes(ax: plt.Axes, title: str, zoom_max: float | None = None) -> None:
    ax.set_xscale("symlog", linthresh=10.0, linscale=1.0, base=10)
    ax.set_yscale("log")
    ax.set_xlim(left=0.0)
    if zoom_max is not None:
        ax.set_xlim(0.0, zoom_max)
    ax.set_xlabel("Generations ago")
    ax.set_ylabel("Effective population size ($N_e$)")
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.25)


def main() -> int:
    args = parse_args()
    report_path = Path(args.report).resolve()
    payload = json.loads(report_path.read_text())
    output_dir = report_path.parent
    pop = payload["population"]

    colors = [
        "black",
        "steelblue",
        "darkorange",
        "forestgreen",
        "purple",
    ]

    series = []
    for idx, run in enumerate(payload["runs"]):
        x, y = load_smcpp_curve(Path(run["csv_path"]))
        series.append((run["label"], x, y, colors[idx % len(colors)]))

    overlay_path = output_dir / f"{pop}_recent_resolution_overlay.png"
    zoom_path = output_dir / f"{pop}_recent_resolution_zoom_0_{int(args.zoom_max)}.png"

    for output_path, zoom_max, title_suffix in [
        (overlay_path, None, "full range"),
        (zoom_path, args.zoom_max, f"zoom 0-{int(args.zoom_max)} generations"),
    ]:
        fig, ax = plt.subplots(figsize=(10, 6))
        for label, x, y, color in series:
            ax.step(x, y, where="post", linewidth=2.2, color=color, label=label)
        configure_axes(ax, f"{pop} recent-resolution sensitivity ({title_suffix})", zoom_max)
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(output_path, dpi=300)
        plt.close(fig)

    manifest = {
        "population": pop,
        "report": str(report_path),
        "overlay_path": str(overlay_path),
        "zoom_path": str(zoom_path),
        "zoom_max": args.zoom_max,
    }
    (output_dir / f"{pop}_recent_resolution_plots.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
