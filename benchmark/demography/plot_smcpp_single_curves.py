#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


DEFAULT_ROOT = Path(__file__).resolve().parent
POPS = ("NCN", "SCN")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot standalone SMC++ demographic curves for each population."
    )
    parser.add_argument(
        "--root",
        default=str(DEFAULT_ROOT),
        help="Benchmark demography root directory. Defaults to benchmark/demography.",
    )
    parser.add_argument(
        "--smcpp-suffix",
        default="_fine_smcpp.csv",
        help=(
            "Suffix for SMC++ CSV files under each population smcpp directory. "
            "For example, '_fine_smcpp.csv' loads NCN_fine_smcpp.csv and SCN_fine_smcpp.csv."
        ),
    )
    parser.add_argument(
        "--output-suffix",
        default="_fine_smcpp_standalone.png",
        help=(
            "Suffix for output PNG files under each population smcpp directory. "
            "For example, '_fine_smcpp_standalone.png' writes NCN_fine_smcpp_standalone.png."
        ),
    )
    return parser.parse_args()


def load_smcpp_curve(smcpp_csv: Path) -> tuple[np.ndarray, np.ndarray]:
    xs: list[float] = []
    ys: list[float] = []
    with smcpp_csv.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            xs.append(float(row["x"]))
            ys.append(float(row["y"]))
    return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)


def plot_single(output_path: Path, pop: str, x: np.ndarray, y: np.ndarray) -> None:
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.step(x, y, where="post", color="steelblue", linewidth=2.5)
    ax.set_xscale("symlog", linthresh=10.0, linscale=1.0, base=10)
    ax.set_yscale("log")
    ax.set_xlim(left=0.0)
    ax.set_title(f"{pop} finer SMC++")
    ax.set_xlabel("Generations ago")
    ax.set_ylabel("Effective population size ($N_e$)")
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    root = Path(args.root).resolve()
    for pop in POPS:
        smcpp_dir = root / pop / "smcpp"
        x, y = load_smcpp_curve(smcpp_dir / f"{pop}{args.smcpp_suffix}")
        plot_single(smcpp_dir / f"{pop}{args.output_suffix}", pop, x, y)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
