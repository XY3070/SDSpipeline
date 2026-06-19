#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT_ROOT = PROJECT_ROOT / "benchmark" / "demography"


def safe_float(text: str) -> float | None:
    try:
        value = float(text)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(value):
        return None
    return value


def load_phlash_history(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    time_grid = []
    median = []
    lower = []
    upper = []
    with path.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            time_grid.append(float(row["time_gen"]))
            median.append(float(row["ne_median"]))
            lower.append(float(row["ne_q025"]))
            upper.append(float(row["ne_q975"]))
    if not time_grid:
        raise ValueError(f"No phlash history rows found in {path}")
    return (
        np.asarray(time_grid, dtype=float),
        np.asarray(median, dtype=float),
        np.asarray(lower, dtype=float),
        np.asarray(upper, dtype=float),
    )


def detect_delimiter(path: Path) -> str:
    sample = path.read_text()[:4096]
    try:
        return csv.Sniffer().sniff(sample).delimiter
    except csv.Error:
        return ","


def choose_column_index(header: list[str] | None, candidate_indices: list[int], role: str) -> int:
    if not candidate_indices:
        raise ValueError(f"No numeric columns available for {role}-axis selection")
    if header:
        preferred_tokens = {
            "x": ("time", "gen", "year", "x"),
            "y": ("ne", "size", "pop", "y"),
        }[role]
        ranked = []
        for index in candidate_indices:
            name = header[index].strip().lower()
            score = sum(token in name for token in preferred_tokens)
            ranked.append((score, -index, index))
        ranked.sort(reverse=True)
        if ranked[0][0] > 0:
            return ranked[0][2]
    return candidate_indices[0 if role == "x" else min(1, len(candidate_indices) - 1)]


def load_smcpp_curve(path: Path) -> tuple[np.ndarray, np.ndarray]:
    delimiter = detect_delimiter(path)
    with path.open() as handle:
        rows = list(csv.reader(handle, delimiter=delimiter))
    if not rows:
        raise ValueError(f"Empty SMC++ CSV: {path}")

    header = None
    data_rows = rows
    numeric_in_first = [index for index, value in enumerate(rows[0]) if safe_float(value) is not None]
    if len(numeric_in_first) < 2:
        header = rows[0]
        data_rows = rows[1:]
    if not data_rows:
        raise ValueError(f"No data rows found in SMC++ CSV: {path}")

    numeric_counts: dict[int, int] = {}
    for row in data_rows:
        for index, value in enumerate(row):
            if safe_float(value) is not None:
                numeric_counts[index] = numeric_counts.get(index, 0) + 1
    candidate_indices = sorted(index for index, count in numeric_counts.items() if count > 0)
    x_index = choose_column_index(header, candidate_indices, "x")
    y_index = choose_column_index(header, [index for index in candidate_indices if index != x_index], "y")

    xs = []
    ys = []
    for row in data_rows:
        if max(x_index, y_index) >= len(row):
            continue
        x_value = safe_float(row[x_index])
        y_value = safe_float(row[y_index])
        if x_value is None or y_value is None or x_value <= 0 or y_value <= 0:
            continue
        xs.append(x_value)
        ys.append(y_value)
    if not xs:
        raise ValueError(f"Could not parse positive numeric SMC++ curve values from {path}")
    return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)


def autodetect_smcpp_csv(pop_dir: Path, pop: str) -> Path:
    smcpp_dir = pop_dir / "smcpp"
    candidates = sorted(smcpp_dir.glob(f"{pop}*_smcpp.png.csv"))
    if not candidates:
        candidates = sorted(smcpp_dir.glob("*.png.csv"))
    if not candidates:
        candidates = sorted(smcpp_dir.glob("*.csv"))
    if not candidates:
        raise FileNotFoundError(f"Could not find an smc++ CSV under {smcpp_dir}")
    return candidates[0]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Overlay phlash and SMC++ demographic benchmark curves for one population."
    )
    parser.add_argument("--pop", required=True, help="Population label, e.g. NCN or SCN.")
    parser.add_argument(
        "--phlash-csv",
        default=None,
        help="phlash history CSV. Defaults to benchmark/demography/<POP>/phlash/<POP>_phlash_history.csv.",
    )
    parser.add_argument(
        "--smcpp-csv",
        default=None,
        help="SMC++ plot CSV. Defaults to auto-detecting under benchmark/demography/<POP>/smcpp.",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output PNG path. Defaults to benchmark/demography/<POP>/plots/<POP>_phlash_vs_smcpp.png.",
    )
    parser.add_argument("--title", default=None, help="Optional plot title override.")
    args = parser.parse_args()
    args.pop = args.pop.upper()
    pop_dir = DEFAULT_OUT_ROOT / args.pop
    if args.phlash_csv is None:
        args.phlash_csv = str(pop_dir / "phlash" / f"{args.pop}_phlash_history.csv")
    if args.smcpp_csv is None:
        args.smcpp_csv = str(autodetect_smcpp_csv(pop_dir, args.pop))
    if args.output is None:
        args.output = str(pop_dir / "plots" / f"{args.pop}_phlash_vs_smcpp.png")
    if args.title is None:
        args.title = f"{args.pop}: phlash vs SMC++"
    return args


def main() -> int:
    args = parse_args()
    phlash_csv = Path(args.phlash_csv).resolve()
    smcpp_csv = Path(args.smcpp_csv).resolve()
    output_path = Path(args.output).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    for path in [phlash_csv, smcpp_csv]:
        if not path.exists():
            raise FileNotFoundError(f"Required input not found: {path}")

    time_grid, median_ne, lower_ne, upper_ne = load_phlash_history(phlash_csv)
    smcpp_x, smcpp_y = load_smcpp_curve(smcpp_csv)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.fill_between(time_grid, lower_ne, upper_ne, color="firebrick", alpha=0.2, label="phlash 95% CI")
    ax.plot(time_grid, median_ne, color="firebrick", linewidth=2.0, label="phlash median")
    ax.plot(smcpp_x, smcpp_y, color="steelblue", linewidth=2.0, label="SMC++")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Generations ago")
    ax.set_ylabel("Effective population size ($N_e$)")
    ax.set_title(args.title)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)

    summary = {
        "population": args.pop,
        "phlash_csv": str(phlash_csv),
        "smcpp_csv": str(smcpp_csv),
        "output": str(output_path),
        "phlash_points": int(time_grid.size),
        "smcpp_points": int(smcpp_x.size),
        "title": args.title,
    }
    output_path.with_suffix(".json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover - CLI error path
        print(f"[Error] {exc}", file=sys.stderr)
        raise
