#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert a Relate .coal file into a simple Ne(t) CSV, manifest, and plot."
    )
    parser.add_argument("--coal", required=True, help="Input .coal file from EstimatePopulationSize.sh.")
    parser.add_argument("--population", required=True, help="Population label for output metadata.")
    parser.add_argument("--output-csv", required=True, help="Output CSV path.")
    parser.add_argument("--output-json", required=True, help="Output manifest JSON path.")
    parser.add_argument("--output-png", required=True, help="Output PNG path.")
    parser.add_argument(
        "--group-index",
        type=int,
        default=0,
        help="Population index to extract from the .coal file. Defaults to 0.",
    )
    parser.add_argument(
        "--years-per-gen",
        type=float,
        default=28.0,
        help="Years per generation used for the derived years column.",
    )
    return parser.parse_args()


def parse_coal(path: Path, group_index: int) -> tuple[list[str], np.ndarray, np.ndarray]:
    lines = [line.strip().split() for line in path.read_text().splitlines() if line.strip()]
    if len(lines) < 3:
        raise ValueError(f"Malformed .coal file: {path}")
    populations = lines[0]
    epochs = np.asarray([float(value) for value in lines[1]], dtype=float)
    data_rows = lines[2:]

    match = None
    for row in data_rows:
        if len(row) != epochs.size + 2:
            continue
        if int(row[0]) == group_index and int(row[1]) == group_index:
            match = row
            break
    if match is None:
        match = data_rows[0]
    rates = np.asarray([float(value) for value in match[2:]], dtype=float)
    if rates.size != epochs.size:
        raise ValueError(f"Epoch/rate length mismatch in {path}")
    return populations, epochs, rates


def write_csv(path: Path, epochs: np.ndarray, rates: np.ndarray, years_per_gen: float) -> np.ndarray:
    ne = 0.5 / rates
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["generation", "years", "coal_rate", "ne"])
        for generation, rate, ne_value in zip(epochs, rates, ne):
            writer.writerow(
                [
                    f"{generation:.10g}",
                    f"{generation * years_per_gen:.10g}",
                    f"{rate:.10g}",
                    f"{ne_value:.10g}",
                ]
            )
    return ne


def plot_ne(path: Path, pop: str, epochs: np.ndarray, ne: np.ndarray) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.step(epochs, ne, where="post", color="darkgreen", linewidth=2.0)
    ax.set_xscale("symlog", linthresh=10.0, linscale=1.0, base=10)
    ax.set_yscale("log")
    ax.set_xlim(left=0.0)
    ax.set_xlabel("Generations ago")
    ax.set_ylabel("Effective population size ($N_e$)")
    ax.set_title(f"{pop}: Relate recent $N_e(t)$")
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def write_manifest(
    path: Path,
    coal_path: Path,
    pop: str,
    populations: list[str],
    epochs: np.ndarray,
    rates: np.ndarray,
    ne: np.ndarray,
    years_per_gen: float,
    group_index: int,
) -> None:
    payload = {
        "coal_path": str(coal_path),
        "population": pop,
        "group_index": group_index,
        "population_labels": populations,
        "years_per_gen": years_per_gen,
        "epochs": int(epochs.size),
        "recent_ne": float(ne[0]),
        "max_ne": float(ne.max()),
        "min_ne": float(ne.min()),
        "min_coal_rate": float(rates.min()),
        "max_coal_rate": float(rates.max()),
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def main() -> int:
    args = parse_args()
    coal_path = Path(args.coal).resolve()
    output_csv = Path(args.output_csv).resolve()
    output_json = Path(args.output_json).resolve()
    output_png = Path(args.output_png).resolve()

    if not coal_path.exists():
        raise FileNotFoundError(f".coal file not found: {coal_path}")

    populations, epochs, rates = parse_coal(coal_path, args.group_index)
    ne = write_csv(output_csv, epochs, rates, args.years_per_gen)
    plot_ne(output_png, args.population, epochs, ne)
    write_manifest(
        output_json,
        coal_path,
        args.population,
        populations,
        epochs,
        rates,
        ne,
        args.years_per_gen,
        args.group_index,
    )

    print(output_csv)
    print(output_json)
    print(output_png)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover - CLI error path
        print(f"[Error] {exc}", file=sys.stderr)
        raise
