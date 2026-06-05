#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert one or more SMC++ step-curve CSVs into backward.py-compatible NPZ arrays."
    )
    parser.add_argument(
        "--curve",
        action="append",
        nargs=3,
        metavar=("POP", "CSV", "STAT"),
        required=True,
        help=(
            "One population curve to include. "
            "Example: --curve NCN <SDS_DEMOGRAPHY_ROOT>/NCN/smcpp/NCN_fine_smcpp.csv median"
        ),
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output NPZ path.",
    )
    return parser.parse_args()


def load_curve(csv_path: Path) -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_csv(csv_path)
    if not {"x", "y"}.issubset(df.columns):
        raise ValueError(f"{csv_path} is missing required columns x,y")
    x = df["x"].astype(float).to_numpy()
    y = df["y"].astype(float).to_numpy()
    if x.ndim != 1 or y.ndim != 1 or len(x) != len(y):
        raise ValueError(f"{csv_path} has malformed x/y arrays")
    order = np.argsort(x, kind="stable")
    return x[order], y[order]


def main() -> None:
    args = parse_args()
    payload: dict[str, np.ndarray] = {}
    for pop, csv_path_text, stat in args.curve:
        csv_path = Path(csv_path_text)
        x, y = load_curve(csv_path)
        payload[f"{pop}_t"] = x
        payload[f"{pop}_{stat}"] = y
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(out_path, **payload)
    print(out_path)


if __name__ == "__main__":
    main()
