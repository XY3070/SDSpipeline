#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path


def load_gamma(path: Path) -> list[tuple[float, float]]:
    rows: list[tuple[float, float]] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            daf_s, shape_s = line.split("\t")
            rows.append((float(daf_s), float(shape_s)))
    return rows


def load_piece_dir(piece_dir: Path) -> list[tuple[float, float]]:
    rows: list[tuple[float, float]] = []
    for piece in sorted(piece_dir.iterdir()):
        if not piece.is_file():
            continue
        if piece.name.startswith("."):
            continue
        rows.extend(load_gamma(piece))
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description="Combine low-frequency gamma pieces with a baseline gamma file.")
    parser.add_argument("--baseline-gfile", required=True)
    parser.add_argument("--low-piece-dir", required=True)
    parser.add_argument("--output-gfile", required=True)
    parser.add_argument("--drop-baseline-below", type=float, default=0.05)
    args = parser.parse_args()

    baseline_rows = [(daf, shape) for daf, shape in load_gamma(Path(args.baseline_gfile)) if daf >= args.drop_baseline_below - 1e-12]
    low_rows = load_piece_dir(Path(args.low_piece_dir))

    merged: dict[float, float] = {}
    for daf, shape in baseline_rows:
        merged[round(daf, 2)] = shape
    for daf, shape in low_rows:
        merged[round(daf, 2)] = shape

    output = Path(args.output_gfile)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as handle:
        for daf in sorted(merged):
            handle.write(f"{daf:.2f}\t{merged[daf]:.16g}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
