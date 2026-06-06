#!/usr/bin/env python3
"""Create a targeted SDS input by filtering only the t_file to BED intervals."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--chr", required=True)
    parser.add_argument("--target-bed", required=True)
    return parser.parse_args()


def normalize_chrom(chrom: str) -> str:
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def load_intervals(path: Path, chrom: str) -> list[tuple[int, int, str]]:
    target = normalize_chrom(chrom)
    intervals = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if normalize_chrom(fields[0]) != target:
                continue
            name = fields[3] if len(fields) > 3 else "target"
            intervals.append((int(fields[1]) + 1, int(fields[2]), name))
    intervals.sort()
    return intervals


def in_intervals(pos: int, intervals: list[tuple[int, int, str]]) -> bool:
    for start, end, _name in intervals:
        if pos < start:
            return False
        if start <= pos <= end:
            return True
    return False


def main() -> int:
    args = parse_args()
    chrom = normalize_chrom(args.chr)
    base = Path(args.base_input_dir)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    intervals = load_intervals(Path(args.target_bed), chrom)
    if not intervals:
        raise SystemExit(f"No target intervals for chr{chrom}")

    for suffix in ["b_file.txt", "s_file.txt", "o_file.txt", "sample_order.txt"]:
        src = base / f"chr{chrom}_{suffix}"
        if src.exists():
            shutil.copy2(src, out / src.name)

    rows_in = 0
    rows_out = 0
    with (base / f"chr{chrom}_t_file.txt").open() as inp, (out / f"chr{chrom}_t_file.txt").open("w") as out_t:
        for raw in inp:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            rows_in += 1
            if in_intervals(int(fields[3]), intervals):
                out_t.write(raw)
                rows_out += 1

    qc_dir = out / "qc"
    qc_dir.mkdir(exist_ok=True)
    with (qc_dir / "target_t_file_summary.tsv").open("w") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"t_rows_in\t{rows_in}\n")
        handle.write(f"t_rows_out\t{rows_out}\n")
        handle.write(f"target_interval_count\t{len(intervals)}\n")
    if rows_out == 0:
        raise SystemExit("Targeted t_file has zero rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
