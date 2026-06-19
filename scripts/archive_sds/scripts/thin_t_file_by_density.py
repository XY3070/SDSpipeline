#!/usr/bin/env python3

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Row:
    pos: int
    text: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Thin an SDS t_file to an approximately uniform per-bin test-SNP density."
    )
    parser.add_argument("input_t_file", help="Path to the source SDS t_file.")
    parser.add_argument("output_t_file", help="Path to the thinned SDS t_file.")
    parser.add_argument("--bin-size", type=int, default=1_000_000, help="Bin size in bp. Default: 1,000,000.")
    parser.add_argument(
        "--target-per-bin",
        type=int,
        default=1_000,
        help="Maximum number of rows to keep per bin. Default: 1,000.",
    )
    return parser.parse_args()


def pick_evenly(rows: list[Row], target: int) -> list[Row]:
    if len(rows) <= target:
        return rows
    if target <= 1:
        return [rows[len(rows) // 2]]

    chosen: list[Row] = []
    step = (len(rows) - 1) / float(target - 1)
    seen: set[int] = set()
    for idx in range(target):
        pick = round(idx * step)
        if pick in seen:
            continue
        chosen.append(rows[pick])
        seen.add(pick)
    if len(chosen) < target:
        for idx, row in enumerate(rows):
            if idx in seen:
                continue
            chosen.append(row)
            if len(chosen) == target:
                break
    return chosen


def iter_bins(path: Path, bin_size: int):
    current_bin: int | None = None
    current_rows: list[Row] = []
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t", 4)
            if len(parts) < 5:
                continue
            pos = int(parts[3])
            bin_id = (pos - 1) // bin_size
            row = Row(pos=pos, text=line)
            if current_bin is None:
                current_bin = bin_id
            if bin_id != current_bin:
                yield current_bin, current_rows
                current_bin = bin_id
                current_rows = [row]
            else:
                current_rows.append(row)
    if current_bin is not None:
        yield current_bin, current_rows


def main() -> int:
    args = parse_args()
    input_path = Path(args.input_t_file).resolve()
    output_path = Path(args.output_t_file).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    kept_total = 0
    bin_count = 0
    with output_path.open("w") as out_handle:
        for _, rows in iter_bins(input_path, args.bin_size):
            bin_count += 1
            kept = pick_evenly(rows, args.target_per_bin)
            for row in kept:
                out_handle.write(row.text)
            kept_total += len(kept)

    print(f"input={input_path}")
    print(f"output={output_path}")
    print(f"bins={bin_count}")
    print(f"target_per_bin={args.target_per_bin}")
    print(f"rows_kept={kept_total}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
