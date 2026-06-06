#!/usr/bin/env python3
"""Summarize raw SDS output in named BED windows."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sds-tsv", required=True)
    parser.add_argument("--windows-bed", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--label", default="")
    return parser.parse_args()


def normalize_chrom(chrom: str) -> str:
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def load_windows(path: Path) -> list[tuple[str, int, int, str]]:
    windows = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            name = fields[3] if len(fields) > 3 else f"{fields[0]}:{fields[1]}-{fields[2]}"
            windows.append((normalize_chrom(fields[0]), int(fields[1]) + 1, int(fields[2]), name))
    return windows


def main() -> int:
    args = parse_args()
    windows = load_windows(Path(args.windows_bed))
    stats = {name: [] for _chrom, _start, _end, name in windows}
    with Path(args.sds_tsv).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            pos = int(row["POS"])
            chrom = None
            if row.get("ID", "").startswith("chr"):
                chrom = normalize_chrom(row["ID"].split(":", 1)[0])
            value = float(row["rSDS"])
            for w_chrom, start, end, name in windows:
                if chrom is not None and chrom != w_chrom:
                    continue
                if start <= pos <= end:
                    stats[name].append(value)

    out = Path(args.output_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["label", "window", "row_count", "mean_rSDS", "mean_abs_rSDS", "max_abs_rSDS", "ge_abs_5", "ge_abs_10"])
        for _chrom, _start, _end, name in windows:
            values = stats[name]
            if values:
                abs_values = [abs(v) for v in values]
                writer.writerow(
                    [
                        args.label,
                        name,
                        len(values),
                        f"{sum(values) / len(values):.6f}",
                        f"{sum(abs_values) / len(abs_values):.6f}",
                        f"{max(abs_values):.6f}",
                        sum(1 for v in abs_values if v >= 5),
                        sum(1 for v in abs_values if v >= 10),
                    ]
                )
            else:
                writer.writerow([args.label, name, 0, "NA", "NA", "NA", 0, 0])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
