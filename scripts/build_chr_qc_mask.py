#!/usr/bin/env python3
"""Build a small BED mask for controlled chr-level SDS QC tests."""

from __future__ import annotations

import argparse
from pathlib import Path


BUILTIN_GRCH38_MASKS = {
    "1": [
        ("chr1_centromere", 122_026_459, 125_184_587),
        # Conservative 1q12 / pericentromeric heterochromatin interval covering
        # the known chr1 125 Mb and 143 Mb artifact peaks in this project.
        ("chr1_1q12_pericentromeric_heterochromatin", 125_000_000, 145_000_000),
    ],
    "6": [
        ("chr6_centromere", 58_830_167, 62_134_745),
    ],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chr", required=True, help="Chromosome number, e.g. 1.")
    parser.add_argument("--output-bed", required=True, help="Output BED path.")
    parser.add_argument(
        "--extra-bed",
        action="append",
        default=[],
        help="Additional BED file(s) to append after built-in intervals.",
    )
    parser.add_argument(
        "--chrom-prefix",
        default="chr",
        help="Chromosome prefix for output intervals. Default: chr.",
    )
    return parser.parse_args()


def normalize_chrom(chrom: str) -> str:
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def main() -> int:
    args = parse_args()
    chrom = normalize_chrom(args.chr)
    out = Path(args.output_bed)
    out.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    out_chrom = f"{args.chrom_prefix}{chrom}" if args.chrom_prefix else chrom
    for name, start_1based, end_1based in BUILTIN_GRCH38_MASKS.get(chrom, []):
        rows.append((out_chrom, start_1based - 1, end_1based, name))

    for bed_path in args.extra_bed:
        with Path(bed_path).open() as handle:
            for raw in handle:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) < 3:
                    raise SystemExit(f"Invalid BED line in {bed_path}: {raw.rstrip()}")
                name = fields[3] if len(fields) > 3 else Path(bed_path).stem
                rows.append((fields[0], int(fields[1]), int(fields[2]), name))

    rows.sort(key=lambda row: (row[0], row[1], row[2], row[3]))
    with out.open("w") as handle:
        for row in rows:
            handle.write("\t".join(map(str, row)) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
