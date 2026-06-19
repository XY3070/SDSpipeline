#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

import pyarrow.feather as feather


def parse_args():
    parser = argparse.ArgumentParser(description="Split a Feather table into fixed-size row chunks.")
    parser.add_argument("input_feather")
    parser.add_argument("output_dir")
    parser.add_argument("--chunk-rows", type=int, required=True, help="Maximum rows per chunk.")
    parser.add_argument("--prefix", default="chunk", help="Chunk filename prefix.")
    parser.add_argument("--summary-csv", help="Optional summary TSV/CSV path.")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.chunk_rows <= 0:
        raise SystemExit("--chunk-rows must be > 0")

    input_path = Path(args.input_feather)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    table = feather.read_table(input_path)
    total_rows = table.num_rows
    chunk_count = 0
    rows_written = 0

    for start in range(0, total_rows, args.chunk_rows):
        end = min(start + args.chunk_rows, total_rows)
        chunk_path = output_dir / f"{args.prefix}_{chunk_count:04d}.feather"
        feather.write_feather(table.slice(start, end - start), chunk_path)
        chunk_count += 1
        rows_written += end - start

    if args.summary_csv:
        summary_path = Path(args.summary_csv)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        with summary_path.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["input_feather", "output_dir", "chunk_rows", "chunk_count", "rows_written"],
            )
            writer.writeheader()
            writer.writerow(
                {
                    "input_feather": str(input_path),
                    "output_dir": str(output_dir),
                    "chunk_rows": args.chunk_rows,
                    "chunk_count": chunk_count,
                    "rows_written": rows_written,
                }
            )


if __name__ == "__main__":
    main()
