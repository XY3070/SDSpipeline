#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

import pyarrow as pa
import pyarrow.feather as feather


def parse_args():
    parser = argparse.ArgumentParser(description="Filter SDS t_file rows to Feather.")
    parser.add_argument("input_tsv")
    parser.add_argument("output_feather")
    parser.add_argument("--start", type=float, required=True)
    parser.add_argument("--end", type=float, required=True)
    parser.add_argument("--summary-csv")
    return parser.parse_args()


def main():
    args = parse_args()
    input_path = Path(args.input_tsv)
    output_path = Path(args.output_feather)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    seen_in_range = False
    with input_path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            try:
                pos = float(parts[3])
            except ValueError:
                continue
            if pos < args.start or pos > args.end:
                if seen_in_range and pos > args.end:
                    break
                continue
            seen_in_range = True
            rows.append(
                {
                    "ID": parts[0],
                    "AA": parts[1],
                    "DA": parts[2],
                    "POS": pos,
                    "GENOTYPES": "\t".join(parts[4:]),
                }
            )

    table = pa.Table.from_pylist(
        rows,
        schema=pa.schema(
            [
                ("ID", pa.string()),
                ("AA", pa.string()),
                ("DA", pa.string()),
                ("POS", pa.float64()),
                ("GENOTYPES", pa.string()),
            ]
        ),
    )
    feather.write_feather(table, output_path)

    if args.summary_csv:
        summary_path = Path(args.summary_csv)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        with summary_path.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["input_tsv", "output_feather", "start", "end", "rows"],
            )
            writer.writeheader()
            writer.writerow(
                {
                    "input_tsv": str(input_path),
                    "output_feather": str(output_path),
                    "start": int(args.start),
                    "end": int(args.end),
                    "rows": len(rows),
                }
            )


if __name__ == "__main__":
    main()
