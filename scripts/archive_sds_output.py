#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

import pyarrow as pa
import pyarrow.csv as pacsv
import pyarrow.parquet as pq


def parse_args():
    parser = argparse.ArgumentParser(description="Archive SDS TSV output to Parquet.")
    parser.add_argument("input_tsv")
    parser.add_argument("output_parquet")
    parser.add_argument("--summary-csv")
    return parser.parse_args()


def main():
    args = parse_args()
    input_path = Path(args.input_tsv)
    output_path = Path(args.output_parquet)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    table = pacsv.read_csv(
        input_path,
        read_options=pacsv.ReadOptions(use_threads=True),
        parse_options=pacsv.ParseOptions(delimiter="\t"),
        convert_options=pacsv.ConvertOptions(
            column_types=pa.schema(
                [
                    ("ID", pa.string()),
                    ("AA", pa.string()),
                    ("DA", pa.string()),
                    ("POS", pa.int64()),
                    ("DAF", pa.float64()),
                    ("nG0", pa.int64()),
                    ("nG1", pa.int64()),
                    ("nG2", pa.int64()),
                    ("rSDS", pa.float64()),
                    ("SuggestedInitPoint", pa.string()),
                ]
            )
        ),
    )
    pq.write_table(table, output_path)

    if args.summary_csv:
        summary_path = Path(args.summary_csv)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        with summary_path.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["input_tsv", "output_parquet", "rows", "columns"],
            )
            writer.writeheader()
            writer.writerow(
                {
                    "input_tsv": str(input_path),
                    "output_parquet": str(output_path),
                    "rows": table.num_rows,
                    "columns": table.num_columns,
                }
            )


if __name__ == "__main__":
    main()
