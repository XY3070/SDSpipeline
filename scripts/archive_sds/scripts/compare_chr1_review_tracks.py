#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Combine chr1 SDS review summaries across tracks."
    )
    parser.add_argument(
        "--track",
        action="append",
        required=True,
        help="TRACK_LABEL=summary.tsv",
    )
    parser.add_argument("--output", required=True, help="Output TSV path")
    return parser.parse_args()


def load_summary(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    args = parse_args()
    rows: list[dict[str, str]] = []
    for spec in args.track:
        if "=" not in spec:
            raise SystemExit(f"Invalid --track spec: {spec}")
        label, path_str = spec.split("=", 1)
        path = Path(path_str)
        if not path.exists():
            raise SystemExit(f"Missing summary file: {path}")
        for row in load_summary(path):
            out = {"track": label}
            out.update(row)
            rows.append(out)

    fieldnames = [
        "track",
        "label",
        "start",
        "end",
        "row_count",
        "max_abs_rSDS",
        "top_id",
        "top_pos",
        "top_rSDS",
        "max_positive_rSDS",
        "min_negative_rSDS",
    ]
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
