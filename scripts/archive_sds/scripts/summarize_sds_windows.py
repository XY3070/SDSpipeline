#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Summarize SDS values in target windows.")
    parser.add_argument("--input", required=True, help="Input chrN.sds.tsv")
    parser.add_argument("--output-prefix", required=True, help="Prefix for .summary.tsv and .top_hits.tsv")
    parser.add_argument(
        "--window",
        action="append",
        default=[],
        help="Window in the form label=start-end",
    )
    return parser.parse_args()


def parse_window(spec):
    label, span = spec.split("=", 1)
    start, end = span.split("-", 1)
    return {"label": label, "start": int(start), "end": int(end)}


def main():
    args = parse_args()
    windows = [parse_window(spec) for spec in args.window]
    if not windows:
        raise SystemExit("At least one --window is required")

    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    summary_path = Path(f"{output_prefix}.summary.tsv")
    hits_path = Path(f"{output_prefix}.top_hits.tsv")

    summary = {
        window["label"]: {
            "label": window["label"],
            "start": window["start"],
            "end": window["end"],
            "row_count": 0,
            "max_abs_rSDS": None,
            "top_id": "",
            "top_pos": "",
            "top_rSDS": None,
            "max_positive_rSDS": None,
            "min_negative_rSDS": None,
        }
        for window in windows
    }
    top_hits = []

    with Path(args.input).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            pos = int(float(row["POS"]))
            rsds = float(row["rSDS"])
            abs_rsds = abs(rsds)
            for window in windows:
                if not (window["start"] <= pos <= window["end"]):
                    continue
                item = summary[window["label"]]
                item["row_count"] += 1
                if item["max_positive_rSDS"] is None or rsds > item["max_positive_rSDS"]:
                    item["max_positive_rSDS"] = rsds
                if item["min_negative_rSDS"] is None or rsds < item["min_negative_rSDS"]:
                    item["min_negative_rSDS"] = rsds
                if item["max_abs_rSDS"] is None or abs_rsds > item["max_abs_rSDS"]:
                    item["max_abs_rSDS"] = abs_rsds
                    item["top_id"] = row["ID"]
                    item["top_pos"] = pos
                    item["top_rSDS"] = rsds
                top_hits.append(
                    {
                        "window": window["label"],
                        "ID": row["ID"],
                        "POS": pos,
                        "rSDS": rsds,
                        "abs_rSDS": abs_rsds,
                        "DAF": row.get("DAF", ""),
                    }
                )

    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
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
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for window in windows:
            writer.writerow(summary[window["label"]])

    top_hits.sort(key=lambda row: (row["window"], -row["abs_rSDS"], row["POS"]))
    with hits_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["window", "ID", "POS", "rSDS", "abs_rSDS", "DAF"],
            delimiter="\t",
        )
        writer.writeheader()
        current_window = None
        kept = 0
        for row in top_hits:
            if row["window"] != current_window:
                current_window = row["window"]
                kept = 0
            if kept >= 10:
                continue
            writer.writerow(row)
            kept += 1


if __name__ == "__main__":
    main()
