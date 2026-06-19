#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare two SDS frequency_bins.tsv tables and write per-bin deltas."
    )
    parser.add_argument("--label-a", required=True, help="Label for baseline track.")
    parser.add_argument("--label-b", required=True, help="Label for comparison track.")
    parser.add_argument("--bins-a", required=True, help="Path to baseline frequency_bins.tsv.")
    parser.add_argument("--bins-b", required=True, help="Path to comparison frequency_bins.tsv.")
    parser.add_argument("--output-tsv", required=True, help="Output TSV path.")
    return parser.parse_args()


def load_rows(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_bin: dict[str, dict[str, str]] = {}
    for row in rows:
        bin_id = (row.get("maf_bin") or row.get("bin_id") or "").strip()
        if not bin_id:
            continue
        by_bin[bin_id] = row
    return by_bin


def as_float(row: dict[str, str], *keys: str) -> float | None:
    for key in keys:
        value = row.get(key, "")
        if value is None or value == "":
            continue
        try:
            return float(value)
        except ValueError:
            continue
    return None


def first_value(row: dict[str, str], *keys: str) -> str:
    for key in keys:
        value = row.get(key, "")
        if value not in (None, ""):
            return value
    return ""


def fmt(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.10g}"


def main() -> int:
    args = parse_args()
    bins_a = load_rows(Path(args.bins_a).resolve())
    bins_b = load_rows(Path(args.bins_b).resolve())
    all_bins = sorted(set(bins_a) | set(bins_b), key=lambda x: (float(x.split("-")[0].strip("[(]")) if "-" in x else 999.0, x))

    out_path = Path(args.output_tsv).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as handle:
        fieldnames = [
            "maf_bin",
            f"{args.label_a}_count",
            f"{args.label_b}_count",
            f"{args.label_a}_mean_rSDS",
            f"{args.label_b}_mean_rSDS",
            "delta_mean_rSDS",
            f"{args.label_a}_sd_rSDS",
            f"{args.label_b}_sd_rSDS",
            "delta_sd_rSDS",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for bin_id in all_bins:
            row_a = bins_a.get(bin_id, {})
            row_b = bins_b.get(bin_id, {})
            count_a = first_value(row_a, "count", "common_variant_count")
            count_b = first_value(row_b, "count", "common_variant_count")
            mean_a = as_float(row_a, "mean_rSDS", "mean", "common_variant_mean")
            mean_b = as_float(row_b, "mean_rSDS", "mean", "common_variant_mean")
            sd_a = as_float(row_a, "sd_rSDS", "sd", "common_variant_sd")
            sd_b = as_float(row_b, "sd_rSDS", "sd", "common_variant_sd")
            delta_mean = None if mean_a is None or mean_b is None else mean_b - mean_a
            delta_sd = None if sd_a is None or sd_b is None else sd_b - sd_a
            writer.writerow(
                {
                    "maf_bin": bin_id,
                    f"{args.label_a}_count": count_a,
                    f"{args.label_b}_count": count_b,
                    f"{args.label_a}_mean_rSDS": fmt(mean_a),
                    f"{args.label_b}_mean_rSDS": fmt(mean_b),
                    "delta_mean_rSDS": fmt(delta_mean),
                    f"{args.label_a}_sd_rSDS": fmt(sd_a),
                    f"{args.label_b}_sd_rSDS": fmt(sd_b),
                    "delta_sd_rSDS": fmt(delta_sd),
                }
            )

    print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
