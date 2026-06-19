#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare multiple SDS diagnostic outputs produced by diagnose_sds_scan.py."
    )
    parser.add_argument(
        "--track",
        action="append",
        required=True,
        help=(
            "Track specification in the form label=/abs/path/prefix where prefix is the "
            "same prefix passed to diagnose_sds_scan.py."
        ),
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Prefix for merged comparison outputs.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="Number of top windows / hits to keep per track in merged tables.",
    )
    return parser.parse_args()


def load_kv(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["key"]: row["value"] for row in reader}


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_track_spec(spec: str) -> tuple[str, Path]:
    if "=" not in spec:
        raise SystemExit(f"Invalid --track value: {spec}")
    label, prefix = spec.split("=", 1)
    label = label.strip()
    prefix = prefix.strip()
    if not label or not prefix:
        raise SystemExit(f"Invalid --track value: {spec}")
    return label, Path(prefix).resolve()


def float_key(row: dict[str, str], key: str) -> float:
    try:
        return float(row.get(key, "") or "nan")
    except ValueError:
        return float("nan")


def int_key(row: dict[str, str], key: str) -> int:
    try:
        return int(float(row.get(key, "") or "0"))
    except ValueError:
        return 0


def write_table(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    args = parse_args()
    output_prefix = Path(args.output_prefix).resolve()
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, object]] = []
    tail_rows: list[dict[str, object]] = []
    peak_rows: list[dict[str, object]] = []
    hit_rows: list[dict[str, object]] = []
    markdown_lines = [
        "# SDS diagnostic comparison",
        "",
        "| Track | Common variants | chr1 >=20 | chr1 >=50 | chr1 >=100 | Top chr1 window |",
        "| --- | ---: | ---: | ---: | ---: | --- |",
    ]

    for spec in args.track:
        label, prefix = parse_track_spec(spec)
        summary = load_kv(prefix.with_name(prefix.name + ".summary.tsv"))
        tail_counts = load_rows(prefix.with_name(prefix.name + ".chrom_tail_counts.tsv"))
        peak_windows = load_rows(prefix.with_name(prefix.name + ".peak_windows.tsv"))
        top_hits = load_rows(prefix.with_name(prefix.name + ".top_hits.tsv"))

        summary_row = {"track": label}
        summary_row.update(summary)
        summary_rows.append(summary_row)

        chr1_row = None
        for row in tail_counts:
            row = dict(row)
            row["track"] = label
            tail_rows.append(row)
            if row.get("chr") == "1":
                chr1_row = row

        peak_windows_sorted = sorted(
            peak_windows,
            key=lambda row: (
                -int_key(row, "hit_count"),
                -float_key(row, "max_neg_log10_p"),
                row.get("chr", ""),
                int_key(row, "window_start"),
            ),
        )
        for rank, row in enumerate(peak_windows_sorted[: args.top_n], start=1):
            merged = dict(row)
            merged["track"] = label
            merged["track_rank"] = rank
            peak_rows.append(merged)

        top_hits_sorted = sorted(
            top_hits,
            key=lambda row: (
                int_key(row, "rank"),
                row.get("chr", ""),
                int_key(row, "pos"),
            ),
        )
        for row in top_hits_sorted[: args.top_n]:
            merged = dict(row)
            merged["track"] = label
            hit_rows.append(merged)

        top_chr1_window = ""
        for row in peak_windows_sorted:
            if row.get("chr") == "1":
                top_chr1_window = (
                    f"{row.get('window_start','')}-{row.get('window_end','')}"
                    f" hits={row.get('hit_count','')} max={row.get('max_neg_log10_p','')}"
                )
                break

        chr1_ge20 = chr1_row.get("count_neglog10_ge_20", "") if chr1_row else ""
        chr1_ge50 = chr1_row.get("count_neglog10_ge_50", "") if chr1_row else ""
        chr1_ge100 = chr1_row.get("count_neglog10_ge_100", "") if chr1_row else ""
        markdown_lines.append(
            f"| {label} | {summary.get('common_variant_count','')} | {chr1_ge20} | {chr1_ge50} | {chr1_ge100} | {top_chr1_window} |"
        )

    write_table(
        output_prefix.with_name(output_prefix.name + ".track_summary.tsv"),
        ["track", "common_variant_count", "neglog10_thresholds", "window_bp", "cluster_threshold"],
        summary_rows,
    )

    tail_fieldnames = ["track"]
    if tail_rows:
        tail_fieldnames.extend(key for key in tail_rows[0].keys() if key != "track")
    write_table(output_prefix.with_name(output_prefix.name + ".chrom_tail_compare.tsv"), tail_fieldnames, tail_rows)

    peak_fieldnames = ["track", "track_rank"]
    if peak_rows:
        peak_fieldnames.extend(key for key in peak_rows[0].keys() if key not in {"track", "track_rank"})
    write_table(output_prefix.with_name(output_prefix.name + ".top_windows.tsv"), peak_fieldnames, peak_rows)

    hit_fieldnames = ["track"]
    if hit_rows:
        hit_fieldnames.extend(key for key in hit_rows[0].keys() if key != "track")
    write_table(output_prefix.with_name(output_prefix.name + ".top_hits.tsv"), hit_fieldnames, hit_rows)

    md_path = output_prefix.with_name(output_prefix.name + ".summary.md")
    md_path.write_text("\n".join(markdown_lines) + "\n")

    print(output_prefix.with_name(output_prefix.name + ".track_summary.tsv"))
    print(output_prefix.with_name(output_prefix.name + ".chrom_tail_compare.tsv"))
    print(output_prefix.with_name(output_prefix.name + ".top_windows.tsv"))
    print(output_prefix.with_name(output_prefix.name + ".top_hits.tsv"))
    print(md_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
