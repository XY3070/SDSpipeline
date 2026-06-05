#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import heapq
import math
import os
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Diagnose whether a genome-wide SDS scan is driven by localized peaks "
            "or by a suspicious excess of extreme common-variant hits."
        )
    )
    parser.add_argument(
        "--input-normalized-tsv",
        required=True,
        help="Path to the normalized.tsv file from postprocess_sds_results.py.",
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Prefix for the TSV/PNG outputs.",
    )
    parser.add_argument(
        "--neglog10-thresholds",
        default="8,20,50,100",
        help="Comma-separated -log10(p) cutoffs used for per-chromosome summaries.",
    )
    parser.add_argument(
        "--cluster-threshold",
        type=float,
        default=20.0,
        help="Minimum -log10(p) used when counting dense peak windows.",
    )
    parser.add_argument(
        "--window-bp",
        type=int,
        default=100000,
        help="Window size in base pairs for dense-peak summaries.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=200,
        help="Number of top common-variant hits to keep.",
    )
    return parser.parse_args()


def parse_thresholds(text: str) -> list[float]:
    thresholds: list[float] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        value = float(part)
        if value <= 0.0:
            raise ValueError("Thresholds must be positive")
        thresholds.append(value)
    if not thresholds:
        raise ValueError("At least one threshold is required")
    return sorted(thresholds)


def normalize_chrom(chrom: str) -> str:
    chrom = chrom.strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    return chrom


def chrom_sort_key(chrom: str) -> tuple[int, str]:
    try:
        return (0, f"{int(chrom):02d}")
    except ValueError:
        return (1, chrom)


def format_float(value: float) -> str:
    if not math.isfinite(value):
        return ""
    return f"{value:.10g}"


def push_top_hit(heap, limit: int, payload: dict[str, object]) -> None:
    key = float(payload["neg_log10_p"])
    item = (key, str(payload["chr"]), int(payload["pos"]), payload)
    if len(heap) < limit:
        heapq.heappush(heap, item)
        return
    if item > heap[0]:
        heapq.heapreplace(heap, item)


def write_summary(summary_path: Path, total_common: int, thresholds: list[float], window_bp: int, cluster_threshold: float) -> None:
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["key", "value"])
        writer.writerow(["common_variant_count", total_common])
        writer.writerow(["neglog10_thresholds", ",".join(format_float(x) for x in thresholds)])
        writer.writerow(["window_bp", window_bp])
        writer.writerow(["cluster_threshold", format_float(cluster_threshold)])


def write_tail_counts(path: Path, thresholds: list[float], chrom_counts, chrom_tail_counts) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        header = ["chr", "common_variant_count"]
        header.extend(f"count_neglog10_ge_{format_float(value)}" for value in thresholds)
        header.extend(f"rate_neglog10_ge_{format_float(value)}" for value in thresholds)
        writer.writerow(header)
        for chrom in sorted(chrom_counts, key=chrom_sort_key):
            total = chrom_counts[chrom]
            row = [chrom, total]
            for threshold in thresholds:
                row.append(chrom_tail_counts[chrom][threshold])
            for threshold in thresholds:
                row.append(format_float(chrom_tail_counts[chrom][threshold] / total if total else math.nan))
            writer.writerow(row)


def write_top_hits(path: Path, heap_items) -> None:
    rows = [item[3] for item in sorted(heap_items, reverse=True)]
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "rank",
                "chr",
                "pos",
                "id",
                "maf",
                "daf",
                "norm_sds",
                "rSDS",
                "p_bothside",
                "neg_log10_p",
            ]
        )
        for rank, row in enumerate(rows, start=1):
            writer.writerow(
                [
                    rank,
                    row["chr"],
                    row["pos"],
                    row["id"],
                    format_float(float(row["maf"])),
                    format_float(float(row["daf"])),
                    format_float(float(row["norm_sds"])),
                    format_float(float(row["rSDS"])),
                    format_float(float(row["p_bothside"])),
                    format_float(float(row["neg_log10_p"])),
                ]
            )


def write_peak_windows(path: Path, peak_windows) -> None:
    rows = []
    for (chrom, window_index), stats in peak_windows.items():
        rows.append(
            {
                "chr": chrom,
                "window_start": window_index,
                "window_end": stats["window_end"],
                "hit_count": stats["hit_count"],
                "max_neg_log10_p": stats["max_neg_log10_p"],
                "top_pos": stats["top_pos"],
                "top_id": stats["top_id"],
            }
        )
    rows.sort(
        key=lambda row: (
            -int(row["hit_count"]),
            -float(row["max_neg_log10_p"]),
            chrom_sort_key(str(row["chr"])),
            int(row["window_start"]),
        )
    )
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "chr",
                "window_start",
                "window_end",
                "hit_count",
                "max_neg_log10_p",
                "top_pos",
                "top_id",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row["chr"],
                    row["window_start"],
                    row["window_end"],
                    row["hit_count"],
                    format_float(float(row["max_neg_log10_p"])),
                    row["top_pos"],
                    row["top_id"],
                ]
            )


def plot_tail_counts(path: Path, thresholds: list[float], chrom_tail_counts) -> None:
    chroms = sorted(chrom_tail_counts, key=chrom_sort_key)
    fig, axes = plt.subplots(len(thresholds), 1, figsize=(14, max(4.5, 2.6 * len(thresholds))), sharex=True)
    if len(thresholds) == 1:
        axes = [axes]

    x = np.arange(len(chroms))
    for ax, threshold in zip(axes, thresholds):
        values = [chrom_tail_counts[chrom][threshold] for chrom in chroms]
        ax.bar(x, values, color="#4E79A7", width=0.8)
        ax.set_ylabel(f">={threshold:g}")
        ax.grid(axis="y", linestyle=":", alpha=0.4)
        ax.set_title(f"Common-variant hits with -log10(p) >= {threshold:g}", fontsize=11)

    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels([f"chr{chrom}" for chrom in chroms], rotation=45, ha="right")
    fig.suptitle("NCN SDS extreme-hit counts by chromosome", fontsize=14, weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=300)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    thresholds = parse_thresholds(args.neglog10_thresholds)
    input_path = Path(args.input_normalized_tsv).resolve()
    output_prefix = Path(args.output_prefix).resolve()
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    chrom_counts: dict[str, int] = defaultdict(int)
    chrom_tail_counts: dict[str, dict[float, int]] = defaultdict(lambda: {value: 0 for value in thresholds})
    peak_windows: dict[tuple[str, int], dict[str, object]] = {}
    top_hits = []
    total_common = 0

    with input_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("is_common_variant") != "1":
                continue
            try:
                p_both = float(row["p_bothside"])
                pos = int(float(row["pos"]))
                maf = float(row["MAF"])
                daf = float(row["DAF"])
                norm_sds = float(row["norm_SDS"])
                rsds = float(row["rSDS"])
            except (KeyError, TypeError, ValueError):
                continue
            if not math.isfinite(p_both) or p_both <= 0.0:
                continue
            chrom = normalize_chrom(row["chr"])
            neg_log10_p = -math.log10(p_both)

            chrom_counts[chrom] += 1
            total_common += 1
            for threshold in thresholds:
                if neg_log10_p >= threshold:
                    chrom_tail_counts[chrom][threshold] += 1

            payload = {
                "chr": chrom,
                "pos": pos,
                "id": row.get("ID", ""),
                "maf": maf,
                "daf": daf,
                "norm_sds": norm_sds,
                "rSDS": rsds,
                "p_bothside": p_both,
                "neg_log10_p": neg_log10_p,
            }
            push_top_hit(top_hits, args.top_n, payload)

            if neg_log10_p >= args.cluster_threshold:
                window_start = ((pos - 1) // args.window_bp) * args.window_bp + 1
                window_end = window_start + args.window_bp - 1
                key = (chrom, window_start)
                stats = peak_windows.setdefault(
                    key,
                    {
                        "window_end": window_end,
                        "hit_count": 0,
                        "max_neg_log10_p": -math.inf,
                        "top_pos": pos,
                        "top_id": row.get("ID", ""),
                    },
                )
                stats["hit_count"] += 1
                if neg_log10_p > float(stats["max_neg_log10_p"]):
                    stats["max_neg_log10_p"] = neg_log10_p
                    stats["top_pos"] = pos
                    stats["top_id"] = row.get("ID", "")

    if total_common == 0:
        raise SystemExit("No common variants were parsed from the normalized TSV.")

    summary_path = output_prefix.with_name(output_prefix.name + ".summary.tsv")
    tail_counts_path = output_prefix.with_name(output_prefix.name + ".chrom_tail_counts.tsv")
    top_hits_path = output_prefix.with_name(output_prefix.name + ".top_hits.tsv")
    peak_windows_path = output_prefix.with_name(output_prefix.name + ".peak_windows.tsv")
    plot_path = output_prefix.with_name(output_prefix.name + ".chrom_tail_counts.png")

    write_summary(summary_path, total_common, thresholds, args.window_bp, args.cluster_threshold)
    write_tail_counts(tail_counts_path, thresholds, chrom_counts, chrom_tail_counts)
    write_top_hits(top_hits_path, top_hits)
    write_peak_windows(peak_windows_path, peak_windows)
    plot_tail_counts(plot_path, thresholds, chrom_tail_counts)

    print(summary_path)
    print(tail_counts_path)
    print(top_hits_path)
    print(peak_windows_path)
    print(plot_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
