#!/usr/bin/env python3

import argparse
import csv
import math
import os
import re
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Summarize how heavily the MAF-binned standardized SDS tails deviate "
            "from a standard normal distribution."
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
        help="Prefix for the PNG/PDF/TSV outputs.",
    )
    parser.add_argument(
        "--thresholds",
        default="3,4,5,6,8",
        help="Comma-separated |Z| cutoffs used for tail exceedance summaries.",
    )
    return parser.parse_args()


def parse_thresholds(raw: str):
    values = []
    for part in raw.split(","):
        part = part.strip()
        if not part:
            continue
        value = float(part)
        if value <= 0:
            raise ValueError("Thresholds must be positive.")
        values.append(value)
    if not values:
        raise ValueError("At least one threshold is required.")
    return sorted(values)


def parse_bin_start(bin_id: str):
    match = re.match(r"^\[([0-9.]+),([0-9.]+)\)$", bin_id)
    if not match:
        return math.inf
    return float(match.group(1))


def two_sided_normal_tail(z_threshold: float):
    return math.erfc(z_threshold / math.sqrt(2.0))


def format_numeric(value):
    if value is None or not math.isfinite(value):
        return ""
    return f"{value:.10g}"


def collect_tail_counts(normalized_tsv: Path, thresholds):
    per_bin = {}
    overall = {"n": 0, "counts": [0] * len(thresholds)}

    with normalized_tsv.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("is_common_variant") != "1":
                continue
            bin_id = row.get("bin_id", "")
            if not bin_id:
                continue
            try:
                z_value = float(row["norm_SDS"])
            except (KeyError, TypeError, ValueError):
                continue
            if not math.isfinite(z_value):
                continue

            stats = per_bin.setdefault(bin_id, {"n": 0, "counts": [0] * len(thresholds)})
            abs_z = abs(z_value)
            stats["n"] += 1
            overall["n"] += 1
            for idx, threshold in enumerate(thresholds):
                if abs_z > threshold:
                    stats["counts"][idx] += 1
                    overall["counts"][idx] += 1
                else:
                    break

    if not per_bin:
        raise SystemExit("No valid common-variant norm_SDS rows were found in the normalized TSV.")

    return per_bin, overall


def build_rows(per_bin, overall, thresholds):
    rows = []

    def emit_rows(bin_id, stats):
        total = stats["n"]
        for idx, threshold in enumerate(thresholds):
            observed_count = stats["counts"][idx]
            observed_rate = observed_count / total
            expected_rate = two_sided_normal_tail(threshold)
            expected_count = total * expected_rate
            enrichment = observed_rate / expected_rate if expected_rate > 0.0 else math.nan
            log10_enrichment = math.log10(enrichment) if enrichment > 0.0 else math.nan
            rows.append(
                {
                    "bin_id": bin_id,
                    "threshold": threshold,
                    "n": total,
                    "observed_count": observed_count,
                    "observed_rate": observed_rate,
                    "expected_count": expected_count,
                    "expected_rate": expected_rate,
                    "enrichment": enrichment,
                    "log10_enrichment": log10_enrichment,
                }
            )

    for bin_id in sorted(per_bin, key=parse_bin_start):
        emit_rows(bin_id, per_bin[bin_id])
    emit_rows("ALL", overall)
    return rows


def write_summary_tsv(rows, output_path: Path):
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "bin_id",
                "threshold",
                "n",
                "observed_count",
                "observed_rate",
                "expected_count",
                "expected_rate",
                "enrichment",
                "log10_enrichment",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row["bin_id"],
                    format_numeric(row["threshold"]),
                    row["n"],
                    row["observed_count"],
                    format_numeric(row["observed_rate"]),
                    format_numeric(row["expected_count"]),
                    format_numeric(row["expected_rate"]),
                    format_numeric(row["enrichment"]),
                    format_numeric(row["log10_enrichment"]),
                ]
            )


def plot_heatmap(rows, thresholds, output_prefix: Path):
    bin_rows = [row for row in rows if row["bin_id"] != "ALL"]
    bin_ids = sorted({row["bin_id"] for row in bin_rows}, key=parse_bin_start)
    threshold_to_idx = {threshold: idx for idx, threshold in enumerate(thresholds)}
    bin_to_idx = {bin_id: idx for idx, bin_id in enumerate(bin_ids)}

    matrix = np.full((len(bin_ids), len(thresholds)), np.nan, dtype=np.float64)
    observed_text = [["" for _ in thresholds] for _ in bin_ids]
    overall = {}

    for row in rows:
        if row["bin_id"] == "ALL":
            overall[row["threshold"]] = row
            continue
        i = bin_to_idx[row["bin_id"]]
        j = threshold_to_idx[row["threshold"]]
        if math.isfinite(row["log10_enrichment"]):
            matrix[i, j] = row["log10_enrichment"]
        observed_text[i][j] = f"{row['observed_count']}/{row['n']}"

    finite_values = matrix[np.isfinite(matrix)]
    vmax = max(1.0, math.ceil(float(finite_values.max()))) if finite_values.size else 1.0
    norm = colors.TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=vmax)
    cmap = plt.get_cmap("coolwarm").copy()
    cmap.set_bad("#f0f0f0")

    fig = plt.figure(figsize=(13, max(12, len(bin_ids) * 0.30)))
    gs = fig.add_gridspec(2, 2, width_ratios=[5.5, 2.5], height_ratios=[18, 2.0], wspace=0.35, hspace=0.14)
    ax = fig.add_subplot(gs[0, 0])
    im = ax.imshow(matrix, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")

    ax.set_xticks(np.arange(len(thresholds)))
    ax.set_xticklabels([f"|Z|>{threshold:g}" for threshold in thresholds], fontsize=11)
    ax.set_yticks(np.arange(len(bin_ids)))
    ax.set_yticklabels(bin_ids, fontsize=8)
    ax.set_xlabel("Tail threshold", fontsize=12)
    ax.set_ylabel("MAF bin", fontsize=12)
    ax.set_title(
        "NCN standardized SDS tail enrichment over a standard normal\n"
        "Color = log10(observed / expected Pr(|Z| > threshold))",
        fontsize=14,
        weight="bold",
    )

    colorbar = fig.colorbar(im, ax=ax, pad=0.02)
    colorbar.set_label("log10(observed / expected)", fontsize=11)

    overall_ax = fig.add_subplot(gs[0, 1])
    overall_obs = np.asarray([overall[threshold]["observed_rate"] for threshold in thresholds], dtype=np.float64)
    overall_exp = np.asarray([overall[threshold]["expected_rate"] for threshold in thresholds], dtype=np.float64)
    overall_ax.plot(thresholds, overall_obs, marker="o", linewidth=2.0, color="#D62728", label="Observed")
    overall_ax.plot(
        thresholds,
        overall_exp,
        marker="o",
        linewidth=2.0,
        linestyle="--",
        color="#4E79A7",
        label="Standard normal",
    )
    overall_ax.set_yscale("log")
    overall_ax.set_xticks(thresholds)
    overall_ax.set_xticklabels([f"|Z|>{threshold:g}" for threshold in thresholds], rotation=25, ha="right")
    overall_ax.set_ylabel("Pr(|Z| > threshold)", fontsize=11)
    overall_ax.set_title("Overall tail rates", fontsize=12, weight="bold")
    overall_ax.grid(axis="y", linestyle=":", linewidth=0.8, alpha=0.6)
    overall_ax.legend(frameon=False, fontsize=10, loc="upper right")

    summary_ax = fig.add_subplot(gs[1, :])
    summary_ax.axis("off")
    overall_lines = []
    for threshold in thresholds:
        row = overall[threshold]
        overall_lines.append(
            f"|Z|>{threshold:g}: {row['enrichment']:.2e}x"
        )
    summary_ax.text(
        0.0,
        0.7,
        "Overall tail rates across all common variants",
        fontsize=11,
        weight="bold",
        transform=summary_ax.transAxes,
    )
    summary_ax.text(
        0.0,
        0.15,
        " | ".join(overall_lines),
        fontsize=9,
        transform=summary_ax.transAxes,
    )

    png_path = Path(f"{output_prefix}.png")
    pdf_path = Path(f"{output_prefix}.pdf")
    fig.subplots_adjust(left=0.11, right=0.96, top=0.94, bottom=0.06)
    fig.savefig(png_path, dpi=400, facecolor="white", bbox_inches="tight")
    fig.savefig(pdf_path, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return png_path, pdf_path


def main():
    args = parse_args()
    thresholds = parse_thresholds(args.thresholds)
    normalized_tsv = Path(args.input_normalized_tsv).resolve()
    output_prefix = Path(args.output_prefix).resolve()
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    per_bin, overall = collect_tail_counts(normalized_tsv, thresholds)
    rows = build_rows(per_bin, overall, thresholds)

    tsv_path = Path(f"{output_prefix}.tsv")
    write_summary_tsv(rows, tsv_path)
    png_path, pdf_path = plot_heatmap(rows, thresholds, output_prefix)

    print(f"input\t{normalized_tsv}")
    print(f"bins\t{len(per_bin)}")
    print(f"thresholds\t{','.join(str(x) for x in thresholds)}")
    print(f"tsv\t{tsv_path}")
    print(f"png\t{png_path}")
    print(f"pdf\t{pdf_path}")


if __name__ == "__main__":
    main()
