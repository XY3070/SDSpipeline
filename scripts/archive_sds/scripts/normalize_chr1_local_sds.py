#!/usr/bin/env python3

import argparse
import csv
import math
from pathlib import Path

import numpy as np
from scipy.stats import norm


def parse_args():
    parser = argparse.ArgumentParser(
        description="Normalize chr1 SDS within local MAF bins and derive pseudo-p values."
    )
    parser.add_argument("--input", required=True, help="Input chr1.sds.tsv")
    parser.add_argument("--output-prefix", required=True, help="Prefix for output tables")
    parser.add_argument("--maf-threshold", type=float, default=0.01)
    parser.add_argument("--maf-bin-width", type=float, default=0.01)
    parser.add_argument(
        "--plot-p-threshold",
        type=float,
        default=None,
        help="Optional p-value threshold to flag points for plotting.",
    )
    return parser.parse_args()


def maf_from_daf(daf: float) -> float:
    return min(daf, 1.0 - daf)


def maf_bin_start(maf: float, maf_threshold: float, maf_bin_width: float):
    if not math.isfinite(maf) or maf <= maf_threshold:
        return None
    raw_idx = int((maf - maf_threshold) / maf_bin_width)
    bin_start = maf_threshold + raw_idx * maf_bin_width
    max_start = max(maf_threshold, 0.5 - maf_bin_width)
    return min(bin_start, max_start)


def maf_bin_id(bin_start: float, maf_bin_width: float) -> str:
    bin_end = min(bin_start + maf_bin_width, 0.5)
    return f"[{bin_start:.2f},{bin_end:.2f})"


def safe_float(text: str) -> float:
    try:
        value = float(text)
    except (TypeError, ValueError):
        return math.nan
    return value


def load_rows(path: Path):
    rows = []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            pos = int(float(row["POS"]))
            daf = safe_float(row["DAF"])
            rsds = safe_float(row["rSDS"])
            maf = maf_from_daf(daf) if math.isfinite(daf) else math.nan
            rows.append(
                {
                    "ID": row["ID"],
                    "AA": row["AA"],
                    "DA": row["DA"],
                    "POS": pos,
                    "DAF": daf,
                    "MAF": maf,
                    "nG0": row["nG0"],
                    "nG1": row["nG1"],
                    "nG2": row["nG2"],
                    "rSDS": rsds,
                    "SuggestedInitPoint": row["SuggestedInitPoint"],
                }
            )
    if not rows:
        raise SystemExit(f"No rows found in {path}")
    return rows


def compute_bin_stats(rows, maf_threshold: float, maf_bin_width: float):
    accumulators = {}
    common_variant_count = 0
    for row in rows:
        maf = row["MAF"]
        rsds = row["rSDS"]
        if not math.isfinite(maf) or not math.isfinite(rsds) or maf <= maf_threshold:
            continue
        bin_start = maf_bin_start(maf, maf_threshold, maf_bin_width)
        if bin_start is None:
            continue
        values = accumulators.setdefault(bin_start, [])
        values.append(rsds)
        common_variant_count += 1

    stats = {}
    for bin_start, values in sorted(accumulators.items()):
        arr = np.asarray(values, dtype=np.float64)
        mean_common = float(np.mean(arr))
        sd_common = float(np.std(arr))
        if not math.isfinite(sd_common) or sd_common <= 0.0:
            continue
        stats[bin_start] = {
            "bin_id": maf_bin_id(bin_start, maf_bin_width),
            "maf_start": bin_start,
            "maf_end": min(bin_start + maf_bin_width, 0.5),
            "common_variant_count": int(arr.size),
            "common_variant_mean": mean_common,
            "common_variant_sd": sd_common,
        }

    if not stats:
        raise SystemExit("No valid chr1 common-variant bins with positive SD were found")
    return stats, common_variant_count


def write_outputs(rows, output_prefix: Path, maf_threshold: float, maf_bin_width: float, plot_p_threshold: float | None):
    bin_stats, common_variant_count = compute_bin_stats(rows, maf_threshold, maf_bin_width)
    normalized_path = Path(f"{output_prefix}.normalized.tsv")
    bins_path = Path(f"{output_prefix}.frequency_bins.tsv")
    summary_path = Path(f"{output_prefix}.summary.tsv")

    common_rows = 0
    normalized_rows = 0
    min_p = 1.0
    max_neglog10 = 0.0
    plotted_rows = 0

    with normalized_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "ID",
                "AA",
                "DA",
                "POS",
                "DAF",
                "MAF",
                "nG0",
                "nG1",
                "nG2",
                "rSDS",
                "SuggestedInitPoint",
                "bin_id",
                "BIN_MEAN",
                "BIN_SD",
                "norm_SDS",
                "p_bothside",
                "neg_log10_p",
                "is_common_variant",
                "passes_plot_filter",
            ]
        )

        for row in rows:
            maf = row["MAF"]
            rsds = row["rSDS"]
            is_common_variant = int(math.isfinite(maf) and maf > maf_threshold)
            if is_common_variant:
                common_rows += 1

            bin_id = ""
            bin_mean = math.nan
            bin_sd = math.nan
            norm_sds = math.nan
            p_both = math.nan
            neg_log10_p = math.nan
            passes_plot = 0

            if is_common_variant and math.isfinite(rsds):
                bin_start = maf_bin_start(maf, maf_threshold, maf_bin_width)
                stats = bin_stats.get(bin_start)
                if stats is not None:
                    bin_id = stats["bin_id"]
                    bin_mean = stats["common_variant_mean"]
                    bin_sd = stats["common_variant_sd"]
                    norm_sds = (rsds - bin_mean) / bin_sd
                    p_left = norm.cdf(norm_sds)
                    p_right = norm.sf(norm_sds)
                    p_both = min(1.0, 2.0 * min(p_left, p_right))
                    neg_log10_p = -math.log10(max(p_both, 1e-300))
                    normalized_rows += 1
                    min_p = min(min_p, p_both)
                    max_neglog10 = max(max_neglog10, neg_log10_p)
                    if plot_p_threshold is None or p_both < plot_p_threshold:
                        passes_plot = 1
                        plotted_rows += 1

            writer.writerow(
                [
                    row["ID"],
                    row["AA"],
                    row["DA"],
                    row["POS"],
                    f"{row['DAF']:.10g}" if math.isfinite(row["DAF"]) else "",
                    f"{maf:.10g}" if math.isfinite(maf) else "",
                    row["nG0"],
                    row["nG1"],
                    row["nG2"],
                    f"{rsds:.10g}" if math.isfinite(rsds) else "",
                    row["SuggestedInitPoint"],
                    bin_id,
                    f"{bin_mean:.10g}" if math.isfinite(bin_mean) else "",
                    f"{bin_sd:.10g}" if math.isfinite(bin_sd) else "",
                    f"{norm_sds:.10g}" if math.isfinite(norm_sds) else "",
                    f"{p_both:.10g}" if math.isfinite(p_both) else "",
                    f"{neg_log10_p:.10g}" if math.isfinite(neg_log10_p) else "",
                    is_common_variant,
                    passes_plot,
                ]
            )

    with bins_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "bin_id",
                "maf_start",
                "maf_end",
                "common_variant_count",
                "common_variant_mean",
                "common_variant_sd",
            ]
        )
        for bin_start in sorted(bin_stats):
            stats = bin_stats[bin_start]
            writer.writerow(
                [
                    stats["bin_id"],
                    f"{stats['maf_start']:.10g}",
                    f"{stats['maf_end']:.10g}",
                    stats["common_variant_count"],
                    f"{stats['common_variant_mean']:.10g}",
                    f"{stats['common_variant_sd']:.10g}",
                ]
            )

    bonferroni_threshold = 0.05 / common_variant_count if common_variant_count > 0 else math.nan
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["key", "value"])
        writer.writerow(["input_rows", len(rows)])
        writer.writerow(["common_variant_rows", common_rows])
        writer.writerow(["normalized_rows", normalized_rows])
        writer.writerow(["plot_rows", plotted_rows])
        writer.writerow(["maf_threshold", f"{maf_threshold:.10g}"])
        writer.writerow(["maf_bin_width", f"{maf_bin_width:.10g}"])
        writer.writerow(["plot_p_threshold", "" if plot_p_threshold is None else f"{plot_p_threshold:.10g}"])
        writer.writerow(["bonferroni_threshold", f"{bonferroni_threshold:.10g}"])
        writer.writerow(["max_neg_log10_p", f"{max_neglog10:.10g}"])
        writer.writerow(["min_p_bothside", f"{min_p:.10g}" if min_p < 1.0 else ""])

    return normalized_path, bins_path, summary_path


def main():
    args = parse_args()
    rows = load_rows(Path(args.input))
    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    write_outputs(
        rows=rows,
        output_prefix=output_prefix,
        maf_threshold=args.maf_threshold,
        maf_bin_width=args.maf_bin_width,
        plot_p_threshold=args.plot_p_threshold,
    )


if __name__ == "__main__":
    main()
