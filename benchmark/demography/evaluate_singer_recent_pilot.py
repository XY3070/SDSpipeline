#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import glob
import json
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import tskit


INDEX_RE = re.compile(r"_(\d+)\.trees$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize SINGER pilot posterior tree sequences with lightweight trace and "
            "windowed mean-TMRCA diagnostics."
        )
    )
    parser.add_argument("--glob", required=True, help="Glob pattern for posterior .trees files.")
    parser.add_argument("--mutation-rate", type=float, required=True, help="Mutation rate.")
    parser.add_argument("--output-prefix", required=True, help="Output prefix for CSV/JSON/PNG.")
    parser.add_argument(
        "--window-size",
        type=float,
        default=1_000_000,
        help="Genomic window size for mean-TMRCA summaries. Default: 1e6.",
    )
    parser.add_argument(
        "--skip-incompatibility",
        action="store_true",
        help="Skip the incompatibility count trace if it is too slow.",
    )
    return parser.parse_args()


def parse_sample_index(path: Path) -> int:
    match = INDEX_RE.search(path.name)
    if not match:
        return -1
    return int(match.group(1))


def count_incompatibility(ts: tskit.TreeSequence) -> int:
    count = 0
    for tree in ts.trees():
        for site in tree.sites():
            num_mutations = len(site.mutations)
            if num_mutations > 2:
                count += 1
            elif num_mutations == 2:
                if site.mutations[0].node != tree.root and site.mutations[1].node != tree.root:
                    count += 1
    return count


def compute_windows(sequence_length: float, window_size: float) -> np.ndarray:
    windows = np.arange(0.0, sequence_length, float(window_size), dtype=float)
    if windows.size == 0 or windows[0] != 0.0:
        windows = np.insert(windows, 0, 0.0)
    if windows[-1] != float(sequence_length):
        windows = np.append(windows, float(sequence_length))
    return windows


def summarize_tree(
    path: Path,
    mutation_rate: float,
    windows: np.ndarray,
    skip_incompatibility: bool,
) -> tuple[dict[str, float], np.ndarray]:
    ts = tskit.load(path)
    site_div = np.asarray(ts.diversity(windows=windows, mode="site"), dtype=float)
    branch_div = np.asarray(ts.diversity(windows=windows, mode="branch"), dtype=float) * mutation_rate
    mean_tmrca_windows = np.asarray(ts.diversity(windows=windows, mode="branch"), dtype=float) / 2.0
    genome_mean_tmrca = float(np.asarray(ts.diversity(mode="branch")).reshape(-1)[0] / 2.0)
    incompatibility = -1 if skip_incompatibility else count_incompatibility(ts)
    metrics = {
        "tree_file": str(path),
        "sample_index": float(parse_sample_index(path)),
        "sequence_length": float(ts.sequence_length),
        "num_trees": float(ts.num_trees),
        "num_sites": float(ts.num_sites),
        "diversity_fit_mse": float(np.mean((site_div - branch_div) ** 2)),
        "incompatibility_sites": float(incompatibility),
        "genome_mean_tmrca": genome_mean_tmrca,
    }
    return metrics, mean_tmrca_windows


def write_metrics_csv(path: Path, rows: list[dict[str, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "tree_file",
        "sample_index",
        "sequence_length",
        "num_trees",
        "num_sites",
        "diversity_fit_mse",
        "incompatibility_sites",
        "genome_mean_tmrca",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_tmrca_csv(path: Path, windows: np.ndarray, tmrca_stack: np.ndarray) -> dict[str, float]:
    median = np.median(tmrca_stack, axis=0)
    lower = np.percentile(tmrca_stack, 2.5, axis=0)
    upper = np.percentile(tmrca_stack, 97.5, axis=0)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["window_start", "window_end", "tmrca_median", "tmrca_q025", "tmrca_q975"])
        for start, end, med, lo, hi in zip(windows[:-1], windows[1:], median, lower, upper):
            writer.writerow(
                [f"{start:.10g}", f"{end:.10g}", f"{med:.10g}", f"{lo:.10g}", f"{hi:.10g}"]
            )
    return {
        "genome_tmrca_median": float(np.median(median)),
        "genome_tmrca_q025": float(np.percentile(median, 2.5)),
        "genome_tmrca_q975": float(np.percentile(median, 97.5)),
    }


def plot_summary(
    path: Path,
    rows: list[dict[str, float]],
    windows: np.ndarray,
    tmrca_stack: np.ndarray,
    skip_incompatibility: bool,
) -> None:
    indices = [row["sample_index"] for row in rows]
    mean_tmrca = [row["genome_mean_tmrca"] for row in rows]
    mse = [row["diversity_fit_mse"] for row in rows]
    incompat = [row["incompatibility_sites"] for row in rows]
    median = np.median(tmrca_stack, axis=0)
    lower = np.percentile(tmrca_stack, 2.5, axis=0)
    upper = np.percentile(tmrca_stack, 97.5, axis=0)
    mids = (windows[:-1] + windows[1:]) / 2.0

    fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=False)

    axes[0].plot(indices, mean_tmrca, color="darkgreen", linewidth=1.8)
    axes[0].set_title("Posterior trace: genome-wide mean TMRCA")
    axes[0].set_xlabel("Posterior sample index")
    axes[0].set_ylabel("Mean TMRCA (generations)")
    axes[0].grid(True, alpha=0.25)

    axes[1].plot(indices, mse, color="firebrick", linewidth=1.6, label="diversity fit MSE")
    axes[1].set_title("Posterior trace: diversity-fit and incompatibility")
    axes[1].set_xlabel("Posterior sample index")
    axes[1].set_ylabel("Diversity-fit MSE")
    axes[1].grid(True, alpha=0.25)
    if not skip_incompatibility:
        twin = axes[1].twinx()
        twin.plot(indices, incompat, color="slateblue", linewidth=1.2, label="incompatibility sites")
        twin.set_ylabel("Incompatibility sites")

    axes[2].fill_between(mids, lower, upper, color="darkgreen", alpha=0.18, label="95% CI")
    axes[2].plot(mids, median, color="darkgreen", linewidth=1.8, label="median")
    axes[2].set_title("Windowed mean TMRCA across posterior")
    axes[2].set_xlabel("Genomic position")
    axes[2].set_ylabel("Mean TMRCA (generations)")
    axes[2].grid(True, alpha=0.25)
    axes[2].legend(frameon=False)

    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    tree_paths = sorted(Path(path).resolve() for path in glob.glob(args.glob))
    if not tree_paths:
        raise FileNotFoundError(f"No tree files matched glob: {args.glob}")

    first_ts = tskit.load(tree_paths[0])
    windows = compute_windows(first_ts.sequence_length, args.window_size)

    rows = []
    tmrca_curves = []
    for path in tree_paths:
        metrics, curve = summarize_tree(path, args.mutation_rate, windows, args.skip_incompatibility)
        rows.append(metrics)
        tmrca_curves.append(curve)

    tmrca_stack = np.asarray(tmrca_curves, dtype=float)
    output_prefix = Path(args.output_prefix).resolve()
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    metrics_csv = output_prefix.with_name(output_prefix.name + "_metrics.csv")
    tmrca_csv = output_prefix.with_name(output_prefix.name + "_tmrca_windows.csv")
    summary_json = output_prefix.with_name(output_prefix.name + "_summary.json")
    plot_png = output_prefix.with_name(output_prefix.name + "_summary.png")

    write_metrics_csv(metrics_csv, rows)
    tmrca_summary = write_tmrca_csv(tmrca_csv, windows, tmrca_stack)
    plot_summary(plot_png, rows, windows, tmrca_stack, args.skip_incompatibility)

    payload = {
        "tree_glob": args.glob,
        "num_tree_files": len(tree_paths),
        "mutation_rate": args.mutation_rate,
        "window_size": args.window_size,
        "skip_incompatibility": bool(args.skip_incompatibility),
        "metrics_csv": str(metrics_csv),
        "tmrca_windows_csv": str(tmrca_csv),
        "summary_png": str(plot_png),
        "mean_tmrca_median": float(np.median([row["genome_mean_tmrca"] for row in rows])),
        "mean_tmrca_q025": float(np.percentile([row["genome_mean_tmrca"] for row in rows], 2.5)),
        "mean_tmrca_q975": float(np.percentile([row["genome_mean_tmrca"] for row in rows], 97.5)),
        "diversity_fit_mse_median": float(np.median([row["diversity_fit_mse"] for row in rows])),
        **tmrca_summary,
    }
    summary_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    print(metrics_csv)
    print(tmrca_csv)
    print(summary_json)
    print(plot_png)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover - CLI error path
        print(f"[Error] {exc}", file=sys.stderr)
        raise
