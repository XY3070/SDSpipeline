#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import pickle
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


DEFAULT_ROOT = Path(__file__).resolve().parent
DEFAULT_MU = 1.25e-8


class _CompatSizeHistory:
    def __new__(cls, t, c):
        obj = object.__new__(cls)
        obj.t = np.asarray(t, dtype=float)
        obj.c = np.asarray(c, dtype=float)
        return obj


class _CompatDemographicModel:
    def __new__(cls, *args):
        obj = object.__new__(cls)
        obj._args = args
        return obj


def _reconstruct_jax_array(reconstruct_func, reconstruct_args, state, attrs):
    arr = reconstruct_func(*reconstruct_args)
    arr.__setstate__(state)
    return np.asarray(arr, dtype=float)


class _CompatPhlashUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module == "phlash.size_history" and name == "SizeHistory":
            return _CompatSizeHistory
        if module == "phlash.size_history" and name == "DemographicModel":
            return _CompatDemographicModel
        if module == "jax._src.array" and name == "_reconstruct_array":
            return _reconstruct_jax_array
        return super().find_class(module, name)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare phlash, fine SMC++, and Relate recent Ne curves on one figure."
    )
    parser.add_argument(
        "--root",
        default=str(DEFAULT_ROOT),
        help="Benchmark demography root directory. Defaults to SDS_DEMOGRAPHY_ROOT.",
    )
    parser.add_argument(
        "--pops",
        default="NCN,SCN",
        help="Comma-separated populations to plot, e.g. NCN or NCN,SCN.",
    )
    parser.add_argument(
        "--mu",
        type=float,
        default=DEFAULT_MU,
        help="Mutation rate used to rescale phlash posterior samples.",
    )
    parser.add_argument(
        "--smcpp-suffix",
        default="_fine_smcpp.csv",
        help="Suffix for SMC++ CSV files under each population smcpp directory.",
    )
    parser.add_argument(
        "--relate-tag",
        default="relate_recent",
        help="Tag used in <POP>_<TAG>_ne.csv under each population relate_recent/popsize directory.",
    )
    parser.add_argument(
        "--output-prefix",
        default="Population_History_Recent_Method_Comparison",
        help="Prefix for output PNG/CSV/manifest files written under the root directory.",
    )
    return parser.parse_args()


def parse_pops(text: str) -> list[str]:
    pops = [item.strip().upper() for item in text.split(",") if item.strip()]
    if not pops:
        raise ValueError("No populations specified")
    return pops


def load_summary(summary_path: Path) -> dict[str, dict[str, list[float]]]:
    return json.loads(summary_path.read_text())


def load_smcpp_curve(smcpp_csv: Path) -> tuple[np.ndarray, np.ndarray]:
    xs: list[float] = []
    ys: list[float] = []
    with smcpp_csv.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            xs.append(float(row["x"]))
            ys.append(float(row["y"]))
    return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)


def load_relate_curve(relate_csv: Path) -> tuple[np.ndarray, np.ndarray]:
    xs: list[float] = []
    ys: list[float] = []
    with relate_csv.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            xs.append(float(row["generation"]))
            ys.append(float(row["ne"]))
    return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)


def load_phlash_posterior(
    pkl_path: Path, target_generation_grid: np.ndarray, mu: float
) -> dict[str, np.ndarray]:
    with pkl_path.open("rb") as handle:
        models = _CompatPhlashUnpickler(handle).load()

    scaled_curves = []
    for model in models:
        eta = model._args[0]
        theta = float(np.asarray(model._args[1]).reshape(-1)[0])
        scale = theta / mu
        scaled_curves.append((1.0 / eta.c) * scale)

    posterior = np.asarray(scaled_curves, dtype=float)
    return {
        "generation": np.asarray(target_generation_grid, dtype=float),
        "median": np.median(posterior, axis=0),
        "lower": np.percentile(posterior, 2.5, axis=0),
        "upper": np.percentile(posterior, 97.5, axis=0),
        "num_models": np.asarray([posterior.shape[0]], dtype=int),
    }


def extend_step_curve(
    x_values: np.ndarray, y_values: np.ndarray, target_max_x: float
) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x_values, dtype=float)
    y = np.asarray(y_values, dtype=float)
    if x.size == 0:
        return x, y
    if target_max_x > x[-1]:
        x = np.append(x, float(target_max_x))
        y = np.append(y, float(y[-1]))
    return x, y


def write_combined_csv(
    output_path: Path,
    pops: list[str],
    phlash_curves: dict[str, dict[str, np.ndarray]],
    smcpp_curves: dict[str, tuple[np.ndarray, np.ndarray]],
    relate_curves: dict[str, tuple[np.ndarray, np.ndarray]],
) -> None:
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["population", "method", "generation", "ne", "ne_lower", "ne_upper"])
        for pop in pops:
            smcpp_x, smcpp_y = smcpp_curves[pop]
            for x, y in zip(smcpp_x, smcpp_y):
                writer.writerow([pop, "smcpp", f"{x:.10g}", f"{y:.10g}", "", ""])
            relate_x, relate_y = relate_curves[pop]
            for x, y in zip(relate_x, relate_y):
                writer.writerow([pop, "relate", f"{x:.10g}", f"{y:.10g}", "", ""])
            phlash = phlash_curves[pop]
            for x, med, lo, hi in zip(
                phlash["generation"], phlash["median"], phlash["lower"], phlash["upper"]
            ):
                writer.writerow(
                    [pop, "phlash", f"{x:.10g}", f"{med:.10g}", f"{lo:.10g}", f"{hi:.10g}"]
                )


def write_manifest(
    output_path: Path,
    root: Path,
    pops: list[str],
    phlash_curves: dict[str, dict[str, np.ndarray]],
    mu: float,
    smcpp_suffix: str,
    relate_tag: str,
) -> None:
    payload: dict[str, object] = {
        "root": str(root),
        "mu": mu,
        "pops": pops,
        "smcpp_suffix": smcpp_suffix,
        "relate_tag": relate_tag,
        "populations": {},
    }
    populations = payload["populations"]
    assert isinstance(populations, dict)
    for pop in pops:
        populations[pop] = {
            "phlash_num_models": int(phlash_curves[pop]["num_models"][0]),
        }
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def plot_comparison(
    output_path: Path,
    pops: list[str],
    phlash_curves: dict[str, dict[str, np.ndarray]],
    smcpp_curves: dict[str, tuple[np.ndarray, np.ndarray]],
    relate_curves: dict[str, tuple[np.ndarray, np.ndarray]],
) -> None:
    fig, axes = plt.subplots(1, len(pops), figsize=(7 * len(pops), 6), sharey=True)
    if len(pops) == 1:
        axes = [axes]

    for ax, pop in zip(axes, pops):
        smcpp_x, smcpp_y = smcpp_curves[pop]
        relate_x, relate_y = relate_curves[pop]
        phlash = phlash_curves[pop]
        phlash_x = np.asarray(phlash["generation"], dtype=float)
        smcpp_x, smcpp_y = extend_step_curve(smcpp_x, smcpp_y, float(phlash_x[-1]))
        relate_x, relate_y = extend_step_curve(relate_x, relate_y, float(phlash_x[-1]))

        ax.step(smcpp_x, smcpp_y, where="post", color="steelblue", linewidth=2.0, label="SMC++ fine")
        ax.step(relate_x, relate_y, where="post", color="darkgreen", linewidth=2.0, label="Relate")
        ax.fill_between(
            phlash_x,
            phlash["lower"],
            phlash["upper"],
            step="post",
            color="firebrick",
            alpha=0.18,
            label="phlash 95% CI",
        )
        ax.step(
            phlash_x,
            phlash["median"],
            where="post",
            color="firebrick",
            linewidth=2.0,
            label="phlash median",
        )
        ax.set_xscale("symlog", linthresh=10.0, linscale=1.0, base=10)
        ax.set_yscale("log")
        ax.set_xlim(left=0.0)
        ax.set_title(pop)
        ax.set_xlabel("Generations ago")
        ax.grid(True, which="both", alpha=0.25)

    axes[0].set_ylabel("Effective population size ($N_e$)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.suptitle("Recent Demography Benchmark: phlash vs SMC++ vs Relate", y=0.98)
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.945),
        ncol=4,
        frameon=False,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    root = Path(args.root).resolve()
    pops = parse_pops(args.pops)
    summary = load_summary(root / "Population_History_Models.json")

    phlash_curves: dict[str, dict[str, np.ndarray]] = {}
    smcpp_curves: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    relate_curves: dict[str, tuple[np.ndarray, np.ndarray]] = {}

    for pop in pops:
        pop_root = root / pop
        phlash_curves[pop] = load_phlash_posterior(
            pop_root / "phlash" / f"{pop}_model_full.pkl",
            np.asarray(summary[pop]["time"], dtype=float),
            args.mu,
        )
        smcpp_curves[pop] = load_smcpp_curve(pop_root / "smcpp" / f"{pop}{args.smcpp_suffix}")
        relate_curves[pop] = load_relate_curve(
            pop_root / "relate_recent" / "popsize" / f"{pop}_{args.relate_tag}_ne.csv"
        )

    output_png = root / f"{args.output_prefix}.png"
    output_csv = root / f"{args.output_prefix}.csv"
    output_manifest = root / f"{args.output_prefix}_manifest.json"

    plot_comparison(output_png, pops, phlash_curves, smcpp_curves, relate_curves)
    write_combined_csv(output_csv, pops, phlash_curves, smcpp_curves, relate_curves)
    write_manifest(
        output_manifest,
        root,
        pops,
        phlash_curves,
        args.mu,
        args.smcpp_suffix,
        args.relate_tag,
    )

    print(output_png)
    print(output_csv)
    print(output_manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
