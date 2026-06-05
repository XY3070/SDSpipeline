#!/usr/bin/env python3

from __future__ import annotations

import argparse
import bisect
import csv
import json
from pathlib import Path


DEFAULT_ROOT = Path(__file__).resolve().parent
DEFAULT_CHECKPOINTS = (10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize recent effective-population-size checkpoints across phlash, "
            "SMC++, and Relate for one or more populations."
        )
    )
    parser.add_argument(
        "--root",
        default=str(DEFAULT_ROOT),
        help="Benchmark demography root directory. Defaults to benchmark/demography.",
    )
    parser.add_argument(
        "--pops",
        default="NCN,SCN",
        help="Comma-separated populations to summarize, e.g. NCN or NCN,SCN.",
    )
    parser.add_argument(
        "--smcpp-suffix",
        default="_fine_smcpp.csv",
        help="Suffix for the SMC++ CSV under <root>/<POP>/smcpp.",
    )
    parser.add_argument(
        "--relate-tag",
        default="relate_recent",
        help="Tag used in <POP>_<TAG>_ne.csv under <root>/<POP>/relate_recent/popsize.",
    )
    parser.add_argument(
        "--checkpoints",
        default=",".join(str(value) for value in DEFAULT_CHECKPOINTS),
        help="Comma-separated generation checkpoints to evaluate.",
    )
    parser.add_argument(
        "--bottleneck-min",
        type=float,
        default=500.0,
        help="Lower bound of the bottleneck search window in generations.",
    )
    parser.add_argument(
        "--bottleneck-max",
        type=float,
        default=5000.0,
        help="Upper bound of the bottleneck search window in generations.",
    )
    parser.add_argument(
        "--output-prefix",
        default="Population_History_Recent_Checkpoints",
        help="Prefix for TSV/JSON outputs written under --root.",
    )
    return parser.parse_args()


def parse_pops(text: str) -> list[str]:
    pops = [item.strip().upper() for item in text.split(",") if item.strip()]
    if not pops:
        raise ValueError("No populations specified")
    return pops


def parse_checkpoints(text: str) -> list[float]:
    values: list[float] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        value = float(part)
        if value < 0.0:
            raise ValueError("Checkpoints must be non-negative")
        values.append(value)
    if not values:
        raise ValueError("At least one checkpoint is required")
    return sorted(values)


def load_phlash_curves(summary_path: Path, pop: str) -> tuple[list[float], list[float]]:
    payload = json.loads(summary_path.read_text())
    if pop not in payload:
        raise KeyError(f"Population {pop} not found in {summary_path}")
    return list(payload[pop]["time"]), list(payload[pop]["ne"])


def load_curve_from_csv(path: Path, x_key: str, y_key: str) -> tuple[list[float], list[float]]:
    xs: list[float] = []
    ys: list[float] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            xs.append(float(row[x_key]))
            ys.append(float(row[y_key]))
    if not xs:
        raise ValueError(f"No rows found in {path}")
    return xs, ys


def step_value(xs: list[float], ys: list[float], generation: float) -> float:
    idx = bisect.bisect_right(xs, generation) - 1
    if idx < 0:
        idx = 0
    return ys[idx]


def bottleneck_summary(
    xs: list[float], ys: list[float], window_start: float, window_end: float
) -> tuple[float, float]:
    candidates = [window_start]
    candidates.extend(x for x in xs if window_start <= x <= window_end)
    values = [(generation, step_value(xs, ys, generation)) for generation in candidates]
    best_generation, best_ne = min(values, key=lambda item: (item[1], item[0]))
    return best_generation, best_ne


def write_checkpoint_tsv(
    path: Path,
    populations: list[str],
    checkpoints: list[float],
    curves: dict[str, dict[str, tuple[list[float], list[float]]]],
) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "population",
                "generation",
                "phlash_ne",
                "smcpp_ne",
                "relate_ne",
                "smcpp_over_phlash",
                "relate_over_phlash",
                "relate_over_smcpp",
            ]
        )
        for pop in populations:
            phlash_x, phlash_y = curves[pop]["phlash"]
            smcpp_x, smcpp_y = curves[pop]["smcpp"]
            relate_x, relate_y = curves[pop]["relate"]
            for generation in checkpoints:
                phlash_ne = step_value(phlash_x, phlash_y, generation)
                smcpp_ne = step_value(smcpp_x, smcpp_y, generation)
                relate_ne = step_value(relate_x, relate_y, generation)
                writer.writerow(
                    [
                        pop,
                        f"{generation:.10g}",
                        f"{phlash_ne:.10g}",
                        f"{smcpp_ne:.10g}",
                        f"{relate_ne:.10g}",
                        f"{(smcpp_ne / phlash_ne):.10g}",
                        f"{(relate_ne / phlash_ne):.10g}",
                        f"{(relate_ne / smcpp_ne):.10g}",
                    ]
                )


def write_bottleneck_tsv(
    path: Path,
    populations: list[str],
    curves: dict[str, dict[str, tuple[list[float], list[float]]]],
    window_start: float,
    window_end: float,
) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "population",
                "method",
                "window_start_generation",
                "window_end_generation",
                "bottleneck_generation",
                "bottleneck_ne",
            ]
        )
        for pop in populations:
            for method in ("phlash", "smcpp", "relate"):
                xs, ys = curves[pop][method]
                generation, ne = bottleneck_summary(xs, ys, window_start, window_end)
                writer.writerow(
                    [
                        pop,
                        method,
                        f"{window_start:.10g}",
                        f"{window_end:.10g}",
                        f"{generation:.10g}",
                        f"{ne:.10g}",
                    ]
                )


def write_manifest(
    path: Path,
    root: Path,
    populations: list[str],
    checkpoints: list[float],
    smcpp_suffix: str,
    relate_tag: str,
    bottleneck_min: float,
    bottleneck_max: float,
) -> None:
    payload = {
        "root": str(root),
        "populations": populations,
        "checkpoints": checkpoints,
        "smcpp_suffix": smcpp_suffix,
        "relate_tag": relate_tag,
        "bottleneck_window_generation": [bottleneck_min, bottleneck_max],
        "sources": {
            "phlash_summary": str(root / "Population_History_Models.json"),
        },
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def main() -> int:
    args = parse_args()
    root = Path(args.root).resolve()
    pops = parse_pops(args.pops)
    checkpoints = parse_checkpoints(args.checkpoints)

    curves: dict[str, dict[str, tuple[list[float], list[float]]]] = {}
    summary_path = root / "Population_History_Models.json"
    for pop in pops:
        pop_root = root / pop
        curves[pop] = {
            "phlash": load_phlash_curves(summary_path, pop),
            "smcpp": load_curve_from_csv(
                pop_root / "smcpp" / f"{pop}{args.smcpp_suffix}",
                "x",
                "y",
            ),
            "relate": load_curve_from_csv(
                pop_root / "relate_recent" / "popsize" / f"{pop}_{args.relate_tag}_ne.csv",
                "generation",
                "ne",
            ),
        }

    output_prefix = root / args.output_prefix
    checkpoint_tsv = output_prefix.with_suffix(".tsv")
    bottleneck_tsv = output_prefix.with_name(output_prefix.name + "_bottleneck.tsv")
    manifest_json = output_prefix.with_name(output_prefix.name + "_manifest.json")

    write_checkpoint_tsv(checkpoint_tsv, pops, checkpoints, curves)
    write_bottleneck_tsv(
        bottleneck_tsv,
        pops,
        curves,
        args.bottleneck_min,
        args.bottleneck_max,
    )
    write_manifest(
        manifest_json,
        root,
        pops,
        checkpoints,
        args.smcpp_suffix,
        args.relate_tag,
        args.bottleneck_min,
        args.bottleneck_max,
    )

    print(checkpoint_tsv)
    print(bottleneck_tsv)
    print(manifest_json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
