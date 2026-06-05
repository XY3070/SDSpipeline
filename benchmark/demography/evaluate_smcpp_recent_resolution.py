#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


DEFAULT_ROOT = Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate recent-time SMC++ sensitivity runs for one population."
    )
    parser.add_argument("--pop", required=True, help="Population label, e.g. NCN or SCN.")
    parser.add_argument(
        "--root",
        default=str(DEFAULT_ROOT),
        help="Benchmark demography root directory. Defaults to SDS_DEMOGRAPHY_ROOT.",
    )
    parser.add_argument(
        "--coarse-base",
        default=None,
        help="Base name for the coarse SMC++ run. Defaults to <POP>.",
    )
    parser.add_argument(
        "--reference-base",
        default=None,
        help="Base name for the current fine SMC++ reference. Defaults to <POP>_fine.",
    )
    parser.add_argument(
        "--candidate-bases",
        nargs="+",
        required=True,
        help="One or more candidate SMC++ base names to compare against the reference.",
    )
    parser.add_argument(
        "--decision-base",
        default=None,
        help="Candidate base to use for pass/fail decisions. Defaults to the first candidate.",
    )
    parser.add_argument(
        "--consistency-base",
        default=None,
        help="Optional second candidate base used for recent-trend consistency checks.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory for report outputs. Defaults to SDS_DEMOGRAPHY_ROOT/<POP>/smcpp/recent_sensitivity.",
    )
    args = parser.parse_args()
    args.pop = args.pop.upper()
    if args.coarse_base is None:
        args.coarse_base = args.pop
    if args.reference_base is None:
        args.reference_base = f"{args.pop}_fine"
    if args.decision_base is None:
        args.decision_base = args.candidate_bases[0]
    if args.consistency_base is None and len(args.candidate_bases) > 1:
        args.consistency_base = args.candidate_bases[1]
    if args.output_dir is None:
        args.output_dir = str(DEFAULT_ROOT / args.pop / "smcpp" / "recent_sensitivity")
    return args


def load_smcpp_curve(path: Path) -> tuple[np.ndarray, np.ndarray]:
    xs: list[float] = []
    ys: list[float] = []
    with path.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            xs.append(float(row["x"]))
            ys.append(float(row["y"]))
    return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)


def find_segments(x: np.ndarray, y: np.ndarray) -> list[dict[str, float]]:
    if x.size == 0:
        return []
    segments: list[dict[str, float]] = []
    start = float(x[0])
    current = float(y[0])
    for idx in range(1, len(x)):
        if not math.isclose(float(y[idx]), current, rel_tol=0.0, abs_tol=1e-12):
            segments.append(
                {
                    "start_generation": start,
                    "end_generation": float(x[idx]),
                    "ne": current,
                }
            )
            start = float(x[idx])
            current = float(y[idx])
    segments.append(
        {
            "start_generation": start,
            "end_generation": float(x[-1]),
            "ne": current,
        }
    )
    return segments


def sign_label(value: float, tol: float = 1e-12) -> str:
    if value > tol:
        return "up"
    if value < -tol:
        return "down"
    return "flat"


def summarize_run(base: str, csv_path: Path) -> dict[str, object]:
    x, y = load_smcpp_curve(csv_path)
    segments = find_segments(x, y)
    change_points = [segment["start_generation"] for segment in segments[1:]]
    first_bin_ne = segments[0]["ne"] if len(segments) >= 1 else None
    second_bin_ne = segments[1]["ne"] if len(segments) >= 2 else None
    third_bin_ne = segments[2]["ne"] if len(segments) >= 3 else None
    ratio = None
    if first_bin_ne is not None and second_bin_ne not in (None, 0):
        ratio = float(first_bin_ne / second_bin_ne)
    three_bin_delta = None
    if first_bin_ne is not None and third_bin_ne is not None:
        three_bin_delta = float(third_bin_ne - first_bin_ne)
    elif first_bin_ne is not None and second_bin_ne is not None:
        three_bin_delta = float(second_bin_ne - first_bin_ne)
    return {
        "base": base,
        "label": base,
        "csv_path": str(csv_path),
        "n_points": int(x.size),
        "n_segments": int(len(segments)),
        "first_change_generation": float(change_points[0]) if change_points else None,
        "first_five_change_points": [float(v) for v in change_points[:5]],
        "first_bin_ne": float(first_bin_ne) if first_bin_ne is not None else None,
        "second_bin_ne": float(second_bin_ne) if second_bin_ne is not None else None,
        "third_bin_ne": float(third_bin_ne) if third_bin_ne is not None else None,
        "first_bin_to_second_bin_ratio": ratio,
        "recent_three_bin_direction": sign_label(three_bin_delta if three_bin_delta is not None else 0.0),
    }


def write_csv(output_path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "base",
        "label",
        "csv_path",
        "n_points",
        "n_segments",
        "first_change_generation",
        "first_bin_ne",
        "second_bin_ne",
        "third_bin_ne",
        "first_bin_to_second_bin_ratio",
        "recent_three_bin_direction",
        "first_five_change_points",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            materialized = dict(row)
            materialized["first_five_change_points"] = ",".join(
                f"{value:.10g}" for value in row["first_five_change_points"]
            )
            writer.writerow(materialized)


def main() -> int:
    args = parse_args()
    root = Path(args.root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    smcpp_dir = root / args.pop / "smcpp"
    bases = [args.coarse_base, args.reference_base, *args.candidate_bases]
    unique_bases: list[str] = []
    for base in bases:
        if base not in unique_bases:
            unique_bases.append(base)

    runs: list[dict[str, object]] = []
    run_by_base: dict[str, dict[str, object]] = {}
    for base in unique_bases:
        csv_path = smcpp_dir / f"{base}_smcpp.csv"
        if not csv_path.exists():
            raise FileNotFoundError(f"Required SMC++ CSV not found: {csv_path}")
        summary = summarize_run(base, csv_path)
        runs.append(summary)
        run_by_base[base] = summary

    coarse = run_by_base[args.coarse_base]
    reference = run_by_base[args.reference_base]
    decision = run_by_base[args.decision_base]
    consistency = run_by_base[args.consistency_base] if args.consistency_base else None

    decision_checks = {
        "criterion_a_first_change_le_50": (
            decision["first_change_generation"] is not None
            and float(decision["first_change_generation"]) <= 50.0
        ),
        "criterion_b_first_second_ratio_ge_0_70": (
            decision["first_bin_to_second_bin_ratio"] is not None
            and float(decision["first_bin_to_second_bin_ratio"]) >= 0.70
        ),
        "criterion_c_recent_trend_consistent": True,
    }
    if consistency is not None:
        direction_a = str(decision["recent_three_bin_direction"])
        direction_b = str(consistency["recent_three_bin_direction"])
        decision_checks["criterion_c_recent_trend_consistent"] = (
            direction_a == direction_b or "flat" in {direction_a, direction_b}
        )

    payload = {
        "population": args.pop,
        "smcpp_dir": str(smcpp_dir),
        "coarse_base": args.coarse_base,
        "reference_base": args.reference_base,
        "candidate_bases": args.candidate_bases,
        "decision_base": args.decision_base,
        "consistency_base": args.consistency_base,
        "runs": runs,
        "comparisons": {
            "coarse_first_change_generation": coarse["first_change_generation"],
            "reference_first_change_generation": reference["first_change_generation"],
            "decision_first_change_generation": decision["first_change_generation"],
        },
        "decision": {
            "checks": decision_checks,
            "decision_base": args.decision_base,
            "passes": all(decision_checks.values()),
        },
    }

    json_path = output_dir / f"{args.pop}_recent_resolution_report.json"
    csv_path = output_dir / f"{args.pop}_recent_resolution_report.csv"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    write_csv(csv_path, runs)
    print(json.dumps(payload, indent=2, sort_keys=True))
    print(json_path)
    print(csv_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
