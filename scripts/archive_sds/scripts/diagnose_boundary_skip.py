#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from compute_SDS import FileReader, PickleCache, SingletonData, iter_test_snps


def parse_scan_window(spec: str) -> tuple[str, int, int]:
    label, span = spec.split("=", 1)
    start_str, end_str = span.split("-", 1)
    return label, int(start_str), int(end_str)


def find_boundary(boundaries: np.ndarray, position: int) -> tuple[int, int] | None:
    for start, end in boundaries:
        if start <= position <= end:
            return int(start), int(end)
    return None


def boundary_stats(
    singletons: SingletonData,
    boundaries: np.ndarray,
    position: int,
    genotypes: np.ndarray,
    skip_fraction: float,
) -> dict[str, object]:
    boundary = find_boundary(boundaries, position)
    valid_mask = genotypes != -1
    daf = None if not np.any(valid_mask) else float(np.mean(genotypes[valid_mask]) / 2.0)
    if boundary is None:
        return {
            "boundary_up": None,
            "boundary_down": None,
            "valid_genotypes": int(np.sum(valid_mask)),
            "daf": daf,
            "upstream_nan_frac": None,
            "downstream_nan_frac": None,
            "skip_boundary": True,
            "status": "no_boundary",
        }

    boundary_up, boundary_down = boundary
    upstream, downstream = singletons.get_intervals(position, boundary_up, boundary_down)
    up_nan = float(np.isnan(upstream).mean())
    dn_nan = float(np.isnan(downstream).mean())
    return {
        "boundary_up": boundary_up,
        "boundary_down": boundary_down,
        "valid_genotypes": int(np.sum(valid_mask)),
        "daf": daf,
        "upstream_nan_frac": up_nan,
        "downstream_nan_frac": dn_nan,
        "skip_boundary": up_nan > skip_fraction or dn_nan > skip_fraction,
        "status": "ok",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Diagnose boundary-based row dropping in compute_SDS.py.")
    parser.add_argument("--s-file", required=True)
    parser.add_argument("--t-file", required=True)
    parser.add_argument("--b-file", required=True)
    parser.add_argument("--max-cols", type=int, default=10000)
    parser.add_argument("--skip-fraction", type=float, default=0.10)
    parser.add_argument("--position", action="append", type=int, default=[])
    parser.add_argument("--scan-window", action="append", default=[])
    parser.add_argument("--output", required=True)
    parser.add_argument("--pickle-cache-dir")
    args = parser.parse_args()

    cache = PickleCache(args.pickle_cache_dir)
    singletons = SingletonData.from_file(args.s_file, args.max_cols, cache)
    boundaries = cache.load_or_create("boundaries", args.b_file, lambda: FileReader.read_matrix(args.b_file))

    target_positions = set(args.position)
    scan_windows = [parse_scan_window(spec) for spec in args.scan_window]
    scan_summaries: dict[str, dict[str, object]] = {
        label: {
            "record_type": "scan_summary",
            "label": label,
            "start": start,
            "end": end,
            "rows_scanned": 0,
            "rows_boundary_ok": 0,
            "first_pos": None,
            "first_boundary_ok_pos": None,
            "best_pos": None,
            "best_max_nan_frac": None,
            "best_upstream_nan_frac": None,
            "best_downstream_nan_frac": None,
        }
        for label, start, end in scan_windows
    }
    position_rows: list[dict[str, object]] = []

    for snp_id, allele_anc, allele_der, position, genotypes in iter_test_snps(args.t_file):
        pos = int(position)

        if pos in target_positions:
            stats = boundary_stats(singletons, boundaries, pos, genotypes, args.skip_fraction)
            position_rows.append(
                {
                    "record_type": "position",
                    "label": snp_id,
                    "position": pos,
                    "allele_anc": allele_anc,
                    "allele_der": allele_der,
                    **stats,
                }
            )

        for label, start, end in scan_windows:
            if pos < start or pos > end:
                continue
            stats = boundary_stats(singletons, boundaries, pos, genotypes, args.skip_fraction)
            summary = scan_summaries[label]
            summary["rows_scanned"] += 1
            if summary["first_pos"] is None:
                summary["first_pos"] = pos
            if not stats["skip_boundary"]:
                summary["rows_boundary_ok"] += 1
                if summary["first_boundary_ok_pos"] is None:
                    summary["first_boundary_ok_pos"] = pos
            if stats["upstream_nan_frac"] is not None and stats["downstream_nan_frac"] is not None:
                max_nan = max(float(stats["upstream_nan_frac"]), float(stats["downstream_nan_frac"]))
                if summary["best_max_nan_frac"] is None or max_nan < summary["best_max_nan_frac"]:
                    summary["best_pos"] = pos
                    summary["best_max_nan_frac"] = max_nan
                    summary["best_upstream_nan_frac"] = stats["upstream_nan_frac"]
                    summary["best_downstream_nan_frac"] = stats["downstream_nan_frac"]

    fieldnames = [
        "record_type",
        "label",
        "position",
        "allele_anc",
        "allele_der",
        "boundary_up",
        "boundary_down",
        "valid_genotypes",
        "daf",
        "upstream_nan_frac",
        "downstream_nan_frac",
        "skip_boundary",
        "status",
        "start",
        "end",
        "rows_scanned",
        "rows_boundary_ok",
        "first_pos",
        "first_boundary_ok_pos",
        "best_pos",
        "best_max_nan_frac",
        "best_upstream_nan_frac",
        "best_downstream_nan_frac",
    ]
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in sorted(position_rows, key=lambda item: item["position"]):
            writer.writerow(row)
        for label, _, _ in scan_windows:
            writer.writerow(scan_summaries[label])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
