#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def summarize(path: Path, maf_cutoff: float) -> dict[str, float]:
    stats = {
        "n_all": 0,
        "n_maf_cutoff": 0,
        "ge20_all": 0,
        "ge40_all": 0,
        "ge20_maf_cutoff": 0,
        "ge40_maf_cutoff": 0,
        "max_all": 0.0,
        "max_maf_cutoff": 0.0,
    }
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                maf = float(row["MAF"])
                if "NEG_LOG10_P_BOTHSIDE" in row and row["NEG_LOG10_P_BOTHSIDE"] not in {"", None}:
                    neglogp = float(row["NEG_LOG10_P_BOTHSIDE"])
                else:
                    p = float(row["p_bothside"])
                    if not math.isfinite(p) or p <= 0.0:
                        neglogp = 300.0
                    else:
                        neglogp = -math.log10(p)
            except Exception:
                continue
            stats["n_all"] += 1
            if neglogp >= 20.0:
                stats["ge20_all"] += 1
            if neglogp >= 40.0:
                stats["ge40_all"] += 1
            stats["max_all"] = max(stats["max_all"], neglogp)
            if maf >= maf_cutoff:
                stats["n_maf_cutoff"] += 1
                if neglogp >= 20.0:
                    stats["ge20_maf_cutoff"] += 1
                if neglogp >= 40.0:
                    stats["ge40_maf_cutoff"] += 1
                stats["max_maf_cutoff"] = max(stats["max_maf_cutoff"], neglogp)
    return stats


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize Manhattan tail behavior before/after a MAF cutoff.")
    parser.add_argument("--input", required=True, help="Normalized TSV path.")
    parser.add_argument("--label", required=True, help="Track label.")
    parser.add_argument("--maf-cutoff", type=float, default=0.05, help="MAF cutoff to summarize.")
    parser.add_argument("--output", default=None, help="Optional TSV output path.")
    args = parser.parse_args()

    stats = summarize(Path(args.input), args.maf_cutoff)
    row = {
        "label": args.label,
        "maf_cutoff": f"{args.maf_cutoff:.4f}",
        **{k: str(v) for k, v in stats.items()},
    }

    if args.output:
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row.keys()), delimiter="\t")
            writer.writeheader()
            writer.writerow(row)
    else:
        print("\t".join(row.keys()))
        print("\t".join(row.values()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
