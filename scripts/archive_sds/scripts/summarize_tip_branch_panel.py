#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize tip-branch dump outputs from the Gravel_CHB Ne(0) panel.")
    parser.add_argument("--manifest", required=True, help="Panel manifest TSV produced by submit_gravel_chb_tip_branch_panel.py")
    parser.add_argument("--output-prefix", required=True, help="Output prefix for summary TSV and Markdown.")
    return parser.parse_args()


def summarize_branch_values(branch_files: list[Path]) -> dict[str, float | int | str]:
    values: list[float] = []
    for path in branch_files:
        with path.open() as handle:
            for line in handle:
                text = line.strip()
                if not text:
                    continue
                values.append(float(text))
    if not values:
        return {
            "n_values": 0,
            "mean": math.nan,
            "sd": math.nan,
            "shape_recomputed": math.nan,
            "q05": math.nan,
            "q25": math.nan,
            "q50": math.nan,
            "q75": math.nan,
            "q95": math.nan,
        }
    arr = np.asarray(values, dtype=float)
    mean = float(arr.mean())
    sd = float(arr.std(ddof=0))
    shape = float((mean * mean) / (sd * sd)) if sd > 0 else math.nan
    return {
        "n_values": int(arr.size),
        "mean": mean,
        "sd": sd,
        "shape_recomputed": shape,
        "q05": float(np.quantile(arr, 0.05)),
        "q25": float(np.quantile(arr, 0.25)),
        "q50": float(np.quantile(arr, 0.50)),
        "q75": float(np.quantile(arr, 0.75)),
        "q95": float(np.quantile(arr, 0.95)),
    }


def main() -> int:
    args = parse_args()
    manifest = Path(args.manifest).resolve()
    output_prefix = Path(args.output_prefix).resolve()
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str]] = []
    with manifest.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append(row)

    summary_rows: list[dict[str, object]] = []
    for row in rows:
        workdir = Path(row["workdir"])
        branch_files = sorted(workdir.glob("tip_branches_*.tab"))
        gamma_path = Path(row["gamma_output"])
        gamma_shape = math.nan
        if gamma_path.exists() and gamma_path.stat().st_size > 0:
            parts = gamma_path.read_text().strip().split("\t")
            if len(parts) >= 2:
                gamma_shape = float(parts[1])
        stats = summarize_branch_values(branch_files)
        summary_rows.append(
            {
                "scenario_label": row["scenario_label"],
                "present_ne": int(row["present_ne"]),
                "daf": float(row["daf"]),
                "sim_reps": int(row["sim_reps"]),
                "tip_branch_files": len(branch_files),
                "gamma_shape_file": gamma_shape,
                **stats,
            }
        )

    tsv_path = output_prefix.with_suffix(".summary.tsv")
    header = [
        "scenario_label",
        "present_ne",
        "daf",
        "sim_reps",
        "tip_branch_files",
        "gamma_shape_file",
        "n_values",
        "mean",
        "sd",
        "shape_recomputed",
        "q05",
        "q25",
        "q50",
        "q75",
        "q95",
    ]
    with tsv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)

    md_path = output_prefix.with_suffix(".summary.md")
    with md_path.open("w") as handle:
        handle.write("# Tip-branch panel summary\n\n")
        handle.write(f"- manifest: `{manifest}`\n")
        handle.write(f"- combinations: `{len(summary_rows)}`\n\n")
        handle.write("| scenario | Ne(0) | DAF | tip_values | mean | sd | shape_from_tips | shape_from_gamma |\n")
        handle.write("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n")
        for row in summary_rows:
            handle.write(
                f"| `{row['scenario_label']}` | {row['present_ne']} | {row['daf']:.2f} | {row['n_values']} | "
                f"{row['mean']:.6f} | {row['sd']:.6f} | {row['shape_recomputed']:.6f} | {row['gamma_shape_file']:.6f} |\n"
            )

    meta = {
        "manifest": str(manifest),
        "summary_tsv": str(tsv_path),
        "summary_md": str(md_path),
        "n_rows": len(summary_rows),
    }
    output_prefix.with_suffix(".summary.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
