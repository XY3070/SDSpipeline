#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_olddefault_map(path: Path) -> dict[str, float]:
    mapping: dict[str, float] = {}
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            daf, shape = line.split("\t")
            mapping[f"{float(daf):.2f}"] = float(shape)
    return mapping


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    parser.add_argument("--olddefault-gfile", default="/data/home/grp-wangyf/xuyuan/sds/g_file.txt")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = Path(args.root).resolve()
    manifest = root / "manifest.tsv"
    out_tsv = root / "summary.tsv"
    out_md = root / "summary.md"
    olddefault = read_olddefault_map(Path(args.olddefault_gfile))

    rows: list[dict[str, str]] = []
    with manifest.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gamma_path = Path(row["gamma_path"])
            exists = gamma_path.is_file() and gamma_path.stat().st_size > 0
            raw_text = gamma_path.read_text().strip() if exists else ""
            shape_txt = raw_text
            if raw_text and "\t" in raw_text:
                parts = raw_text.split("\t")
                shape_txt = parts[-1]
            rows.append(
                {
                    "label": row["label"],
                    "makefile_source": row["makefile_source"],
                    "daf": row["daf"],
                    "present_ne": row["present_ne"],
                    "sim_reps": row["sim_reps"],
                    "gamma_path": str(gamma_path),
                    "shape": shape_txt,
                    "olddefault_shape": f"{olddefault.get(row['daf'], float('nan'))}",
                    "status": "ok" if exists else "missing",
                }
            )

    header = [
        "label",
        "makefile_source",
        "daf",
        "present_ne",
        "sim_reps",
        "status",
        "shape",
        "olddefault_shape",
        "gamma_path",
    ]
    with out_tsv.open("w") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(row[col] for col in header) + "\n")

    lines = [
        f"Root: `{root}`",
        "",
        "| Label | Makefile | DAF | Present `Ne(0)` | Reps | Status | Shape | Olddefault shape |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in rows:
        lines.append(
            f"| `{row['label']}` | `{row['makefile_source']}` | `{row['daf']}` | "
            f"`{row['present_ne']}` | `{row['sim_reps']}` | `{row['status']}` | `{row['shape']}` | `{row['olddefault_shape']}` |"
        )
    out_md.write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
