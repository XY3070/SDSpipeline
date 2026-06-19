#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a parent-like merged s_file from split-cohort s_files by removing singleton positions observed in the opposite cohort VCF."
    )
    parser.add_argument("--ncn-s-file", required=True)
    parser.add_argument("--scn-s-file", required=True)
    parser.add_argument("--ncn-mask-pos", required=True)
    parser.add_argument("--scn-mask-pos", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def load_mask(path: Path) -> set[str]:
    mask: set[str] = set()
    with path.open() as handle:
        for line in handle:
            pos = line.strip()
            if pos:
                mask.add(pos)
    return mask


def filter_row(line: str, mask: set[str]) -> str:
    line = line.rstrip("\n")
    if not line or line == "NA":
        return "NA"
    kept = [tok for tok in line.split("\t") if tok and tok not in mask]
    return "\t".join(kept) if kept else "NA"


def main() -> None:
    args = parse_args()
    ncn_mask = load_mask(Path(args.ncn_mask_pos))
    scn_mask = load_mask(Path(args.scn_mask_pos))

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("w") as out_handle:
        with Path(args.ncn_s_file).open() as handle:
            for line in handle:
                out_handle.write(filter_row(line, ncn_mask) + "\n")
        with Path(args.scn_s_file).open() as handle:
            for line in handle:
                out_handle.write(filter_row(line, scn_mask) + "\n")


if __name__ == "__main__":
    main()
