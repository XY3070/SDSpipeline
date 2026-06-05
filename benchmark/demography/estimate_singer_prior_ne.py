#!/usr/bin/env python3

from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Estimate a rough diploid Ne prior for SINGER from phased no-missing VCFs "
            "using variant-only pairwise diversity over variant span."
        )
    )
    parser.add_argument("--population", required=True, help="Population label for metadata.")
    parser.add_argument("--mu", type=float, default=1.25e-8, help="Mutation rate.")
    parser.add_argument("--output", required=True, help="Output JSON path.")
    parser.add_argument("vcfs", nargs="+", help="Input phased clean VCFs.")
    return parser.parse_args()


def summarize_vcf(path: Path) -> dict[str, float]:
    first_pos = None
    last_pos = None
    total_pi = 0.0
    variant_count = 0
    hap_count = None

    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            pos = int(fields[1])
            alt_count = 0
            local_haps = 0
            for sample_field in fields[9:]:
                gt = sample_field.split(":", 1)[0]
                alleles = gt.replace("|", "/").split("/")
                if any(allele == "." for allele in alleles):
                    raise ValueError(f"Missing genotype found in supposedly clean VCF: {path}")
                for allele in alleles:
                    if allele not in {"0", "1"}:
                        raise ValueError(f"Non-biallelic allele {allele!r} found in {path}")
                    alt_count += int(allele)
                    local_haps += 1
            if hap_count is None:
                hap_count = local_haps
            elif hap_count != local_haps:
                raise ValueError(f"Inconsistent haplotype count across variants in {path}")

            p = alt_count / local_haps
            total_pi += 2.0 * p * (1.0 - p) * (local_haps / (local_haps - 1.0))
            variant_count += 1
            if first_pos is None:
                first_pos = pos
            last_pos = pos

    if variant_count == 0 or first_pos is None or last_pos is None or hap_count is None:
        raise ValueError(f"No usable variants found in {path}")

    return {
        "path": str(path),
        "variant_count": float(variant_count),
        "first_pos": float(first_pos),
        "last_pos": float(last_pos),
        "span_bp": float(last_pos - first_pos + 1),
        "hap_count": float(hap_count),
        "pi_sum": float(total_pi),
    }


def main() -> int:
    args = parse_args()
    vcf_paths = [Path(item).resolve() for item in args.vcfs]
    for path in vcf_paths:
        if not path.exists():
            raise FileNotFoundError(f"VCF not found: {path}")

    per_vcf = [summarize_vcf(path) for path in vcf_paths]
    total_span = sum(item["span_bp"] for item in per_vcf)
    total_pi = sum(item["pi_sum"] for item in per_vcf)
    if total_span <= 0:
        raise ValueError("Total variant span must be positive")
    pi_estimate = total_pi / total_span
    ne_diploid = pi_estimate / (4.0 * args.mu)

    payload = {
        "population": args.population,
        "mu": args.mu,
        "vcf_count": len(vcf_paths),
        "pi_estimate_variant_span": pi_estimate,
        "ne_diploid_estimate": ne_diploid,
        "note": (
            "Approximate pi from variant-only phased clean VCFs, using per-chromosome "
            "variant span as denominator. Intended only as a rough SINGER prior."
        ),
        "per_vcf": per_vcf,
    }

    output_path = Path(args.output).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover - CLI error path
        print(f"[Error] {exc}", file=sys.stderr)
        raise
