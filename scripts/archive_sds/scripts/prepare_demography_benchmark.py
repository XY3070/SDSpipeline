#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import random
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SDS_ENV_PREFIX = Path("/data/home/grp-wangyf/intern/miniforge3/envs/sds")
DEFAULT_BCFTOOLS = SDS_ENV_PREFIX / "bin" / "bcftools"
DEFAULT_TABIX = SDS_ENV_PREFIX / "bin" / "tabix"
DEFAULT_SEED = 20260423


def parse_chromosomes(text: str) -> list[int]:
    values: list[int] = []
    for token in text.split(","):
        item = token.strip()
        if not item:
            continue
        if "-" in item:
            start_text, end_text = item.split("-", 1)
            start = int(start_text)
            end = int(end_text)
            if end < start:
                raise ValueError(f"Invalid chromosome range: {item}")
            values.extend(range(start, end + 1))
        else:
            values.append(int(item))
    if not values:
        raise ValueError("No chromosomes specified")
    return sorted(dict.fromkeys(values))


def run_checked(cmd: list[str]) -> None:
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def run_capture(cmd: list[str]) -> str:
    result = subprocess.run(cmd, check=True, text=True, capture_output=True)
    return result.stdout


def read_nonempty_lines(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def write_lines(path: Path, lines: list[str]) -> None:
    path.write_text("".join(f"{line}\n" for line in lines))


def bcftools_query_samples(bcftools: Path, vcf_path: Path) -> list[str]:
    output = run_capture([str(bcftools), "query", "-l", str(vcf_path)])
    samples = [line.strip() for line in output.splitlines() if line.strip()]
    if not samples:
        raise ValueError(f"No samples found in VCF: {vcf_path}")
    return samples


def select_subset(
    pool_samples: list[str],
    vcf_samples: list[str],
    subset_size: int,
    seed: int,
) -> tuple[list[str], list[str]]:
    pool_set = set(pool_samples)
    eligible_sorted = sorted(pool_set.intersection(vcf_samples))
    if len(eligible_sorted) < subset_size:
        raise ValueError(
            f"Requested {subset_size} samples but only found {len(eligible_sorted)} "
            "eligible samples in the VCF header intersection"
        )
    rng = random.Random(seed)
    chosen = set(rng.sample(eligible_sorted, subset_size))
    subset_in_vcf_order = [sample for sample in vcf_samples if sample in chosen]
    if len(subset_in_vcf_order) != subset_size:
        raise AssertionError("Subset size changed after reordering to VCF order")
    return eligible_sorted, subset_in_vcf_order


def choose_distinguished_pairs(
    subset_samples: list[str],
    num_pairs: int,
    seed: int,
    pair_mode: str,
) -> list[tuple[str, str]]:
    rng = random.Random(seed + 1)
    stable_subset = sorted(subset_samples)
    if pair_mode == "self":
        if len(stable_subset) < num_pairs:
            raise ValueError(
                f"Need at least {num_pairs} samples for self-pairs, got {len(stable_subset)}"
            )
        chosen = rng.sample(stable_subset, num_pairs)
        return [(sample, sample) for sample in chosen]
    if len(stable_subset) < 2 * num_pairs:
        raise ValueError(
            f"Need at least {2 * num_pairs} samples for cross-pairs, got {len(stable_subset)}"
        )
    chosen = rng.sample(stable_subset, 2 * num_pairs)
    return [(chosen[2 * i], chosen[2 * i + 1]) for i in range(num_pairs)]


def materialize_subset_vcfs(
    *,
    pop: str,
    vcf_root: Path,
    chromosomes: list[int],
    subset_file: Path,
    subset_size: int,
    subset_dir: Path,
    bcftools: Path,
    tabix: Path,
    threads: int,
    force: bool,
) -> list[str]:
    outputs: list[str] = []
    for chrom in chromosomes:
        input_vcf = vcf_root / f"UKBQC_{pop}_chr{chrom}.vcf.gz"
        if not input_vcf.exists():
            raise FileNotFoundError(f"VCF not found: {input_vcf}")
        output_vcf = subset_dir / f"UKBQC_{pop}_subset{subset_size}_chr{chrom}.vcf.gz"
        outputs.append(str(output_vcf))
        if output_vcf.exists() and not force:
            continue
        cmd = [
            str(bcftools),
            "view",
            "--threads",
            str(threads),
            "-S",
            str(subset_file),
            "-m2",
            "-M2",
            "-v",
            "snps",
            "-O",
            "z",
            "-o",
            str(output_vcf),
            str(input_vcf),
        ]
        run_checked(cmd)
        run_checked([str(tabix), "-f", "-p", "vcf", str(output_vcf)])
    return outputs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare fixed sample subsets and optional per-chromosome subset VCFs for "
            "the NCN/SCN phlash vs SMC++ demographic benchmark."
        )
    )
    parser.add_argument("--pop", required=True, help="Population label, e.g. NCN or SCN.")
    parser.add_argument(
        "--sample-pool",
        default=None,
        help="Text file containing cohort sample IDs. Defaults to sds/data/<POP>.txt.",
    )
    parser.add_argument(
        "--vcf-root",
        default=None,
        help="Directory containing cohort VCFs. Defaults to sds/data/vcf/<POP>.",
    )
    parser.add_argument(
        "--out-root",
        default=None,
        help="Benchmark output root. Defaults to benchmark/demography under the project root.",
    )
    parser.add_argument("--subset-size", type=int, default=100, help="Subset size.")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="Random seed.")
    parser.add_argument(
        "--num-distinguished-pairs",
        type=int,
        default=4,
        help="Number of SMC++ distinguished pairs to write.",
    )
    parser.add_argument(
        "--pair-mode",
        choices=("self", "cross"),
        default="self",
        help=(
            "Whether each distinguished pair should reuse the same sample twice "
            "(recommended for unphased VCF input) or pair different samples."
        ),
    )
    parser.add_argument(
        "--chromosomes",
        default="1-22",
        help="Chromosome list/ranges to materialize, e.g. 1-22 or 1,2,22.",
    )
    parser.add_argument(
        "--materialize-subset-vcfs",
        action="store_true",
        help="Also write subset SNP VCFs for each requested chromosome.",
    )
    parser.add_argument(
        "--bcftools",
        default=str(DEFAULT_BCFTOOLS),
        help="Path to bcftools. Defaults to the sds environment bcftools.",
    )
    parser.add_argument(
        "--tabix",
        default=str(DEFAULT_TABIX),
        help="Path to tabix. Defaults to the sds environment tabix.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Threads passed to bcftools view when materializing subset VCFs.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing subset VCFs and metadata files.",
    )
    args = parser.parse_args()
    args.pop = args.pop.upper()
    if args.sample_pool is None:
        args.sample_pool = str(PROJECT_ROOT / "sds" / "data" / f"{args.pop}.txt")
    if args.vcf_root is None:
        args.vcf_root = str(PROJECT_ROOT / "sds" / "data" / "vcf" / args.pop)
    if args.out_root is None:
        args.out_root = str(PROJECT_ROOT / "benchmark" / "demography")
    args.chromosomes = parse_chromosomes(args.chromosomes)
    return args


def main() -> int:
    args = parse_args()
    sample_pool = Path(args.sample_pool).resolve()
    vcf_root = Path(args.vcf_root).resolve()
    out_root = Path(args.out_root).resolve()
    bcftools = Path(args.bcftools).resolve()
    tabix = Path(args.tabix).resolve()
    first_vcf = vcf_root / f"UKBQC_{args.pop}_chr{args.chromosomes[0]}.vcf.gz"

    for path in [sample_pool, vcf_root, bcftools, tabix, first_vcf]:
        if not path.exists():
            raise FileNotFoundError(f"Required path not found: {path}")

    pop_dir = out_root / args.pop
    subset_dir = pop_dir / "subset_vcf"
    smcpp_dir = pop_dir / "smcpp"
    phlash_dir = pop_dir / "phlash"
    plots_dir = pop_dir / "plots"
    for directory in [pop_dir, subset_dir, smcpp_dir, phlash_dir, plots_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    pool_samples = read_nonempty_lines(sample_pool)
    vcf_samples = bcftools_query_samples(bcftools, first_vcf)
    eligible_samples, subset_samples = select_subset(
        pool_samples=pool_samples,
        vcf_samples=vcf_samples,
        subset_size=args.subset_size,
        seed=args.seed,
    )
    distinguished_pairs = choose_distinguished_pairs(
        subset_samples=subset_samples,
        num_pairs=args.num_distinguished_pairs,
        seed=args.seed,
        pair_mode=args.pair_mode,
    )

    subset_file = pop_dir / f"subset_{args.subset_size}.samples.txt"
    eligible_file = pop_dir / "sample_pool.intersection.txt"
    pair_file = pop_dir / "smcpp_distinguished_pairs.tsv"
    manifest_file = pop_dir / "benchmark_manifest.json"

    write_lines(subset_file, subset_samples)
    write_lines(eligible_file, eligible_samples)
    pair_lines = ["pair_id\tsample_id_1\tsample_id_2"]
    for index, (sample_a, sample_b) in enumerate(distinguished_pairs, start=1):
        pair_lines.append(f"pair{index}\t{sample_a}\t{sample_b}")
    write_lines(pair_file, pair_lines)

    subset_vcfs: list[str] = []
    if args.materialize_subset_vcfs:
        subset_vcfs = materialize_subset_vcfs(
            pop=args.pop,
            vcf_root=vcf_root,
            chromosomes=args.chromosomes,
            subset_file=subset_file,
            subset_size=args.subset_size,
            subset_dir=subset_dir,
            bcftools=bcftools,
            tabix=tabix,
            threads=args.threads,
            force=args.force,
        )

    manifest = {
        "population": args.pop,
        "seed": args.seed,
        "subset_size": args.subset_size,
        "pair_mode": args.pair_mode,
        "num_distinguished_pairs": args.num_distinguished_pairs,
        "chromosomes": args.chromosomes,
        "sample_pool": str(sample_pool),
        "vcf_root": str(vcf_root),
        "eligible_count": len(eligible_samples),
        "subset_samples_file": str(subset_file),
        "eligible_samples_file": str(eligible_file),
        "smcpp_distinguished_pairs_file": str(pair_file),
        "subset_vcf_dir": str(subset_dir),
        "subset_vcfs": subset_vcfs,
        "smcpp_output_dir": str(smcpp_dir),
        "phlash_output_dir": str(phlash_dir),
        "plots_dir": str(plots_dir),
    }
    manifest_file.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    print(
        json.dumps(
            {
                "population": args.pop,
                "eligible_count": len(eligible_samples),
                "subset_size": len(subset_samples),
                "subset_samples_file": str(subset_file),
                "distinguished_pairs_file": str(pair_file),
                "subset_vcfs_materialized": len(subset_vcfs),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover - CLI error path
        print(f"[Error] {exc}", file=sys.stderr)
        raise
