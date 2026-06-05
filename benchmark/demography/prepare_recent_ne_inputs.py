#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_WORKSPACE_ROOT = Path(os.environ.get("SDS_WORKSPACE_ROOT", PROJECT_ROOT.parent / "SDSworkspace")).resolve()
DEFAULT_RESULTS_ROOT = Path(os.environ.get("SDS_RESULTS_ROOT", DEFAULT_WORKSPACE_ROOT / "results")).resolve()
DEFAULT_EXTERNAL_ROOT = Path(os.environ.get("SDS_EXTERNAL_ROOT", DEFAULT_WORKSPACE_ROOT / "external")).resolve()
DEFAULT_ROOT = Path(os.environ.get("SDS_DEMOGRAPHY_ROOT", DEFAULT_RESULTS_ROOT / "production" / "demography")).resolve()
DEFAULT_RELATE_DIR = Path(os.environ.get("SDS_RELATE_DIR", DEFAULT_EXTERNAL_ROOT / "relate")).resolve()
DEFAULT_SHARED_RELATE_ROOT = DEFAULT_ROOT / "relate_shared"
DEFAULT_RELATE_DOWNLOAD_SCRIPT = PROJECT_ROOT / "scripts" / "run_relate_download_refs.sh"
DEFAULT_RELATE_PHASED_VCF_ROOT = Path(
    os.environ.get("SDS_VCF_ROOT", DEFAULT_WORKSPACE_ROOT / "input" / "raw" / "vcf")
).resolve()


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


def run_checked_capture_lines(cmd: list[str]) -> list[str]:
    print("+", " ".join(cmd), flush=True)
    result = subprocess.run(cmd, check=True, text=True, capture_output=True)
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def read_nonempty_lines(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def strip_known_suffix(path: Path, suffix: str) -> Path:
    if not path.name.endswith(suffix):
        raise ValueError(f"Expected {path} to end with {suffix}")
    return path.with_name(path.name[: -len(suffix)])


def materialize_reference(path: Path, force: bool) -> Path:
    if path.suffix != ".gz":
        return path
    out_path = Path(str(path)[: -3])
    if out_path.exists() and not force:
        return out_path
    with subprocess.Popen(["gzip", "-cd", str(path)], stdout=subprocess.PIPE) as proc:
        assert proc.stdout is not None
        out_path.write_bytes(proc.stdout.read())
        ret = proc.wait()
    if ret != 0:
        raise RuntimeError(f"Failed to materialize reference: {path}")
    return out_path


def ensure_file_removed(path: Path) -> None:
    if path.exists() or path.is_symlink():
        path.unlink()


def ensure_nonempty_targets(paths: list[Path]) -> bool:
    return all(path.exists() and path.stat().st_size > 0 for path in paths)


def targets_outdated(targets: list[Path], sources: list[Path]) -> bool:
    if not targets or not sources:
        return False
    newest_source = max(path.stat().st_mtime for path in sources)
    oldest_target = min(path.stat().st_mtime for path in targets)
    return newest_source > oldest_target


def normalize_annot_file(path: Path) -> None:
    if not path.exists():
        return
    lines = [line.rstrip("\n").rstrip(";") + "\n" for line in path.read_text().splitlines()]
    path.write_text("".join(lines))


def ensure_relate_refs(shared_root: Path, download_script: Path, force: bool) -> Path:
    ref_root = shared_root / "refs"
    sentinel = ref_root / "hg38_ancestor_chr1.fa.gz"
    if not sentinel.exists():
        shared_root.mkdir(parents=True, exist_ok=True)
        run_checked([str(download_script), "--out-root", str(shared_root)])
    if not sentinel.exists():
        raise FileNotFoundError(f"Relate reference sentinel missing after download: {sentinel}")
    return ref_root


def build_poplabels(subset_samples: list[str], pop: str, out_path: Path) -> None:
    lines = ["sample\tpopulation\tgroup\tsex"]
    for sample in subset_samples:
        lines.append(f"{sample}\t{pop}\t{pop}\t0")
    out_path.write_text("".join(f"{line}\n" for line in lines))


def stage_phased_subset_vcf(source_vcf: Path, subset_file: Path, target_vcf: Path, force: bool) -> Path:
    target_tbi = Path(f"{target_vcf}.tbi")
    target_vcf.parent.mkdir(parents=True, exist_ok=True)
    if force:
        ensure_file_removed(target_vcf)
        ensure_file_removed(target_tbi)
    if not target_vcf.exists():
        run_checked(
            [
                "bcftools",
                "view",
                "-S",
                str(subset_file),
                "-Oz",
                "-o",
                str(target_vcf),
                str(source_vcf),
            ]
        )
    if force or not target_tbi.exists():
        run_checked(["tabix", "-f", "-p", "vcf", str(target_vcf)])
    return target_vcf


def stage_symlink(src: Path, dst: Path, force: bool) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.is_symlink() or dst.exists():
        if not force:
            return
        dst.unlink()
    dst.symlink_to(src)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare reusable benchmark inputs for recent-time Ne sensitivity runs with "
            "Relate and SINGER."
        )
    )
    parser.add_argument("--pop", required=True, help="Population label, e.g. NCN or SCN.")
    parser.add_argument(
        "--root",
        default=str(DEFAULT_ROOT),
        help="Benchmark demography root. Defaults to benchmark/demography.",
    )
    parser.add_argument(
        "--chromosomes",
        default="1-22",
        help="Chromosome list/ranges to stage, e.g. 1-22 or 20-22.",
    )
    parser.add_argument(
        "--relate-dir",
        default=str(DEFAULT_RELATE_DIR),
        help="Path to the local Relate checkout.",
    )
    parser.add_argument(
        "--shared-relate-root",
        default=str(DEFAULT_SHARED_RELATE_ROOT),
        help="Shared Relate reference root, used with run_relate_download_refs.sh.",
    )
    parser.add_argument(
        "--download-relate-refs-script",
        default=str(DEFAULT_RELATE_DOWNLOAD_SCRIPT),
        help="Path to run_relate_download_refs.sh.",
    )
    parser.add_argument(
        "--relate-phased-vcf-dir",
        default=None,
        help=(
            "Directory containing phased per-chromosome VCFs named "
            "UKBQC_<POP>_chr<CHR>.phased.vcf.gz. Defaults to SDS_VCF_ROOT/<POP>/shapeit5."
        ),
    )
    parser.add_argument(
        "--skip-relate",
        action="store_true",
        help="Skip preparing Relate raw/prepared inputs.",
    )
    parser.add_argument(
        "--skip-singer",
        action="store_true",
        help="Skip staging SINGER input symlinks.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing prepared files and symlinks.",
    )
    args = parser.parse_args()
    args.pop = args.pop.upper()
    args.chromosomes = parse_chromosomes(args.chromosomes)
    return args


def main() -> int:
    args = parse_args()

    root = Path(args.root).resolve()
    pop_root = root / args.pop
    subset_file = pop_root / "subset_100.samples.txt"
    subset_vcf_dir = pop_root / "subset_vcf"
    phlash_clean_dir = pop_root / "phlash_clean"
    relate_dir = Path(args.relate_dir).resolve()
    shared_relate_root = Path(args.shared_relate_root).resolve()
    download_script = Path(args.download_relate_refs_script).resolve()
    if args.relate_phased_vcf_dir is None:
        relate_phased_vcf_dir = (DEFAULT_RELATE_PHASED_VCF_ROOT / args.pop / "shapeit5").resolve()
    else:
        relate_phased_vcf_dir = Path(args.relate_phased_vcf_dir).resolve()

    for path in [subset_file, subset_vcf_dir, phlash_clean_dir, relate_dir, download_script]:
        if not path.exists():
            raise FileNotFoundError(f"Required path not found: {path}")
    subset_samples = read_nonempty_lines(subset_file)
    if not subset_samples:
        raise ValueError(f"No subset samples found in {subset_file}")
    subset_size = len(subset_samples)
    subset_base = f"UKBQC_{args.pop}_subset{subset_size}"

    ref_root: Path | None = None

    relate_root = pop_root / "relate_recent"
    singer_root = pop_root / "singer_recent"
    poplabels_path = relate_root / f"{args.pop}.poplabels"
    if not args.skip_relate:
        relate_root.mkdir(parents=True, exist_ok=True)

    manifest: dict[str, object] = {
        "population": args.pop,
        "subset_samples_file": str(subset_file),
        "subset_size": subset_size,
        "chromosomes": args.chromosomes,
        "skip_relate": bool(args.skip_relate),
        "skip_singer": bool(args.skip_singer),
        "shared_relate_root": str(shared_relate_root),
        "relate_root": str(relate_root),
        "relate_phased_vcf_dir": str(relate_phased_vcf_dir),
        "singer_root": str(singer_root),
        "chromosome_outputs": {},
    }
    if not args.skip_relate:
        manifest["relate_poplabels"] = str(poplabels_path)

    chrom_outputs = manifest["chromosome_outputs"]
    assert isinstance(chrom_outputs, dict)
    poplabels_written = False

    for chrom in args.chromosomes:
        subset_vcf = subset_vcf_dir / f"{subset_base}_chr{chrom}.vcf.gz"
        clean_vcf = phlash_clean_dir / f"{subset_base}_chr{chrom}.phlash.vcf.gz"
        if not subset_vcf.exists():
            raise FileNotFoundError(f"Subset VCF not found: {subset_vcf}")
        if not clean_vcf.exists():
            raise FileNotFoundError(f"phlash-clean VCF not found: {clean_vcf}")

        entry: dict[str, object] = {
            "subset_vcf": str(subset_vcf),
            "phlash_clean_vcf": str(clean_vcf),
        }

        if not args.skip_relate:
            phased_source_vcf = relate_phased_vcf_dir / f"UKBQC_{args.pop}_chr{chrom}.phased.vcf.gz"
            raw_dir = relate_root / "raw"
            prepared_dir = relate_root / "prepared"
            phased_subset_dir = relate_root / "phased_subset_vcf"
            raw_dir.mkdir(parents=True, exist_ok=True)
            prepared_dir.mkdir(parents=True, exist_ok=True)

            raw_prefix = raw_dir / f"{subset_base}_chr{chrom}"
            prepared_prefix = prepared_dir / f"{subset_base}_chr{chrom}"
            raw_targets = [raw_prefix.with_suffix(".haps"), raw_prefix.with_suffix(".sample")]
            prepared_targets = [
                prepared_prefix.with_suffix(".haps.gz"),
                prepared_prefix.with_suffix(".sample.gz"),
                prepared_prefix.with_suffix(".dist.gz"),
                prepared_prefix.with_suffix(".annot"),
            ]
            phased_subset_vcf = phased_subset_dir / f"{subset_base}_chr{chrom}.phased.vcf.gz"
            needs_relate_prep = (
                args.force
                or not ensure_nonempty_targets(prepared_targets)
                or not poplabels_path.exists()
            )

            ancestor_path = shared_relate_root / "refs" / f"hg38_ancestor_chr{chrom}.fa"
            mask_path = shared_relate_root / "refs" / f"hg38_mask_chr{chrom}.fa"

            if needs_relate_prep:
                if not relate_phased_vcf_dir.exists():
                    raise FileNotFoundError(
                        f"Required phased Relate VCF directory not found: {relate_phased_vcf_dir}"
                    )
                if not phased_source_vcf.exists():
                    raise FileNotFoundError(f"Phased Relate source VCF not found: {phased_source_vcf}")
                if ref_root is None:
                    ref_root = ensure_relate_refs(shared_relate_root, download_script, args.force)

                ancestor_path = materialize_reference(ref_root / f"hg38_ancestor_chr{chrom}.fa.gz", args.force)
                mask_path = materialize_reference(ref_root / f"hg38_mask_chr{chrom}.fa.gz", args.force)

                phased_subset_vcf = stage_phased_subset_vcf(
                    phased_source_vcf,
                    subset_file,
                    phased_subset_vcf,
                    args.force,
                )
                if not poplabels_written and (args.force or not poplabels_path.exists()):
                    build_poplabels(
                        run_checked_capture_lines(["bcftools", "query", "-l", str(phased_subset_vcf)]),
                        args.pop,
                        poplabels_path,
                    )
                    poplabels_written = True

            else:
                if not poplabels_path.exists():
                    raise FileNotFoundError(f"Relate poplabels file not found: {poplabels_path}")

            raw_needs_rebuild = (
                needs_relate_prep
                and (
                    args.force
                or not ensure_nonempty_targets(raw_targets)
                or targets_outdated(raw_targets, [phased_subset_vcf])
                )
            )
            if raw_needs_rebuild:
                for path in raw_targets:
                    ensure_file_removed(path)
                run_checked(
                    [
                        str(relate_dir / "bin" / "RelateFileFormats"),
                        "--mode",
                        "ConvertFromVcf",
                        "--haps",
                        str(raw_prefix.with_suffix(".haps")),
                        "--sample",
                        str(raw_prefix.with_suffix(".sample")),
                        "-i",
                        str(strip_known_suffix(phased_subset_vcf, ".vcf.gz")),
                    ]
                )

            prepared_needs_rebuild = (
                needs_relate_prep
                and (
                    args.force
                or not ensure_nonempty_targets(prepared_targets)
                or targets_outdated(
                    prepared_targets,
                    raw_targets + [poplabels_path, ancestor_path, mask_path],
                )
                )
            )
            if prepared_needs_rebuild:
                for suffix in [
                    ".haps.gz",
                    ".sample.gz",
                    ".dist.gz",
                    ".annot",
                    "_biall.haps",
                    "_ancest.haps",
                    "_filtered.haps",
                    "_filtered.dist",
                ]:
                    ensure_file_removed(prepared_prefix.with_name(prepared_prefix.name + suffix))
                run_checked(
                    [
                        "env",
                        "-u",
                        "SHELLOPTS",
                        "bash",
                        str(relate_dir / "scripts" / "PrepareInputFiles" / "PrepareInputFiles.sh"),
                        "--haps",
                        str(raw_prefix.with_suffix(".haps")),
                        "--sample",
                        str(raw_prefix.with_suffix(".sample")),
                        "--ancestor",
                        str(ancestor_path),
                        "--mask",
                        str(mask_path),
                        "--poplabels",
                        str(poplabels_path),
                        "-o",
                        str(prepared_prefix),
                    ]
                )
            normalize_annot_file(prepared_prefix.with_suffix(".annot"))

            entry["relate_raw_prefix"] = str(raw_prefix)
            entry["relate_prepared_prefix"] = str(prepared_prefix)
            entry["relate_ancestor_fasta"] = str(ancestor_path)
            entry["relate_mask_fasta"] = str(mask_path)
            entry["relate_phased_source_vcf"] = str(phased_source_vcf)
            entry["relate_phased_subset_vcf"] = str(phased_subset_vcf)

        if not args.skip_singer:
            singer_input_dir = singer_root / "input"
            singer_input_dir.mkdir(parents=True, exist_ok=True)
            singer_vcf = singer_input_dir / clean_vcf.name
            singer_tbi = singer_input_dir / f"{clean_vcf.name}.tbi"
            stage_symlink(clean_vcf, singer_vcf, args.force)
            if clean_vcf.with_suffix(clean_vcf.suffix + ".tbi").exists():
                stage_symlink(clean_vcf.with_suffix(clean_vcf.suffix + ".tbi"), singer_tbi, args.force)
            entry["singer_input_prefix"] = str(strip_known_suffix(singer_vcf, ".vcf.gz"))

        chrom_outputs[f"chr{chrom}"] = entry

    manifest_path = pop_root / "recent_ne_inputs_manifest.json"
    write_json(manifest_path, manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover - CLI error path
        print(f"[Error] {exc}", file=sys.stderr)
        raise
