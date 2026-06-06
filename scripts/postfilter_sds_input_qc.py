#!/usr/bin/env python3
"""Apply controlled QC post-filters to an existing SDS b/s/t/o input set."""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import statistics
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class VariantMetric:
    pos: int
    qual: float
    dp_avg: float | None
    exc_het: float | None
    f_missing: float
    ac: int | None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--chr", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--sample-list", required=True)
    parser.add_argument("--mask-bed", default=None)
    parser.add_argument("--target-bed", default=None, help="Optional BED intervals for the output t_file only.")
    parser.add_argument("--apply-mask", action="store_true")
    parser.add_argument("--apply-variant-qc", action="store_true")
    parser.add_argument("--apply-singleton-density-qc", action="store_true")
    parser.add_argument("--apply-excess-sample-qc", action="store_true")
    parser.add_argument("--observability-mode", choices=["copy", "callrate"], default="copy")
    parser.add_argument("--call-rate-min", type=float, default=0.90)
    parser.add_argument("--qual-min", type=float, default=56.0)
    parser.add_argument("--dp-low-factor", type=float, default=0.5)
    parser.add_argument("--dp-high-factor", type=float, default=2.0)
    parser.add_argument("--hwe-min", type=float, default=1e-6)
    parser.add_argument("--singleton-window-size", type=int, default=20_000)
    parser.add_argument("--singleton-window-step", type=int, default=10_000)
    parser.add_argument("--singleton-density-z", type=float, default=4.0)
    parser.add_argument("--sample-singleton-z", type=float, default=4.0)
    parser.add_argument("--metrics-tsv", default=None)
    return parser.parse_args()


def normalize_chrom(chrom: str) -> str:
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def load_sample_order(base_dir: Path, chrom: str, sample_list: Path) -> list[str]:
    sample_order = base_dir / f"chr{chrom}_sample_order.txt"
    source = sample_order if sample_order.exists() else sample_list
    samples = []
    with source.open() as handle:
        for raw in handle:
            line = raw.strip()
            if line and not line.startswith("#"):
                samples.append(line.split()[0])
    if not samples:
        raise SystemExit(f"No samples found in {source}")
    return samples


def load_mask_positions(mask_bed: Path | None, chrom: str) -> list[tuple[int, int, str]]:
    if mask_bed is None:
        return []
    target = normalize_chrom(chrom)
    intervals = []
    with mask_bed.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise SystemExit(f"Invalid BED line: {raw.rstrip()}")
            bed_chrom = normalize_chrom(fields[0])
            if bed_chrom != target:
                continue
            name = fields[3] if len(fields) > 3 else "mask"
            intervals.append((int(fields[1]) + 1, int(fields[2]), name))
    intervals.sort()
    return intervals


def in_intervals(pos: int, intervals: list[tuple[int, int, str]]) -> bool:
    for start, end, _name in intervals:
        if pos < start:
            return False
        if start <= pos <= end:
            return True
    return False


def run_bcftools_metrics(vcf: Path, sample_list: Path, chrom: str, out_tsv: Path) -> None:
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    target_chroms = [f"chr{normalize_chrom(chrom)}", normalize_chrom(chrom)]
    region = None
    header = subprocess.run(
        ["bcftools", "view", "-h", str(vcf)],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    for candidate in target_chroms:
        if f"##contig=<ID={candidate}" in header:
            region = candidate
            break
    region = region or target_chroms[0]

    command = (
        f"bcftools view -r {region} -S {sample_list} -f PASS -m2 -M2 -v snps --force-samples -Ou {vcf} | "
        "bcftools +fill-tags -Ou -- -t AC,F_MISSING | "
        r"bcftools query -f '%POS\t%QUAL\t%INFO/DP_AVG\t%INFO/ExcHet\t%INFO/F_MISSING\t%INFO/AC\n'"
    )
    with out_tsv.open("w") as out:
        subprocess.run(["bash", "-lc", command], check=True, stdout=out)


def parse_float(value: str) -> float | None:
    if value in {".", ""}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def load_metrics(metrics_tsv: Path) -> tuple[dict[int, VariantMetric], float]:
    metrics: dict[int, VariantMetric] = {}
    dp_values = []
    with metrics_tsv.open() as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            pos = int(fields[0])
            qual = parse_float(fields[1])
            dp_avg = parse_float(fields[2])
            exc_het = parse_float(fields[3])
            f_missing = parse_float(fields[4])
            ac_value = parse_float(fields[5].split(",")[0])
            if qual is None or f_missing is None:
                continue
            ac = int(ac_value) if ac_value is not None else None
            if dp_avg is not None and dp_avg > 0:
                dp_values.append(dp_avg)
            metrics[pos] = VariantMetric(pos, qual, dp_avg, exc_het, f_missing, ac)
    if not metrics:
        raise SystemExit(f"No variant metrics loaded from {metrics_tsv}")
    dp_mean = statistics.mean(dp_values) if dp_values else math.nan
    return metrics, dp_mean


def passes_variant_qc(metric: VariantMetric | None, args: argparse.Namespace, dp_mean: float) -> bool:
    if metric is None:
        return False
    if metric.f_missing > (1.0 - args.call_rate_min):
        return False
    if metric.qual < args.qual_min:
        return False
    if metric.exc_het is not None and metric.exc_het <= args.hwe_min:
        return False
    if metric.dp_avg is not None and not math.isnan(dp_mean):
        if metric.dp_avg < args.dp_low_factor * dp_mean:
            return False
        if metric.dp_avg > args.dp_high_factor * dp_mean:
            return False
    return True


def read_s_file(path: Path, sample_names: list[str]) -> list[list[int]]:
    rows = []
    with path.open() as handle:
        for idx, raw in enumerate(handle):
            fields = raw.strip().split()
            if fields and idx < len(sample_names) and fields[0] == sample_names[idx]:
                fields = fields[1:]
            if not fields or fields == ["NA"]:
                rows.append([])
            else:
                rows.append([int(value) for value in fields if value != "NA"])
    return rows


def write_s_file(path: Path, rows: list[list[int]], keep_samples: list[int]) -> None:
    with path.open("w") as handle:
        for idx in keep_samples:
            positions = rows[idx]
            handle.write("NA\n" if not positions else "\t".join(map(str, positions)) + "\n")


def density_outlier_positions(
    singleton_rows: list[list[int]],
    window_size: int,
    window_step: int,
    z: float,
) -> tuple[set[int], list[tuple[int, int, int, bool]]]:
    positions = sorted({pos for row in singleton_rows for pos in row})
    if not positions:
        return set(), []
    min_pos = (positions[0] // window_step) * window_step + 1
    max_pos = positions[-1]
    counts = []
    left = 0
    right = 0
    n = len(positions)
    start = min_pos
    while start <= max_pos:
        end = start + window_size - 1
        while left < n and positions[left] < start:
            left += 1
        while right < n and positions[right] <= end:
            right += 1
        counts.append((start, end, right - left))
        start += window_step
    values = [count for _start, _end, count in counts]
    mean = statistics.mean(values)
    sd = statistics.pstdev(values) if len(values) > 1 else 0.0
    threshold = mean + z * sd
    outlier_windows = [(s, e, c, c > threshold) for s, e, c in counts]
    outlier_positions = set()
    for start, end, count, is_outlier in outlier_windows:
        if not is_outlier:
            continue
        for pos in positions:
            if pos < start:
                continue
            if pos > end:
                break
            outlier_positions.add(pos)
    return outlier_positions, outlier_windows


def filter_singletons(
    rows: list[list[int]],
    metrics: dict[int, VariantMetric],
    dp_mean: float,
    intervals: list[tuple[int, int, str]],
    args: argparse.Namespace,
) -> tuple[list[list[int]], set[int], list[tuple[int, int, int, bool]]]:
    filtered = []
    for row in rows:
        new_row = []
        for pos in row:
            if args.apply_mask and in_intervals(pos, intervals):
                continue
            if args.apply_variant_qc and not passes_variant_qc(metrics.get(pos), args, dp_mean):
                continue
            metric = metrics.get(pos)
            if args.apply_variant_qc and (metric is None or metric.ac != 1):
                continue
            new_row.append(pos)
        filtered.append(new_row)
    outlier_positions: set[int] = set()
    windows: list[tuple[int, int, int, bool]] = []
    if args.apply_singleton_density_qc:
        outlier_positions, windows = density_outlier_positions(
            filtered,
            args.singleton_window_size,
            args.singleton_window_step,
            args.singleton_density_z,
        )
        filtered = [[pos for pos in row if pos not in outlier_positions] for row in filtered]
    return filtered, outlier_positions, windows


def select_keep_samples(rows: list[list[int]], z: float, apply: bool) -> tuple[list[int], list[tuple[int, int, bool]], float]:
    counts = [len(row) for row in rows]
    if not counts:
        return [], [], math.nan
    mean = statistics.mean(counts)
    sd = statistics.pstdev(counts) if len(counts) > 1 else 0.0
    threshold = mean + z * sd
    keep = []
    detail = []
    for idx, count in enumerate(counts):
        drop = apply and count > threshold
        detail.append((idx, count, drop))
        if not drop:
            keep.append(idx)
    return keep, detail, threshold


def filter_t_file(
    in_path: Path,
    out_path: Path,
    keep_samples: list[int],
    metrics: dict[int, VariantMetric],
    dp_mean: float,
    intervals: list[tuple[int, int, str]],
    target_intervals: list[tuple[int, int, str]],
    args: argparse.Namespace,
) -> tuple[int, int, list[int], list[int]]:
    rows_in = 0
    rows_out = 0
    total_calls = [0 for _ in keep_samples]
    missing_calls = [0 for _ in keep_samples]
    with in_path.open() as inp, out_path.open("w") as out:
        for raw in inp:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            rows_in += 1
            pos = int(fields[3])
            if target_intervals and not in_intervals(pos, target_intervals):
                continue
            if args.apply_mask and in_intervals(pos, intervals):
                continue
            if args.apply_variant_qc and not passes_variant_qc(metrics.get(pos), args, dp_mean):
                continue
            genotypes = fields[4:]
            kept_genotypes = [genotypes[idx] for idx in keep_samples]
            for i, gt in enumerate(kept_genotypes):
                total_calls[i] += 1
                if gt == "NA":
                    missing_calls[i] += 1
            out.write("\t".join(fields[:4] + kept_genotypes) + "\n")
            rows_out += 1
    return rows_in, rows_out, total_calls, missing_calls


def write_vector(path: Path, values: list[float]) -> None:
    with path.open("w") as handle:
        handle.write("\t".join(f"{value:.8f}" for value in values) + "\n")


def main() -> int:
    args = parse_args()
    chrom = normalize_chrom(args.chr)
    base = Path(args.base_input_dir)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    qc_dir = out / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)

    metrics_tsv = Path(args.metrics_tsv) if args.metrics_tsv else qc_dir / f"chr{chrom}.variant_metrics.tsv"
    metrics: dict[int, VariantMetric] = {}
    dp_mean = math.nan
    if args.apply_variant_qc:
        if not metrics_tsv.exists():
            run_bcftools_metrics(Path(args.vcf), Path(args.sample_list), chrom, metrics_tsv)
        metrics, dp_mean = load_metrics(metrics_tsv)

    intervals = load_mask_positions(Path(args.mask_bed) if args.mask_bed else None, chrom)
    target_intervals = load_mask_positions(Path(args.target_bed) if args.target_bed else None, chrom)
    samples = load_sample_order(base, chrom, Path(args.sample_list))
    s_rows = read_s_file(base / f"chr{chrom}_s_file.txt", samples)
    if len(s_rows) != len(samples):
        raise SystemExit(f"s_file row count {len(s_rows)} != sample count {len(samples)}")

    filtered_s, outlier_positions, density_windows = filter_singletons(s_rows, metrics, dp_mean, intervals, args)
    keep_samples, sample_detail, sample_threshold = select_keep_samples(
        filtered_s, args.sample_singleton_z, args.apply_excess_sample_qc
    )
    kept_sample_names = [samples[idx] for idx in keep_samples]

    shutil.copy2(base / f"chr{chrom}_b_file.txt", out / f"chr{chrom}_b_file.txt")
    write_s_file(out / f"chr{chrom}_s_file.txt", filtered_s, keep_samples)
    with (out / f"chr{chrom}_sample_order.txt").open("w") as handle:
        for sample in kept_sample_names:
            handle.write(sample + "\n")

    rows_in, rows_out, total_calls, missing_calls = filter_t_file(
        base / f"chr{chrom}_t_file.txt",
        out / f"chr{chrom}_t_file.txt",
        keep_samples,
        metrics,
        dp_mean,
        intervals,
        target_intervals,
        args,
    )
    if args.observability_mode == "copy":
        old_values = []
        with (base / f"chr{chrom}_o_file.txt").open() as handle:
            old_values = handle.read().strip().split()
        values = [float(old_values[idx]) for idx in keep_samples]
    else:
        values = [
            1.0 - (missing / total if total > 0 else 1.0)
            for total, missing in zip(total_calls, missing_calls)
        ]
    write_vector(out / f"chr{chrom}_o_file.txt", values)

    with (qc_dir / "singleton_density_windows.tsv").open("w") as handle:
        handle.write("start\tend\tsingleton_count\tis_outlier\n")
        for start, end, count, is_outlier in density_windows:
            handle.write(f"{start}\t{end}\t{count}\t{int(is_outlier)}\n")

    with (qc_dir / "singleton_sample_counts.tsv").open("w") as handle:
        handle.write("sample_id\tsingleton_count\tdropped_excess_singleton\n")
        for idx, count, dropped in sample_detail:
            handle.write(f"{samples[idx]}\t{count}\t{int(dropped)}\n")

    with (qc_dir / "input_qc_summary.tsv").open("w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerow(["dp_mean", f"{dp_mean:.6f}" if not math.isnan(dp_mean) else "NA"])
        writer.writerow(["mask_intervals", len(intervals)])
        writer.writerow(["target_intervals", len(target_intervals)])
        writer.writerow(["variant_metrics_rows", len(metrics)])
        writer.writerow(["t_rows_in", rows_in])
        writer.writerow(["t_rows_out", rows_out])
        writer.writerow(["samples_in", len(samples)])
        writer.writerow(["samples_out", len(keep_samples)])
        writer.writerow(["sample_singleton_drop_threshold", f"{sample_threshold:.6f}"])
        writer.writerow(["singleton_positions_density_dropped", len(outlier_positions)])
        writer.writerow(["apply_mask", int(args.apply_mask)])
        writer.writerow(["apply_variant_qc", int(args.apply_variant_qc)])
        writer.writerow(["apply_singleton_density_qc", int(args.apply_singleton_density_qc)])
        writer.writerow(["apply_excess_sample_qc", int(args.apply_excess_sample_qc)])
        writer.writerow(["observability_mode", args.observability_mode])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
