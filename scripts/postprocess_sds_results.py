#!/usr/bin/env python3

import argparse
import csv
import math
import os
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import numpy as np
import pyarrow as pa
import pyarrow.csv as pacsv
from scipy.stats import norm

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_WORKSPACE_ROOT = Path(os.environ.get("SDS_WORKSPACE_ROOT", PIPELINE_ROOT.parent / "SDSworkspace"))
DEFAULT_SDS_OUTPUT_ROOT = Path(
    os.environ.get("SDS_SDS_OUTPUT_ROOT", DEFAULT_WORKSPACE_ROOT / "results" / "production" / "sds_output")
)

SPLIT_BIN_LOW_WIDTH = 0.05
SPLIT_BIN_HIGH_WIDTH = 0.1
SPLIT_BIN_TRANSITION = 0.1

RESULT_SCHEMA = pa.schema(
    [
        ("ID", pa.string()),
        ("AA", pa.string()),
        ("DA", pa.string()),
        ("POS", pa.int64()),
        ("DAF", pa.float64()),
        ("nG0", pa.int64()),
        ("nG1", pa.int64()),
        ("nG2", pa.int64()),
        ("rSDS", pa.float64()),
        ("SuggestedInitPoint", pa.string()),
    ]
)


@dataclass(frozen=True)
class RegionTopHit:
    snp_id: str = ""
    af: float = math.nan
    maf: float = math.nan
    sds: float = math.nan
    pval: float = math.inf


class RegionSummary:
    __slots__ = ("chrom", "start", "end", "significant", "top_hit")

    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.significant = 0
        self.top_hit = RegionTopHit()


@dataclass(frozen=True)
class BinStats:
    bin_id: str
    maf_start: float
    maf_end: float
    mean_common: float
    sd_common: float
    common_variant_count: int


@dataclass(frozen=True)
class NormalizationStats:
    bin_stats: dict
    common_variant_count: int
    total_snvs: int
    chrom_max_pos: dict
    chrom_counts: dict


@dataclass(frozen=True)
class PlotTsvStats:
    plot_point_count: int
    plot_max_neg_log10_p: float
    chrom_centers: list
    genome_size: int
    bonferroni_threshold_y: float


def normalize_chrom_label(chrom):
    value = str(chrom).strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    return value


def parse_args():
    parser = argparse.ArgumentParser(
        description="MAF-binned common-variant normalization and Manhattan plot for SDS results."
    )
    parser.add_argument(
        "--input-dir",
        default=str(DEFAULT_SDS_OUTPUT_ROOT / "NCN"),
        help="Directory containing chr*_p.sds.tsv and chr*_q.sds.tsv files.",
    )
    parser.add_argument(
        "--output-prefix",
        default=None,
        help="Prefix for outputs. Defaults to <input-dir>/<POP>.",
    )
    parser.add_argument(
        "--pop",
        default="NCN",
        help="Population label used in output file names and plot title.",
    )
    parser.add_argument(
        "--maf-threshold",
        type=float,
        default=0.01,
        help="Minor-allele-frequency threshold used to define common variants.",
    )
    parser.add_argument(
        "--region-gap",
        type=int,
        default=200_000,
        help="Maximum gap in base pairs to merge consecutive significant SNPs into one region.",
    )
    parser.add_argument(
        "--plot-points-tsv",
        default=None,
        help=(
            "Existing Manhattan plot-points TSV to render from directly. "
            "When set, the script will not regenerate or overwrite that TSV."
        ),
    )
    parser.add_argument(
        "--bonferroni-threshold",
        type=float,
        default=None,
        help=(
            "Optional Bonferroni p-value threshold for the horizontal line. "
            "Useful with --plot-points-tsv when reusing an existing TSV."
        ),
    )
    parser.add_argument(
        "--plot-p-threshold",
        type=float,
        default=None,
        help=(
            "Optional p-value cutoff for writing Manhattan plot points. "
            "When omitted, all common variants with finite p-values are plotted."
        ),
    )
    parser.add_argument(
        "--exclude-regions-tsv",
        default=None,
        help=(
            "Optional TSV of regions to exclude before normalization and plotting. "
            "Expected columns: chrom start end."
        ),
    )
    return parser.parse_args()


def discover_input_files(input_dir: Path):
    files = []
    for path in input_dir.glob("chr*_*.sds.tsv"):
        stem = path.name.replace(".sds.tsv", "")
        if not (stem.endswith("_p") or stem.endswith("_q")):
            continue
        chrom = stem.split("_", 1)[0][3:]
        chrom_num = chrom_key(chrom)
        files.append((chrom_num, chrom, path))
    files.sort(key=lambda item: (item[0], item[2].name))
    return files


def chrom_key(chrom: str):
    if chrom.isdigit():
        return (0, int(chrom))
    if chrom == "X":
        return (1, 23)
    if chrom == "Y":
        return (1, 24)
    if chrom in {"MT", "M"}:
        return (1, 25)
    return (2, chrom)


def open_reader(path: Path):
    return pacsv.open_csv(
        path,
        read_options=pacsv.ReadOptions(block_size=1 << 20, use_threads=True),
        parse_options=pacsv.ParseOptions(delimiter="\t"),
        convert_options=pacsv.ConvertOptions(column_types=RESULT_SCHEMA),
    )


def load_exclude_regions(path: Path):
    regions = {}
    with path.open() as handle:
        for line_num, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if line_num == 1 and fields[0].lower() in {"chrom", "chr", "chromosome"}:
                continue
            if len(fields) < 3:
                raise SystemExit(f"Invalid exclude-regions line {line_num}: expected at least 3 tab-separated fields")
            chrom = normalize_chrom_label(fields[0])
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as exc:
                raise SystemExit(f"Invalid exclude-regions coordinates on line {line_num}: {raw_line.rstrip()}") from exc
            if start > end:
                raise SystemExit(f"Invalid exclude-regions interval on line {line_num}: start > end")
            regions.setdefault(chrom, []).append((start, end))
    for chrom in regions:
        regions[chrom].sort()
    return regions


def build_exclude_mask(chrom: str, positions: np.ndarray, exclude_regions: dict):
    spans = exclude_regions.get(normalize_chrom_label(chrom))
    if not spans or positions.size == 0:
        return np.zeros(positions.shape, dtype=bool)
    mask = np.zeros(positions.shape, dtype=bool)
    for start, end in spans:
        mask |= (positions >= start) & (positions <= end)
    return mask


def compute_maf(daf: np.ndarray):
    return np.minimum(daf, 1.0 - daf)


def safe_float(value):
    if value is None:
        return math.nan
    return float(value)


def safe_int(value):
    if value is None:
        return 0
    return int(value)


def maf_bin_start(maf: float, maf_threshold: float):
    if not math.isfinite(maf) or maf < maf_threshold:
        return None
    if maf < SPLIT_BIN_TRANSITION:
        bw = SPLIT_BIN_LOW_WIDTH
        raw_idx = int((maf - maf_threshold) / bw)
        return maf_threshold + raw_idx * bw
    else:
        bw = SPLIT_BIN_HIGH_WIDTH
        raw_idx = int((maf - SPLIT_BIN_TRANSITION) / bw)
        bin_start = SPLIT_BIN_TRANSITION + raw_idx * bw
        max_start = max(SPLIT_BIN_TRANSITION, 0.5 - bw)
        return min(bin_start, max_start)


def _split_bin_width_for_start(bin_start: float) -> float:
    """Return the bin width for a given bin_start in split-bin mode."""
    if bin_start < SPLIT_BIN_TRANSITION:
        return SPLIT_BIN_LOW_WIDTH
    return SPLIT_BIN_HIGH_WIDTH


def maf_bin_id(bin_start: float):
    bw = _split_bin_width_for_start(bin_start)
    if bin_start < SPLIT_BIN_TRANSITION:
        bin_end = min(bin_start + bw, SPLIT_BIN_TRANSITION)
    else:
        bin_end = min(bin_start + bw, 0.5)
    return f"[{bin_start:.2f},{bin_end:.2f})"


def first_pass(files, maf_threshold: float, exclude_regions: dict):
    bin_accumulators = {}
    common_variant_count = 0
    chrom_max_pos = {}
    chrom_counts = {}
    total_snvs = 0

    for _, chrom, path in files:
        reader = open_reader(path)
        for batch in reader:
            data = batch.to_pydict()
            rows = len(data["ID"])
            if rows == 0:
                continue

            rsds = np.asarray(data["rSDS"], dtype=np.float64)
            daf = np.asarray(data["DAF"], dtype=np.float64)
            pos = np.asarray(data["POS"], dtype=np.int64)
            exclude_mask = build_exclude_mask(chrom, pos, exclude_regions)

            valid_rsds = np.isfinite(rsds) & ~exclude_mask
            batch_total = int(np.count_nonzero(valid_rsds))
            if batch_total > 0:
                total_snvs += batch_total
                chrom_counts[chrom] = chrom_counts.get(chrom, 0) + batch_total

            if pos.size > 0:
                batch_max = int(pos.max())
                chrom_max_pos[chrom] = max(chrom_max_pos.get(chrom, 0), batch_max)

            valid_common = valid_rsds & np.isfinite(daf)
            if not np.any(valid_common):
                continue

            maf = compute_maf(daf[valid_common])
            common_mask = np.isfinite(maf) & (maf > maf_threshold)
            if not np.any(common_mask):
                continue

            common_rsds = rsds[valid_common][common_mask]
            common_maf = maf[common_mask]
            common_variant_count += common_rsds.size

            for rsds_value, maf_value in zip(common_rsds, common_maf):
                bin_start = maf_bin_start(float(maf_value), maf_threshold)
                if bin_start is None:
                    continue
                acc = bin_accumulators.setdefault(bin_start, [0.0, 0.0, 0])
                acc[0] += float(rsds_value)
                acc[1] += float(rsds_value * rsds_value)
                acc[2] += 1

    if total_snvs == 0:
        raise SystemExit("No valid SNP rows were found in the input files")
    if common_variant_count == 0:
        raise SystemExit(f"No common variants passed MAF > {maf_threshold}")

    bin_stats = {}
    for bin_start, (sum_rsds, sumsq_rsds, count) in sorted(bin_accumulators.items()):
        mean_common = sum_rsds / count
        variance_common = sumsq_rsds / count - mean_common * mean_common
        variance_common = max(variance_common, 0.0)
        sd_common = math.sqrt(variance_common)
        if not math.isfinite(sd_common) or sd_common <= 0.0:
            raise SystemExit(
                f"Common-variant SDS standard deviation is zero or non-finite in MAF bin {maf_bin_id(bin_start)}"
            )
        effective_bw = _split_bin_width_for_start(bin_start)
        bin_stats[bin_start] = BinStats(
            bin_id=maf_bin_id(bin_start),
            maf_start=bin_start,
            maf_end=min(bin_start + effective_bw, 0.5),
            mean_common=mean_common,
            sd_common=sd_common,
            common_variant_count=count,
        )

    return NormalizationStats(
        bin_stats=bin_stats,
        common_variant_count=common_variant_count,
        total_snvs=total_snvs,
        chrom_max_pos=chrom_max_pos,
        chrom_counts=chrom_counts,
    )


def chromosome_offsets(chrom_max_pos):
    ordered = sorted(chrom_max_pos.items(), key=lambda item: chrom_key(item[0]))
    offsets = {}
    centers = []
    running = 0
    for chrom, max_pos in ordered:
        offsets[chrom] = running
        centers.append((chrom, running + max_pos / 2.0))
        running += max_pos
    return offsets, centers, running


def second_pass(
    files,
    output_prefix: Path,
    maf_threshold: float,
    bin_stats: dict,
    bonferroni_threshold: float,
    chrom_offsets,
    region_gap: int,
    plot_p_threshold: float | None,
    exclude_regions: dict,
):
    normalized_tsv = output_prefix.with_name(output_prefix.name + ".normalized.tsv")
    stats_tsv = output_prefix.with_name(output_prefix.name + ".frequency_bins.tsv")
    plot_tsv = output_prefix.with_name(output_prefix.name + ".manhattan_points.tsv")
    regions_tsv = output_prefix.with_name(output_prefix.name + ".significant_regions.tsv")

    significant_hits = []
    plot_point_count = 0
    plot_max_neg_log10_p = 0.0

    with normalized_tsv.open("w", newline="") as norm_handle, plot_tsv.open(
        "w", newline=""
    ) as plot_handle:
        norm_writer = csv.writer(norm_handle, delimiter="\t")
        plot_writer = csv.writer(plot_handle, delimiter="\t")

        norm_writer.writerow(
            [
                "bin_id",
                "chr",
                "pos",
                "alt",
                "AA",
                "DA",
                "POS",
                "DAF",
                "MAF",
                "nG0",
                "nG1",
                "nG2",
                "rSDS",
                "SuggestedInitPoint",
                "BIN_MEAN",
                "BIN_SD",
                "norm_SDS",
                "pleft",
                "pright",
                "p_bothside",
                "ID",
                "is_common_variant",
            ]
        )
        plot_writer.writerow(["plot_x", "neg_log10_p", "chr_index", "is_significant"])

        chrom_to_index = {
            chrom: idx + 1
            for idx, (chrom, _) in enumerate(
                sorted(chrom_offsets.items(), key=lambda item: chrom_key(item[0]))
            )
        }

        for _, chrom, path in files:
            reader = open_reader(path)
            offset = chrom_offsets[chrom]
            chr_idx = chrom_to_index[chrom]
            for batch in reader:
                data = batch.to_pydict()
                rows = len(data["ID"])
                if rows == 0:
                    continue

                pos = np.asarray(data["POS"], dtype=np.int64)
                keep_mask = ~build_exclude_mask(chrom, pos, exclude_regions)
                if not np.any(keep_mask):
                    continue
                if not np.all(keep_mask):
                    keep_idx = np.flatnonzero(keep_mask)
                    data = {key: [values[i] for i in keep_idx] for key, values in data.items()}
                    rows = len(keep_idx)
                    pos = pos[keep_mask]

                daf = np.asarray(data["DAF"], dtype=np.float64)
                rsds = np.asarray(data["rSDS"], dtype=np.float64)

                maf = np.full(rows, np.nan, dtype=np.float64)
                valid_daf = np.isfinite(daf)
                maf[valid_daf] = compute_maf(daf[valid_daf])

                bin_ids = [""] * rows
                bin_means = np.full(rows, np.nan, dtype=np.float64)
                bin_sds = np.full(rows, np.nan, dtype=np.float64)
                for i in range(rows):
                    bin_start = maf_bin_start(float(maf[i]), maf_threshold)
                    if bin_start is None:
                        continue
                    stats = bin_stats.get(bin_start)
                    if stats is None:
                        continue
                    bin_ids[i] = stats.bin_id
                    bin_means[i] = stats.mean_common
                    bin_sds[i] = stats.sd_common

                norm_sds = np.full(rows, np.nan, dtype=np.float64)
                valid_norm = np.isfinite(rsds) & np.isfinite(bin_means) & np.isfinite(bin_sds) & (bin_sds > 0.0)
                norm_sds[valid_norm] = (rsds[valid_norm] - bin_means[valid_norm]) / bin_sds[valid_norm]

                pleft = np.full(rows, np.nan, dtype=np.float64)
                pright = np.full(rows, np.nan, dtype=np.float64)
                p_both = np.full(rows, np.nan, dtype=np.float64)
                pleft[valid_norm] = norm.cdf(norm_sds[valid_norm])
                pright[valid_norm] = norm.sf(norm_sds[valid_norm])
                p_both[valid_norm] = 2.0 * np.minimum(pleft[valid_norm], pright[valid_norm])
                p_both = np.minimum(p_both, 1.0)

                is_common = np.isfinite(maf) & (maf > maf_threshold)
                plot_daf_ok = np.isfinite(daf) & (daf > maf_threshold) & (daf < (1.0 - maf_threshold))
                plot_mask = is_common & plot_daf_ok & np.isfinite(p_both)
                if plot_p_threshold is not None:
                    plot_mask &= p_both < plot_p_threshold

                for i in range(rows):
                    norm_writer.writerow(
                        [
                            bin_ids[i],
                            chrom,
                            safe_int(pos[i]),
                            data["DA"][i],
                            data["AA"][i],
                            data["DA"][i],
                            safe_int(pos[i]),
                            safe_float(daf[i]),
                            "" if not np.isfinite(maf[i]) else f"{maf[i]:.10g}",
                            safe_int(data["nG0"][i]),
                            safe_int(data["nG1"][i]),
                            safe_int(data["nG2"][i]),
                            safe_float(rsds[i]),
                            data["SuggestedInitPoint"][i],
                            "" if not np.isfinite(bin_means[i]) else f"{bin_means[i]:.10g}",
                            "" if not np.isfinite(bin_sds[i]) else f"{bin_sds[i]:.10g}",
                            "" if not np.isfinite(norm_sds[i]) else f"{norm_sds[i]:.10g}",
                            "" if not np.isfinite(pleft[i]) else f"{pleft[i]:.10g}",
                            "" if not np.isfinite(pright[i]) else f"{pright[i]:.10g}",
                            "" if not np.isfinite(p_both[i]) else f"{p_both[i]:.10g}",
                            data["ID"][i],
                            int(is_common[i]),
                        ]
                    )

                if not np.any(plot_mask):
                    continue

                plot_positions = pos[plot_mask]
                plot_pvals = p_both[plot_mask]
                plot_daf = daf[plot_mask]
                plot_maf = maf[plot_mask]
                plot_norm = norm_sds[plot_mask]
                plot_ids = np.asarray(data["ID"], dtype=object)[plot_mask]

                neg_log10_p = -np.log10(np.maximum(plot_pvals, 1e-300))
                is_significant = plot_pvals < bonferroni_threshold
                plot_x = offset + plot_positions

                for x, y, sig in zip(plot_x, neg_log10_p, is_significant):
                    plot_writer.writerow([int(x), f"{y:.10g}", chr_idx, int(sig)])

                plot_point_count += int(plot_x.size)
                plot_max_neg_log10_p = max(plot_max_neg_log10_p, float(np.max(neg_log10_p)))

                for plot_pos, plot_id, daf_val, maf_val, norm_val, pval_val, sig in zip(
                    plot_positions,
                    plot_ids,
                    plot_daf,
                    plot_maf,
                    plot_norm,
                    plot_pvals,
                    is_significant,
                ):
                    if sig:
                        significant_hits.append(
                            (
                                chrom_key(chrom),
                                chrom,
                                int(plot_pos),
                                str(plot_id),
                                float(daf_val),
                                float(maf_val),
                                float(norm_val),
                                float(pval_val),
                            )
                        )

    with stats_tsv.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "bin_id",
                "maf_start",
                "maf_end",
                "common_variant_count",
                "common_variant_mean",
                "common_variant_sd",
            ]
        )
        for bin_start in sorted(bin_stats):
            stats = bin_stats[bin_start]
            writer.writerow(
                [
                    stats.bin_id,
                    f"{stats.maf_start:.10g}",
                    f"{stats.maf_end:.10g}",
                    stats.common_variant_count,
                    f"{stats.mean_common:.10g}",
                    f"{stats.sd_common:.10g}",
                ]
            )

    with regions_tsv.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "Region",
                "Sig_SNVs",
                "Top significant SNV ID",
                "AF",
                "MAF",
                "SDS",
                "Genes",
                "BonferroniThreshold",
                "TopP",
            ]
        )

        clusters = []
        current = None
        for _, chrom, pos, snp_id, af, maf, sds, pval in sorted(significant_hits):
            if current is None or chrom != current.chrom or pos - current.end > region_gap:
                current = RegionSummary(chrom, pos, pos)
                clusters.append(current)
            current.end = pos
            current.significant += 1
            if pval < current.top_hit.pval:
                current.top_hit = RegionTopHit(
                    snp_id=snp_id,
                    af=af,
                    maf=maf,
                    sds=sds,
                    pval=pval,
                )

        for summary in clusters:
            top_hit = summary.top_hit
            writer.writerow(
                [
                    f"chr{summary.chrom}:{summary.start}-{summary.end}",
                    summary.significant,
                    top_hit.snp_id,
                    "" if not np.isfinite(top_hit.af) else f"{top_hit.af:.10g}",
                    "" if not np.isfinite(top_hit.maf) else f"{top_hit.maf:.10g}",
                    "" if not np.isfinite(top_hit.sds) else f"{top_hit.sds:.10g}",
                    "",
                    f"{bonferroni_threshold:.10g}",
                    "" if not np.isfinite(top_hit.pval) else f"{top_hit.pval:.10g}",
                ]
            )

    return normalized_tsv, stats_tsv, plot_tsv, regions_tsv, plot_point_count, plot_max_neg_log10_p


def summarize_plot_tsv(plot_tsv: Path):
    plot_point_count = 0
    plot_max_neg_log10_p = 0.0
    chrom_min = {}
    chrom_max = {}
    min_sig_y = math.inf
    with plot_tsv.open("r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) < 4:
                continue
            plot_point_count += 1
            x = int(float(row[0]))
            y = float(row[1])
            chr_idx = int(row[2])
            is_significant = row[3] == "1"
            if y > plot_max_neg_log10_p:
                plot_max_neg_log10_p = y
            chrom_min[chr_idx] = x if chr_idx not in chrom_min else min(chrom_min[chr_idx], x)
            chrom_max[chr_idx] = x if chr_idx not in chrom_max else max(chrom_max[chr_idx], x)
            if is_significant and y < min_sig_y:
                min_sig_y = y

    chrom_centers = []
    genome_size = 0
    for chr_idx in sorted(chrom_min):
        chrom_centers.append((str(chr_idx), (chrom_min[chr_idx] + chrom_max[chr_idx]) / 2.0))
        genome_size = max(genome_size, chrom_max[chr_idx])

    return PlotTsvStats(
        plot_point_count=plot_point_count,
        plot_max_neg_log10_p=plot_max_neg_log10_p,
        chrom_centers=chrom_centers,
        genome_size=genome_size,
        bonferroni_threshold_y=min_sig_y if math.isfinite(min_sig_y) else 0.0,
    )


def main():
    args = parse_args()
    if not (0.0 < args.maf_threshold < 0.5):
        raise SystemExit("--maf-threshold must be between 0 and 0.5")
    if args.plot_p_threshold is not None and not (0.0 < args.plot_p_threshold <= 1.0):
        raise SystemExit("--plot-p-threshold must be between 0 and 1")
    if args.plot_points_tsv and args.exclude_regions_tsv:
        raise SystemExit("--exclude-regions-tsv cannot be combined with --plot-points-tsv")

    input_dir = Path(args.input_dir).resolve()
    output_prefix = (
        Path(args.output_prefix).resolve()
        if args.output_prefix
        else input_dir / args.pop
    )
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    exclude_regions_tsv = None if args.exclude_regions_tsv is None else Path(args.exclude_regions_tsv).resolve()
    exclude_regions = {} if exclude_regions_tsv is None else load_exclude_regions(exclude_regions_tsv)

    files = []
    stats = None
    normalized_tsv = None
    stats_tsv = None
    regions_tsv = None
    if args.plot_points_tsv:
        plot_tsv = Path(args.plot_points_tsv).resolve()
        if not plot_tsv.exists():
            raise SystemExit(f"--plot-points-tsv does not exist: {plot_tsv}")
        plot_stats = summarize_plot_tsv(plot_tsv)
        plot_point_count = plot_stats.plot_point_count
        plot_max_neg_log10_p = plot_stats.plot_max_neg_log10_p
        chrom_centers = plot_stats.chrom_centers
        genome_size = plot_stats.genome_size
        if args.bonferroni_threshold is not None:
            bonferroni_threshold = args.bonferroni_threshold
        elif plot_stats.bonferroni_threshold_y > 0.0:
            bonferroni_threshold = 10.0 ** (-plot_stats.bonferroni_threshold_y)
        else:
            bonferroni_threshold = 0.0
    else:
        files = discover_input_files(input_dir)
        if not files:
            raise SystemExit(f"No chr*_p/q.sds.tsv files found under {input_dir}")
        stats = first_pass(files, args.maf_threshold, exclude_regions)
        bonferroni_threshold = 0.05 / stats.common_variant_count
        offsets, chrom_centers, genome_size = chromosome_offsets(stats.chrom_max_pos)
        (
            normalized_tsv,
            stats_tsv,
            plot_tsv,
            regions_tsv,
            plot_point_count,
            plot_max_neg_log10_p,
        ) = second_pass(
            files=files,
            output_prefix=output_prefix,
            maf_threshold=args.maf_threshold,
            bin_stats=stats.bin_stats,
            bonferroni_threshold=bonferroni_threshold,
            chrom_offsets=offsets,
            region_gap=args.region_gap,
            plot_p_threshold=args.plot_p_threshold,
            exclude_regions=exclude_regions,
        )
    plot_path = None

    summary_path = output_prefix.with_name(output_prefix.name + ".postprocess_summary.tsv")
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["key", "value"])
        writer.writerow(["input_dir", str(input_dir)])
        writer.writerow(["population", args.pop])
        writer.writerow(["summary_mode", "plot_only_reuse" if args.plot_points_tsv else "full_recompute"])
        writer.writerow(["input_files", "" if stats is None else len(files)])
        writer.writerow(["total_snvs", "" if stats is None else stats.total_snvs])
        writer.writerow(["common_variant_count", "" if stats is None else stats.common_variant_count])
        writer.writerow(["maf_threshold", f"{args.maf_threshold:.10g}"])
        writer.writerow(["bin_mode", "split"])
        writer.writerow(["split_bin_low_width", f"{SPLIT_BIN_LOW_WIDTH:.10g}"])
        writer.writerow(["split_bin_high_width", f"{SPLIT_BIN_HIGH_WIDTH:.10g}"])
        writer.writerow(["split_bin_transition", f"{SPLIT_BIN_TRANSITION:.10g}"])
        writer.writerow(["plot_daf_min", f"{args.maf_threshold:.10g}"])
        writer.writerow(["plot_daf_max", f"{1.0 - args.maf_threshold:.10g}"])
        writer.writerow(["plot_maf_threshold", f"{args.maf_threshold:.10g}"])
        writer.writerow(["exclude_regions_tsv", "" if exclude_regions_tsv is None else str(exclude_regions_tsv)])
        writer.writerow(["plot_p_threshold", "" if args.plot_p_threshold is None else f"{args.plot_p_threshold:.10g}"])
        writer.writerow(
            [
                "plot_neg_log10_threshold",
                "" if args.plot_p_threshold is None else f"{-math.log10(args.plot_p_threshold):.10g}",
            ]
        )
        writer.writerow(["normalization_method", "maf_binned_common_variant"])
        writer.writerow(["normalization_bins", "" if stats is None else len(stats.bin_stats)])
        writer.writerow(["bonferroni_threshold", f"{bonferroni_threshold:.10g}"])
        writer.writerow(["normalized_tsv", "" if normalized_tsv is None else str(normalized_tsv)])
        writer.writerow(["frequency_bins_tsv", "" if stats_tsv is None else str(stats_tsv)])
        writer.writerow(["regions_tsv", "" if regions_tsv is None else str(regions_tsv)])
        writer.writerow(["plot_points_tsv", str(plot_tsv)])
        writer.writerow(["plot_points_written", plot_point_count])
        writer.writerow(["manhattan_plot", "" if plot_path is None else str(plot_path)])
        writer.writerow(["plot_backend", "none"])
        writer.writerow(
            [
                "summary_warning",
                "core normalization stats not recomputed in this run" if args.plot_points_tsv else "",
            ]
        )
        writer.writerow(
            [
                "plot_points_source",
                "existing_tsv" if args.plot_points_tsv else "regenerated_from_input",
            ]
        )
        if stats is not None:
            for chrom, count in sorted(stats.chrom_counts.items(), key=lambda item: chrom_key(item[0])):
                writer.writerow([f"chr{chrom}_snvs", count])

    print(f"normalized_tsv\t{'' if normalized_tsv is None else normalized_tsv}")
    print(f"frequency_bins_tsv\t{'' if stats_tsv is None else stats_tsv}")
    print(f"regions_tsv\t{'' if regions_tsv is None else regions_tsv}")
    print(f"plot_points_tsv\t{plot_tsv}")
    print(f"manhattan_plot\t{'' if plot_path is None else plot_path}")
    print(f"exclude_regions_tsv\t{'' if exclude_regions_tsv is None else exclude_regions_tsv}")
    print(f"summary_tsv\t{summary_path}")


if __name__ == "__main__":
    main()
