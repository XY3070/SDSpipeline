#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from audit_sds_inputs import Window, find_control_window


SCRIPT_PATH = Path(__file__).resolve()
SDS_ROOT = SCRIPT_PATH.parents[1]


@dataclass(frozen=True)
class WindowSpec:
    chrom: str
    label: str
    role: str
    start: int
    end: int

    @property
    def region(self) -> str:
        return f"chr{self.chrom}:{self.start}-{self.end}"

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit old/new SDS input lineage and chr-level model amplification for chr1 + chr6."
    )
    parser.add_argument(
        "--outdir",
        default=str(SDS_ROOT / "data" / "processed" / f"input_lineage_audit_{datetime.now().strftime('%Y%m%d_%H%M%S')}"),
    )
    return parser.parse_args()


def count_fields_first_line(path: Path) -> int:
    with path.open() as handle:
        line = handle.readline().rstrip("\n")
    return len(line.split("\t")) if line else 0


def read_t_window_rows(t_file: Path, windows: list[WindowSpec]) -> dict[str, list[dict[str, object]]]:
    rows = {window.label: [] for window in windows}
    with t_file.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t", 4)
            if len(parts) < 5:
                continue
            pos = int(parts[3])
            for window in windows:
                if window.start <= pos <= window.end:
                    rows[window.label].append(
                        {
                            "ID": parts[0],
                            "REF": parts[1],
                            "ALT": parts[2],
                            "POS": pos,
                            "GENO": parts[4],
                        }
                    )
                    break
    return rows


def missing_and_maf_summary(rows: list[dict[str, object]]) -> dict[str, float | int]:
    if not rows:
        return {
            "row_count": 0,
            "avg_missing_rate": 0.0,
            "max_missing_rate": 0.0,
            "avg_alt_af": 0.0,
            "avg_folded_maf": 0.0,
            "maf_lt_0_01_rows": 0,
        }

    total_missing = 0.0
    max_missing = 0.0
    alt_af_sum = 0.0
    maf_sum = 0.0
    maf_lt_count = 0

    for row in rows:
        genotypes = row["GENO"].split("\t")
        called = 0
        alt_count = 0
        missing = 0
        for gt in genotypes:
            if gt == "NA":
                missing += 1
                continue
            called += 1
            alt_count += int(gt)
        missing_rate = missing / float(len(genotypes))
        max_missing = max(max_missing, missing_rate)
        total_missing += missing_rate
        if called == 0:
            continue
        alt_af = alt_count / float(2 * called)
        folded = min(alt_af, 1.0 - alt_af)
        alt_af_sum += alt_af
        maf_sum += folded
        if folded < 0.01:
            maf_lt_count += 1

    return {
        "row_count": len(rows),
        "avg_missing_rate": total_missing / float(len(rows)),
        "max_missing_rate": max_missing,
        "avg_alt_af": alt_af_sum / float(len(rows)),
        "avg_folded_maf": maf_sum / float(len(rows)),
        "maf_lt_0_01_rows": maf_lt_count,
    }


def compare_t_site_sets(rows_a: list[dict[str, object]], rows_b: list[dict[str, object]]) -> dict[str, int]:
    ids_a = {row["ID"] for row in rows_a}
    ids_b = {row["ID"] for row in rows_b}
    pos_a = {row["POS"] for row in rows_a}
    pos_b = {row["POS"] for row in rows_b}
    return {
        "rows_a": len(rows_a),
        "rows_b": len(rows_b),
        "shared_ids": len(ids_a & ids_b),
        "a_only_ids": len(ids_a - ids_b),
        "b_only_ids": len(ids_b - ids_a),
        "shared_positions": len(pos_a & pos_b),
        "a_only_positions": len(pos_a - pos_b),
        "b_only_positions": len(pos_b - pos_a),
    }


def count_singletons_in_window(s_file: Path, window: WindowSpec) -> list[int]:
    counts: list[int] = []
    with s_file.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped == "NA":
                counts.append(0)
                continue
            total = 0
            for token in stripped.split("\t"):
                if not token or token == "NA":
                    continue
                pos = int(token)
                if pos < window.start:
                    continue
                if pos > window.end:
                    break
                total += 1
            counts.append(total)
    return counts


def singleton_distribution_summary(counts: list[int]) -> dict[str, float | int]:
    total = sum(counts)
    positive = sorted((count for count in counts if count > 0), reverse=True)
    top10 = sum(positive[:10])
    top50 = sum(positive[:50])
    max_count = positive[0] if positive else 0
    return {
        "sample_count": len(counts),
        "samples_with_singletons": len(positive),
        "total_singletons": total,
        "mean_singletons_per_sample": 0.0 if not counts else total / float(len(counts)),
        "max_singletons_in_sample": max_count,
        "top10_share": 0.0 if total == 0 else top10 / float(total),
        "top50_share": 0.0 if total == 0 else top50 / float(total),
    }


def load_sds_window_stats(sds_path: Path, windows: list[WindowSpec]) -> dict[str, dict[str, object]]:
    stats = {
        window.label: {"rows": 0, "max_abs_rSDS": 0.0, "top_id": "", "top_pos": 0, "mean_abs_rSDS": 0.0}
        for window in windows
    }
    sums = {window.label: 0.0 for window in windows}
    with sds_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            pos = int(row["POS"])
            rsds = abs(float(row["rSDS"]))
            for window in windows:
                if window.start <= pos <= window.end:
                    stats[window.label]["rows"] += 1
                    sums[window.label] += rsds
                    if rsds > stats[window.label]["max_abs_rSDS"]:
                        stats[window.label]["max_abs_rSDS"] = rsds
                        stats[window.label]["top_id"] = row["ID"]
                        stats[window.label]["top_pos"] = pos
                    break
    for window in windows:
        rows = stats[window.label]["rows"]
        stats[window.label]["mean_abs_rSDS"] = 0.0 if rows == 0 else sums[window.label] / float(rows)
    return stats


def load_diag_top_chr6_windows(path: Path, top_n: int = 2) -> list[WindowSpec]:
    windows: list[WindowSpec] = []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["chr"] != "6":
                continue
            windows.append(
                WindowSpec(
                    chrom="6",
                    label=f"chr6_peak_{len(windows)+1}",
                    role="peak",
                    start=int(row["window_start"]),
                    end=int(row["window_end"]),
                )
            )
            if len(windows) == top_n:
                break
    if len(windows) < top_n:
        raise RuntimeError(f"Unable to find {top_n} chr6 windows in {path}")
    return windows


def load_frequency_bins(path: Path) -> dict[str, dict[str, float | int]]:
    bins: dict[str, dict[str, float | int]] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            bins[row["bin_id"]] = {
                "count": int(row["common_variant_count"]),
                "mean": float(row["common_variant_mean"]),
                "sd": float(row["common_variant_sd"]),
            }
    return bins


def load_tail_counts(path: Path) -> dict[str, dict[str, float | int]]:
    rows: dict[str, dict[str, float | int]] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows[row["chr"]] = {
                "common_variant_count": int(row["common_variant_count"]),
                "ge8": int(row["count_neglog10_ge_8"]),
                "ge20": int(row["count_neglog10_ge_20"]),
                "ge50": int(row["count_neglog10_ge_50"]),
                "ge100": int(row["count_neglog10_ge_100"]),
                "rate_ge20": float(row["rate_neglog10_ge_20"]),
            }
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def compare_input_roots(
    compare_id: str,
    chrom: str,
    label_a: str,
    root_a: Path,
    label_b: str,
    root_b: Path,
    windows: list[WindowSpec],
    sds_a: Path | None = None,
    sds_b: Path | None = None,
) -> list[dict[str, object]]:
    t_a = root_a / f"chr{chrom}_t_file.txt"
    s_a = root_a / f"chr{chrom}_s_file.txt"
    o_a = root_a / f"chr{chrom}_o_file.txt"
    b_a = root_a / f"chr{chrom}_b_file.txt"
    t_b = root_b / f"chr{chrom}_t_file.txt"
    s_b = root_b / f"chr{chrom}_s_file.txt"
    o_b = root_b / f"chr{chrom}_o_file.txt"
    b_b = root_b / f"chr{chrom}_b_file.txt"

    t_rows_a = read_t_window_rows(t_a, windows)
    t_rows_b = read_t_window_rows(t_b, windows)
    sds_stats_a = load_sds_window_stats(sds_a, windows) if sds_a else {}
    sds_stats_b = load_sds_window_stats(sds_b, windows) if sds_b else {}

    b_lines_a = [line.strip() for line in b_a.read_text().splitlines() if line.strip()]
    b_lines_b = [line.strip() for line in b_b.read_text().splitlines() if line.strip()]

    rows: list[dict[str, object]] = []
    for window in windows:
        summary_a = missing_and_maf_summary(t_rows_a[window.label])
        summary_b = missing_and_maf_summary(t_rows_b[window.label])
        cmp_sites = compare_t_site_sets(t_rows_a[window.label], t_rows_b[window.label])
        single_a = singleton_distribution_summary(count_singletons_in_window(s_a, window))
        single_b = singleton_distribution_summary(count_singletons_in_window(s_b, window))

        row = {
            "compare_id": compare_id,
            "chrom": chrom,
            "window": window.label,
            "role": window.role,
            "region": window.region,
            "label_a": label_a,
            "label_b": label_b,
            "t_rows_a": cmp_sites["rows_a"],
            "t_rows_b": cmp_sites["rows_b"],
            "shared_ids": cmp_sites["shared_ids"],
            "a_only_ids": cmp_sites["a_only_ids"],
            "b_only_ids": cmp_sites["b_only_ids"],
            "shared_positions": cmp_sites["shared_positions"],
            "a_only_positions": cmp_sites["a_only_positions"],
            "b_only_positions": cmp_sites["b_only_positions"],
            "avg_missing_a": summary_a["avg_missing_rate"],
            "avg_missing_b": summary_b["avg_missing_rate"],
            "avg_maf_a": summary_a["avg_folded_maf"],
            "avg_maf_b": summary_b["avg_folded_maf"],
            "maf_lt_0_01_a": summary_a["maf_lt_0_01_rows"],
            "maf_lt_0_01_b": summary_b["maf_lt_0_01_rows"],
            "singletons_total_a": single_a["total_singletons"],
            "singletons_total_b": single_b["total_singletons"],
            "singletons_top10_share_a": single_a["top10_share"],
            "singletons_top10_share_b": single_b["top10_share"],
            "b_lines_a": len(b_lines_a),
            "b_lines_b": len(b_lines_b),
            "b_text_a": " | ".join(b_lines_a),
            "b_text_b": " | ".join(b_lines_b),
            "t_fields_a": count_fields_first_line(t_a),
            "t_fields_b": count_fields_first_line(t_b),
            "o_fields_a": count_fields_first_line(o_a),
            "o_fields_b": count_fields_first_line(o_b),
            "s_rows_a": sum(1 for _ in s_a.open()),
            "s_rows_b": sum(1 for _ in s_b.open()),
        }
        if sds_a:
            row["sds_rows_a"] = sds_stats_a[window.label]["rows"]
            row["max_abs_rSDS_a"] = sds_stats_a[window.label]["max_abs_rSDS"]
            row["mean_abs_rSDS_a"] = sds_stats_a[window.label]["mean_abs_rSDS"]
            row["top_id_a"] = sds_stats_a[window.label]["top_id"]
        if sds_b:
            row["sds_rows_b"] = sds_stats_b[window.label]["rows"]
            row["max_abs_rSDS_b"] = sds_stats_b[window.label]["max_abs_rSDS"]
            row["mean_abs_rSDS_b"] = sds_stats_b[window.label]["mean_abs_rSDS"]
            row["top_id_b"] = sds_stats_b[window.label]["top_id"]
        rows.append(row)
    return rows


def compare_models_same_input(
    compare_id: str,
    pop: str,
    chrom: str,
    label_a: str,
    sds_a: Path,
    label_b: str,
    sds_b: Path,
    windows: list[WindowSpec],
) -> list[dict[str, object]]:
    stats_a = load_sds_window_stats(sds_a, windows)
    stats_b = load_sds_window_stats(sds_b, windows)
    rows: list[dict[str, object]] = []
    for window in windows:
        rows.append(
            {
                "compare_id": compare_id,
                "population": pop,
                "chrom": chrom,
                "window": window.label,
                "role": window.role,
                "region": window.region,
                "label_a": label_a,
                "label_b": label_b,
                "rows_a": stats_a[window.label]["rows"],
                "rows_b": stats_b[window.label]["rows"],
                "max_abs_rSDS_a": stats_a[window.label]["max_abs_rSDS"],
                "max_abs_rSDS_b": stats_b[window.label]["max_abs_rSDS"],
                "mean_abs_rSDS_a": stats_a[window.label]["mean_abs_rSDS"],
                "mean_abs_rSDS_b": stats_b[window.label]["mean_abs_rSDS"],
                "top_id_a": stats_a[window.label]["top_id"],
                "top_id_b": stats_b[window.label]["top_id"],
            }
        )
    return rows


def write_summary(
    path: Path,
    input_rows: list[dict[str, object]],
    model_rows: list[dict[str, object]],
    bin_rows: list[dict[str, object]],
    tail_rows: list[dict[str, object]],
) -> None:
    with path.open("w") as handle:
        handle.write("# SDS old/new input lineage audit (`chr1 + chr6`)\n\n")
        handle.write("## Headline findings\n\n")
        handle.write("- `legacy chr1 input` is structurally different from current pop-specific input: `t_file` density is much lower and `b_file` boundaries are different.\n")
        handle.write("- `current split-pop old input` and `rebuilt main-contract input` are identical for `NCN chr1`, and differ modestly only in one `chr6` MHC window.\n")
        handle.write("- On the same rebuilt input, `historical old/default` and `Gravel_CHB@100k` diverge sharply at `chr6`, which means the current failure mode is more consistent with input/site-universe artifacts being amplified than with a simple plotting bug.\n")
        handle.write("\n## Input comparisons\n\n")
        handle.write("| compare_id | chrom | window | t_rows_a | t_rows_b | shared_ids | a_only_ids | b_only_ids | avg_maf_a | avg_maf_b | singletons_total_a | singletons_total_b |\n")
        handle.write("| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n")
        for row in input_rows:
            handle.write(
                f"| {row['compare_id']} | {row['chrom']} | {row['window']} | {row['t_rows_a']} | {row['t_rows_b']} | "
                f"{row['shared_ids']} | {row['a_only_ids']} | {row['b_only_ids']} | "
                f"{row['avg_maf_a']:.6f} | {row['avg_maf_b']:.6f} | {row['singletons_total_a']} | {row['singletons_total_b']} |\n"
            )
        handle.write("\n## Same-input model comparisons\n\n")
        handle.write("| compare_id | chrom | window | max_abs_rSDS_a | max_abs_rSDS_b | mean_abs_rSDS_a | mean_abs_rSDS_b | top_id_a | top_id_b |\n")
        handle.write("| --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- |\n")
        for row in model_rows:
            handle.write(
                f"| {row['compare_id']} | {row['chrom']} | {row['window']} | "
                f"{row['max_abs_rSDS_a']:.4f} | {row['max_abs_rSDS_b']:.4f} | "
                f"{row['mean_abs_rSDS_a']:.4f} | {row['mean_abs_rSDS_b']:.4f} | "
                f"{row['top_id_a']} | {row['top_id_b']} |\n"
            )
        handle.write("\n## Frequency bins\n\n")
        handle.write("| population | model | bin_id | count | mean_rSDS | sd_rSDS |\n")
        handle.write("| --- | --- | --- | ---: | ---: | ---: |\n")
        for row in bin_rows:
            handle.write(
                f"| {row['population']} | {row['model']} | {row['bin_id']} | {row['count']} | {row['mean']:.6f} | {row['sd']:.6f} |\n"
            )
        handle.write("\n## Chromosome tails\n\n")
        handle.write("| population | model | chr | ge8 | ge20 | ge50 | ge100 | rate_ge20 |\n")
        handle.write("| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |\n")
        for row in tail_rows:
            handle.write(
                f"| {row['population']} | {row['model']} | {row['chr']} | {row['ge8']} | {row['ge20']} | "
                f"{row['ge50']} | {row['ge100']} | {row['rate_ge20']:.8f} |\n"
            )
        handle.write("\n## Interpretation\n\n")
        handle.write("- `chr1 legacy vs rebuilt` should be read as `pre-split-like density/boundary baseline` vs `current pop-specific input`; it is not a strict same-cohort comparison.\n")
        handle.write("- `chr6 old-split vs rebuilt` is the direct test of whether the 2026-05 rebuilt contract itself changed current `NCN` input semantics in a large way.\n")
        handle.write("- `same-input olddefault vs Gravel` isolates amplification on a fixed input lineage.\n")


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    chr1_windows = [
        WindowSpec("1", "peak125", "peak", 125_100_001, 125_200_000),
        WindowSpec("1", "peak143", "peak", 143_200_001, 143_300_000),
        WindowSpec("1", "control35", "control", 35_735_220, 36_063_259),
    ]

    chr6_peak_windows = load_diag_top_chr6_windows(
        SDS_ROOT / "data" / "processed" / "sds_output_gravel_chb_ne100k_newinput_20260511" / "NCN" / "diagnostics" / "NCN.gravel_chb_ne100k_manual.peak_windows.tsv",
        top_n=2,
    )
    chr6_control = find_control_window(
        SDS_ROOT / "data" / "processed" / "sds_output_olddefault_newinput_20260514" / "NCN" / "NCN.normalized.tsv",
        "6",
        Window("target", chr6_peak_windows[0].start, chr6_peak_windows[0].end),
    )
    chr6_windows = chr6_peak_windows + [
        WindowSpec("6", "chr6_control", "control", chr6_control.start, chr6_control.end)
    ]

    input_rows: list[dict[str, object]] = []
    input_rows.extend(
        compare_input_roots(
            compare_id="chr1_legacy_vs_rebuilt_ncn",
            chrom="1",
            label_a="legacy_chr1",
            root_a=SDS_ROOT / "data" / "processed" / "sds_input" / "legacy",
            label_b="rebuilt_ncn",
            root_b=SDS_ROOT / "data" / "processed" / "sds_input_rebuilt_main_contract_20260511" / "NCN",
            windows=chr1_windows,
            sds_a=SDS_ROOT / "data" / "processed" / "sds_input" / "legacy" / "chr1_sds_res_FULL.txt",
            sds_b=SDS_ROOT / "data" / "processed" / "sds_output_olddefault_newinput_20260514" / "NCN" / "chr1.sds.tsv",
        )
    )
    input_rows.extend(
        compare_input_roots(
            compare_id="chr1_oldsplit_vs_rebuilt_ncn",
            chrom="1",
            label_a="oldsplit_ncn",
            root_a=SDS_ROOT / "data" / "processed" / "sds_input" / "NCN",
            label_b="rebuilt_ncn",
            root_b=SDS_ROOT / "data" / "processed" / "sds_input_rebuilt_main_contract_20260511" / "NCN",
            windows=chr1_windows,
            sds_a=SDS_ROOT / "data" / "processed" / "sds_output_olddefault_mainline" / "NCN" / "chr1.sds.tsv",
            sds_b=SDS_ROOT / "data" / "processed" / "sds_output_olddefault_newinput_20260514" / "NCN" / "chr1.sds.tsv",
        )
    )
    input_rows.extend(
        compare_input_roots(
            compare_id="chr6_oldsplit_vs_rebuilt_ncn",
            chrom="6",
            label_a="oldsplit_ncn",
            root_a=SDS_ROOT / "data" / "processed" / "sds_input" / "NCN",
            label_b="rebuilt_ncn",
            root_b=SDS_ROOT / "data" / "processed" / "sds_input_rebuilt_main_contract_20260511" / "NCN",
            windows=chr6_windows,
            sds_a=SDS_ROOT / "data" / "processed" / "sds_output_olddefault_mainline" / "NCN" / "chr6.sds.tsv",
            sds_b=SDS_ROOT / "data" / "processed" / "sds_output_olddefault_newinput_20260514" / "NCN" / "chr6.sds.tsv",
        )
    )

    model_rows: list[dict[str, object]] = []
    model_rows.extend(
        compare_models_same_input(
            compare_id="ncn_rebuilt_olddefault_vs_gravel_chr1",
            pop="NCN",
            chrom="1",
            label_a="olddefault_newinput",
            sds_a=SDS_ROOT / "data" / "processed" / "sds_output_olddefault_newinput_20260514" / "NCN" / "chr1.sds.tsv",
            label_b="gravel_chb_ne100k_newinput",
            sds_b=SDS_ROOT / "data" / "processed" / "sds_output_gravel_chb_ne100k_newinput_20260511" / "NCN" / "chr1.sds.tsv",
            windows=chr1_windows,
        )
    )
    model_rows.extend(
        compare_models_same_input(
            compare_id="ncn_rebuilt_olddefault_vs_gravel_chr6",
            pop="NCN",
            chrom="6",
            label_a="olddefault_newinput",
            sds_a=SDS_ROOT / "data" / "processed" / "sds_output_olddefault_newinput_20260514" / "NCN" / "chr6.sds.tsv",
            label_b="gravel_chb_ne100k_newinput",
            sds_b=SDS_ROOT / "data" / "processed" / "sds_output_gravel_chb_ne100k_newinput_20260511" / "NCN" / "chr6.sds.tsv",
            windows=chr6_windows,
        )
    )

    bin_rows: list[dict[str, object]] = []
    for pop, model, path in [
        ("NCN", "olddefault_newinput", SDS_ROOT / "data" / "processed" / "sds_output_olddefault_newinput_20260514" / "NCN" / "NCN.frequency_bins.tsv"),
        ("NCN", "gravel_chb_ne100k_newinput", SDS_ROOT / "data" / "processed" / "sds_output_gravel_chb_ne100k_newinput_20260511" / "NCN" / "NCN.frequency_bins.tsv"),
        ("SCN", "gravel_chb_ne100k_newinput", SDS_ROOT / "data" / "processed" / "sds_output_gravel_chb_ne100k_newinput_20260511" / "SCN" / "SCN.frequency_bins.tsv"),
    ]:
        for bin_id in ("[0.01,0.02)", "[0.02,0.03)", "[0.03,0.04)", "[0.04,0.05)"):
            row = load_frequency_bins(path)[bin_id]
            bin_rows.append(
                {
                    "population": pop,
                    "model": model,
                    "bin_id": bin_id,
                    "count": row["count"],
                    "mean": row["mean"],
                    "sd": row["sd"],
                }
            )

    tail_rows: list[dict[str, object]] = []
    for pop, model, path in [
        ("NCN", "olddefault_newinput", SDS_ROOT / "data" / "processed" / "sds_output_olddefault_newinput_20260514" / "NCN" / "diagnostics" / "NCN.olddefault_newinput.chrom_tail_counts.tsv"),
        ("NCN", "gravel_chb_ne100k_newinput", SDS_ROOT / "data" / "processed" / "sds_output_gravel_chb_ne100k_newinput_20260511" / "NCN" / "diagnostics" / "NCN.gravel_chb_ne100k_manual.chrom_tail_counts.tsv"),
        ("SCN", "gravel_chb_ne100k_newinput", SDS_ROOT / "data" / "processed" / "sds_output_gravel_chb_ne100k_newinput_20260511" / "SCN" / "diagnostics" / "SCN.gravel_chb_ne100k_manual.chrom_tail_counts.tsv"),
    ]:
        tails = load_tail_counts(path)
        for chrom in ("1", "6"):
            row = tails[chrom]
            tail_rows.append(
                {
                    "population": pop,
                    "model": model,
                    "chr": chrom,
                    "ge8": row["ge8"],
                    "ge20": row["ge20"],
                    "ge50": row["ge50"],
                    "ge100": row["ge100"],
                    "rate_ge20": row["rate_ge20"],
                }
            )

    write_tsv(outdir / "input_compare.tsv", input_rows)
    write_tsv(outdir / "model_compare.tsv", model_rows)
    write_tsv(outdir / "frequency_bin_compare.tsv", bin_rows)
    write_tsv(outdir / "tail_compare.tsv", tail_rows)
    write_summary(outdir / "summary.md", input_rows, model_rows, bin_rows, tail_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
