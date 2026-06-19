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
DEFAULT_CHROM = "1"


@dataclass(frozen=True)
class WindowSpec:
    label: str
    role: str
    start: int
    end: int

    @property
    def region(self) -> str:
        return f"chr{DEFAULT_CHROM}:{self.start}-{self.end}"

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare chr1 legacy pre-split SDS inputs against current NCN inputs.")
    parser.add_argument("--chrom", default=DEFAULT_CHROM)
    parser.add_argument(
        "--current-input-root",
        default=str(SDS_ROOT / "data" / "processed" / "sds_input" / "NCN"),
    )
    parser.add_argument(
        "--legacy-input-root",
        default=str(SDS_ROOT / "data" / "processed" / "sds_input" / "legacy"),
    )
    parser.add_argument(
        "--normalized-table",
        default=str(
            SDS_ROOT / "data" / "processed" / "sds_output_olddefault_mainline" / "NCN" / "NCN.normalized.tsv"
        ),
    )
    parser.add_argument(
        "--top-windows",
        default=str(
            SDS_ROOT / "data" / "processed" / "sds_diagnostic_comparison" / "NCN_gamma_tracks.top_windows.tsv"
        ),
    )
    parser.add_argument(
        "--legacy-sds",
        default=str(SDS_ROOT / "data" / "processed" / "sds_input" / "legacy" / "chr1_sds_res_FULL.txt"),
    )
    parser.add_argument(
        "--current-sds",
        default=str(SDS_ROOT / "data" / "processed" / "sds_output_olddefault_mainline" / "NCN" / "chr1.sds.tsv"),
    )
    parser.add_argument(
        "--outdir",
        default=str(SDS_ROOT / "tmp" / f"chr1_legacy_vs_current_{datetime.now().strftime('%Y%m%d_%H%M%S')}"),
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
            genotype_text = parts[4]
            for window in windows:
                if window.start <= pos <= window.end:
                    rows[window.label].append(
                        {
                            "ID": parts[0],
                            "REF": parts[1],
                            "ALT": parts[2],
                            "POS": pos,
                            "GENO": genotype_text,
                        }
                    )
                    break
    return rows


def missing_and_maf_summary(rows: list[dict[str, object]]) -> dict[str, object]:
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


def compare_t_site_sets(current_rows: list[dict[str, object]], legacy_rows: list[dict[str, object]]) -> dict[str, object]:
    current_ids = {row["ID"] for row in current_rows}
    legacy_ids = {row["ID"] for row in legacy_rows}
    current_pos = {row["POS"] for row in current_rows}
    legacy_pos = {row["POS"] for row in legacy_rows}
    shared_ids = current_ids & legacy_ids
    return {
        "current_rows": len(current_rows),
        "legacy_rows": len(legacy_rows),
        "shared_ids": len(shared_ids),
        "current_only_ids": len(current_ids - legacy_ids),
        "legacy_only_ids": len(legacy_ids - current_ids),
        "shared_positions": len(current_pos & legacy_pos),
        "current_only_positions": len(current_pos - legacy_pos),
        "legacy_only_positions": len(legacy_pos - current_pos),
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


def singleton_distribution_summary(counts: list[int]) -> dict[str, object]:
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


def load_peak_window(top_windows_path: Path, rank: str, label: str, role: str) -> WindowSpec:
    with top_windows_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("track") != "olddefault":
                continue
            if row.get("chr") != DEFAULT_CHROM:
                continue
            if row.get("track_rank") != rank:
                continue
            return WindowSpec(label, role, int(row["window_start"]), int(row["window_end"]))
    raise RuntimeError(f"Unable to find chr{DEFAULT_CHROM} rank {rank} window in {top_windows_path}")


def pick_matched_control(normalized_table: Path, target: WindowSpec, exclude: list[WindowSpec]) -> WindowSpec:
    candidate = find_control_window(normalized_table, DEFAULT_CHROM, Window("target", target.start, target.end))
    overlaps = any(not (candidate.end < item.start or candidate.start > item.end) for item in exclude)
    if not overlaps:
        return WindowSpec("C2_matched", "matched_control", candidate.start, candidate.end)

    rows: list[tuple[int, float]] = []
    with normalized_table.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("chr") != DEFAULT_CHROM:
                continue
            if row.get("is_common_variant") not in {"1", "1.0", "true", "True"}:
                continue
            pos = int(row["POS"])
            p_text = row.get("p_bothside", "").strip()
            if not p_text:
                continue
            p_value = float(p_text)
            neglog10 = float("inf") if p_value <= 0.0 else -math.log10(p_value)
            rows.append((pos, neglog10))
    rows.sort()
    best: tuple[float, int, int] | None = None
    right = 0
    max_queue: list[tuple[float, int]] = []
    window_len = target.length
    for left in range(0, len(rows), 20):
        start_pos = rows[left][0]
        end_pos = start_pos + window_len - 1
        while right < len(rows) and rows[right][0] <= end_pos:
            neglog10 = rows[right][1]
            while max_queue and max_queue[-1][0] <= neglog10:
                max_queue.pop()
            max_queue.append((neglog10, right))
            right += 1
        while max_queue and max_queue[0][1] < left:
            max_queue.pop(0)
        if right - left < 150:
            continue
        if any(not (end_pos < item.start or start_pos > item.end) for item in exclude):
            continue
        max_signal = max_queue[0][0] if max_queue else float("inf")
        candidate_tuple = (max_signal, start_pos, end_pos)
        if best is None or candidate_tuple < best:
            best = candidate_tuple
    if best is None:
        raise RuntimeError("Unable to choose matched control window")
    return WindowSpec("C2_matched", "matched_control", best[1], best[2])


def load_sds_window_stats(sds_path: Path, windows: list[WindowSpec]) -> dict[str, dict[str, object]]:
    stats = {window.label: {"rows": 0, "max_abs_rSDS": 0.0, "top_id": "", "top_pos": 0} for window in windows}
    with sds_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            pos = int(row["POS"])
            rsds = abs(float(row["rSDS"]))
            for window in windows:
                if window.start <= pos <= window.end:
                    stats[window.label]["rows"] += 1
                    if rsds > stats[window.label]["max_abs_rSDS"]:
                        stats[window.label]["max_abs_rSDS"] = rsds
                        stats[window.label]["top_id"] = row["ID"]
                        stats[window.label]["top_pos"] = pos
                    break
    return stats


def write_table(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, metadata: dict[str, object], windows: list[WindowSpec], rows: list[dict[str, object]]) -> None:
    with path.open("w") as handle:
        handle.write("# chr1 Legacy vs Current Input Audit Summary\n\n")
        handle.write("## Metadata\n\n")
        for key, value in metadata.items():
            handle.write(f"- {key}: `{value}`\n")
        handle.write("\n## Window Manifest\n\n")
        for window in windows:
            handle.write(f"- `{window.label}` `{window.role}`: `{window.region}`\n")
        handle.write("\n## Window Results\n\n")
        handle.write("| window | role | current_t_rows | legacy_t_rows | shared_ids | current_only_ids | legacy_only_ids | current_avg_missing | legacy_avg_missing | current_avg_maf | legacy_avg_maf | current_total_singletons | legacy_total_singletons | current_top10_share | legacy_top10_share | current_max_abs_rSDS | legacy_max_abs_rSDS |\n")
        handle.write("| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n")
        for row in rows:
            handle.write(
                f"| {row['window']} | {row['role']} | {row['current_t_rows']} | {row['legacy_t_rows']} | "
                f"{row['shared_t_ids']} | {row['current_only_t_ids']} | {row['legacy_only_t_ids']} | "
                f"{row['current_avg_missing_rate']:.6f} | {row['legacy_avg_missing_rate']:.6f} | "
                f"{row['current_avg_folded_maf']:.6f} | {row['legacy_avg_folded_maf']:.6f} | "
                f"{row['current_total_singletons']} | {row['legacy_total_singletons']} | "
                f"{row['current_top10_share']:.4f} | {row['legacy_top10_share']:.4f} | "
                f"{row['current_max_abs_rSDS']:.4f} | {row['legacy_max_abs_rSDS']:.4f} |\n"
            )
        handle.write("\n## Interpretation Guide\n\n")
        handle.write("- `s_file` is expected to differ because singleton status depends on cohort size.\n")
        handle.write("- The main audit targets are whether peak windows show abnormal singleton concentration, unusual missingness, or shifted site composition relative to controls.\n")


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    current_input_root = Path(args.current_input_root).resolve()
    legacy_input_root = Path(args.legacy_input_root).resolve()
    normalized_table = Path(args.normalized_table).resolve()
    top_windows = Path(args.top_windows).resolve()
    current_sds = Path(args.current_sds).resolve()
    legacy_sds = Path(args.legacy_sds).resolve()

    current_t = current_input_root / f"chr{args.chrom}_t_file.txt"
    current_s = current_input_root / f"chr{args.chrom}_s_file.txt"
    current_o = current_input_root / f"chr{args.chrom}_o_file.txt"
    current_b = current_input_root / f"chr{args.chrom}_b_file.txt"

    legacy_t = legacy_input_root / f"chr{args.chrom}_t_file.txt"
    legacy_s = legacy_input_root / f"chr{args.chrom}_s_file.txt"
    legacy_o = legacy_input_root / f"chr{args.chrom}_o_file.txt"
    legacy_b = legacy_input_root / f"chr{args.chrom}_b_file.txt"

    w1 = load_peak_window(top_windows, "1", "W1_peak", "peak")
    w2 = load_peak_window(top_windows, "2", "W2_secondary", "secondary_peak")
    c1 = WindowSpec("C1_quiet", "control", 35735220, 36063259)
    c2 = pick_matched_control(normalized_table, w2, [w1, w2, c1])
    windows = [w1, w2, c1, c2]

    current_t_rows = read_t_window_rows(current_t, windows)
    legacy_t_rows = read_t_window_rows(legacy_t, windows)
    current_sds_stats = load_sds_window_stats(current_sds, windows)
    legacy_sds_stats = load_sds_window_stats(legacy_sds, windows)

    table_rows: list[dict[str, object]] = []
    for window in windows:
        current_t_summary = missing_and_maf_summary(current_t_rows[window.label])
        legacy_t_summary = missing_and_maf_summary(legacy_t_rows[window.label])
        t_cmp = compare_t_site_sets(current_t_rows[window.label], legacy_t_rows[window.label])

        current_singletons = singleton_distribution_summary(count_singletons_in_window(current_s, window))
        legacy_singletons = singleton_distribution_summary(count_singletons_in_window(legacy_s, window))

        table_rows.append(
            {
                "window": window.label,
                "role": window.role,
                "region": window.region,
                "current_t_rows": t_cmp["current_rows"],
                "legacy_t_rows": t_cmp["legacy_rows"],
                "shared_t_ids": t_cmp["shared_ids"],
                "current_only_t_ids": t_cmp["current_only_ids"],
                "legacy_only_t_ids": t_cmp["legacy_only_ids"],
                "shared_positions": t_cmp["shared_positions"],
                "current_only_positions": t_cmp["current_only_positions"],
                "legacy_only_positions": t_cmp["legacy_only_positions"],
                "current_avg_missing_rate": current_t_summary["avg_missing_rate"],
                "legacy_avg_missing_rate": legacy_t_summary["avg_missing_rate"],
                "current_max_missing_rate": current_t_summary["max_missing_rate"],
                "legacy_max_missing_rate": legacy_t_summary["max_missing_rate"],
                "current_avg_alt_af": current_t_summary["avg_alt_af"],
                "legacy_avg_alt_af": legacy_t_summary["avg_alt_af"],
                "current_avg_folded_maf": current_t_summary["avg_folded_maf"],
                "legacy_avg_folded_maf": legacy_t_summary["avg_folded_maf"],
                "current_maf_lt_0_01_rows": current_t_summary["maf_lt_0_01_rows"],
                "legacy_maf_lt_0_01_rows": legacy_t_summary["maf_lt_0_01_rows"],
                "current_sample_count": current_singletons["sample_count"],
                "legacy_sample_count": legacy_singletons["sample_count"],
                "current_samples_with_singletons": current_singletons["samples_with_singletons"],
                "legacy_samples_with_singletons": legacy_singletons["samples_with_singletons"],
                "current_total_singletons": current_singletons["total_singletons"],
                "legacy_total_singletons": legacy_singletons["total_singletons"],
                "current_mean_singletons_per_sample": current_singletons["mean_singletons_per_sample"],
                "legacy_mean_singletons_per_sample": legacy_singletons["mean_singletons_per_sample"],
                "current_max_singletons_in_sample": current_singletons["max_singletons_in_sample"],
                "legacy_max_singletons_in_sample": legacy_singletons["max_singletons_in_sample"],
                "current_top10_share": current_singletons["top10_share"],
                "legacy_top10_share": legacy_singletons["top10_share"],
                "current_top50_share": current_singletons["top50_share"],
                "legacy_top50_share": legacy_singletons["top50_share"],
                "current_sds_rows": current_sds_stats[window.label]["rows"],
                "legacy_sds_rows": legacy_sds_stats[window.label]["rows"],
                "current_max_abs_rSDS": current_sds_stats[window.label]["max_abs_rSDS"],
                "legacy_max_abs_rSDS": legacy_sds_stats[window.label]["max_abs_rSDS"],
                "current_top_id": current_sds_stats[window.label]["top_id"],
                "legacy_top_id": legacy_sds_stats[window.label]["top_id"],
            }
        )

    metadata = {
        "current_t_fields": count_fields_first_line(current_t),
        "legacy_t_fields": count_fields_first_line(legacy_t),
        "current_o_fields": count_fields_first_line(current_o),
        "legacy_o_fields": count_fields_first_line(legacy_o),
        "current_b_lines": sum(1 for _ in current_b.open()),
        "legacy_b_lines": sum(1 for _ in legacy_b.open()),
        "current_s_rows": sum(1 for _ in current_s.open()),
        "legacy_s_rows": sum(1 for _ in legacy_s.open()),
        "current_sds_rows_total": sum(1 for _ in current_sds.open()) - 1,
        "legacy_sds_rows_total": sum(1 for _ in legacy_sds.open()) - 1,
    }

    write_table(outdir / "window_summary.tsv", table_rows)
    write_summary(outdir / "summary.md", metadata, windows, table_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
