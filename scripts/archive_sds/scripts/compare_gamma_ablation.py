#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.stats import norm, spearmanr


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_POP = "NCN"
DEFAULT_CHROM = "1"
DEFAULT_START = 143205259
DEFAULT_END = 143533298
DEFAULT_OLD_G = str(REPO_ROOT / "g_file.txt")
DEFAULT_NEW_G = "/data/home/grp-wangyf/xuyuan/ms/scripts/sds_input.gamma_shapes.NCN.final"
DEFAULT_NORMALIZED = str(REPO_ROOT / "tmp/ncn_postprocess_binned/NCN_binned.normalized.tsv")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare old vs phlash gamma effects with the current SDS compute script.")
    parser.add_argument("--pop", default=DEFAULT_POP)
    parser.add_argument("--chrom", default=DEFAULT_CHROM)
    parser.add_argument("--start", type=int, default=DEFAULT_START)
    parser.add_argument("--end", type=int, default=DEFAULT_END)
    parser.add_argument("--old-g-file", default=DEFAULT_OLD_G)
    parser.add_argument("--new-g-file", default=DEFAULT_NEW_G)
    parser.add_argument("--normalized-table", default=DEFAULT_NORMALIZED)
    parser.add_argument("--init", default="0.00001")
    parser.add_argument("--s-file-ncol", default="20000")
    parser.add_argument(
        "--outdir",
        default=str(REPO_ROOT / "tmp" / f"gamma_ablation_{datetime.now().strftime('%Y%m%d_%H%M%S')}"),
    )
    return parser.parse_args()


def safe_float(text: str | None) -> float | None:
    if text is None:
        return None
    value = text.strip()
    if not value or value.upper() == "NA":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def safe_int(text: str | None) -> int | None:
    value = safe_float(text)
    if value is None:
        return None
    return int(value)


def format_metric(value: float | None, digits: int = 6) -> str:
    if value is None or not math.isfinite(value):
        return "NA"
    return f"{value:.{digits}f}"


def pearson_corr(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 2:
        return None
    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    if np.std(x) == 0.0 or np.std(y) == 0.0:
        return None
    return float(np.corrcoef(x, y)[0, 1])


def spearman_corr(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 2:
        return None
    corr, _ = spearmanr(xs, ys)
    if corr is None or not math.isfinite(float(corr)):
        return None
    return float(corr)


def normal_two_sided_metrics(z_value: float | None) -> tuple[float | None, float | None]:
    if z_value is None or not math.isfinite(z_value):
        return None, None
    log_p = math.log(2.0) + float(norm.logsf(abs(z_value)))
    neg_log10_p = -log_p / math.log(10.0)
    if log_p < -745:
        p_value = 0.0
    else:
        p_value = math.exp(log_p)
    return p_value, neg_log10_p


def filter_t_region(source_t: Path, dest_t: Path, start: int, end: int) -> int:
    rows = 0
    with source_t.open() as src, dest_t.open("w") as dst:
        for line in src:
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            pos = int(parts[3])
            if start <= pos <= end:
                dst.write(line)
                rows += 1
    return rows


def run_compute(
    s_file: Path,
    t_file: Path,
    o_file: Path,
    b_file: Path,
    g_file: Path,
    init: str,
    s_file_ncol: str,
    output_tsv: Path,
    summary_csv: Path,
    cache_dir: Path,
) -> None:
    cmd = [
        sys.executable,
        str(REPO_ROOT / "scripts" / "compute_SDS.py"),
        str(s_file),
        str(t_file),
        str(o_file),
        str(b_file),
        str(g_file),
        init,
        s_file_ncol,
        "--output",
        str(output_tsv),
        "--summary-csv",
        str(summary_csv),
        "--pickle-cache-dir",
        str(cache_dir),
    ]
    subprocess.run(cmd, check=True, cwd=REPO_ROOT)


def load_sds_rows(path: Path) -> dict[str, dict[str, object]]:
    rows: dict[str, dict[str, object]] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            snp_id = row.get("ID", "")
            if not snp_id:
                continue
            rows[snp_id] = {
                "ID": snp_id,
                "POS": safe_int(row.get("POS")),
                "DAF": safe_float(row.get("DAF")),
                "nG0": safe_int(row.get("nG0")),
                "nG1": safe_int(row.get("nG1")),
                "nG2": safe_int(row.get("nG2")),
                "rSDS": safe_float(row.get("rSDS")),
                "SuggestedInitPoint": row.get("SuggestedInitPoint", ""),
            }
    return rows


def load_current_reference(pop: str, chrom: str, start: int, end: int) -> dict[str, dict[str, object]]:
    out_root = REPO_ROOT / "data/processed/sds_output" / pop
    candidates = [
        out_root / f"chr{chrom}_p.sds.tsv",
        out_root / f"chr{chrom}_q.sds.tsv",
        out_root / f"chr{chrom}.sds.tsv",
    ]
    rows: dict[str, dict[str, object]] = {}
    for path in candidates:
        if not path.exists():
            continue
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                pos = safe_int(row.get("POS"))
                snp_id = row.get("ID", "")
                if not snp_id or pos is None:
                    continue
                if start <= pos <= end:
                    rows[snp_id] = {
                        "ID": snp_id,
                        "POS": pos,
                        "rSDS": safe_float(row.get("rSDS")),
                    }
    return rows


def load_normalized_rows(path: Path, chrom: str, start: int, end: int) -> dict[str, dict[str, object]]:
    rows: dict[str, dict[str, object]] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("chr") != chrom:
                continue
            if row.get("is_common_variant") not in {"1", "1.0", "true", "True"}:
                continue
            pos = safe_int(row.get("POS"))
            snp_id = row.get("ID", "")
            if pos is None or not snp_id or not (start <= pos <= end):
                continue
            current_norm = safe_float(row.get("norm_SDS"))
            current_p, current_neg = normal_two_sided_metrics(current_norm)
            rows[snp_id] = {
                "ID": snp_id,
                "POS": pos,
                "BIN_MEAN": safe_float(row.get("BIN_MEAN")),
                "BIN_SD": safe_float(row.get("BIN_SD")),
                "current_norm_SDS": current_norm,
                "current_NEG_LOG10_P": current_neg,
                "current_p_bothside": current_p,
            }
    return rows


def determine_verdict(raw_corr: float | None, retained_fraction: float | None) -> tuple[str, str]:
    if raw_corr is not None and retained_fraction is not None:
        if raw_corr >= 0.95 and retained_fraction >= 0.8:
            return (
                "gamma-unlikely",
                "Replacing the phlash gamma with the old g_file barely changes the chr1 peak, so gamma/history is unlikely to be the main explanation.",
            )
        if raw_corr < 0.8 or retained_fraction < 0.5:
            return (
                "gamma-plausible",
                "Replacing the gamma file substantially alters the chr1 peak, so population-history / gamma remains a plausible explanation.",
            )
    return (
        "ambiguous",
        "The old-vs-new gamma comparison changes the peak partially, but not enough for a clean attribution.",
    )


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: ("" if value is None else value) for key, value in row.items()})


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    s_file = REPO_ROOT / "data/processed/sds_input" / args.pop / f"chr{args.chrom}_s_file.txt"
    t_file = REPO_ROOT / "data/processed/sds_input" / args.pop / f"chr{args.chrom}_t_file.txt"
    o_file = REPO_ROOT / "data/processed/sds_input" / args.pop / f"chr{args.chrom}_o_file.txt"
    b_file = REPO_ROOT / "data/processed/sds_input" / args.pop / f"chr{args.chrom}_b_file.txt"
    old_g_file = Path(args.old_g_file)
    new_g_file = Path(args.new_g_file)
    normalized_table = Path(args.normalized_table)

    region_t = outdir / "chr1_peak.t.tsv"
    old_out = outdir / "old_gamma.sds.tsv"
    new_out = outdir / "new_gamma.sds.tsv"
    current_ref_out = outdir / "current_reference.sds.tsv"
    joined_out = outdir / "gamma_joined.tsv"
    summary_out = outdir / "summary.md"

    rows_in_region = filter_t_region(t_file, region_t, args.start, args.end)
    if rows_in_region == 0:
        raise RuntimeError(f"No t_file rows in region chr{args.chrom}:{args.start}-{args.end}")

    run_compute(
        s_file=s_file,
        t_file=region_t,
        o_file=o_file,
        b_file=b_file,
        g_file=old_g_file,
        init=args.init,
        s_file_ncol=args.s_file_ncol,
        output_tsv=old_out,
        summary_csv=outdir / "old_gamma.summary.csv",
        cache_dir=outdir / "cache_old",
    )
    run_compute(
        s_file=s_file,
        t_file=region_t,
        o_file=o_file,
        b_file=b_file,
        g_file=new_g_file,
        init=args.init,
        s_file_ncol=args.s_file_ncol,
        output_tsv=new_out,
        summary_csv=outdir / "new_gamma.summary.csv",
        cache_dir=outdir / "cache_new",
    )

    old_rows = load_sds_rows(old_out)
    new_rows = load_sds_rows(new_out)
    current_ref = load_current_reference(args.pop, args.chrom, args.start, args.end)
    normalized_rows = load_normalized_rows(normalized_table, args.chrom, args.start, args.end)

    ref_rows = []
    for snp_id, row in current_ref.items():
        ref_rows.append({"ID": snp_id, "POS": row["POS"], "current_reference_rSDS": row["rSDS"]})
    if ref_rows:
        write_tsv(current_ref_out, ["ID", "POS", "current_reference_rSDS"], sorted(ref_rows, key=lambda row: (row["POS"], row["ID"])))

    overlap_ids = sorted(
        set(old_rows) & set(new_rows) & set(normalized_rows) & set(current_ref),
        key=lambda snp_id: (normalized_rows[snp_id]["POS"], snp_id),
    )
    if not overlap_ids:
        raise RuntimeError("No overlapping rows between old/new gamma outputs and current reference")

    joined_rows: list[dict[str, object]] = []
    old_raw: list[float] = []
    new_raw: list[float] = []
    old_proj: list[float] = []
    new_proj: list[float] = []
    ref_delta_values: list[float] = []
    for snp_id in overlap_ids:
        old_row = old_rows[snp_id]
        new_row = new_rows[snp_id]
        norm_row = normalized_rows[snp_id]
        ref_row = current_ref[snp_id]
        bin_mean = norm_row["BIN_MEAN"]
        bin_sd = norm_row["BIN_SD"]
        old_proj_z = None
        new_proj_z = None
        if (
            bin_mean is not None
            and bin_sd is not None
            and bin_sd > 0.0
            and old_row["rSDS"] is not None
            and new_row["rSDS"] is not None
        ):
            old_proj_z = (old_row["rSDS"] - bin_mean) / bin_sd
            new_proj_z = (new_row["rSDS"] - bin_mean) / bin_sd
        old_p, old_neg = normal_two_sided_metrics(old_proj_z)
        new_p, new_neg = normal_two_sided_metrics(new_proj_z)
        ref_delta = None
        if new_row["rSDS"] is not None and ref_row["rSDS"] is not None:
            ref_delta = new_row["rSDS"] - ref_row["rSDS"]
            ref_delta_values.append(abs(ref_delta))
        if old_row["rSDS"] is not None and new_row["rSDS"] is not None:
            old_raw.append(old_row["rSDS"])
            new_raw.append(new_row["rSDS"])
        if old_proj_z is not None and new_proj_z is not None:
            old_proj.append(old_proj_z)
            new_proj.append(new_proj_z)
        joined_rows.append(
            {
                "ID": snp_id,
                "POS": norm_row["POS"],
                "old_gamma_rSDS": old_row["rSDS"],
                "new_gamma_rSDS": new_row["rSDS"],
                "current_reference_rSDS": ref_row["rSDS"],
                "new_vs_reference_delta_rSDS": ref_delta,
                "BIN_MEAN_current": bin_mean,
                "BIN_SD_current": bin_sd,
                "old_gamma_projected_norm_SDS": old_proj_z,
                "new_gamma_projected_norm_SDS": new_proj_z,
                "old_gamma_NEG_LOG10_P_projected": old_neg,
                "new_gamma_NEG_LOG10_P_projected": new_neg,
                "current_NEG_LOG10_P_reference": norm_row["current_NEG_LOG10_P"],
                "old_gamma_p_projected": old_p,
                "new_gamma_p_projected": new_p,
            }
        )

    joined_rows.sort(
        key=lambda row: (
            -(row["current_NEG_LOG10_P_reference"] if row["current_NEG_LOG10_P_reference"] is not None else -1.0),
            row["POS"],
            row["ID"],
        )
    )
    write_tsv(
        joined_out,
        [
            "ID",
            "POS",
            "old_gamma_rSDS",
            "new_gamma_rSDS",
            "current_reference_rSDS",
            "new_vs_reference_delta_rSDS",
            "BIN_MEAN_current",
            "BIN_SD_current",
            "old_gamma_projected_norm_SDS",
            "new_gamma_projected_norm_SDS",
            "old_gamma_NEG_LOG10_P_projected",
            "new_gamma_NEG_LOG10_P_projected",
            "current_NEG_LOG10_P_reference",
            "old_gamma_p_projected",
            "new_gamma_p_projected",
        ],
        joined_rows,
    )

    raw_pearson = pearson_corr(old_raw, new_raw)
    raw_spearman = spearman_corr(old_raw, new_raw)
    proj_pearson = pearson_corr(old_proj, new_proj)
    proj_spearman = spearman_corr(old_proj, new_proj)
    current_extreme = [
        row for row in joined_rows if row["current_NEG_LOG10_P_reference"] is not None and row["current_NEG_LOG10_P_reference"] >= 20.0
    ]
    retained = [
        row
        for row in current_extreme
        if row["old_gamma_NEG_LOG10_P_projected"] is not None
        and row["old_gamma_NEG_LOG10_P_projected"] >= 20.0
        and row["old_gamma_projected_norm_SDS"] is not None
        and row["new_gamma_projected_norm_SDS"] is not None
        and math.copysign(1.0, row["old_gamma_projected_norm_SDS"]) == math.copysign(1.0, row["new_gamma_projected_norm_SDS"])
    ]
    retained_fraction = len(retained) / len(current_extreme) if current_extreme else None
    verdict, interpretation = determine_verdict(raw_pearson, retained_fraction)

    target_snp = "chr1:143273067:C:T"
    target = next((row for row in joined_rows if row["ID"] == target_snp), None)
    max_ref_delta = max(ref_delta_values) if ref_delta_values else None

    lines = [
        "# Gamma Ablation Summary",
        "",
        f"- Population: `{args.pop}`",
        f"- Region: `chr{args.chrom}:{args.start}-{args.end}`",
        f"- Old gamma: `{old_g_file}`",
        f"- New gamma: `{new_g_file}`",
        f"- Verdict: `{verdict}`",
        f"- Interpretation: {interpretation}",
        "",
        "## Counts",
        "",
        f"- Filtered region rows: `{rows_in_region}`",
        f"- Old gamma rows: `{len(old_rows)}`",
        f"- New gamma rows: `{len(new_rows)}`",
        f"- Joined common-variant rows: `{len(joined_rows)}`",
        "",
        "## Correlations",
        "",
        f"- Raw `rSDS` Pearson: `{format_metric(raw_pearson)}`",
        f"- Raw `rSDS` Spearman: `{format_metric(raw_spearman)}`",
        f"- Projected standardized SDS Pearson: `{format_metric(proj_pearson)}`",
        f"- Projected standardized SDS Spearman: `{format_metric(proj_spearman)}`",
        f"- New gamma vs current-reference max abs delta `rSDS`: `{format_metric(max_ref_delta)}`",
        "",
        "## Extreme Peak Retention",
        "",
        f"- Current extreme SNPs (`NEG_LOG10_P >= 20`): `{len(current_extreme)}`",
        f"- Old gamma retains extreme signal with same sign: `{len(retained)}`",
        f"- Retained fraction: `{format_metric(retained_fraction)}`",
        "",
        "## Target SNP",
        "",
    ]
    if target is not None:
        lines.extend(
            [
                f"- Target SNP: `{target_snp}`",
                f"- Old gamma `rSDS`: `{format_metric(target['old_gamma_rSDS'], 4)}`",
                f"- New gamma `rSDS`: `{format_metric(target['new_gamma_rSDS'], 4)}`",
                f"- Old gamma projected `-log10(p)`: `{format_metric(target['old_gamma_NEG_LOG10_P_projected'], 4)}`",
                f"- New gamma projected `-log10(p)`: `{format_metric(target['new_gamma_NEG_LOG10_P_projected'], 4)}`",
            ]
        )
    else:
        lines.append(f"- Target SNP `{target_snp}` not present in joined rows.")

    lines.extend(
        [
            "",
            "## Top 20 Current Extreme SNPs",
            "",
            "| ID | POS | Old gamma rSDS | New gamma rSDS | Old gamma -log10(p) | New gamma -log10(p) |",
            "| --- | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in joined_rows[:20]:
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["ID"]),
                    str(row["POS"]),
                    format_metric(row["old_gamma_rSDS"], 4),
                    format_metric(row["new_gamma_rSDS"], 4),
                    format_metric(row["old_gamma_NEG_LOG10_P_projected"], 4),
                    format_metric(row["new_gamma_NEG_LOG10_P_projected"], 4),
                ]
            )
            + " |"
        )

    summary_out.write_text("\n".join(lines) + "\n")
    print(outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
