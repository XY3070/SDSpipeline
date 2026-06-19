#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
import subprocess
from collections import Counter
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PEAK_START = 143205259
DEFAULT_PEAK_END = 143533298
DEFAULT_TARGET_SNP = "chr1:143273067:C:T"


@dataclass(frozen=True)
class Window:
    label: str
    start: int
    end: int

    @property
    def region_str(self) -> str:
        return f"{self.start}-{self.end}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit SDS input alignment against VCF/sample order.")
    parser.add_argument("--chrom", default="1")
    parser.add_argument("--peak-start", type=int, default=DEFAULT_PEAK_START)
    parser.add_argument("--peak-end", type=int, default=DEFAULT_PEAK_END)
    parser.add_argument("--target-snp", default=DEFAULT_TARGET_SNP)
    parser.add_argument(
        "--normalized-table",
        default=str(REPO_ROOT / "tmp/ncn_postprocess_binned/NCN_binned.normalized.tsv"),
    )
    parser.add_argument(
        "--outdir",
        default=str(REPO_ROOT / "tmp" / f"input_order_audit_{datetime.now().strftime('%Y%m%d_%H%M%S')}"),
    )
    parser.add_argument("--common-per-window", type=int, default=10)
    parser.add_argument("--singleton-per-window", type=int, default=10)
    return parser.parse_args()


def clean_samples(path: Path) -> list[str]:
    samples: list[str] = []
    with path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            samples.append(stripped.split()[0])
    return samples


def run_cmd(cmd: list[str]) -> str:
    res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return res.stdout


def run_bcf_pipe(view_cmd: list[str], query_cmd: list[str]) -> str:
    p1 = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert p1.stdout is not None
    p2 = subprocess.run(query_cmd, stdin=p1.stdout, capture_output=True, text=True)
    p1.stdout.close()
    stderr1 = p1.stderr.read().decode("utf-8", errors="replace")
    rc1 = p1.wait()
    stderr2 = p2.stderr or ""
    allowed_view_codes = {0}
    if p2.returncode == 0:
        allowed_view_codes.add(-13)
    if rc1 not in allowed_view_codes or p2.returncode != 0:
        raise RuntimeError(
            "bcftools pipeline failed:\n"
            + "VIEW: "
            + " ".join(view_cmd)
            + f"\nrc1={rc1}\n{stderr1}\nQUERY: "
            + " ".join(query_cmd)
            + f"\nrc2={p2.returncode}\n{stderr2}"
        )
    return p2.stdout


def gt_to_sds_code(gt: str) -> str:
    if gt in {"0/0", "0|0"}:
        return "0"
    if gt in {"0/1", "1/0", "0|1", "1|0"}:
        return "1"
    if gt in {"1/1", "1|1"}:
        return "2"
    return "NA"


def normal_two_sided_neglog10(z_value: float) -> float:
    p_value = math.erfc(abs(z_value) / math.sqrt(2.0))
    if p_value <= 0.0:
        return float("inf")
    return -math.log10(p_value)


def read_t_rows_in_window(t_file: Path, start: int, end: int) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with t_file.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            pos = int(parts[3])
            if start <= pos <= end:
                rows.append(
                    {
                        "ID": parts[0],
                        "REF": parts[1],
                        "ALT": parts[2],
                        "POS": pos,
                        "gts": parts[4:],
                    }
                )
    rows.sort(key=lambda row: (row["POS"], row["ID"]))
    return rows


def estimate_daf(codes: list[str]) -> float | None:
    vals = [int(code) for code in codes if code != "NA"]
    if not vals:
        return None
    return sum(vals) / (2.0 * len(vals))


def pick_evenly(items: list[dict[str, object]], n: int, preferred_id: str | None = None) -> list[dict[str, object]]:
    if not items:
        return []
    picked: list[dict[str, object]] = []
    used_ids: set[str] = set()
    if preferred_id:
        preferred = next((item for item in items if item["ID"] == preferred_id), None)
        if preferred is not None:
            picked.append(preferred)
            used_ids.add(preferred["ID"])
    if len(items) <= n:
        for item in items:
            if item["ID"] not in used_ids:
                picked.append(item)
        return picked[:n]
    if n <= len(picked):
        return picked[:n]

    remaining = n - len(picked)
    step = (len(items) - 1) / max(remaining - 1, 1)
    for idx in range(remaining):
        item = items[round(idx * step)]
        if item["ID"] in used_ids:
            continue
        picked.append(item)
        used_ids.add(item["ID"])
        if len(picked) == n:
            break
    if len(picked) < n:
        for item in items:
            if item["ID"] in used_ids:
                continue
            picked.append(item)
            if len(picked) == n:
                break
    return picked


def query_sample_order(vcf_file: Path, sample_list: Path) -> list[str]:
    view_cmd = ["bcftools", "view", "-S", str(sample_list), "-Ou", str(vcf_file)]
    query_cmd = ["bcftools", "query", "-l"]
    output = run_bcf_pipe(view_cmd, query_cmd)
    return [line.strip() for line in output.splitlines() if line.strip()]


def query_vcf_rows(vcf_file: Path, sample_list: Path, chrom: str, start: int, end: int) -> dict[str, dict[str, object]]:
    regions = [f"chr{chrom}:{start}-{end}", f"{chrom}:{start}-{end}"]
    seen_regions: set[str] = set()
    for region in regions:
        if region in seen_regions:
            continue
        seen_regions.add(region)
        view_cmd = [
            "bcftools",
            "view",
            "-r",
            region,
            "-S",
            str(sample_list),
            "-m2",
            "-M2",
            "-v",
            "snps",
            "--force-samples",
            "-Ou",
            str(vcf_file),
        ]
        query_cmd = ["bcftools", "query", "-f", "%ID\t%REF\t%ALT\t%POS[\t%GT]\n"]
        output = run_bcf_pipe(view_cmd, query_cmd)
        if not output.strip():
            continue
        rows: dict[str, dict[str, object]] = {}
        for line in output.splitlines():
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            rows[parts[0]] = {
                "ID": parts[0],
                "REF": parts[1],
                "ALT": parts[2],
                "POS": int(parts[3]),
                "gt_strings": parts[4:],
                "gts": [gt_to_sds_code(gt) for gt in parts[4:]],
            }
        if rows:
            return rows
    raise RuntimeError(f"No VCF rows found for chr{chrom}:{start}-{end} in {vcf_file}")


def find_control_window(normalized_table: Path, chrom: str, peak_window: Window) -> Window:
    window_len = peak_window.end - peak_window.start + 1
    rows: list[tuple[int, float]] = []
    with normalized_table.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("chr") != chrom:
                continue
            if row.get("is_common_variant") not in {"1", "1.0", "true", "True"}:
                continue
            pos = int(row["POS"])
            p_text = row.get("p_bothside", "").strip()
            if p_text:
                p_value = float(p_text)
                neglog10 = float("inf") if p_value <= 0.0 else -math.log10(p_value)
            else:
                neglog10 = normal_two_sided_neglog10(float(row["norm_SDS"]))
            rows.append((pos, neglog10))
    rows.sort()
    if not rows:
        raise RuntimeError(f"No common variants found in normalized table {normalized_table}")

    best: tuple[float, int, int] | None = None
    right = 0
    max_queue: list[tuple[float, int]] = []
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
        count = right - left
        if count < 150:
            continue
        if not (end_pos < peak_window.start or start_pos > peak_window.end):
            continue
        max_signal = max_queue[0][0] if max_queue else float("inf")
        candidate = (max_signal, start_pos, end_pos)
        if best is None or candidate < best:
            best = candidate
    if best is None:
        fallback_start = peak_window.end + 5_000_000
        return Window("control", fallback_start, fallback_start + window_len - 1)
    return Window("control", best[1], best[2])


def build_sfile_carriers(s_file: Path, samples: list[str], target_positions: set[int]) -> dict[int, list[str]]:
    carriers: dict[int, list[str]] = defaultdict(list)
    with s_file.open() as handle:
        for sample, line in zip(samples, handle):
            parts = line.strip().split("\t")
            if not parts or parts == ["NA"]:
                continue
            for token in parts:
                if token == "NA":
                    continue
                pos = int(token)
                if pos in target_positions:
                    carriers[pos].append(sample)
    return carriers


def compare_vectors(expected: list[str], observed: list[str], samples: list[str]) -> tuple[int, str]:
    mismatch_count = 0
    first_mismatch = ""
    for sample, left, right in zip(samples, expected, observed):
        if left != right:
            mismatch_count += 1
            if not first_mismatch:
                first_mismatch = f"{sample}:{left}!={right}"
    if len(expected) != len(observed):
        mismatch_count += abs(len(expected) - len(observed))
        if not first_mismatch:
            first_mismatch = f"length:{len(expected)}!={len(observed)}"
    return mismatch_count, first_mismatch


def load_saved_order(order_file: Path) -> list[str]:
    return [line.strip() for line in order_file.read_text().splitlines() if line.strip()]


def audit_population(
    pop: str,
    chrom: str,
    windows: list[Window],
    target_snp: str,
    common_per_window: int,
    singleton_per_window: int,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    sample_list = REPO_ROOT / "data" / f"{pop}.txt"
    vcf_file = REPO_ROOT / "data" / "vcf" / pop / f"UKBQC_{pop}_chr{chrom}.vcf.gz"
    input_dir = REPO_ROOT / "data" / "processed" / "sds_input" / pop
    s_file = input_dir / f"chr{chrom}_s_file.txt"
    o_file = input_dir / f"chr{chrom}_o_file.txt"
    samples = clean_samples(sample_list)
    vcf_samples = [line.strip() for line in run_cmd(["bcftools", "query", "-l", str(vcf_file)]).splitlines() if line.strip()]

    results: list[dict[str, object]] = []

    missing_in_vcf = sorted(set(samples) - set(vcf_samples))
    extra_in_vcf = sorted(set(vcf_samples) - set(samples))
    results.append(
        {
            "population": pop,
            "chromosome": chrom,
            "window": "all",
            "check_type": "sample_set",
            "locus_id": "",
            "status": "pass" if not missing_in_vcf and not extra_in_vcf else "fail",
            "mismatch_count": len(missing_in_vcf) + len(extra_in_vcf),
            "first_mismatch": (missing_in_vcf[:1] or extra_in_vcf[:1] or [""])[0],
            "detail": f"missing_in_vcf={len(missing_in_vcf)} extra_in_vcf={len(extra_in_vcf)}",
        }
    )

    subset_order = query_sample_order(vcf_file, sample_list)
    mismatch_count, first_mismatch = compare_vectors(samples, subset_order, samples)
    results.append(
        {
            "population": pop,
            "chromosome": chrom,
            "window": "all",
            "check_type": "sample_order_bcftools",
            "locus_id": "",
            "status": "pass" if mismatch_count == 0 else "fail",
            "mismatch_count": mismatch_count,
            "first_mismatch": first_mismatch,
            "detail": f"sample_count={len(samples)}",
        }
    )

    order_passes = 0
    for chr_num in range(1, 23):
        order_file = input_dir / "tmp" / f"chr{chr_num}_order_s.txt"
        if not order_file.exists():
            continue
        saved_order = load_saved_order(order_file)
        mm_count, first_mm = compare_vectors(samples, saved_order, samples)
        status = "pass" if mm_count == 0 else "fail"
        if status == "pass":
            order_passes += 1
        results.append(
            {
                "population": pop,
                "chromosome": str(chr_num),
                "window": "all",
                "check_type": "saved_order_file",
                "locus_id": "",
                "status": status,
                "mismatch_count": mm_count,
                "first_mismatch": first_mm,
                "detail": str(order_file),
            }
        )

    for chr_num in range(1, 23):
        ofile = input_dir / f"chr{chr_num}_o_file.txt"
        if not ofile.exists():
            continue
        vector = ofile.read_text().strip().split("\t")
        bad_values = [value for value in vector if value != "1"]
        status = "pass" if len(vector) == len(samples) and not bad_values else "fail"
        results.append(
            {
                "population": pop,
                "chromosome": str(chr_num),
                "window": "all",
                "check_type": "o_file_vector",
                "locus_id": "",
                "status": status,
                "mismatch_count": abs(len(vector) - len(samples)) + len(bad_values),
                "first_mismatch": bad_values[0] if bad_values else "",
                "detail": f"length={len(vector)} expected={len(samples)}",
            }
        )

    summary: dict[str, object] = {
        "population": pop,
        "sample_count": len(samples),
        "saved_order_passes": order_passes,
        "windows": {},
    }

    for window in windows:
        t_file = input_dir / f"chr{chrom}_t_file.txt"
        t_rows = read_t_rows_in_window(t_file, window.start, window.end)
        vcf_rows = query_vcf_rows(vcf_file, sample_list, chrom, window.start, window.end)

        pos_counts = Counter(int(row["POS"]) for row in t_rows)
        common_candidates = [row for row in t_rows if (daf := estimate_daf(row["gts"])) is not None and 0.05 <= daf <= 0.95]
        singleton_candidates: list[dict[str, object]] = []
        for row in t_rows:
            vcf_row = vcf_rows.get(str(row["ID"]))
            if vcf_row is None:
                continue
            codes = list(vcf_row["gts"])
            missing_fraction = codes.count("NA") / max(len(codes), 1)
            if sum(int(gt) for gt in codes if gt != "NA") != 1:
                continue
            if "2" in codes:
                continue
            if missing_fraction > 0.005:
                continue
            if pos_counts[int(row["POS"])] != 1:
                continue
            singleton_candidates.append(row)
        common_loci = pick_evenly(common_candidates, common_per_window, preferred_id=target_snp if pop == "NCN" and window.label == "peak" else None)
        singleton_loci = pick_evenly(singleton_candidates, singleton_per_window)

        target_positions = {int(row["POS"]) for row in singleton_loci}
        sfile_carriers = build_sfile_carriers(s_file, samples, target_positions)

        for row in common_loci:
            vcf_row = vcf_rows.get(str(row["ID"]))
            if vcf_row is None:
                results.append(
                    {
                        "population": pop,
                        "chromosome": chrom,
                        "window": window.label,
                        "check_type": "t_genotype_alignment",
                        "locus_id": row["ID"],
                        "status": "fail",
                        "mismatch_count": len(samples),
                        "first_mismatch": "missing_in_vcf_query",
                        "detail": f"pos={row['POS']}",
                    }
                )
                continue
            mm_count, first_mm = compare_vectors(row["gts"], vcf_row["gts"], samples)
            results.append(
                {
                    "population": pop,
                    "chromosome": chrom,
                    "window": window.label,
                    "check_type": "t_genotype_alignment",
                    "locus_id": row["ID"],
                    "status": "pass" if mm_count == 0 else "fail",
                    "mismatch_count": mm_count,
                    "first_mismatch": first_mm,
                    "detail": f"pos={row['POS']}",
                }
            )

        for row in singleton_loci:
            vcf_row = vcf_rows.get(str(row["ID"]))
            if vcf_row is None:
                results.append(
                    {
                        "population": pop,
                        "chromosome": chrom,
                        "window": window.label,
                        "check_type": "s_singleton_alignment",
                        "locus_id": row["ID"],
                        "status": "fail",
                        "mismatch_count": 1,
                        "first_mismatch": "missing_in_vcf_query",
                        "detail": f"pos={row['POS']}",
                    }
                )
                continue
            carriers = [sample for sample, gt in zip(samples, vcf_row["gts"]) if gt == "1"]
            sfile_hits = sfile_carriers.get(int(row["POS"]), [])
            status = "pass" if carriers == sfile_hits and len(carriers) == 1 else "fail"
            results.append(
                {
                    "population": pop,
                    "chromosome": chrom,
                    "window": window.label,
                    "check_type": "s_singleton_alignment",
                    "locus_id": row["ID"],
                    "status": status,
                    "mismatch_count": 0 if status == "pass" else abs(len(carriers) - len(sfile_hits)) + 1,
                    "first_mismatch": carriers[0] if carriers else (sfile_hits[0] if sfile_hits else ""),
                    "detail": f"vcf_carriers={','.join(carriers)} sfile_hits={','.join(sfile_hits)} pos={row['POS']}",
                }
            )

        summary["windows"][window.label] = {
            "start": window.start,
            "end": window.end,
            "t_rows": len(t_rows),
            "common_candidates": len(common_candidates),
            "singleton_candidates": len(singleton_candidates),
            "common_audited": len(common_loci),
            "singleton_audited": len(singleton_loci),
        }

    return results, summary


def write_tsv(path: Path, rows: Iterable[dict[str, object]]) -> None:
    rows = list(rows)
    if not rows:
        return
    fieldnames = [
        "population",
        "chromosome",
        "window",
        "check_type",
        "locus_id",
        "status",
        "mismatch_count",
        "first_mismatch",
        "detail",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_summary(path: Path, peak_window: Window, control_window: Window, summaries: list[dict[str, object]], rows: list[dict[str, object]]) -> None:
    failing = [row for row in rows if row["status"] != "pass"]
    conclusion = (
        "No direct evidence of sample/genotype order mismatch was found in the audited SDS inputs."
        if not failing
        else "At least one direct input-alignment check failed; inspect the failing rows below before blaming gamma/history."
    )
    lines = [
        "# SDS Input Audit Summary",
        "",
        f"- Peak window: `chr1:{peak_window.start}-{peak_window.end}`",
        f"- Control window: `chr1:{control_window.start}-{control_window.end}`",
        f"- Overall conclusion: {conclusion}",
        f"- Total checks: `{len(rows)}`",
        f"- Failing checks: `{len(failing)}`",
        "",
        "## Population Summary",
        "",
    ]
    for summary in summaries:
        lines.append(f"### {summary['population']}")
        lines.append("")
        lines.append(f"- Sample count: `{summary['sample_count']}`")
        lines.append(f"- Saved `chr*_order_s.txt` passes: `{summary['saved_order_passes']}`")
        for label, info in summary["windows"].items():
            lines.append(
                f"- {label}: `chr1:{info['start']}-{info['end']}`, "
                f"`t_rows={info['t_rows']}`, `common_audited={info['common_audited']}`, `singleton_audited={info['singleton_audited']}`"
            )
        lines.append("")

    lines.extend(
        [
            "## Failing Checks",
            "",
        ]
    )
    if failing:
        lines.append("| Population | Chr | Window | Check | Locus | Mismatches | First mismatch | Detail |")
        lines.append("| --- | ---: | --- | --- | --- | ---: | --- | --- |")
        for row in failing[:50]:
            lines.append(
                f"| {row['population']} | {row['chromosome']} | {row['window']} | {row['check_type']} | {row['locus_id']} | "
                f"{row['mismatch_count']} | {row['first_mismatch']} | {row['detail']} |"
            )
    else:
        lines.append("No failing checks.")
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    peak_window = Window("peak", args.peak_start, args.peak_end)
    control_window = find_control_window(Path(args.normalized_table), args.chrom, peak_window)
    windows = [peak_window, control_window]

    all_rows: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []
    for pop in ("NCN", "SCN"):
        rows, summary = audit_population(
            pop=pop,
            chrom=args.chrom,
            windows=windows,
            target_snp=args.target_snp,
            common_per_window=args.common_per_window,
            singleton_per_window=args.singleton_per_window,
        )
        all_rows.extend(rows)
        summaries.append(summary)

    write_tsv(outdir / "audit_checks.tsv", all_rows)
    write_summary(outdir / "summary.md", peak_window, control_window, summaries, all_rows)
    print(outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
