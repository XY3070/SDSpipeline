#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
MS_ROOT = Path("/data/home/grp-wangyf/xuyuan/ms")
PHLASH_ROOT = Path("/data/home/grp-wangyf/xuyuan/phlash")
PHLASH_PYTHON = Path("/data/home/grp-wangyf/intern/miniforge3/envs/phlash/bin/python")
SDS_PYTHON = Path("/data/home/grp-wangyf/intern/miniforge3/envs/sds/bin/python")
MS_MAKE_DIR = MS_ROOT / "scripts"
MS_BINARY = MS_ROOT / "msdir" / "ms"
BACKWARD_SCRIPT = MS_ROOT / "scripts" / "backward.py"
TARGET_GENERATIONS = (0, 1, 10, 50, 100, 500, 1000, 5000, 10000)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Scan gamma-shape / Ne sensitivity for a single SDS SNP using phlash posterior percentiles."
    )
    parser.add_argument("--pop", default="NCN", help="Population label. Defaults to NCN.")
    parser.add_argument("--chrom", default="4", help="Chromosome containing the target SNP.")
    parser.add_argument("--snp-id", default="chr4:99317841:T:C", help="Exact SNP ID in the t_file / SDS output.")
    parser.add_argument("--pos", type=int, default=99317841, help="Target SNP position.")
    parser.add_argument(
        "--input-root",
        default=None,
        help="Input root containing s/o/b/t files. Defaults to data/processed/sds_input/<POP>.",
    )
    parser.add_argument(
        "--output-root",
        default=None,
        help="Processed SDS output root. Defaults to data/processed/sds_output/<POP>.",
    )
    parser.add_argument(
        "--normalized-table",
        default=None,
        help="Existing normalized SDS table. Defaults to data/processed/sds_output/<POP>/<POP>.normalized.tsv.",
    )
    parser.add_argument(
        "--current-g-file",
        default=None,
        help="Current canonical gamma-shape file. Defaults to /data/home/.../ms/scripts/sds_input.gamma_shapes.<POP>.final.",
    )
    parser.add_argument(
        "--phlash-pickle",
        default=None,
        help="phlash posterior pickle. Defaults to /data/home/.../phlash/results/<POP>/<POP>_model_full.pkl.",
    )
    parser.add_argument(
        "--percentiles",
        default="2.5,25,50,75,97.5",
        help="Comma-separated posterior percentiles to scan.",
    )
    parser.add_argument(
        "--constant-ne",
        type=float,
        default=None,
        help="Optional constant diploid Ne scenario to add, e.g. 152000.",
    )
    parser.add_argument(
        "--sim-reps",
        type=int,
        default=1000,
        help="Number of neutral replicates per DAF when building each minimal gamma file.",
    )
    parser.add_argument(
        "--init",
        default="0.0001",
        help="Initial SDS optimizer scale passed to compute_SDS.py. Defaults to 1e-4 for this target.",
    )
    parser.add_argument("--s-file-ncol", default="20000", help="Maximum singleton columns per individual.")
    parser.add_argument(
        "--outdir",
        default=str(REPO_ROOT / "tmp" / f"gamma_single_snp_{datetime.now().strftime('%Y%m%d_%H%M%S')}"),
        help="Output directory for intermediate files and summary tables.",
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Only prepare target metadata and phlash percentile scenario inputs, then exit.",
    )
    parser.add_argument(
        "--no-reuse-existing-gamma",
        action="store_true",
        help="Force regeneration even when a gamma fragment file already exists and is valid.",
    )
    parser.add_argument(
        "--allow-partial",
        action="store_true",
        help="Summarize completed scenarios only; do not regenerate missing or empty gamma fragments.",
    )
    parser.add_argument("--article-label", default=None, help="Optional external cohort/article label for comparison.")
    parser.add_argument("--article-af", type=float, default=None, help="Optional external cohort allele frequency.")
    parser.add_argument("--article-sds", type=float, default=None, help="Optional external cohort SDS statistic.")
    parser.add_argument("--phlash-python", default=str(PHLASH_PYTHON), help="Python executable with phlash installed.")
    parser.add_argument("--sds-python", default=str(SDS_PYTHON), help="Python executable for compute_SDS.py.")
    parser.add_argument("--ms-make-dir", default=str(MS_MAKE_DIR), help="Directory containing the MS Makefile.")
    parser.add_argument("--ms-binary", default=str(MS_BINARY), help="Path to the ms binary.")
    parser.add_argument("--backward-script", default=str(BACKWARD_SCRIPT), help="Path to backward.py.")
    return parser.parse_args()


def default_input_root(pop: str) -> Path:
    return REPO_ROOT / "data" / "processed" / "sds_input" / pop


def default_output_root(pop: str) -> Path:
    return REPO_ROOT / "data" / "processed" / "sds_output" / pop


def default_normalized_table(pop: str) -> Path:
    return REPO_ROOT / "data" / "processed" / "sds_output" / pop / f"{pop}.normalized.tsv"


def default_current_g_file(pop: str) -> Path:
    return MS_ROOT / "scripts" / f"sds_input.gamma_shapes.{pop}.final"


def default_phlash_pickle(pop: str) -> Path:
    return PHLASH_ROOT / "results" / pop / f"{pop}_model_full.pkl"


def parse_percentiles(text: str) -> list[float]:
    percentiles = []
    for item in text.split(","):
        stripped = item.strip()
        if not stripped:
            continue
        value = float(stripped)
        if value < 0.0 or value > 100.0:
            raise ValueError(f"Percentile out of range: {value}")
        percentiles.append(value)
    if not percentiles:
        raise ValueError("No percentiles specified")
    return percentiles


def percentile_key(percentile: float) -> str:
    return f"p_{str(percentile).replace('.', '_')}"


def percentile_label(percentile: float) -> str:
    return f"q{percentile:g}"


def percentile_file_label(percentile: float) -> str:
    return percentile_label(percentile).replace(".", "p")


def format_float(value: float, digits: int = 6) -> str:
    if value is None or not math.isfinite(value):
        return "NA"
    return f"{value:.{digits}f}"


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


def run_checked(cmd: list[str], cwd: Path | None = None) -> None:
    subprocess.run(cmd, check=True, cwd=str(cwd) if cwd is not None else None)


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def load_target_normalized_row(normalized_table: Path, snp_id: str) -> dict[str, str] | None:
    with normalized_table.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("ID") == snp_id:
                return row
    return None


def load_target_raw_sds_row(output_root: Path, chrom: str, snp_id: str) -> tuple[dict[str, str] | None, Path | None]:
    candidates = [
        output_root / f"chr{chrom}_p.sds.tsv",
        output_root / f"chr{chrom}_q.sds.tsv",
        output_root / f"chr{chrom}.sds.tsv",
    ]
    for path in candidates:
        if not path.exists():
            continue
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                if row.get("ID") == snp_id:
                    return row, path
    return None, None


def load_target_t_row(t_file: Path, snp_id: str, pos: int) -> tuple[str, dict[str, object]]:
    matches: list[tuple[str, dict[str, object]]] = []
    with t_file.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            current_id = parts[0]
            current_pos = int(float(parts[3]))
            if current_id != snp_id or current_pos != pos:
                continue
            genotypes = []
            for token in parts[4:]:
                genotypes.append(-1 if token == "NA" else int(token))
            valid = [value for value in genotypes if value >= 0]
            if not valid:
                raise RuntimeError(f"Target SNP has no valid genotypes: {snp_id}")
            n0 = sum(value == 0 for value in valid)
            n1 = sum(value == 1 for value in valid)
            n2 = sum(value == 2 for value in valid)
            sample_count = len(valid)
            daf = (n1 + 2 * n2) / (2.0 * sample_count)
            matches.append(
                (
                    raw_line if raw_line.endswith("\n") else raw_line + "\n",
                    {
                        "AA": parts[1],
                        "DA": parts[2],
                        "POS": current_pos,
                        "nG0": n0,
                        "nG1": n1,
                        "nG2": n2,
                        "sample_count": sample_count,
                        "daf_exact": daf,
                        "daf_complement_exact": 1.0 - daf,
                    },
                )
            )
    if not matches:
        raise RuntimeError(f"Target SNP not found in t_file: {snp_id}")
    if len(matches) > 1:
        raise RuntimeError(f"Target SNP matched multiple rows in t_file: {snp_id}")
    return matches[0]


def write_target_t_file(path: Path, target_line: str) -> None:
    path.write_text(target_line)


def extract_percentile_curves(
    phlash_python: Path,
    phlash_pickle: Path,
    percentiles: list[float],
    out_npz: Path,
) -> dict[float, np.ndarray]:
    helper = """
import pickle
import sys
import numpy as np
from scipy.interpolate import interp1d

MU = 1.25e-8

def pct_key(value):
    return "p_" + str(value).replace(".", "_")

pkl_path, out_path, percentiles_text = sys.argv[1:4]
percentiles = [float(item) for item in percentiles_text.split(",") if item.strip()]

with open(pkl_path, "rb") as handle:
    samples = pickle.load(handle, encoding="latin1")

target_t = np.concatenate([
    np.arange(0, 5001, 1, dtype=float),
    np.exp(np.linspace(np.log(5005.0), np.log(300000.0), 1000)),
])

all_ne = []
for model in samples:
    rescaled = model.rescale(MU)
    abs_t = np.asarray(rescaled.eta.t, dtype=float)
    abs_ne = 1.0 / np.asarray(rescaled.eta.c, dtype=float)
    step = interp1d(
        abs_t,
        abs_ne,
        kind="zero",
        bounds_error=False,
        fill_value=(abs_ne[0], abs_ne[-1]),
    )
    all_ne.append(step(target_t))

all_ne = np.asarray(all_ne, dtype=float)
payload = {"t": target_t}
for pct in percentiles:
    payload[pct_key(pct)] = np.percentile(all_ne, pct, axis=0)

np.savez(out_path, **payload)
"""
    run_checked(
        [
            str(phlash_python),
            "-c",
            helper,
            str(phlash_pickle),
            str(out_npz),
            ",".join(str(value) for value in percentiles),
        ],
        cwd=REPO_ROOT,
    )
    loaded = np.load(out_npz)
    curves = {}
    for percentile in percentiles:
        curves[percentile] = np.asarray(loaded[percentile_key(percentile)], dtype=float)
    curves[0.0] = np.asarray(loaded["t"], dtype=float)
    return curves


def generation_index(t_grid: np.ndarray, generation: int) -> int:
    idx = np.searchsorted(t_grid, generation, side="right") - 1
    return max(0, min(idx, len(t_grid) - 1))


def scenario_curve_summary(t_grid: np.ndarray, ne_curve: np.ndarray) -> dict[str, float]:
    summary = {"present_diploid_pop_size": float(ne_curve[0])}
    for generation in TARGET_GENERATIONS:
        idx = generation_index(t_grid, generation)
        summary[f"ne_t{generation}"] = float(ne_curve[idx])
    return summary


def write_scenario_npz(path: Path, pop: str, t_grid: np.ndarray, ne_curve: np.ndarray) -> None:
    np.savez(path, **{f"{pop}_t": t_grid, f"{pop}_median": ne_curve})


def constant_ne_label(value: float) -> str:
    if math.isclose(value, round(value), rel_tol=0.0, abs_tol=1e-9):
        return f"constNe{int(round(value))}"
    return f"constNe{str(value).replace('.', 'p')}"


def build_constant_ne_curve(t_grid: np.ndarray, diploid_ne: float) -> np.ndarray:
    return np.full_like(t_grid, float(diploid_ne), dtype=float)


def run_make_single_daf(
    make_dir: Path,
    ms_binary: Path,
    backward_script: Path,
    pop: str,
    scenario_npz: Path,
    present_diploid_pop_size: int,
    sim_reps: int,
    daf: float,
    gamma_prefix: Path,
    workdir: Path,
) -> Path:
    daf_text = repr(float(daf))
    output_path = Path(f"{gamma_prefix}.{daf_text}")
    cmd = [
        "make",
        "-C",
        str(make_dir),
        "sim_single_daf",
        f"DAF={daf_text}",
        f"POP_MODEL={pop}",
        "NE_STAT=median",
        f"NPZ_PATH={scenario_npz}",
        f"SIM_NUM_REPLICATIONS={sim_reps}",
        f"PRESENT_DIPLOID_POPULATION_SIZE={present_diploid_pop_size}",
        f"GAMMA_PREFIX={gamma_prefix}",
        f"WORKDIR={workdir}",
        f"MS={ms_binary}",
        f"SIMUPOP_BACKWARD={backward_script} -m {pop} --npz_path {scenario_npz} --ne_stat median",
    ]
    run_checked(cmd, cwd=make_dir)
    if not output_path.exists():
        raise RuntimeError(f"Expected gamma output missing after make run: {output_path}")
    if output_path.stat().st_size == 0:
        raise RuntimeError(f"Expected non-empty gamma output after make run: {output_path}")
    return output_path


def load_single_gamma_value(path: Path) -> tuple[float, float]:
    rows = []
    with path.open() as handle:
        for line in handle:
            values = line.split()
            if len(values) < 2:
                continue
            rows.append((float(values[0]), float(values[1])))
    if not rows:
        raise RuntimeError(
            f"Gamma output is empty: {path}. Increase --sim-reps or inspect the upstream make run."
        )
    if len(rows) != 1:
        raise RuntimeError(f"Expected exactly one gamma row in {path}, found {len(rows)}")
    return rows[0]


def try_load_single_gamma_value(path: Path) -> tuple[float, float] | None:
    if not path.exists() or path.stat().st_size == 0:
        return None
    try:
        return load_single_gamma_value(path)
    except RuntimeError:
        return None


def write_minimal_g_file(path: Path, gamma_rows: list[tuple[float, float]]) -> None:
    ordered = sorted(gamma_rows, key=lambda row: row[0])
    lines = [f"{frequency:.16g}\t{shape:.16g}" for frequency, shape in ordered]
    path.write_text("\n".join(lines) + "\n")


def load_gamma_points(path: Path) -> list[tuple[float, float]]:
    rows = []
    with path.open() as handle:
        for line in handle:
            values = line.split()
            if len(values) < 2:
                continue
            rows.append((float(values[0]), float(values[1])))
    if not rows:
        raise RuntimeError(f"No gamma rows found in {path}")
    rows.sort(key=lambda row: row[0])
    return rows


def interpolate_gamma_shape(gamma_rows: list[tuple[float, float]], frequency: float) -> float:
    if frequency <= gamma_rows[0][0]:
        return gamma_rows[0][1]
    if frequency >= gamma_rows[-1][0]:
        return gamma_rows[-1][1]
    for index in range(1, len(gamma_rows)):
        left_x, left_y = gamma_rows[index - 1]
        right_x, right_y = gamma_rows[index]
        if frequency <= right_x:
            if right_x == left_x:
                return right_y
            frac = (frequency - left_x) / (right_x - left_x)
            return left_y + (right_y - left_y) * frac
    return gamma_rows[-1][1]


def run_compute(
    sds_python: Path,
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
    run_checked(
        [
            str(sds_python),
            str(REPO_ROOT / "scripts" / "compute_SDS.py"),
            str(s_file),
            str(t_file),
            str(o_file),
            str(b_file),
            str(g_file),
            str(init),
            str(s_file_ncol),
            "--output",
            str(output_tsv),
            "--summary-csv",
            str(summary_csv),
            "--pickle-cache-dir",
            str(cache_dir),
            "--progress-every",
            "1",
        ],
        cwd=REPO_ROOT,
    )


def load_single_sds_row(path: Path) -> dict[str, object]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if len(rows) != 1:
        raise RuntimeError(f"Expected exactly one SDS output row in {path}, found {len(rows)}")
    row = rows[0]
    return {
        "ID": row["ID"],
        "AA": row["AA"],
        "DA": row["DA"],
        "POS": safe_int(row["POS"]),
        "DAF": safe_float(row["DAF"]),
        "nG0": safe_int(row["nG0"]),
        "nG1": safe_int(row["nG1"]),
        "nG2": safe_int(row["nG2"]),
        "rSDS": safe_float(row["rSDS"]),
        "SuggestedInitPoint": row["SuggestedInitPoint"],
    }


def normal_two_sided_metrics(z_value: float | None) -> tuple[float | None, float | None]:
    if z_value is None or not math.isfinite(z_value):
        return None, None
    p_value = math.erfc(abs(z_value) / math.sqrt(2.0))
    neg_log10_p = math.inf if p_value == 0.0 else -math.log10(p_value)
    return p_value, neg_log10_p


def classify_sensitivity(control_rsds: float, posterior_rsds: list[float]) -> str:
    if not posterior_rsds:
        return "unknown"
    denominator = abs(control_rsds)
    if denominator < 1e-12:
        return "undefined"
    max_abs_delta = max(abs(value - control_rsds) for value in posterior_rsds)
    ratio = max_abs_delta / denominator
    if ratio < 0.10:
        return "small"
    if ratio < 0.50:
        return "moderate"
    return "large"


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: ("" if row.get(key) is None else row.get(key)) for key in fieldnames})


def write_scenario_manifest(path: Path, scenario_inputs: list[dict[str, object]]) -> None:
    rows = []
    for scenario in scenario_inputs:
        rows.append(
            {
                "label": scenario["label"],
                "file_label": scenario["file_label"],
                "source": scenario["source"],
                "percentile": scenario["percentile"],
                "scenario_dir": str(scenario["dir"]),
                "npz": str(scenario["npz"]),
                "present_diploid_pop_size": scenario["summary"]["present_diploid_pop_size"],
            }
        )
    write_tsv(
        path,
        ["label", "file_label", "source", "percentile", "scenario_dir", "npz", "present_diploid_pop_size"],
        rows,
    )


def write_summary_md(
    path: Path,
    args: argparse.Namespace,
    target_metadata: dict[str, object],
    summary_rows: list[dict[str, object]],
) -> None:
    control_row = next(row for row in summary_rows if row["scenario"] == "current_gfile")
    posterior_rows = [row for row in summary_rows if row["scenario"] != "current_gfile"]
    completed_rows = [row for row in posterior_rows if row.get("status") == "ok" and row["rSDS"] is not None]
    incomplete_rows = [row for row in posterior_rows if row.get("status") != "ok" or row["rSDS"] is None]
    posterior_rsds = [float(row["rSDS"]) for row in completed_rows if row["rSDS"] is not None]
    min_rsds = min(posterior_rsds) if posterior_rsds else math.nan
    max_rsds = max(posterior_rsds) if posterior_rsds else math.nan
    rsds_range = max_rsds - min_rsds if posterior_rsds else math.nan
    sign_flip = any(float(row["rSDS"]) * float(control_row["rSDS"]) < 0.0 for row in completed_rows if row["rSDS"] is not None)
    classification = classify_sensitivity(float(control_row["rSDS"]), posterior_rsds)

    lines = [
        "# Single-SNP Gamma Sensitivity Summary",
        "",
        f"- Population: `{args.pop}`",
        f"- Target SNP: `{args.snp_id}`",
        f"- Position: `chr{args.chrom}:{args.pos}`",
        f"- Exact DAF used for gamma generation: `{target_metadata['daf_exact']:.16f}`",
        f"- Complementary frequency used by SDS: `{target_metadata['daf_complement_exact']:.16f}`",
        f"- Existing processed raw `rSDS`: `{format_float(safe_float(None if target_metadata.get('reference_raw_rsds') is None else str(target_metadata['reference_raw_rsds'])), 4)}`",
        f"- Existing processed raw source: `{target_metadata.get('reference_raw_source') or 'NA'}`",
        f"- Existing normalized-table `norm_SDS`: `{format_float(safe_float(None if target_metadata.get('reference_norm_sds') is None else str(target_metadata['reference_norm_sds'])), 4)}`",
        f"- Control re-run with current g_file: `{format_float(float(control_row['rSDS']), 4)}`",
        f"- Completed non-control scenarios: `{len(completed_rows)}`",
        f"- Incomplete non-control scenarios: `{len(incomplete_rows)}`",
        f"- Non-control min/max/range `rSDS`: `{format_float(min_rsds, 4)}` / `{format_float(max_rsds, 4)}` / `{format_float(rsds_range, 4)}`",
        f"- Sign flip versus control: `{'yes' if sign_flip else 'no'}`",
        f"- Heuristic sensitivity classification: `{classification}`",
    ]
    if target_metadata.get("article_label") or target_metadata.get("article_sds") is not None:
        lines.append(
            f"- External comparison ({target_metadata.get('article_label') or 'external cohort'}): "
            f"`AF={format_float(safe_float(None if target_metadata.get('article_af') is None else str(target_metadata['article_af'])), 4)}`, "
            f"`SDS={format_float(safe_float(None if target_metadata.get('article_sds') is None else str(target_metadata['article_sds'])), 4)}`"
        )
    if target_metadata.get("reference_common_mean") is None or target_metadata.get("reference_common_sd") is None:
        lines.append("- Projected normalization diagnostics are unavailable because this SNP is absent from the normalized common-variant reference table.")
    lines.extend(
        [
            "",
            "## Scenario Table",
            "",
            "| Scenario | Status | Percentile | rSDS | Delta vs current g_file | Delta vs q50 | gamma(1-DAF) | gamma(DAF) | Projected norm_SDS | Projected -log10(p) |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in summary_rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["scenario"]),
                    str(row.get("status", "NA")),
                    "" if row["percentile"] is None else format_float(float(row["percentile"]), 1),
                    format_float(safe_float(None if row["rSDS"] is None else str(row["rSDS"])), 4),
                    format_float(safe_float(None if row["delta_vs_current_gfile"] is None else str(row["delta_vs_current_gfile"])), 4),
                    format_float(safe_float(None if row["delta_vs_q50"] is None else str(row["delta_vs_q50"])), 4),
                    format_float(safe_float(None if row["gamma_shape_complement"] is None else str(row["gamma_shape_complement"])), 6),
                    format_float(safe_float(None if row["gamma_shape_daf"] is None else str(row["gamma_shape_daf"])), 6),
                    format_float(row["projected_norm_sds"], 4),
                    format_float(row["projected_neg_log10_p"], 4),
                ]
            )
            + " |"
        )
    if incomplete_rows:
        lines.extend(
            [
                "",
                "## Incomplete Scenarios",
                "",
            ]
        )
        for row in incomplete_rows:
            note = row.get("note") or "missing gamma fragment"
            lines.append(f"- `{row['scenario']}`: {note}")
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- `current_gfile` is a control re-run on an exact one-row t_file; when an existing processed raw row is available, it should be close but may differ slightly because the local one-row rerun rebuilds the likelihood inputs.",
            "- `q50` is the posterior median curve recomputed only at the two exact frequencies used by this SNP, so small deltas versus `current_gfile` can still reflect simulation noise.",
            "- Projected normalization columns stay `NA` when this SNP is missing from the normalized common-variant reference table.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()

    input_root = (Path(args.input_root) if args.input_root else default_input_root(args.pop)).resolve()
    output_root = (Path(args.output_root) if args.output_root else default_output_root(args.pop)).resolve()
    normalized_table = (Path(args.normalized_table) if args.normalized_table else default_normalized_table(args.pop)).resolve()
    current_g_file = (Path(args.current_g_file) if args.current_g_file else default_current_g_file(args.pop)).resolve()
    phlash_pickle = (Path(args.phlash_pickle) if args.phlash_pickle else default_phlash_pickle(args.pop)).resolve()
    percentiles = parse_percentiles(args.percentiles)
    has_q50_requested = any(math.isclose(percentile, 50.0, rel_tol=0.0, abs_tol=1e-9) for percentile in percentiles)

    phlash_python = Path(args.phlash_python).resolve()
    sds_python = Path(args.sds_python).resolve()
    ms_make_dir = Path(args.ms_make_dir).resolve()
    ms_binary = Path(args.ms_binary).resolve()
    backward_script = Path(args.backward_script).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    s_file = input_root / f"chr{args.chrom}_s_file.txt"
    t_file = input_root / f"chr{args.chrom}_t_file.txt"
    o_file = input_root / f"chr{args.chrom}_o_file.txt"
    b_file = input_root / f"chr{args.chrom}_b_file.txt"
    for path in [s_file, t_file, o_file, b_file, current_g_file, phlash_pickle, phlash_python, sds_python, ms_binary, backward_script]:
        if not path.exists():
            raise RuntimeError(f"Required path does not exist: {path}")
    if args.normalized_table and not normalized_table.exists():
        raise RuntimeError(f"Requested normalized table does not exist: {normalized_table}")
    if args.output_root and not output_root.exists():
        raise RuntimeError(f"Requested output root does not exist: {output_root}")

    normalized_row = load_target_normalized_row(normalized_table, args.snp_id) if normalized_table.exists() else None
    reference_raw_row, reference_raw_path = load_target_raw_sds_row(output_root, args.chrom, args.snp_id) if output_root.exists() else (None, None)
    target_line, target_info = load_target_t_row(t_file, args.snp_id, args.pos)

    target_t_file = outdir / "target_snp.t.tsv"
    write_target_t_file(target_t_file, target_line)

    reference_raw_rsds = safe_float(reference_raw_row.get("rSDS")) if reference_raw_row is not None else None
    if reference_raw_rsds is None and normalized_row is not None:
        reference_raw_rsds = safe_float(normalized_row.get("rSDS"))
    reference_raw_source = str(reference_raw_path) if reference_raw_path is not None else (str(normalized_table) if normalized_row is not None else None)
    reference_norm_sds = safe_float(normalized_row.get("norm_SDS")) if normalized_row is not None else None
    common_mean = safe_float(normalized_row.get("COMMON_MEAN")) if normalized_row is not None else None
    common_sd = safe_float(normalized_row.get("COMMON_SD")) if normalized_row is not None else None

    target_metadata = {
        "snp_id": args.snp_id,
        "chrom": args.chrom,
        "pos": args.pos,
        "AA": target_info["AA"],
        "DA": target_info["DA"],
        "nG0": target_info["nG0"],
        "nG1": target_info["nG1"],
        "nG2": target_info["nG2"],
        "sample_count": target_info["sample_count"],
        "daf_exact": target_info["daf_exact"],
        "daf_complement_exact": target_info["daf_complement_exact"],
        "reference_raw_rsds": reference_raw_rsds,
        "reference_raw_source": reference_raw_source,
        "reference_raw_daf": safe_float(reference_raw_row.get("DAF")) if reference_raw_row is not None else None,
        "reference_norm_sds": reference_norm_sds,
        "reference_norm_source": str(normalized_table) if normalized_row is not None else None,
        "reference_common_mean": common_mean,
        "reference_common_sd": common_sd,
        "article_label": args.article_label,
        "article_af": args.article_af,
        "article_sds": args.article_sds,
    }
    write_json(outdir / "target_metadata.json", target_metadata)

    posterior_npz = outdir / "phlash_percentiles.npz"
    print(f"[scan] extracting phlash posterior percentiles -> {posterior_npz}", flush=True)
    curves = extract_percentile_curves(phlash_python, phlash_pickle, percentiles, posterior_npz)
    t_grid = curves.pop(0.0)

    scenario_inputs: list[dict[str, object]] = []
    for percentile in percentiles:
        label = percentile_label(percentile)
        file_label = percentile_file_label(percentile)
        scenario_dir = outdir / file_label
        scenario_dir.mkdir(parents=True, exist_ok=True)
        ne_curve = curves[percentile]
        scenario_summary = scenario_curve_summary(t_grid, ne_curve)
        scenario_npz = scenario_dir / f"{args.pop}_{file_label}.npz"
        write_scenario_npz(scenario_npz, args.pop, t_grid, ne_curve)
        write_json(scenario_dir / "curve_summary.json", scenario_summary)
        scenario_inputs.append(
            {
                "label": label,
                "file_label": file_label,
                "source": "phlash_percentile",
                "percentile": percentile,
                "dir": scenario_dir,
                "summary": scenario_summary,
                "npz": scenario_npz,
            }
        )

    if args.constant_ne is not None:
        label = constant_ne_label(args.constant_ne)
        file_label = label
        scenario_dir = outdir / file_label
        scenario_dir.mkdir(parents=True, exist_ok=True)
        ne_curve = build_constant_ne_curve(t_grid, args.constant_ne)
        scenario_summary = scenario_curve_summary(t_grid, ne_curve)
        scenario_npz = scenario_dir / f"{args.pop}_{file_label}.npz"
        write_scenario_npz(scenario_npz, args.pop, t_grid, ne_curve)
        write_json(scenario_dir / "curve_summary.json", scenario_summary)
        scenario_inputs.append(
            {
                "label": label,
                "file_label": file_label,
                "source": "constant_ne",
                "percentile": None,
                "dir": scenario_dir,
                "summary": scenario_summary,
                "npz": scenario_npz,
            }
        )

    write_scenario_manifest(outdir / "scenario_manifest.tsv", scenario_inputs)

    if args.prepare_only:
        print(f"[scan] prepared scenario inputs only -> {outdir}", flush=True)
        print(outdir)
        return 0

    current_gamma_rows = load_gamma_points(current_g_file)
    current_gamma_complement = interpolate_gamma_shape(current_gamma_rows, float(target_info["daf_complement_exact"]))
    current_gamma_daf = interpolate_gamma_shape(current_gamma_rows, float(target_info["daf_exact"]))

    cache_dir = outdir / "cache"
    summary_rows: list[dict[str, object]] = []
    q50_rsds = None
    q50_label = None

    current_dir = outdir / "current_gfile"
    current_dir.mkdir(parents=True, exist_ok=True)
    current_output = current_dir / "current_gfile.sds.tsv"
    print(f"[scan] running control compute with current g_file -> {current_output}", flush=True)
    run_compute(
        sds_python=sds_python,
        s_file=s_file,
        t_file=target_t_file,
        o_file=o_file,
        b_file=b_file,
        g_file=current_g_file,
        init=args.init,
        s_file_ncol=args.s_file_ncol,
        output_tsv=current_output,
        summary_csv=current_dir / "compute.summary.csv",
        cache_dir=cache_dir,
    )
    current_result = load_single_sds_row(current_output)
    projected_norm = None
    if common_mean is not None and common_sd is not None and common_sd > 0.0:
        projected_norm = (float(current_result["rSDS"]) - common_mean) / common_sd
    projected_p, projected_neg = normal_two_sided_metrics(projected_norm)
    summary_rows.append(
        {
            "scenario": "current_gfile",
            "source": "current_gfile",
            "percentile": None,
            "g_file": str(current_g_file),
            "output_tsv": str(current_output),
            "gamma_shape_complement": current_gamma_complement,
            "gamma_shape_daf": current_gamma_daf,
            "rSDS": current_result["rSDS"],
            "projected_norm_sds": projected_norm,
            "projected_p_bothside": projected_p,
            "projected_neg_log10_p": projected_neg,
            "suggested_init_point": current_result["SuggestedInitPoint"],
            "delta_vs_existing_output": None if reference_raw_rsds is None else float(current_result["rSDS"]) - reference_raw_rsds,
            "reference_raw_rsds": reference_raw_rsds,
            "reference_raw_source": reference_raw_source,
            "reference_norm_sds": reference_norm_sds,
            "reference_common_mean": common_mean,
            "reference_common_sd": common_sd,
            "delta_vs_current_gfile": 0.0,
            "abs_delta_vs_current_gfile": 0.0,
            "delta_vs_q50": None,
            "abs_delta_vs_q50": None,
            "status": "ok",
            "note": "",
        }
    )

    reuse_existing_gamma = not args.no_reuse_existing_gamma

    for scenario in scenario_inputs:
        label = scenario["label"]
        file_label = scenario["file_label"]
        source = scenario["source"]
        percentile = scenario["percentile"]
        scenario_dir = scenario["dir"]
        scenario_summary = scenario["summary"]
        scenario_npz = scenario["npz"]
        print(f"[scan] scenario {label}: present Ne={int(round(float(scenario_summary['present_diploid_pop_size'])))}", flush=True)

        gamma_rows = []
        gamma_by_label: dict[str, tuple[float, float] | None] = {}
        missing_notes: list[str] = []
        for freq_label, frequency in [
            ("complement", float(target_info["daf_complement_exact"])),
            ("daf", float(target_info["daf_exact"])),
        ]:
            gamma_prefix = scenario_dir / f"{file_label}_{freq_label}.gamma"
            workdir = scenario_dir / f"{file_label}_{freq_label}.work"
            gamma_piece = Path(f"{gamma_prefix}.{repr(float(frequency))}")
            existing = try_load_single_gamma_value(gamma_piece) if reuse_existing_gamma else None
            if existing is not None:
                print(
                    f"[scan] scenario {label}: reusing gamma for {freq_label} freq={frequency:.16f}",
                    flush=True,
                )
                gamma_rows.append(existing)
                gamma_by_label[freq_label] = existing
                continue
            if args.allow_partial:
                note = f"{freq_label} gamma missing or empty: {gamma_piece}"
                print(f"[scan] scenario {label}: {note}; marking scenario incomplete", flush=True)
                missing_notes.append(note)
                gamma_by_label[freq_label] = None
                continue
            print(
                f"[scan] scenario {label}: building gamma for {freq_label} freq={frequency:.16f}",
                flush=True,
            )
            gamma_piece = run_make_single_daf(
                make_dir=ms_make_dir,
                ms_binary=ms_binary,
                backward_script=backward_script,
                pop=args.pop,
                scenario_npz=scenario_npz,
                present_diploid_pop_size=int(round(float(scenario_summary["present_diploid_pop_size"]))),
                sim_reps=args.sim_reps,
                daf=frequency,
                gamma_prefix=gamma_prefix,
                workdir=workdir,
            )
            gamma_value = load_single_gamma_value(gamma_piece)
            gamma_rows.append(gamma_value)
            gamma_by_label[freq_label] = gamma_value

        if missing_notes:
            row = {
                "scenario": label,
                "source": source,
                "percentile": percentile,
                "g_file": str(scenario_dir / f"{file_label}.minimal_g.tsv"),
                "output_tsv": str(scenario_dir / f"{file_label}.sds.tsv"),
                "gamma_shape_complement": None if gamma_by_label.get("complement") is None else gamma_by_label["complement"][1],
                "gamma_shape_daf": None if gamma_by_label.get("daf") is None else gamma_by_label["daf"][1],
                "rSDS": None,
                "projected_norm_sds": None,
                "projected_p_bothside": None,
                "projected_neg_log10_p": None,
                "suggested_init_point": None,
                "delta_vs_existing_output": None,
                "reference_raw_rsds": reference_raw_rsds,
                "reference_raw_source": reference_raw_source,
                "reference_norm_sds": reference_norm_sds,
                "reference_common_mean": common_mean,
                "reference_common_sd": common_sd,
                "status": "incomplete",
                "note": "; ".join(missing_notes),
            }
            row.update(scenario_summary)
            summary_rows.append(row)
            continue

        minimal_g_file = scenario_dir / f"{file_label}.minimal_g.tsv"
        write_minimal_g_file(minimal_g_file, gamma_rows)

        output_tsv = scenario_dir / f"{file_label}.sds.tsv"
        print(f"[scan] scenario {label}: running compute -> {output_tsv}", flush=True)
        run_compute(
            sds_python=sds_python,
            s_file=s_file,
            t_file=target_t_file,
            o_file=o_file,
            b_file=b_file,
            g_file=minimal_g_file,
            init=args.init,
            s_file_ncol=args.s_file_ncol,
            output_tsv=output_tsv,
            summary_csv=scenario_dir / "compute.summary.csv",
            cache_dir=cache_dir,
        )
        result = load_single_sds_row(output_tsv)

        projected_norm = None
        if common_mean is not None and common_sd is not None and common_sd > 0.0:
            projected_norm = (float(result["rSDS"]) - common_mean) / common_sd
        projected_p, projected_neg = normal_two_sided_metrics(projected_norm)

        row = {
            "scenario": label,
            "source": source,
            "percentile": percentile,
            "g_file": str(minimal_g_file),
            "output_tsv": str(output_tsv),
            "gamma_shape_complement": interpolate_gamma_shape(gamma_rows, float(target_info["daf_complement_exact"])),
            "gamma_shape_daf": interpolate_gamma_shape(gamma_rows, float(target_info["daf_exact"])),
            "rSDS": result["rSDS"],
            "projected_norm_sds": projected_norm,
            "projected_p_bothside": projected_p,
            "projected_neg_log10_p": projected_neg,
            "suggested_init_point": result["SuggestedInitPoint"],
            "delta_vs_existing_output": None if reference_raw_rsds is None else float(result["rSDS"]) - reference_raw_rsds,
            "reference_raw_rsds": reference_raw_rsds,
            "reference_raw_source": reference_raw_source,
            "reference_norm_sds": reference_norm_sds,
            "reference_common_mean": common_mean,
            "reference_common_sd": common_sd,
            "status": "ok",
            "note": "",
        }
        row.update(scenario_summary)
        summary_rows.append(row)

        if source == "phlash_percentile" and percentile is not None and math.isclose(percentile, 50.0, rel_tol=0.0, abs_tol=1e-9):
            q50_rsds = float(result["rSDS"])
            q50_label = label

    control_rsds = float(summary_rows[0]["rSDS"])
    if has_q50_requested and q50_rsds is None and not args.allow_partial:
        raise RuntimeError("The requested percentile grid does not include 50, so delta_vs_q50 cannot be computed")

    for row in summary_rows:
        row["delta_vs_current_gfile"] = None if row["rSDS"] is None else float(row["rSDS"]) - control_rsds
        row["abs_delta_vs_current_gfile"] = None if row["delta_vs_current_gfile"] is None else abs(float(row["delta_vs_current_gfile"]))
        row["delta_vs_q50"] = None if row["rSDS"] is None or q50_rsds is None else float(row["rSDS"]) - q50_rsds
        row["abs_delta_vs_q50"] = None if row["delta_vs_q50"] is None else abs(float(row["delta_vs_q50"]))
        row["q50_label"] = q50_label
        row["target_snp_id"] = args.snp_id
        row["target_pos"] = args.pos
        row["target_daf_exact"] = target_info["daf_exact"]
        row["target_daf_complement_exact"] = target_info["daf_complement_exact"]

    fieldnames = [
        "scenario",
        "source",
        "percentile",
        "target_snp_id",
        "target_pos",
        "target_daf_exact",
        "target_daf_complement_exact",
        "g_file",
        "output_tsv",
        "gamma_shape_complement",
        "gamma_shape_daf",
        "rSDS",
        "delta_vs_current_gfile",
        "abs_delta_vs_current_gfile",
        "delta_vs_q50",
        "abs_delta_vs_q50",
        "delta_vs_existing_output",
        "projected_norm_sds",
        "projected_p_bothside",
        "projected_neg_log10_p",
        "suggested_init_point",
        "reference_raw_rsds",
        "reference_raw_source",
        "reference_norm_sds",
        "reference_common_mean",
        "reference_common_sd",
        "present_diploid_pop_size",
        "ne_t0",
        "ne_t1",
        "ne_t10",
        "ne_t50",
        "ne_t100",
        "ne_t500",
        "ne_t1000",
        "ne_t5000",
        "ne_t10000",
        "q50_label",
        "status",
        "note",
    ]
    write_tsv(outdir / "summary.tsv", fieldnames, summary_rows)
    write_summary_md(outdir / "summary.md", args, target_metadata, summary_rows)

    print(outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
