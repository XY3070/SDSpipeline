#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a concise final report from single-SNP gamma sensitivity scan outputs."
    )
    parser.add_argument("--outdir", required=True, help="Output directory produced by scan_single_snp_gamma_sensitivity.py")
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


def format_float(value: float | None, digits: int = 4) -> str:
    if value is None or not math.isfinite(value):
        return "NA"
    return f"{value:.{digits}f}"


def gamma_ratio(row: dict[str, object]) -> float | None:
    left = row.get("gamma_shape_complement")
    right = row.get("gamma_shape_daf")
    if left is None or right is None:
        return None
    left_f = float(left)
    right_f = float(right)
    if not math.isfinite(left_f) or not math.isfinite(right_f) or abs(right_f) < 1e-300:
        return None
    return left_f / right_f


def gamma_geomean(row: dict[str, object]) -> float | None:
    left = row.get("gamma_shape_complement")
    right = row.get("gamma_shape_daf")
    if left is None or right is None:
        return None
    left_f = float(left)
    right_f = float(right)
    if left_f <= 0.0 or right_f <= 0.0:
        return None
    return math.sqrt(left_f * right_f)


def classify_sensitivity(control_rsds: float, posterior_rsds: list[float]) -> tuple[str, float]:
    if not posterior_rsds:
        return "unknown", math.nan
    denom = abs(control_rsds)
    if denom < 1e-12:
        return "undefined", math.inf
    max_abs_delta = max(abs(value - control_rsds) for value in posterior_rsds)
    ratio = max_abs_delta / denom
    if ratio < 0.10:
        return "small", ratio
    if ratio < 0.50:
        return "moderate", ratio
    return "large", ratio


def monotonic_direction(xs: list[float], ys: list[float]) -> str:
    if len(xs) < 2:
        return "unknown"
    diffs = []
    for left, right in zip(ys[:-1], ys[1:]):
        delta = right - left
        if abs(delta) < 1e-12:
            continue
        diffs.append(delta)
    if not diffs:
        return "flat"
    if all(delta > 0 for delta in diffs):
        return "increasing"
    if all(delta < 0 for delta in diffs):
        return "decreasing"
    return "non-monotonic"


def load_summary_rows(path: Path) -> list[dict[str, object]]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = []
        for row in reader:
            rows.append(
                {
                    "scenario": row["scenario"],
                    "source": row.get("source") or "",
                    "percentile": safe_float(row.get("percentile")),
                    "rSDS": safe_float(row.get("rSDS")),
                    "delta_vs_current_gfile": safe_float(row.get("delta_vs_current_gfile")),
                    "delta_vs_q50": safe_float(row.get("delta_vs_q50")),
                    "gamma_shape_complement": safe_float(row.get("gamma_shape_complement")),
                    "gamma_shape_daf": safe_float(row.get("gamma_shape_daf")),
                    "projected_norm_sds": safe_float(row.get("projected_norm_sds")),
                    "projected_neg_log10_p": safe_float(row.get("projected_neg_log10_p")),
                    "status": row.get("status") or "ok",
                    "note": row.get("note") or "",
                }
            )
    if not rows:
        raise RuntimeError(f"No rows found in summary TSV: {path}")
    return rows


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir).resolve()
    summary_tsv = outdir / "summary.tsv"
    target_json = outdir / "target_metadata.json"
    report_md = outdir / "final_report.md"
    report_json = outdir / "final_report.json"

    if not summary_tsv.exists():
        raise RuntimeError(f"Missing summary.tsv: {summary_tsv}")
    if not target_json.exists():
        raise RuntimeError(f"Missing target metadata: {target_json}")

    rows = load_summary_rows(summary_tsv)
    with target_json.open() as handle:
        target = json.load(handle)

    control = next(row for row in rows if row["scenario"] == "current_gfile")
    scenario_rows = [row for row in rows if row["scenario"] != "current_gfile"]
    scenario_rows.sort(key=lambda row: (row["percentile"] is None, row["percentile"], row["scenario"]))
    completed_rows = [row for row in scenario_rows if row["status"] == "ok" and row["rSDS"] is not None]
    incomplete_rows = [row for row in scenario_rows if row["status"] != "ok" or row["rSDS"] is None]
    if not completed_rows:
        raise RuntimeError(f"No completed non-control scenarios in summary TSV: {summary_tsv}")

    percentile_rows = [row for row in completed_rows if row["percentile"] is not None]
    constant_rows = [row for row in completed_rows if row["source"] == "constant_ne"]
    q50 = next((row for row in percentile_rows if abs(float(row["percentile"]) - 50.0) < 1e-9), None)
    q25 = next((row for row in completed_rows if row["scenario"] == "q25"), None)
    q97p5 = next((row for row in completed_rows if row["scenario"] == "q97.5"), None)
    const_row = constant_rows[0] if constant_rows else None

    completed_rsds = [float(row["rSDS"]) for row in completed_rows if row["rSDS"] is not None]
    control_rsds = float(control["rSDS"])
    min_row = min(completed_rows, key=lambda row: float(row["rSDS"]))
    max_row = max(completed_rows, key=lambda row: float(row["rSDS"]))
    max_abs_delta_row = max(completed_rows, key=lambda row: abs(float(row["delta_vs_current_gfile"])))
    sign_flip = any(float(row["rSDS"]) * control_rsds < 0.0 for row in completed_rows)
    monotonic = monotonic_direction(
        [float(row["percentile"]) for row in percentile_rows],
        [float(row["rSDS"]) for row in percentile_rows],
    )
    sensitivity_label, relative_shift = classify_sensitivity(control_rsds, completed_rsds)
    q25_ratio = None if q25 is None else gamma_ratio(q25)
    q50_ratio = None if q50 is None else gamma_ratio(q50)
    q97p5_ratio = None if q97p5 is None else gamma_ratio(q97p5)
    q25_geomean = None if q25 is None else gamma_geomean(q25)
    q50_geomean = None if q50 is None else gamma_geomean(q50)
    q97p5_geomean = None if q97p5 is None else gamma_geomean(q97p5)
    q97p5_vs_q50_scale = None
    if q97p5_geomean is not None and q50_geomean is not None and q50_geomean > 0.0:
        q97p5_vs_q50_scale = q97p5_geomean / q50_geomean
    q25_sign_flip = False if q25 is None or q50 is None else float(q25["rSDS"]) * float(q50["rSDS"]) < 0.0
    q97p5_pulls_back = False if q97p5 is None or q50 is None else abs(float(q97p5["rSDS"])) < abs(float(q50["rSDS"]))

    article_af = target.get("article_af")
    article_sds = target.get("article_sds")
    article_label = target.get("article_label")
    article_same_sign = None
    article_magnitude_ratio = None
    if article_sds is not None and math.isfinite(float(article_sds)) and not math.isclose(float(article_sds), 0.0, rel_tol=0.0, abs_tol=1e-12):
        article_same_sign = control_rsds * float(article_sds) > 0.0
        article_magnitude_ratio = abs(control_rsds) / abs(float(article_sds))

    metrics = {
        "target_snp": target["snp_id"],
        "control_rsds": control_rsds,
        "reference_raw_rsds": target.get("reference_raw_rsds"),
        "reference_raw_source": target.get("reference_raw_source"),
        "reference_norm_sds": target.get("reference_norm_sds"),
        "scenario_min_rsds": float(min_row["rSDS"]),
        "scenario_min_scenario": min_row["scenario"],
        "scenario_max_rsds": float(max_row["rSDS"]),
        "scenario_max_scenario": max_row["scenario"],
        "scenario_range_rsds": float(max_row["rSDS"]) - float(min_row["rSDS"]),
        "max_abs_delta_vs_control": abs(float(max_abs_delta_row["delta_vs_current_gfile"])),
        "max_abs_delta_scenario": max_abs_delta_row["scenario"],
        "sign_flip_vs_control": sign_flip,
        "percentile_monotonicity": monotonic,
        "sensitivity_label": sensitivity_label,
        "relative_max_shift_vs_control": relative_shift,
        "q50_rsds": None if q50 is None else float(q50["rSDS"]),
        "completed_noncontrol_scenarios": [str(row["scenario"]) for row in completed_rows],
        "completed_percentile_scenarios": [str(row["scenario"]) for row in percentile_rows],
        "incomplete_noncontrol_scenarios": [str(row["scenario"]) for row in incomplete_rows],
        "q25_gamma_ratio": q25_ratio,
        "q50_gamma_ratio": q50_ratio,
        "q97p5_gamma_ratio": q97p5_ratio,
        "q25_gamma_geomean": q25_geomean,
        "q50_gamma_geomean": q50_geomean,
        "q97p5_gamma_geomean": q97p5_geomean,
        "q97p5_vs_q50_gamma_geomean_scale_factor": q97p5_vs_q50_scale,
        "article_label": article_label,
        "article_af": article_af,
        "article_sds": article_sds,
        "article_same_sign_as_control": article_same_sign,
        "article_control_magnitude_ratio": article_magnitude_ratio,
    }
    report_json.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")

    lines = [
        "# Final Gamma Sensitivity Report",
        "",
        f"- Target SNP: `{target['snp_id']}`",
        f"- Counts: `nG0={target['nG0']}`, `nG1={target['nG1']}`, `nG2={target['nG2']}`",
        f"- Exact DAF / 1-DAF: `{target['daf_exact']:.16f}` / `{target['daf_complement_exact']:.16f}`",
        f"- Existing pipeline baseline: raw `rSDS={format_float(target.get('reference_raw_rsds'))}`, `norm_SDS={format_float(target.get('reference_norm_sds'))}`",
        f"- Existing raw baseline source: `{target.get('reference_raw_source') or 'NA'}`",
        f"- Control re-run with current g_file: raw `rSDS={format_float(control_rsds)}`",
        f"- Completed non-control scenarios: `{len(completed_rows)}`",
        f"- Incomplete non-control scenarios: `{len(incomplete_rows)}`",
        "",
        "## Main Conclusion",
        "",
        f"- Across all completed non-control scenarios, raw `rSDS` spans `{format_float(metrics['scenario_min_rsds'])}` to `{format_float(metrics['scenario_max_rsds'])}`, a range of `{format_float(metrics['scenario_range_rsds'])}`.",
        f"- The largest absolute shift versus the current g_file control is `{format_float(metrics['max_abs_delta_vs_control'])}` at `{metrics['max_abs_delta_scenario']}`.",
        f"- Relative to the control, the maximum shift is classified as `{sensitivity_label}` (`{format_float(relative_shift, 3)}` of `|control rSDS|`).",
        f"- Sign flip versus the control: `{'yes' if sign_flip else 'no'}`.",
        f"- Percentile-to-`rSDS` trend: `{monotonic}`.",
    ]
    if q50 is not None:
        lines.append(f"- Posterior median (`q50`) raw `rSDS`: `{format_float(q50['rSDS'])}`.")
    if const_row is not None:
        lines.append(f"- Constant-Ne scenario `{const_row['scenario']}` gives raw `rSDS={format_float(const_row['rSDS'])}`.")
    if incomplete_rows:
        lines.append(f"- Incomplete scenarios excluded from range calculations: `{', '.join(str(row['scenario']) for row in incomplete_rows)}`.")

    if article_label or article_af is not None or article_sds is not None:
        lines.extend(
            [
                "",
                "## Cross-Cohort Context",
                "",
                f"- External record ({article_label or 'external cohort'}): `AF={format_float(article_af)}`, `SDS={format_float(article_sds)}`",
                f"- NCN control re-run: `AF={format_float(target.get('daf_exact'), 4)}`, `rSDS={format_float(control_rsds)}`",
                f"- Direction agreement between NCN control and external SDS: `{'yes' if article_same_sign else 'no' if article_same_sign is not None else 'NA'}`",
                f"- Magnitude ratio `|NCN control rSDS| / |external SDS|`: `{format_float(article_magnitude_ratio, 3)}`",
                "- Cross-cohort equality is not expected here because cohort composition, phasing, singleton structure, and demographic calibration differ; the robustness question is whether the NCN estimate keeps a stable direction and moderate range under gamma changes.",
            ]
        )

    lines.extend(
        [
            "",
            "## Why The Sign Changes",
            "",
            "The key quantities used by `compute_SDS.py` are:",
            "",
            "```text",
            "A1 = gamma(1 - DAF)",
            "A2 = gamma(DAF)",
            "rSDS = log(E1_hat) - log(E2_hat)",
            "B1 = A1 / E1",
            "B2 = A2 / E2",
            "",
            "l0 ~ 2*A1*(log(B1) - mean(log(d0 + B1))) + log(A1) - 2*mean(log(d0 + B1))",
            "l2 ~ 2*A2*(log(B2) - mean(log(d2 + B2))) + log(A2) - 2*mean(log(d2 + B2))",
            "l1 ~ mixed term containing both A1 and A2",
            "```",
            "",
            "So gamma shape is not a final rescaling step. `A1` and `A2` directly reshape the likelihood surface for `(log(E1), log(E2))`.",
            "The sign of `rSDS` changes whenever the maximum-likelihood point moves across the diagonal `log(E1) = log(E2)`.",
        ]
    )

    if q25 is not None and q50 is not None:
        lines.extend(
            [
                "",
                "For the `q25` versus `q50` comparison:",
                f"- `q50`: `A1={format_float(q50['gamma_shape_complement'], 6)}`, `A2={format_float(q50['gamma_shape_daf'], 6)}`, `A1/A2={format_float(q50_ratio, 3)}`, `rSDS={format_float(q50['rSDS'])}`",
                f"- `q25`: `A1={format_float(q25['gamma_shape_complement'], 6)}`, `A2={format_float(q25['gamma_shape_daf'], 6)}`, `A1/A2={format_float(q25_ratio, 3)}`, `rSDS={format_float(q25['rSDS'])}`",
            ]
        )
        if q25_sign_flip:
            lines.extend(
                [
                    "The ratio `A1/A2` expands enough from the `q50` regime to the `q25` regime that the complement-frequency and DAF-frequency parts of the likelihood pull the optimizer across the `log(E1)=log(E2)` diagonal.",
                    "That is why `rSDS` flips sign.",
                ]
            )
        else:
            lines.extend(
                [
                    "The ratio `A1/A2` changes between `q50` and `q25`, so the likelihood surface is reshaped, but not enough to move the optimum across the `log(E1)=log(E2)` diagonal.",
                    "That is why this SNP can move in magnitude under gamma changes while still keeping the same sign.",
                ]
            )

    if q97p5 is not None and q50 is not None:
        lines.extend(
            [
                "",
                "## Why `q97.5` Pulls The Value Back" if q97p5_pulls_back else "## Why `q97.5` Shifts The Value",
                "",
                f"- `q97.5`: `A1={format_float(q97p5['gamma_shape_complement'], 6)}`, `A2={format_float(q97p5['gamma_shape_daf'], 6)}`, `A1/A2={format_float(q97p5_ratio, 3)}`, `sqrt(A1*A2)={format_float(q97p5_geomean, 6)}`, `rSDS={format_float(q97p5['rSDS'])}`",
                f"- `q50`: `A1/A2={format_float(q50_ratio, 3)}`, `sqrt(A1*A2)={format_float(q50_geomean, 6)}`",
                f"- Shared gamma scale change: `sqrt(A1*A2)` is larger by `~{format_float(q97p5_vs_q50_scale, 1)}x` at `q97.5` than at `q50`.",
            ]
        )
        if q97p5_pulls_back:
            lines.extend(
                [
                    "At `q97.5`, the directional asymmetry is mostly gone because `A1/A2` moves back close to `1`, so the likelihood is no longer strongly tilted toward either side of the `log(E1)=log(E2)` diagonal.",
                    "At the same time, both gamma shapes increase in absolute scale, so the `A1`- and `A2`-weighted parts of the likelihood become steeper on both axes.",
                    "That combination pulls the optimizer back toward `log(E1) ≈ log(E2)`, so `rSDS` moves closer to zero.",
                ]
            )
        else:
            lines.extend(
                [
                    "At `q97.5`, the asymmetry term `A1/A2` changes again, while the shared scale `sqrt(A1*A2)` also changes.",
                    "That reweights the two sides of the likelihood and moves the optimum to a different point in `(log(E1), log(E2))` space.",
                    "For this SNP, the shift changes the magnitude of `rSDS` without forcing it back toward zero.",
                ]
            )

    lines.extend(
        [
            "",
            "## Scenario Table",
            "",
            "| Scenario | Status | Percentile | rSDS | Delta vs control | Delta vs q50 | gamma(1-DAF) | gamma(DAF) | Projected norm_SDS | Projected -log10(p) |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in [control] + scenario_rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["scenario"]),
                    str(row["status"]),
                    "" if row["percentile"] is None else format_float(row["percentile"], 1),
                    format_float(row["rSDS"]),
                    format_float(row["delta_vs_current_gfile"]),
                    format_float(row["delta_vs_q50"]),
                    format_float(row["gamma_shape_complement"], 6),
                    format_float(row["gamma_shape_daf"], 6),
                    format_float(row["projected_norm_sds"]),
                    format_float(row["projected_neg_log10_p"]),
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
            note = row["note"] if row["note"] else "missing gamma fragment"
            lines.append(f"- `{row['scenario']}`: {note}")

    report_md.write_text("\n".join(lines) + "\n")
    print(report_md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
