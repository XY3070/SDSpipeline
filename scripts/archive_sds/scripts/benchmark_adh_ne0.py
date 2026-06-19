#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np

from scan_single_snp_gamma_sensitivity import (
    BACKWARD_SCRIPT,
    MS_BINARY,
    MS_MAKE_DIR,
    PHLASH_PYTHON,
    REPO_ROOT,
    SDS_PYTHON,
    default_input_root,
    default_normalized_table,
    default_phlash_pickle,
    extract_percentile_curves,
    format_float,
    load_single_gamma_value,
    normal_two_sided_metrics,
    parse_percentiles,
    percentile_file_label,
    percentile_label,
    run_compute,
    run_checked,
    run_make_single_daf,
    safe_float,
    safe_int,
    scenario_curve_summary,
    try_load_single_gamma_value,
    write_json,
    write_minimal_g_file,
    write_scenario_npz,
    write_tsv,
)


RUN_GAMMA_CHUNK_SCRIPT = REPO_ROOT / "scripts" / "run_single_snp_gamma_chunk.sh"
AGGREGATE_GAMMA_CHUNKS_SCRIPT = REPO_ROOT / "scripts" / "aggregate_single_snp_gamma_chunks.sh"


@dataclass
class RegionDefinition:
    region_key: str
    chrom: str
    start: int
    end: int
    top_pos: int
    top_snv: str
    top_snp_id: str
    article_af: float | None
    article_sds: float | None
    genes: str | None
    source_index: int
    source_path: str


@dataclass
class SentinelRecord:
    role: str
    snp_id: str
    chrom: str
    pos: int
    aa: str
    da: str
    daf: float
    maf: float | None
    ng0: int | None
    ng1: int | None
    ng2: int | None
    baseline_rsds: float | None
    baseline_norm_sds: float | None
    common_mean: float | None
    common_sd: float | None
    is_common_variant: bool
    distance_to_top_snv_bp: int
    t_line: str | None = None

    @property
    def daf_complement(self) -> float:
        return 1.0 - self.daf


@dataclass
class ScenarioRecord:
    label: str
    file_label: str
    scenario_type: str
    source: str
    percentile: float | None
    target_ne0: float | None
    scale_factor: float | None
    interval_eligible: bool
    smoke_only: bool
    scenario_dir: Path
    scenario_npz: Path
    summary: dict[str, float]
    status: str = "prepared"
    note: str = ""


@dataclass
class PieceFragment:
    scenario_label: str
    purpose: str
    frequency: float
    sim_reps: int
    piece_path: Path
    sidecar_path: Path
    gamma_value: tuple[float, float] | None
    status: str
    source_label: str | None = None


@dataclass
class PreparedBenchmark:
    outdir: Path
    region: RegionDefinition
    sentinels: list[SentinelRecord]
    scenarios: list[ScenarioRecord]
    q50_ne0: float
    panel_t_file: Path
    s_file: Path
    t_file: Path
    o_file: Path
    b_file: Path
    scenario_aliases: dict[str, str]
    aliases_by_canonical: dict[str, list[str]]
    benchmark_context: dict[str, object]
    unique_frequencies: list[float]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Benchmark one positive-control region across low-side Ne(0) scenarios with resumable gamma pieces."
    )
    parser.add_argument("--pop", default="NCN", help="Population label. Defaults to NCN.")
    parser.add_argument(
        "--chrom",
        default=None,
        help="Optional chromosome override/validation. Defaults to the selected region chromosome.",
    )
    parser.add_argument(
        "--input-root",
        default=None,
        help="Input root containing s/o/b/t files. Defaults to data/processed/sds_input/<POP>.",
    )
    parser.add_argument(
        "--normalized-table",
        default=None,
        help="Normalized SDS table. Defaults to data/processed/sds_output/<POP>/<POP>.normalized.tsv.",
    )
    parser.add_argument(
        "--phlash-pickle",
        default=None,
        help="phlash posterior pickle. Defaults to /data/home/.../phlash/results/<POP>/<POP>_model_full.pkl.",
    )
    parser.add_argument(
        "--region-source-file",
        default=None,
        help="Region catalog used to pick the benchmark locus. Supports legacy one-row files and significant-region tables.",
    )
    parser.add_argument(
        "--region-key",
        default=None,
        help="Explicit region to benchmark, e.g. chr12:110887118-113506244. If omitted, auto-pick from the source table.",
    )
    parser.add_argument(
        "--auto-pick-rank",
        type=int,
        default=1,
        help="1-based rank after sorting remaining candidate regions by descending |SDS|.",
    )
    parser.add_argument(
        "--exclude-regions",
        default="chr4:98693689-99710113",
        help="Comma-separated region keys to exclude when auto-picking. Defaults to the ADH region.",
    )
    parser.add_argument(
        "--posterior-context-percentiles",
        default="2.5,25,50",
        help="Comma-separated posterior percentiles to include as context scenarios.",
    )
    parser.add_argument(
        "--target-ne0-grid",
        default="100000,152000,250000,500000,750000,1000000,q50",
        help="Comma-separated scaled-q50 target Ne(0) grid. Supports the token q50.",
    )
    parser.add_argument(
        "--diagnostic-ne0-grid",
        default="2500000,5000000",
        help="Comma-separated high-side diagnostic scaled-q50 Ne(0) grid.",
    )
    parser.add_argument(
        "--skip-high-diagnostics",
        action="store_true",
        help="Skip the high-side diagnostic scaled-q50 scenarios.",
    )
    parser.add_argument(
        "--skip-q97p5-smoke",
        action="store_true",
        help="Skip the q97.5 posterior smoke-only numerical stability check.",
    )
    parser.add_argument(
        "--scenario-include",
        default=None,
        help="Optional comma-separated scenario labels to keep, e.g. scaled_ne0_100000.",
    )
    parser.add_argument(
        "--gamma-generation-mode",
        choices=["serial", "chunked"],
        default="serial",
        help="Gamma generation mode for missing pieces. Defaults to serial.",
    )
    parser.add_argument(
        "--gamma-chunk-size",
        type=int,
        default=100,
        help="Replicates per chunk when --gamma-generation-mode=chunked.",
    )
    parser.add_argument(
        "--sim-reps",
        type=int,
        default=1000,
        help="Number of neutral replicates per DAF for full gamma generation.",
    )
    parser.add_argument(
        "--smoke-sim-reps",
        type=int,
        default=1,
        help="Number of neutral replicates per DAF for smoke-only gamma checks.",
    )
    parser.add_argument(
        "--spacing-bp",
        type=int,
        default=20000,
        help="Minimum spacing for non-anchor peak sentinels.",
    )
    parser.add_argument(
        "--large-ne-smoke-threshold",
        type=float,
        default=1.0e7,
        help="Run a smoke gamma pass before full compute when present Ne(0) reaches this threshold.",
    )
    parser.add_argument(
        "--init",
        default="0.0001",
        help="Initial optimizer scale passed to compute_SDS.py.",
    )
    parser.add_argument("--s-file-ncol", default="20000", help="Maximum singleton columns per individual.")
    parser.add_argument(
        "--outdir",
        default=None,
        help="Output directory for manifests, scenario files, and summaries. Defaults to a region-specific tmp path.",
    )
    parser.add_argument(
        "--resume-outdir",
        default=None,
        help="Reuse or finalize an existing prepared outdir instead of creating a fresh timestamped directory.",
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Only prepare the region sentinel panel and scenario inputs, then exit.",
    )
    parser.add_argument(
        "--piece-worker",
        action="store_true",
        help="Generate exactly one gamma piece for one scenario/purpose/frequency triple, then exit.",
    )
    parser.add_argument(
        "--chunk-worker",
        action="store_true",
        help="Generate exactly one gamma chunk for one scenario/purpose/frequency/range tuple, then exit.",
    )
    parser.add_argument(
        "--aggregate-worker",
        action="store_true",
        help="Aggregate one chunked gamma piece and write canonical cache metadata, then exit.",
    )
    parser.add_argument(
        "--finalize-only",
        action="store_true",
        help="Read a prepared outdir plus completed gamma pieces and rebuild the final summaries.",
    )
    parser.add_argument(
        "--no-reuse-existing-gamma",
        action="store_true",
        help="Force regeneration even when cached gamma fragments are present and validated.",
    )
    parser.add_argument("--scenario-label", default=None, help="Scenario label used by --piece-worker.")
    parser.add_argument("--gamma-purpose", choices=["full", "smoke"], default=None, help="Gamma purpose used by --piece-worker.")
    parser.add_argument("--frequency", type=float, default=None, help="Gamma frequency used by --piece-worker.")
    parser.add_argument("--chunk-index", type=int, default=None, help="1-based chunk index used by --chunk-worker.")
    parser.add_argument("--chunk-start-rep", type=int, default=None, help="Inclusive replication start used by --chunk-worker.")
    parser.add_argument("--chunk-end-rep", type=int, default=None, help="Inclusive replication end used by --chunk-worker.")
    parser.add_argument("--phlash-python", default=str(PHLASH_PYTHON), help="Python executable with phlash installed.")
    parser.add_argument("--sds-python", default=str(SDS_PYTHON), help="Python executable for compute_SDS.py.")
    parser.add_argument("--ms-make-dir", default=str(MS_MAKE_DIR), help="Directory containing the MS Makefile.")
    parser.add_argument("--ms-binary", default=str(MS_BINARY), help="Path to the ms binary.")
    parser.add_argument("--backward-script", default=str(BACKWARD_SCRIPT), help="Path to backward.py.")
    return parser.parse_args(argv)


def normalize_chromosome(text: str) -> str:
    value = text.strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    return value


def parse_csv_labels(text: str | None) -> list[str]:
    if text is None:
        return []
    labels = [item.strip() for item in text.split(",") if item.strip()]
    seen: set[str] = set()
    ordered: list[str] = []
    for label in labels:
        if label in seen:
            continue
        ordered.append(label)
        seen.add(label)
    return ordered


def region_key(chrom: str, start: int, end: int) -> str:
    return f"chr{normalize_chromosome(chrom)}:{int(start)}-{int(end)}"


def parse_region_key(text: str) -> tuple[str, int, int]:
    cleaned = text.strip()
    chrom_part, interval_text = cleaned.split(":")
    start_text, end_text = interval_text.split("-")
    return normalize_chromosome(chrom_part), int(start_text), int(end_text)


def canonicalize_top_snp_id(text: str) -> str:
    cleaned = text.strip()
    if cleaned.lower().startswith("chr") and ":" in cleaned:
        return cleaned
    parts = cleaned.split("_")
    if len(parts) != 4:
        raise RuntimeError(f"Unexpected top SNV format: {text}")
    chrom, pos, aa, da = parts
    return f"chr{normalize_chromosome(chrom)}:{int(pos)}:{aa}:{da}"


def maybe_get(row: dict[str, str], *keys: str) -> str | None:
    for key in keys:
        value = row.get(key)
        if value is not None and str(value).strip():
            return str(value).strip()
    return None


def load_region_catalog(path: Path) -> list[RegionDefinition]:
    rows: list[RegionDefinition] = []
    with path.open() as probe:
        first_line = probe.readline()
    delimiter = "\t" if "\t" in first_line else " "
    with path.open() as handle:
        if delimiter == "\t":
            reader = csv.DictReader(handle, delimiter="\t")
        else:
            reader = csv.DictReader(handle, delimiter=" ", skipinitialspace=True)
        for index, row in enumerate(reader, start=1):
            region_text = maybe_get(row, "Region", "region")
            top_snv = maybe_get(row, "Top significant SNV ID", "top_SNV")
            if region_text is None or top_snv is None:
                continue
            chrom, start, end = parse_region_key(region_text)
            rows.append(
                RegionDefinition(
                    region_key=region_key(chrom, start, end),
                    chrom=chrom,
                    start=start,
                    end=end,
                    top_pos=int(canonicalize_top_snp_id(top_snv).split(":")[1]),
                    top_snv=top_snv,
                    top_snp_id=canonicalize_top_snp_id(top_snv),
                    article_af=safe_float(maybe_get(row, "AF")),
                    article_sds=safe_float(maybe_get(row, "SDS")),
                    genes=maybe_get(row, "Genes"),
                    source_index=index,
                    source_path=str(path),
                )
            )
    if not rows:
        raise RuntimeError(f"No region rows found in {path}")
    return rows


def parse_exclude_regions(text: str) -> set[str]:
    excluded: set[str] = set()
    for item in text.split(","):
        stripped = item.strip()
        if not stripped:
            continue
        chrom, start, end = parse_region_key(stripped)
        excluded.add(region_key(chrom, start, end))
    return excluded


def select_region_definition(
    regions: list[RegionDefinition],
    region_key_text: str | None,
    auto_pick_rank: int,
    excluded_regions: set[str],
) -> RegionDefinition:
    if auto_pick_rank < 1:
        raise RuntimeError("--auto-pick-rank must be >= 1")
    if region_key_text is not None:
        chrom, start, end = parse_region_key(region_key_text)
        requested = region_key(chrom, start, end)
        for region in regions:
            if region.region_key == requested:
                return region
        raise RuntimeError(f"Requested region {requested} was not found in the source catalog")

    candidates = [region for region in regions if region.region_key not in excluded_regions and region.article_sds is not None]
    if not candidates:
        if len(regions) == 1:
            return regions[0]
        raise RuntimeError("No candidate regions remain after exclusions")
    candidates.sort(key=lambda region: (-abs(float(region.article_sds)), region.source_index))
    if auto_pick_rank > len(candidates):
        raise RuntimeError(f"--auto-pick-rank {auto_pick_rank} exceeds the {len(candidates)} available candidate regions")
    return candidates[auto_pick_rank - 1]


def parse_ne0_grid(text: str, q50_ne0: float) -> list[float]:
    values: list[float] = []
    for item in text.split(","):
        token = item.strip().lower()
        if not token:
            continue
        if token == "q50":
            values.append(float(q50_ne0))
        else:
            values.append(float(token))
    deduped: list[float] = []
    for value in values:
        if any(math.isclose(value, existing, rel_tol=0.0, abs_tol=1e-9) for existing in deduped):
            continue
        deduped.append(value)
    return deduped


def scenario_ne0_label(value: float) -> str:
    return f"{int(round(float(value)))}"


def frequency_key(value: float) -> str:
    return repr(float(value))


def frequency_file_label(value: float) -> str:
    return (
        frequency_key(value)
        .replace(".", "p")
        .replace("-", "m")
        .replace("+", "")
    )


SIGNATURE_CACHE: dict[tuple[str, str], dict[str, object]] = {}


def file_signature(path: Path) -> dict[str, object]:
    resolved = path.resolve()
    cache_key = ("file", str(resolved))
    cached = SIGNATURE_CACHE.get(cache_key)
    if cached is not None:
        return cached
    stat = resolved.stat()
    digest = hashlib.sha256()
    with resolved.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    payload = {
        "size": stat.st_size,
        "sha256": digest.hexdigest(),
    }
    SIGNATURE_CACHE[cache_key] = payload
    return payload


def scenario_npz_signature(path: Path, pop: str) -> dict[str, object]:
    resolved = path.resolve()
    cache_key = ("npz", f"{resolved}:{pop}")
    cached = SIGNATURE_CACHE.get(cache_key)
    if cached is not None:
        return cached
    with np.load(resolved) as loaded:
        t_grid = np.asarray(loaded[f"{pop}_t"], dtype=float)
        ne_curve = np.asarray(loaded[f"{pop}_median"], dtype=float)
    digest = hashlib.sha256()
    digest.update(t_grid.tobytes(order="C"))
    digest.update(ne_curve.tobytes(order="C"))
    payload = {
        "point_count": int(len(t_grid)),
        "sha256": digest.hexdigest(),
    }
    SIGNATURE_CACHE[cache_key] = payload
    return payload


def bool_from_text(text: str | None) -> bool:
    if text is None:
        return False
    return text.strip() in {"1", "true", "TRUE", "True"}


def load_region_rows(normalized_table: Path, chrom: str, start: int, end: int) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with normalized_table.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            row_chrom = normalize_chromosome(row.get("chr", ""))
            if row_chrom != chrom:
                continue
            pos = safe_int(row.get("POS"))
            if pos is None or pos < start or pos > end:
                continue
            rows.append(row)
    return rows


def build_sentinel_record(row: dict[str, str], role: str, top_pos: int) -> SentinelRecord:
    pos = safe_int(row.get("POS"))
    if pos is None:
        raise RuntimeError(f"Region row missing POS: {row}")
    snp_id = row.get("ID")
    if not snp_id:
        raise RuntimeError(f"Region row missing ID at position {pos}")
    chrom = normalize_chromosome(row.get("chr", ""))
    return SentinelRecord(
        role=role,
        snp_id=snp_id,
        chrom=chrom,
        pos=pos,
        aa=row.get("AA", ""),
        da=row.get("DA", ""),
        daf=float(row["DAF"]),
        maf=safe_float(row.get("MAF")),
        ng0=safe_int(row.get("nG0")),
        ng1=safe_int(row.get("nG1")),
        ng2=safe_int(row.get("nG2")),
        baseline_rsds=safe_float(row.get("rSDS")),
        baseline_norm_sds=safe_float(row.get("norm_SDS")),
        common_mean=safe_float(row.get("COMMON_MEAN")),
        common_sd=safe_float(row.get("COMMON_SD")),
        is_common_variant=bool_from_text(row.get("is_common_variant")),
        distance_to_top_snv_bp=abs(pos - top_pos),
    )


def sentinel_role_order(role: str) -> tuple[int, int]:
    if role == "top_exact":
        return (0, 0)
    if role == "top_proxy":
        return (0, 1)
    if role.startswith("peak_"):
        return (1, int(role.split("_", 1)[1]))
    return (2, 0)


def has_t_row(t_file: Path, snp_id: str, pos: int) -> bool:
    with t_file.open() as handle:
        for raw_line in handle:
            if not raw_line:
                continue
            parts = raw_line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            if parts[0] == snp_id and int(float(parts[3])) == pos:
                return True
    return False


def find_region_row_by_id(region_rows: list[dict[str, str]], snp_id: str) -> dict[str, str] | None:
    for row in region_rows:
        if row.get("ID") == snp_id:
            return row
    return None


def select_sentinel_panel(
    region_rows: list[dict[str, str]],
    region: RegionDefinition,
    spacing_bp: int,
    use_exact_anchor: bool,
) -> list[SentinelRecord]:
    common_rows = [row for row in region_rows if bool_from_text(row.get("is_common_variant"))]
    if len(common_rows) < 4:
        raise RuntimeError(
            f"Expected at least 4 common variants in {region.region_key}, found {len(common_rows)}"
        )

    exact_row = find_region_row_by_id(region_rows, region.top_snp_id)
    if use_exact_anchor and exact_row is not None:
        anchor_row = exact_row
        anchor_role = "top_exact"
    else:
        anchor_row = min(
            common_rows,
            key=lambda row: (
                abs(int(row["POS"]) - region.top_pos),
                -abs(float(row["norm_SDS"])),
                int(row["POS"]),
            ),
        )
        anchor_role = "top_proxy"

    selected: list[dict[str, str]] = [anchor_row]
    selected_ids = {anchor_row["ID"]}

    ranked_rows = sorted(
        common_rows,
        key=lambda row: (
            -abs(float(row["norm_SDS"])),
            abs(int(row["POS"]) - region.top_pos),
            int(row["POS"]),
        ),
    )

    for row in ranked_rows:
        if row["ID"] in selected_ids:
            continue
        pos = int(row["POS"])
        if any(abs(pos - int(existing["POS"])) < spacing_bp for existing in selected):
            continue
        selected.append(row)
        selected_ids.add(row["ID"])
        if len(selected) == 5:
            break

    if len(selected) < 5:
        for row in ranked_rows:
            if row["ID"] in selected_ids:
                continue
            selected.append(row)
            selected_ids.add(row["ID"])
            if len(selected) == 5:
                break

    if len(selected) != 5:
        raise RuntimeError(f"Failed to select 5 sentinels in {region.region_key}, got {len(selected)}")

    sentinels = [build_sentinel_record(anchor_row, anchor_role, region.top_pos)]
    peak_rows = [row for row in selected if row["ID"] != anchor_row["ID"]]
    peak_rows.sort(key=lambda row: (-abs(float(row["norm_SDS"])), int(row["POS"])))
    for index, row in enumerate(peak_rows, start=1):
        sentinels.append(build_sentinel_record(row, f"peak_{index}", region.top_pos))
    sentinels.sort(key=lambda item: sentinel_role_order(item.role))
    return sentinels


def load_t_lines_for_sentinels(t_file: Path, sentinels: Iterable[SentinelRecord]) -> dict[str, str]:
    expected = {sentinel.snp_id: sentinel.pos for sentinel in sentinels}
    found: dict[str, str] = {}
    with t_file.open() as handle:
        for raw_line in handle:
            if len(found) == len(expected):
                break
            line = raw_line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            snp_id = parts[0]
            if snp_id not in expected:
                continue
            pos = int(float(parts[3]))
            if pos != expected[snp_id]:
                continue
            found[snp_id] = raw_line if raw_line.endswith("\n") else raw_line + "\n"
    missing = sorted(set(expected) - set(found))
    if missing:
        raise RuntimeError(f"Sentinel IDs missing from t_file {t_file}: {', '.join(missing)}")
    return found


def write_panel_t_file(path: Path, sentinels: list[SentinelRecord]) -> None:
    ordered = sorted(sentinels, key=lambda item: item.pos)
    path.write_text("".join(sentinel.t_line for sentinel in ordered if sentinel.t_line is not None))


def load_sds_rows(path: Path) -> dict[str, dict[str, object]]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = {}
        for row in reader:
            rows[row["ID"]] = {
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
    return rows


def load_json_if_exists(path: Path) -> dict[str, object] | None:
    if not path.exists():
        return None
    with path.open() as handle:
        return json.load(handle)


def percentile_npz_key(percentile: float) -> str:
    return f"p_{str(percentile).replace('.', '_')}"


def q50_equivalent_ne0(value: float, q50_ne0: float) -> bool:
    abs_tol = max(1e-6, abs(float(value)) * 1e-12, abs(float(q50_ne0)) * 1e-12)
    return math.isclose(float(value), float(q50_ne0), rel_tol=0.0, abs_tol=abs_tol)


def load_existing_percentile_curves(path: Path, percentiles: list[float]) -> dict[float, np.ndarray]:
    with np.load(path) as loaded:
        if "t" not in loaded:
            raise KeyError("t")
        curves = {0.0: np.asarray(loaded["t"], dtype=float)}
        for percentile in percentiles:
            key = percentile_npz_key(percentile)
            if key not in loaded:
                raise KeyError(key)
            curves[percentile] = np.asarray(loaded[key], dtype=float)
    return curves


def load_sentinel_manifest(path: Path) -> list[SentinelRecord]:
    sentinels: list[SentinelRecord] = []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sentinels.append(
                SentinelRecord(
                    role=row["role"],
                    snp_id=row["snp_id"],
                    chrom=normalize_chromosome(row["chrom"]),
                    pos=int(row["pos"]),
                    aa=row["aa"],
                    da=row["da"],
                    daf=float(row["daf"]),
                    maf=safe_float(row.get("maf")),
                    ng0=safe_int(row.get("ng0")),
                    ng1=safe_int(row.get("ng1")),
                    ng2=safe_int(row.get("ng2")),
                    baseline_rsds=safe_float(row.get("baseline_rsds")),
                    baseline_norm_sds=safe_float(row.get("baseline_norm_sds")),
                    common_mean=safe_float(row.get("common_mean")),
                    common_sd=safe_float(row.get("common_sd")),
                    is_common_variant=bool_from_text(row.get("is_common_variant")),
                    distance_to_top_snv_bp=int(row["distance_to_top_snv_bp"]),
                )
            )
    if not sentinels:
        raise RuntimeError(f"No sentinel rows found in {path}")
    sentinels.sort(key=lambda item: sentinel_role_order(item.role))
    return sentinels


def sentinel_manifest_rows(sentinels: list[SentinelRecord]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for sentinel in sentinels:
        payload = asdict(sentinel)
        payload["t_line"] = None
        rows.append(payload)
    return rows


def scenario_manifest_rows(scenarios: list[ScenarioRecord]) -> list[dict[str, object]]:
    return [
        {
            "label": scenario.label,
            "file_label": scenario.file_label,
            "scenario_type": scenario.scenario_type,
            "source": scenario.source,
            "percentile": scenario.percentile,
            "target_ne0": scenario.target_ne0,
            "scale_factor": scenario.scale_factor,
            "interval_eligible": int(scenario.interval_eligible),
            "smoke_only": int(scenario.smoke_only),
            "scenario_dir": str(scenario.scenario_dir),
            "scenario_npz": str(scenario.scenario_npz),
            "present_diploid_pop_size": scenario.summary["present_diploid_pop_size"],
            "status": scenario.status,
            "note": scenario.note,
        }
        for scenario in scenarios
    ]


def write_scenario_manifest(path: Path, scenarios: list[ScenarioRecord]) -> None:
    write_tsv(
        path,
        [
            "label",
            "file_label",
            "scenario_type",
            "source",
            "percentile",
            "target_ne0",
            "scale_factor",
            "interval_eligible",
            "smoke_only",
            "scenario_dir",
            "scenario_npz",
            "present_diploid_pop_size",
            "status",
            "note",
        ],
        scenario_manifest_rows(scenarios),
    )


def purpose_sim_reps(purpose: str, args: argparse.Namespace) -> int:
    if purpose == "smoke":
        return int(args.smoke_sim_reps)
    if purpose == "full":
        return int(args.sim_reps)
    raise RuntimeError(f"Unsupported gamma purpose: {purpose}")


def gamma_prefix_for(scenario_dir: Path, file_label: str, purpose: str) -> Path:
    return scenario_dir / f"{purpose}_gamma" / f"{file_label}.{purpose}"


def piece_sidecar_path(piece_path: Path) -> Path:
    return Path(f"{piece_path}.meta.json")


def piece_path_for(scenario_dir: Path, file_label: str, purpose: str, frequency: float) -> Path:
    return Path(f"{gamma_prefix_for(scenario_dir, file_label, purpose)}.{frequency_key(frequency)}")


def chunk_root_for(scenario_dir: Path, purpose: str, frequency: float) -> Path:
    return scenario_dir / f"{purpose}_chunks" / frequency_file_label(frequency)


def chunk_workdir_for(chunk_root: Path, chunk_index: int) -> Path:
    return chunk_root / f"chunk_{int(chunk_index):03d}"


def chunk_ranges(sim_reps: int, chunk_size: int) -> list[tuple[int, int, int]]:
    if sim_reps <= 0:
        raise RuntimeError("sim_reps must be positive")
    if chunk_size <= 0:
        raise RuntimeError("chunk_size must be positive")
    ranges: list[tuple[int, int, int]] = []
    start = 1
    index = 1
    while start <= sim_reps:
        end = min(sim_reps, start + chunk_size - 1)
        ranges.append((index, start, end))
        start = end + 1
        index += 1
    return ranges


def chunk_result_file_count(workdir: Path, frequency: float) -> int:
    pattern = f"res_{frequency_key(frequency)}_*.tab"
    return sum(1 for path in workdir.glob(pattern) if path.is_file() and path.stat().st_size > 0)


def chunk_is_complete(workdir: Path, frequency: float, start_rep: int, end_rep: int) -> bool:
    expected = int(end_rep) - int(start_rep) + 1
    if expected <= 0:
        raise RuntimeError("chunk replicate range must be increasing")
    return chunk_result_file_count(workdir, frequency) == expected


def filter_scenarios(
    scenarios: list[ScenarioRecord],
    scenario_aliases: dict[str, str],
    requested_labels: list[str],
) -> tuple[list[ScenarioRecord], dict[str, str], dict[str, list[str]]]:
    if not requested_labels:
        aliases_by_canonical: dict[str, list[str]] = {}
        for alias_label, canonical_label in scenario_aliases.items():
            aliases_by_canonical.setdefault(canonical_label, []).append(alias_label)
        for labels in aliases_by_canonical.values():
            labels.sort()
        return scenarios, scenario_aliases, aliases_by_canonical

    existing_labels = {scenario.label for scenario in scenarios}
    canonical_requested: set[str] = set()
    unknown: list[str] = []
    for label in requested_labels:
        if label in existing_labels:
            canonical_requested.add(label)
        elif label in scenario_aliases:
            canonical_requested.add(scenario_aliases[label])
        else:
            unknown.append(label)
    if unknown:
        raise RuntimeError(f"Unknown scenario labels: {', '.join(unknown)}")

    filtered_scenarios = [scenario for scenario in scenarios if scenario.label in canonical_requested]
    filtered_aliases = {
        alias_label: canonical_label
        for alias_label, canonical_label in scenario_aliases.items()
        if canonical_label in canonical_requested
    }
    aliases_by_canonical: dict[str, list[str]] = {}
    for alias_label, canonical_label in filtered_aliases.items():
        aliases_by_canonical.setdefault(canonical_label, []).append(alias_label)
    for labels in aliases_by_canonical.values():
        labels.sort()
    return filtered_scenarios, filtered_aliases, aliases_by_canonical


def piece_metadata_payload(
    args: argparse.Namespace,
    scenario_label: str,
    purpose: str,
    frequency: float,
    sim_reps: int,
    scenario_npz: Path,
) -> dict[str, object]:
    return {
        "metadata_version": 1,
        "scenario_label": scenario_label,
        "purpose": purpose,
        "frequency": float(frequency),
        "sim_reps": int(sim_reps),
        "scenario_npz_signature": scenario_npz_signature(scenario_npz, args.pop),
        "ms_binary_signature": file_signature(Path(args.ms_binary)),
        "backward_script_signature": file_signature(Path(args.backward_script)),
    }


def find_scenario(scenarios: list[ScenarioRecord], label: str) -> ScenarioRecord | None:
    for scenario in scenarios:
        if scenario.label == label:
            return scenario
    return None


def resolve_frequency_from_panel(panel_frequencies: Iterable[float], requested_frequency: float) -> float:
    requested_value = float(requested_frequency)
    for value in panel_frequencies:
        if math.isclose(float(value), requested_value, rel_tol=0.0, abs_tol=1e-9):
            return float(value)
    available = ", ".join(sorted(frequency_key(value) for value in panel_frequencies))
    raise RuntimeError(f"Frequency {requested_frequency} is not in the panel set: {available}")


def scenario_requires_smoke_prepass(scenario: ScenarioRecord, args: argparse.Namespace) -> bool:
    if scenario.smoke_only:
        return True
    return float(scenario.summary["present_diploid_pop_size"]) >= float(args.large_ne_smoke_threshold)


def candidate_piece_specs(
    prepared: PreparedBenchmark,
    args: argparse.Namespace,
    scenario: ScenarioRecord,
    purpose: str,
    frequency: float,
) -> list[dict[str, object]]:
    specs: list[dict[str, object]] = [
        {
            "piece_path": piece_path_for(scenario.scenario_dir, scenario.file_label, purpose, frequency),
            "scenario_npz": scenario.scenario_npz,
            "source_label": scenario.label,
        }
    ]
    for alias_label in prepared.aliases_by_canonical.get(scenario.label, []):
        alias_dir = prepared.outdir / alias_label
        alias_npz = alias_dir / f"{args.pop}_{alias_label}.npz"
        if not alias_npz.exists():
            continue
        specs.append(
            {
                "piece_path": piece_path_for(alias_dir, alias_label, purpose, frequency),
                "scenario_npz": alias_npz,
                "source_label": alias_label,
            }
        )
    return specs


def discover_piece_fragment(
    prepared: PreparedBenchmark,
    scenario: ScenarioRecord,
    purpose: str,
    frequency: float,
    args: argparse.Namespace,
    *,
    allow_legacy_backfill: bool,
) -> PieceFragment:
    sim_reps = purpose_sim_reps(purpose, args)
    for spec in candidate_piece_specs(prepared, args, scenario, purpose, frequency):
        piece_path = Path(spec["piece_path"])
        sidecar_path = piece_sidecar_path(piece_path)
        expected_metadata = piece_metadata_payload(
            args=args,
            scenario_label=scenario.label,
            purpose=purpose,
            frequency=frequency,
            sim_reps=sim_reps,
            scenario_npz=Path(spec["scenario_npz"]),
        )
        if piece_path.exists() and piece_path.stat().st_size > 0 and sidecar_path.exists():
            if load_json_if_exists(sidecar_path) == expected_metadata:
                gamma_value = try_load_single_gamma_value(piece_path)
                if gamma_value is not None:
                    return PieceFragment(
                        scenario_label=scenario.label,
                        purpose=purpose,
                        frequency=frequency,
                        sim_reps=sim_reps,
                        piece_path=piece_path,
                        sidecar_path=sidecar_path,
                        gamma_value=gamma_value,
                        status="reused" if spec["source_label"] == scenario.label else "reused_alias",
                        source_label=str(spec["source_label"]),
                    )
        if piece_path.exists() and piece_path.stat().st_size > 0 and not sidecar_path.exists() and allow_legacy_backfill:
            gamma_value = try_load_single_gamma_value(piece_path)
            if gamma_value is None:
                continue
            write_json(sidecar_path, expected_metadata)
            return PieceFragment(
                scenario_label=scenario.label,
                purpose=purpose,
                frequency=frequency,
                sim_reps=sim_reps,
                piece_path=piece_path,
                sidecar_path=sidecar_path,
                gamma_value=gamma_value,
                status="legacy_backfilled" if spec["source_label"] == scenario.label else "legacy_alias_backfilled",
                source_label=str(spec["source_label"]),
            )
    canonical_piece = piece_path_for(scenario.scenario_dir, scenario.file_label, purpose, frequency)
    return PieceFragment(
        scenario_label=scenario.label,
        purpose=purpose,
        frequency=frequency,
        sim_reps=sim_reps,
        piece_path=canonical_piece,
        sidecar_path=piece_sidecar_path(canonical_piece),
        gamma_value=None,
        status="missing",
        source_label=scenario.label,
    )


def generate_piece_fragment(
    prepared: PreparedBenchmark,
    scenario: ScenarioRecord,
    purpose: str,
    frequency: float,
    args: argparse.Namespace,
    *,
    reuse_existing_gamma: bool,
) -> PieceFragment:
    if reuse_existing_gamma:
        reusable = discover_piece_fragment(
            prepared=prepared,
            scenario=scenario,
            purpose=purpose,
            frequency=frequency,
            args=args,
            allow_legacy_backfill=True,
        )
        if reusable.gamma_value is not None:
            return reusable

    sim_reps = purpose_sim_reps(purpose, args)
    purpose_dir = scenario.scenario_dir / f"{purpose}_gamma"
    purpose_dir.mkdir(parents=True, exist_ok=True)
    gamma_prefix = purpose_dir / f"{scenario.file_label}.{purpose}"
    workdir = scenario.scenario_dir / f"{purpose}_work" / frequency_file_label(frequency)
    workdir.mkdir(parents=True, exist_ok=True)
    piece_path = run_make_single_daf(
        make_dir=Path(args.ms_make_dir).resolve(),
        ms_binary=Path(args.ms_binary).resolve(),
        backward_script=Path(args.backward_script).resolve(),
        pop=args.pop,
        scenario_npz=scenario.scenario_npz,
        present_diploid_pop_size=int(round(float(scenario.summary["present_diploid_pop_size"]))),
        sim_reps=sim_reps,
        daf=frequency,
        gamma_prefix=gamma_prefix,
        workdir=workdir,
    )
    gamma_value = load_single_gamma_value(piece_path)
    sidecar_path = piece_sidecar_path(piece_path)
    write_json(
        sidecar_path,
        piece_metadata_payload(
            args=args,
            scenario_label=scenario.label,
            purpose=purpose,
            frequency=frequency,
            sim_reps=sim_reps,
            scenario_npz=scenario.scenario_npz,
        ),
    )
    return PieceFragment(
        scenario_label=scenario.label,
        purpose=purpose,
        frequency=frequency,
        sim_reps=sim_reps,
        piece_path=piece_path,
        sidecar_path=sidecar_path,
        gamma_value=gamma_value,
        status="generated",
        source_label=scenario.label,
    )


def aggregate_piece_fragment(
    prepared: PreparedBenchmark,
    scenario: ScenarioRecord,
    purpose: str,
    frequency: float,
    args: argparse.Namespace,
    *,
    reuse_existing_gamma: bool,
) -> PieceFragment:
    if reuse_existing_gamma:
        reusable = discover_piece_fragment(
            prepared=prepared,
            scenario=scenario,
            purpose=purpose,
            frequency=frequency,
            args=args,
            allow_legacy_backfill=True,
        )
        if reusable.gamma_value is not None:
            return reusable

    sim_reps = purpose_sim_reps(purpose, args)
    ranges = chunk_ranges(sim_reps, int(args.gamma_chunk_size))
    chunk_root = chunk_root_for(scenario.scenario_dir, purpose, frequency)
    missing_chunks = [
        f"chunk_{chunk_index:03d}"
        for chunk_index, start_rep, end_rep in ranges
        if not chunk_is_complete(chunk_workdir_for(chunk_root, chunk_index), frequency, start_rep, end_rep)
    ]
    if missing_chunks:
        raise RuntimeError(
            f"Cannot aggregate {scenario.label} {purpose} {frequency_key(frequency)}; missing chunks: "
            + ", ".join(missing_chunks)
        )

    gamma_prefix = gamma_prefix_for(scenario.scenario_dir, scenario.file_label, purpose)
    piece_path = piece_path_for(scenario.scenario_dir, scenario.file_label, purpose, frequency)
    piece_path.parent.mkdir(parents=True, exist_ok=True)
    run_checked(
        [
            "bash",
            str(AGGREGATE_GAMMA_CHUNKS_SCRIPT),
            "--daf",
            frequency_key(frequency),
            "--chunk-root",
            str(chunk_root),
            "--gamma-prefix",
            str(gamma_prefix),
        ],
        cwd=REPO_ROOT,
    )
    gamma_value = load_single_gamma_value(piece_path)
    sidecar_path = piece_sidecar_path(piece_path)
    write_json(
        sidecar_path,
        piece_metadata_payload(
            args=args,
            scenario_label=scenario.label,
            purpose=purpose,
            frequency=frequency,
            sim_reps=sim_reps,
            scenario_npz=scenario.scenario_npz,
        ),
    )
    return PieceFragment(
        scenario_label=scenario.label,
        purpose=purpose,
        frequency=frequency,
        sim_reps=sim_reps,
        piece_path=piece_path,
        sidecar_path=sidecar_path,
        gamma_value=gamma_value,
        status="generated_chunked",
        source_label=scenario.label,
    )


def collect_gamma_fragments(
    prepared: PreparedBenchmark,
    scenario: ScenarioRecord,
    frequencies: list[float],
    purpose: str,
    args: argparse.Namespace,
    *,
    generate_missing: bool,
    reuse_existing_gamma: bool,
    allow_legacy_backfill: bool,
    write_manifest: bool,
) -> tuple[dict[str, tuple[float, float]], list[dict[str, object]], list[float]]:
    gamma_values: dict[str, tuple[float, float]] = {}
    manifest_rows: list[dict[str, object]] = []
    missing_frequencies: list[float] = []

    for frequency in frequencies:
        if generate_missing:
            fragment = generate_piece_fragment(
                prepared=prepared,
                scenario=scenario,
                purpose=purpose,
                frequency=frequency,
                args=args,
                reuse_existing_gamma=reuse_existing_gamma,
            )
        else:
            fragment = discover_piece_fragment(
                prepared=prepared,
                scenario=scenario,
                purpose=purpose,
                frequency=frequency,
                args=args,
                allow_legacy_backfill=allow_legacy_backfill,
            )

        output_frequency: float | None = None
        shape: float | None = None
        if fragment.gamma_value is not None:
            gamma_values[frequency_key(frequency)] = fragment.gamma_value
            output_frequency = float(fragment.gamma_value[0])
            shape = float(fragment.gamma_value[1])
        else:
            missing_frequencies.append(frequency)

        manifest_rows.append(
            {
                "scenario": scenario.label,
                "purpose": purpose,
                "frequency": output_frequency if output_frequency is not None else float(frequency),
                "shape": shape,
                "sim_reps": fragment.sim_reps,
                "cache_status": fragment.status,
                "piece_path": str(fragment.piece_path),
                "source_label": fragment.source_label,
            }
        )

    if write_manifest:
        write_tsv(
            scenario.scenario_dir / f"{purpose}_gamma_manifest.tsv",
            ["scenario", "purpose", "frequency", "shape", "sim_reps", "cache_status", "piece_path", "source_label"],
            manifest_rows,
        )
    return gamma_values, manifest_rows, missing_frequencies


def build_panel_summary_rows(
    scenario: ScenarioRecord,
    sentinels: list[SentinelRecord],
    result_rows: dict[str, dict[str, object]],
    gamma_values: dict[str, tuple[float, float]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for sentinel in sentinels:
        result = result_rows.get(sentinel.snp_id)
        if result is None:
            raise RuntimeError(f"Scenario {scenario.label} missing sentinel result for {sentinel.snp_id}")
        rsds = safe_float(None if result["rSDS"] is None else str(result["rSDS"]))
        projected_norm = None
        if rsds is not None and sentinel.common_mean is not None and sentinel.common_sd is not None and sentinel.common_sd > 0.0:
            projected_norm = (rsds - sentinel.common_mean) / sentinel.common_sd
        projected_p, projected_neg = normal_two_sided_metrics(projected_norm)
        same_sign = None
        if rsds is not None and sentinel.baseline_rsds is not None:
            same_sign = (rsds == 0.0 and sentinel.baseline_rsds == 0.0) or (rsds * sentinel.baseline_rsds > 0.0)
        rows.append(
            {
                "scenario": scenario.label,
                "scenario_type": scenario.scenario_type,
                "source": scenario.source,
                "percentile": scenario.percentile,
                "target_ne0": scenario.target_ne0,
                "scale_factor": scenario.scale_factor,
                "interval_eligible": int(scenario.interval_eligible),
                "smoke_only": int(scenario.smoke_only),
                "status": scenario.status,
                "sentinel_role": sentinel.role,
                "snp_id": sentinel.snp_id,
                "chrom": sentinel.chrom,
                "pos": sentinel.pos,
                "distance_to_top_snv_bp": sentinel.distance_to_top_snv_bp,
                "daf": sentinel.daf,
                "daf_complement": sentinel.daf_complement,
                "baseline_rSDS": sentinel.baseline_rsds,
                "baseline_norm_SDS": sentinel.baseline_norm_sds,
                "scenario_rSDS": rsds,
                "delta_vs_baseline_rSDS": None if rsds is None or sentinel.baseline_rsds is None else rsds - sentinel.baseline_rsds,
                "projected_norm_SDS": projected_norm,
                "projected_p_bothside": projected_p,
                "projected_neg_log10_p": projected_neg,
                "same_sign_vs_baseline": None if same_sign is None else int(same_sign),
                "gamma_shape_daf": gamma_values[frequency_key(sentinel.daf)][1],
                "gamma_shape_complement": gamma_values[frequency_key(sentinel.daf_complement)][1],
                "suggested_init_point": result["SuggestedInitPoint"],
                "present_diploid_pop_size": scenario.summary["present_diploid_pop_size"],
            }
        )
    return rows


def build_region_summary_row(
    scenario: ScenarioRecord,
    panel_rows: list[dict[str, object]],
    sentinels: list[SentinelRecord],
) -> dict[str, object]:
    baseline_max = max(
        abs(float(sentinel.baseline_norm_sds))
        for sentinel in sentinels
        if sentinel.baseline_norm_sds is not None and math.isfinite(float(sentinel.baseline_norm_sds))
    )
    projected_values = [
        abs(float(row["projected_norm_SDS"]))
        for row in panel_rows
        if row["projected_norm_SDS"] is not None and math.isfinite(float(row["projected_norm_SDS"]))
    ]
    projected_max = max(projected_values) if projected_values else None
    anchor_row = next(row for row in panel_rows if row["sentinel_role"] in {"top_exact", "top_proxy"})
    same_sign_count = sum(int(row["same_sign_vs_baseline"]) for row in panel_rows if row["same_sign_vs_baseline"] is not None)
    retained_peak_fraction = None if projected_max is None or baseline_max <= 0.0 else projected_max / baseline_max
    region_pass = bool(
        anchor_row["same_sign_vs_baseline"] == 1
        and same_sign_count >= 3
        and retained_peak_fraction is not None
        and retained_peak_fraction >= 0.7
    )
    return {
        "scenario": scenario.label,
        "scenario_type": scenario.scenario_type,
        "source": scenario.source,
        "percentile": scenario.percentile,
        "target_ne0": scenario.target_ne0,
        "scale_factor": scenario.scale_factor,
        "interval_eligible": int(scenario.interval_eligible),
        "smoke_only": int(scenario.smoke_only),
        "status": scenario.status,
        "anchor_same_sign": anchor_row["same_sign_vs_baseline"],
        "same_sign_count": same_sign_count,
        "sentinel_count": len(sentinels),
        "baseline_panel_max_abs_norm": baseline_max,
        "scenario_panel_max_abs_projected_norm": projected_max,
        "retained_peak_fraction": retained_peak_fraction,
        "region_pass": int(region_pass),
        "note": scenario.note,
        "present_diploid_pop_size": scenario.summary["present_diploid_pop_size"],
    }


def build_incomplete_region_summary_row(
    scenario: ScenarioRecord,
    sentinels: list[SentinelRecord],
    note: str,
) -> dict[str, object]:
    return {
        "scenario": scenario.label,
        "scenario_type": scenario.scenario_type,
        "source": scenario.source,
        "percentile": scenario.percentile,
        "target_ne0": scenario.target_ne0,
        "scale_factor": scenario.scale_factor,
        "interval_eligible": int(scenario.interval_eligible),
        "smoke_only": int(scenario.smoke_only),
        "status": scenario.status,
        "anchor_same_sign": None,
        "same_sign_count": None,
        "sentinel_count": len(sentinels),
        "baseline_panel_max_abs_norm": None,
        "scenario_panel_max_abs_projected_norm": None,
        "retained_peak_fraction": None,
        "region_pass": 0,
        "note": note,
        "present_diploid_pop_size": scenario.summary["present_diploid_pop_size"],
    }


def choose_recommended_interval(region_rows: list[dict[str, object]]) -> dict[str, object] | None:
    eligible = [
        row for row in region_rows
        if int(row["interval_eligible"]) == 1 and row["status"] == "ok" and int(row["region_pass"]) == 1
    ]
    if not eligible:
        return None
    ordered_all = sorted(
        [
            row for row in region_rows
            if int(row["interval_eligible"]) == 1
            and row["target_ne0"] is not None
        ],
        key=lambda row: float(row["target_ne0"]),
    )

    runs: list[list[dict[str, object]]] = []
    current: list[dict[str, object]] = []
    for row in ordered_all:
        if row["status"] == "ok" and int(row["region_pass"]) == 1:
            current.append(row)
        else:
            if current:
                runs.append(current)
                current = []
    if current:
        runs.append(current)
    if not runs:
        return None

    def run_priority(run: list[dict[str, object]]) -> tuple[int, int, float]:
        start = float(run[0]["target_ne0"])
        end = float(run[-1]["target_ne0"])
        contains_152k = 1 if start <= 152000.0 <= end else 0
        return (len(run), contains_152k, -start)

    best = max(runs, key=run_priority)
    return {
        "start_ne0": float(best[0]["target_ne0"]),
        "end_ne0": float(best[-1]["target_ne0"]),
        "scenario_count": len(best),
        "scenario_labels": [row["scenario"] for row in best],
    }


def write_final_report(
    path: Path,
    args: argparse.Namespace,
    region: RegionDefinition,
    sentinels: list[SentinelRecord],
    scenarios: list[ScenarioRecord],
    region_rows: list[dict[str, object]],
    recommended_interval: dict[str, object] | None,
    q50_ne0: float,
) -> None:
    scenario_lookup = {row["scenario"]: row for row in region_rows}
    lines = [
        "# Region Ne(0) Benchmark",
        "",
        f"- Population: `{args.pop}`",
        f"- Region: `chr{region.chrom}:{region.start}-{region.end}`",
        f"- Top SNV anchor: `{region.top_snv}` -> `{region.top_snp_id}`",
        f"- Posterior median present Ne(0): `{format_float(q50_ne0, 3)}`",
        f"- Sentinel count: `{len(sentinels)}`",
        f"- Region source file: `{region.source_path}`",
    ]
    if region.article_af is not None or region.article_sds is not None:
        lines.append(
            f"- Source row summary: `AF={format_float(region.article_af, 4)}`, `SDS={format_float(region.article_sds, 4)}`"
        )
    if region.genes:
        lines.append(f"- Genes: `{region.genes}`")
    if recommended_interval is None:
        lines.append("- Recommended low-side Ne(0) interval: `none passed the region rule`")
    else:
        lines.append(
            f"- Recommended low-side Ne(0) interval: "
            f"`{format_float(recommended_interval['start_ne0'], 0)} - {format_float(recommended_interval['end_ne0'], 0)}`"
        )
    lines.extend(
        [
            "",
            "## Sentinels",
            "",
            "| Role | SNP | Position | Dist to Top SNV | DAF | Baseline rSDS | Baseline norm_SDS |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for sentinel in sentinels:
        lines.append(
            "| "
            + " | ".join(
                [
                    sentinel.role,
                    sentinel.snp_id,
                    str(sentinel.pos),
                    str(sentinel.distance_to_top_snv_bp),
                    format_float(sentinel.daf, 4),
                    format_float(sentinel.baseline_rsds, 4),
                    format_float(sentinel.baseline_norm_sds, 4),
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Region Summary",
            "",
            "| Scenario | Type | Target Ne(0) | Status | Anchor same sign | Same-sign count | Peak retained | Pass |",
            "| --- | --- | ---: | --- | ---: | ---: | ---: | ---: |",
        ]
    )
    for scenario in scenarios:
        row = scenario_lookup.get(scenario.label)
        if row is None:
            lines.append(
                "| "
                + " | ".join(
                    [
                        scenario.label,
                        scenario.scenario_type,
                        "" if scenario.target_ne0 is None else format_float(scenario.target_ne0, 0),
                        scenario.status,
                        "",
                        "",
                        "",
                        "",
                    ]
                )
                + " |"
            )
            continue
        lines.append(
            "| "
            + " | ".join(
                [
                    scenario.label,
                    scenario.scenario_type,
                    "" if row["target_ne0"] is None else format_float(float(row["target_ne0"]), 0),
                    str(row["status"]),
                    "" if row["anchor_same_sign"] is None else str(row["anchor_same_sign"]),
                    str(row["same_sign_count"]),
                    format_float(safe_float(None if row["retained_peak_fraction"] is None else str(row["retained_peak_fraction"])), 4),
                    str(row["region_pass"]),
                ]
            )
            + " |"
        )

    smoke_rows = [row for row in region_rows if int(row["smoke_only"]) == 1]
    if smoke_rows:
        lines.extend(
            [
                "",
                "## Smoke Diagnostics",
                "",
            ]
        )
        for row in smoke_rows:
            note = row["note"] or "completed"
            lines.append(f"- `{row['scenario']}`: `{row['status']}` ({note})")

    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- `region_pass` requires the top-anchor sentinel to keep the baseline sign, at least 3 of 5 sentinels to keep sign, and the panel max projected |norm_SDS| to retain at least 70% of baseline.",
            "- `target-ne0-grid` scenarios scale the full q50 `Ne(t)` curve, they do not replace demography with a bare Ne(0).",
            "- `q97.5` is smoke-only here because high-Ne ms runs can still expose numerical precision issues even after the recent branch-length formatting fix.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def region_slug(region: RegionDefinition) -> str:
    return f"chr{region.chrom}_{region.start}_{region.end}"


def default_outdir_for_region(region: RegionDefinition) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return REPO_ROOT / "tmp" / f"region_ne0_benchmark_{region_slug(region)}_{timestamp}"


def resolve_default_region_source_file(entrypoint_name: str) -> Path:
    if entrypoint_name == "benchmark_adh_ne0.py":
        return REPO_ROOT / "Nyuwa_ADH_top_SNV.txt"
    return REPO_ROOT / "Nyuwa_significant_regions.tsv"


def ensure_scenario_artifacts(
    scenario_dir: Path,
    file_label: str,
    pop: str,
    t_grid: np.ndarray,
    ne_curve: np.ndarray,
) -> tuple[Path, dict[str, float]]:
    scenario_dir.mkdir(parents=True, exist_ok=True)
    scenario_npz = scenario_dir / f"{pop}_{file_label}.npz"
    if not scenario_npz.exists():
        write_scenario_npz(scenario_npz, pop, t_grid, ne_curve)
    summary_path = scenario_dir / "curve_summary.json"
    summary = load_json_if_exists(summary_path)
    if summary is None:
        payload = scenario_curve_summary(t_grid, ne_curve)
        write_json(summary_path, payload)
        summary = payload
    return scenario_npz, {str(key): float(value) for key, value in summary.items()}


def build_scenarios(
    outdir: Path,
    args: argparse.Namespace,
    t_grid: np.ndarray,
    curves: dict[float, np.ndarray],
    posterior_context_percentiles: list[float],
    target_ne0_grid: list[float],
    diagnostic_ne0_grid: list[float],
    include_q97p5_smoke: bool,
    q50_ne0: float,
) -> tuple[list[ScenarioRecord], dict[str, str]]:
    scenarios: list[ScenarioRecord] = []
    scenario_aliases: dict[str, str] = {}
    q50_curve = np.asarray(curves[50.0], dtype=float)

    def append_scenario(
        *,
        label: str,
        file_label: str,
        scenario_type: str,
        source: str,
        percentile: float | None,
        target_ne0: float | None,
        scale_factor: float | None,
        interval_eligible: bool,
        smoke_only: bool,
        ne_curve: np.ndarray,
    ) -> None:
        scenario_dir = outdir / file_label
        scenario_npz, summary = ensure_scenario_artifacts(
            scenario_dir=scenario_dir,
            file_label=file_label,
            pop=args.pop,
            t_grid=t_grid,
            ne_curve=ne_curve,
        )
        scenarios.append(
            ScenarioRecord(
                label=label,
                file_label=file_label,
                scenario_type=scenario_type,
                source=source,
                percentile=percentile,
                target_ne0=target_ne0,
                scale_factor=scale_factor,
                interval_eligible=interval_eligible,
                smoke_only=smoke_only,
                scenario_dir=scenario_dir,
                scenario_npz=scenario_npz,
                summary=summary,
            )
        )

    for percentile in sorted(set(posterior_context_percentiles + [50.0])):
        file_stub = percentile_file_label(percentile)
        append_scenario(
            label=f"{percentile_label(percentile)}_curve",
            file_label=f"{file_stub}_curve",
            scenario_type="posterior_context",
            source="phlash_percentile",
            percentile=percentile,
            target_ne0=float(curves[percentile][0]),
            scale_factor=None,
            interval_eligible=False,
            smoke_only=False,
            ne_curve=np.asarray(curves[percentile], dtype=float),
        )

    for target_ne0 in target_ne0_grid:
        label = f"scaled_ne0_{scenario_ne0_label(target_ne0)}"
        if q50_equivalent_ne0(target_ne0, q50_ne0):
            if label != "q50_curve":
                scenario_aliases[label] = "q50_curve"
            continue
        scale_factor = float(target_ne0) / q50_ne0
        append_scenario(
            label=label,
            file_label=label,
            scenario_type="scaled_ne0_main",
            source="scaled_q50",
            percentile=None,
            target_ne0=float(target_ne0),
            scale_factor=scale_factor,
            interval_eligible=True,
            smoke_only=False,
            ne_curve=q50_curve * scale_factor,
        )

    for target_ne0 in diagnostic_ne0_grid:
        scale_factor = float(target_ne0) / q50_ne0
        label = f"scaled_ne0_diag_{scenario_ne0_label(target_ne0)}"
        append_scenario(
            label=label,
            file_label=label,
            scenario_type="scaled_ne0_diagnostic",
            source="scaled_q50",
            percentile=None,
            target_ne0=float(target_ne0),
            scale_factor=scale_factor,
            interval_eligible=False,
            smoke_only=False,
            ne_curve=q50_curve * scale_factor,
        )

    if include_q97p5_smoke:
        append_scenario(
            label="q97p5_smoke",
            file_label="q97p5_smoke",
            scenario_type="posterior_smoke",
            source="phlash_percentile",
            percentile=97.5,
            target_ne0=float(curves[97.5][0]),
            scale_factor=None,
            interval_eligible=False,
            smoke_only=True,
            ne_curve=np.asarray(curves[97.5], dtype=float),
        )

    return scenarios, scenario_aliases


def validate_args(args: argparse.Namespace) -> None:
    mode_count = sum(
        bool(flag)
        for flag in [
            args.prepare_only,
            args.piece_worker,
            args.chunk_worker,
            args.aggregate_worker,
            args.finalize_only,
        ]
    )
    if mode_count > 1:
        raise RuntimeError(
            "Use at most one of --prepare-only, --piece-worker, --chunk-worker, --aggregate-worker, or --finalize-only"
        )
    if args.resume_outdir and args.outdir:
        if Path(args.resume_outdir).resolve() != Path(args.outdir).resolve():
            raise RuntimeError("--resume-outdir and --outdir must resolve to the same directory when both are provided")
    if int(args.gamma_chunk_size) <= 0:
        raise RuntimeError("--gamma-chunk-size must be positive")
    if args.piece_worker:
        if args.scenario_label is None or args.gamma_purpose is None or args.frequency is None:
            raise RuntimeError("--piece-worker requires --scenario-label, --gamma-purpose, and --frequency")
    elif args.chunk_worker:
        if (
            args.scenario_label is None
            or args.gamma_purpose is None
            or args.frequency is None
            or args.chunk_index is None
            or args.chunk_start_rep is None
            or args.chunk_end_rep is None
        ):
            raise RuntimeError(
                "--chunk-worker requires --scenario-label, --gamma-purpose, --frequency, --chunk-index, --chunk-start-rep, and --chunk-end-rep"
            )
    elif args.aggregate_worker:
        if args.scenario_label is None or args.gamma_purpose is None or args.frequency is None:
            raise RuntimeError("--aggregate-worker requires --scenario-label, --gamma-purpose, and --frequency")
    elif any(
        value is not None
        for value in [
            args.scenario_label,
            args.gamma_purpose,
            args.frequency,
            args.chunk_index,
            args.chunk_start_rep,
            args.chunk_end_rep,
        ]
    ):
        raise RuntimeError(
            "--scenario-label/--gamma-purpose/--frequency and chunk range flags are only valid with worker modes"
        )
    if (args.piece_worker or args.chunk_worker or args.aggregate_worker or args.finalize_only) and not (
        args.resume_outdir or args.outdir
    ):
        raise RuntimeError(
            "--piece-worker, --chunk-worker, --aggregate-worker, and --finalize-only require --resume-outdir or --outdir"
        )


def prepare_benchmark(args: argparse.Namespace) -> PreparedBenchmark:
    validate_args(args)
    entrypoint_name = Path(sys.argv[0]).name
    selected_outdir = Path(args.resume_outdir or args.outdir).resolve() if (args.resume_outdir or args.outdir) else None
    existing_context = None
    if selected_outdir is not None:
        existing_context = load_json_if_exists(selected_outdir / "benchmark_context.json")

    if existing_context is not None:
        region = RegionDefinition(**existing_context["region"])
        selected_chrom = normalize_chromosome(args.chrom) if args.chrom else region.chrom
        if selected_chrom != region.chrom:
            raise RuntimeError(
                f"--chrom {args.chrom} does not match selected region chromosome {region.chrom} from {region.region_key}"
            )
        outdir = selected_outdir
        region_source_file = Path(
            existing_context.get("region_source_file") or resolve_default_region_source_file(entrypoint_name)
        ).resolve()
        posterior_context_percentiles = [float(value) for value in existing_context.get("posterior_context_percentiles", [])]
        target_ne0_grid = [float(value) for value in existing_context.get("target_ne0_grid", [])]
        diagnostic_ne0_grid = [float(value) for value in existing_context.get("diagnostic_ne0_grid", [])]
        include_q97p5_smoke = (outdir / "q97p5_smoke").exists() or not args.skip_q97p5_smoke
        excluded_regions = set(str(item) for item in existing_context.get("excluded_regions", []))
    else:
        region_source_file = (
            Path(args.region_source_file).resolve()
            if args.region_source_file
            else resolve_default_region_source_file(entrypoint_name).resolve()
        )
        if not region_source_file.exists():
            raise RuntimeError(f"Required path does not exist: {region_source_file}")
        regions = load_region_catalog(region_source_file)
        excluded_regions = parse_exclude_regions(args.exclude_regions)
        region = select_region_definition(regions, args.region_key, args.auto_pick_rank, excluded_regions)
        selected_chrom = normalize_chromosome(args.chrom) if args.chrom else region.chrom
        if selected_chrom != region.chrom:
            raise RuntimeError(
                f"--chrom {args.chrom} does not match selected region chromosome {region.chrom} from {region.region_key}"
            )
        outdir = selected_outdir if selected_outdir is not None else default_outdir_for_region(region)
        posterior_context_percentiles = parse_percentiles(args.posterior_context_percentiles)
        target_ne0_grid = []
        diagnostic_ne0_grid = []
        include_q97p5_smoke = not args.skip_q97p5_smoke

    outdir.mkdir(parents=True, exist_ok=True)
    input_root = (Path(args.input_root) if args.input_root else default_input_root(args.pop)).resolve()
    normalized_table = (Path(args.normalized_table) if args.normalized_table else default_normalized_table(args.pop)).resolve()
    phlash_pickle = (Path(args.phlash_pickle) if args.phlash_pickle else default_phlash_pickle(args.pop)).resolve()
    s_file = input_root / f"chr{selected_chrom}_s_file.txt"
    t_file = input_root / f"chr{selected_chrom}_t_file.txt"
    o_file = input_root / f"chr{selected_chrom}_o_file.txt"
    b_file = input_root / f"chr{selected_chrom}_b_file.txt"

    for path in [Path(args.ms_make_dir), Path(args.ms_binary), Path(args.backward_script)]:
        if not path.exists():
            raise RuntimeError(f"Required path does not exist: {path}")
    for path in [s_file, t_file, o_file, b_file]:
        if not path.exists():
            raise RuntimeError(f"Required path does not exist: {path}")

    needed_percentiles = sorted(set(posterior_context_percentiles + [50.0] + ([97.5] if include_q97p5_smoke else [])))
    posterior_npz = outdir / "phlash_percentiles.npz"
    if posterior_npz.exists():
        try:
            curves = load_existing_percentile_curves(posterior_npz, needed_percentiles)
        except KeyError:
            if not phlash_pickle.exists():
                raise RuntimeError(f"Required path does not exist: {phlash_pickle}")
            if not Path(args.phlash_python).exists():
                raise RuntimeError(f"Required path does not exist: {args.phlash_python}")
            curves = extract_percentile_curves(Path(args.phlash_python).resolve(), phlash_pickle, needed_percentiles, posterior_npz)
    else:
        if not phlash_pickle.exists():
            raise RuntimeError(f"Required path does not exist: {phlash_pickle}")
        if not Path(args.phlash_python).exists():
            raise RuntimeError(f"Required path does not exist: {args.phlash_python}")
        curves = extract_percentile_curves(Path(args.phlash_python).resolve(), phlash_pickle, needed_percentiles, posterior_npz)
    t_grid = np.asarray(curves[0.0], dtype=float)
    q50_curve = np.asarray(curves[50.0], dtype=float)
    q50_ne0 = float(q50_curve[0])

    if existing_context is None:
        target_ne0_grid = parse_ne0_grid(args.target_ne0_grid, q50_ne0)
        diagnostic_ne0_grid = [] if args.skip_high_diagnostics else parse_ne0_grid(args.diagnostic_ne0_grid, q50_ne0)

    panel_t_file = outdir / "region_panel.t.tsv"
    sentinel_manifest_path = outdir / "sentinel_manifest.tsv"
    if sentinel_manifest_path.exists():
        sentinels = load_sentinel_manifest(sentinel_manifest_path)
    else:
        if not normalized_table.exists():
            raise RuntimeError(f"Required path does not exist: {normalized_table}")
        region_rows = load_region_rows(normalized_table, region.chrom, region.start, region.end)
        exact_row = find_region_row_by_id(region_rows, region.top_snp_id)
        use_exact_anchor = (
            exact_row is not None
            and safe_int(exact_row.get("POS")) is not None
            and has_t_row(t_file, region.top_snp_id, int(exact_row["POS"]))
        )
        sentinels = select_sentinel_panel(region_rows, region, args.spacing_bp, use_exact_anchor)
        t_lines = load_t_lines_for_sentinels(t_file, sentinels)
        for sentinel in sentinels:
            sentinel.t_line = t_lines[sentinel.snp_id]
        write_tsv(
            sentinel_manifest_path,
            [
                "role",
                "snp_id",
                "chrom",
                "pos",
                "aa",
                "da",
                "daf",
                "maf",
                "ng0",
                "ng1",
                "ng2",
                "baseline_rsds",
                "baseline_norm_sds",
                "common_mean",
                "common_sd",
                "is_common_variant",
                "distance_to_top_snv_bp",
            ],
            sentinel_manifest_rows(sentinels),
        )
        write_panel_t_file(panel_t_file, sentinels)

    if not panel_t_file.exists():
        t_lines = load_t_lines_for_sentinels(t_file, sentinels)
        for sentinel in sentinels:
            sentinel.t_line = t_lines[sentinel.snp_id]
        write_panel_t_file(panel_t_file, sentinels)

    scenarios, scenario_aliases = build_scenarios(
        outdir=outdir,
        args=args,
        t_grid=t_grid,
        curves=curves,
        posterior_context_percentiles=posterior_context_percentiles,
        target_ne0_grid=target_ne0_grid,
        diagnostic_ne0_grid=diagnostic_ne0_grid,
        include_q97p5_smoke=include_q97p5_smoke,
        q50_ne0=q50_ne0,
    )
    requested_scenarios = parse_csv_labels(args.scenario_include)
    scenarios, scenario_aliases, aliases_by_canonical = filter_scenarios(
        scenarios=scenarios,
        scenario_aliases=scenario_aliases,
        requested_labels=requested_scenarios,
    )

    benchmark_context_payload = {
        "region": asdict(region),
        "region_source_file": str(region_source_file),
        "region_key_request": args.region_key,
        "auto_pick_rank": args.auto_pick_rank,
        "excluded_regions": sorted(excluded_regions),
        "posterior_context_percentiles": posterior_context_percentiles,
        "target_ne0_grid": target_ne0_grid,
        "diagnostic_ne0_grid": diagnostic_ne0_grid,
        "q50_ne0": q50_ne0,
        "sentinel_ids": [sentinel.snp_id for sentinel in sentinels],
        "anchor_mode": next((sentinel.role for sentinel in sentinels if sentinel.role in {"top_exact", "top_proxy"}), "unknown"),
        "scenario_aliases": scenario_aliases,
        "scenario_include": requested_scenarios,
    }
    if existing_context is None:
        write_json(outdir / "benchmark_context.json", benchmark_context_payload)
    write_json(outdir / "scenario_aliases.json", {"scenario_aliases": scenario_aliases})
    write_scenario_manifest(outdir / "scenario_manifest.tsv", scenarios)

    unique_frequencies = sorted(
        {float(sentinel.daf) for sentinel in sentinels}
        | {float(sentinel.daf_complement) for sentinel in sentinels}
    )
    return PreparedBenchmark(
        outdir=outdir,
        region=region,
        sentinels=sentinels,
        scenarios=scenarios,
        q50_ne0=q50_ne0,
        panel_t_file=panel_t_file,
        s_file=s_file,
        t_file=t_file,
        o_file=o_file,
        b_file=b_file,
        scenario_aliases=scenario_aliases,
        aliases_by_canonical=aliases_by_canonical,
        benchmark_context=existing_context or benchmark_context_payload,
        unique_frequencies=unique_frequencies,
    )


def finalize_benchmark(prepared: PreparedBenchmark, args: argparse.Namespace) -> None:
    cache_dir = prepared.outdir / "compute_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    panel_summary_rows: list[dict[str, object]] = []
    region_summary_rows: list[dict[str, object]] = []

    for scenario in prepared.scenarios:
        try:
            if scenario.smoke_only:
                _, _, missing = collect_gamma_fragments(
                    prepared=prepared,
                    scenario=scenario,
                    frequencies=prepared.unique_frequencies,
                    purpose="smoke",
                    args=args,
                    generate_missing=False,
                    reuse_existing_gamma=True,
                    allow_legacy_backfill=True,
                    write_manifest=True,
                )
                if missing:
                    raise RuntimeError(
                        "missing smoke gamma pieces: "
                        + ", ".join(frequency_key(value) for value in missing)
                    )
                scenario.status = "smoke_only_ok"
                scenario.note = "gamma smoke completed"
                region_summary_rows.append(build_incomplete_region_summary_row(scenario, prepared.sentinels, scenario.note))
                continue

            if scenario_requires_smoke_prepass(scenario, args):
                _, _, missing_smoke = collect_gamma_fragments(
                    prepared=prepared,
                    scenario=scenario,
                    frequencies=prepared.unique_frequencies,
                    purpose="smoke",
                    args=args,
                    generate_missing=False,
                    reuse_existing_gamma=True,
                    allow_legacy_backfill=True,
                    write_manifest=True,
                )
                if missing_smoke:
                    raise RuntimeError(
                        "missing smoke gamma pieces: "
                        + ", ".join(frequency_key(value) for value in missing_smoke)
                    )

            gamma_values, _, missing_full = collect_gamma_fragments(
                prepared=prepared,
                scenario=scenario,
                frequencies=prepared.unique_frequencies,
                purpose="full",
                args=args,
                generate_missing=False,
                reuse_existing_gamma=True,
                allow_legacy_backfill=True,
                write_manifest=True,
            )
            if missing_full:
                raise RuntimeError(
                    "missing full gamma pieces: "
                    + ", ".join(frequency_key(value) for value in missing_full)
                )

            minimal_g_file = scenario.scenario_dir / f"{scenario.file_label}.minimal_g.tsv"
            write_minimal_g_file(minimal_g_file, list(gamma_values.values()))
            output_tsv = scenario.scenario_dir / f"{scenario.file_label}.panel.sds.tsv"
            run_compute(
                sds_python=Path(args.sds_python).resolve(),
                s_file=prepared.s_file,
                t_file=prepared.panel_t_file,
                o_file=prepared.o_file,
                b_file=prepared.b_file,
                g_file=minimal_g_file,
                init=args.init,
                s_file_ncol=args.s_file_ncol,
                output_tsv=output_tsv,
                summary_csv=scenario.scenario_dir / "compute.summary.csv",
                cache_dir=cache_dir,
            )
            result_rows = load_sds_rows(output_tsv)
            scenario.status = "ok"
            scenario.note = ""
            scenario_panel_rows = build_panel_summary_rows(scenario, prepared.sentinels, result_rows, gamma_values)
            panel_summary_rows.extend(scenario_panel_rows)
            region_summary_rows.append(build_region_summary_row(scenario, scenario_panel_rows, prepared.sentinels))
        except Exception as exc:
            scenario.status = "numerically_unstable" if scenario.smoke_only else "error"
            scenario.note = str(exc)
            region_summary_rows.append(build_incomplete_region_summary_row(scenario, prepared.sentinels, scenario.note))

    write_tsv(
        prepared.outdir / "panel_summary.tsv",
        [
            "scenario",
            "scenario_type",
            "source",
            "percentile",
            "target_ne0",
            "scale_factor",
            "interval_eligible",
            "smoke_only",
            "status",
            "sentinel_role",
            "snp_id",
            "chrom",
            "pos",
            "distance_to_top_snv_bp",
            "daf",
            "daf_complement",
            "baseline_rSDS",
            "baseline_norm_SDS",
            "scenario_rSDS",
            "delta_vs_baseline_rSDS",
            "projected_norm_SDS",
            "projected_p_bothside",
            "projected_neg_log10_p",
            "same_sign_vs_baseline",
            "gamma_shape_daf",
            "gamma_shape_complement",
            "suggested_init_point",
            "present_diploid_pop_size",
        ],
        panel_summary_rows,
    )
    summary_fields = [
        "scenario",
        "scenario_type",
        "source",
        "percentile",
        "target_ne0",
        "scale_factor",
        "interval_eligible",
        "smoke_only",
        "status",
        "anchor_same_sign",
        "same_sign_count",
        "sentinel_count",
        "baseline_panel_max_abs_norm",
        "scenario_panel_max_abs_projected_norm",
        "retained_peak_fraction",
        "region_pass",
        "note",
        "present_diploid_pop_size",
    ]
    write_tsv(prepared.outdir / "region_summary.tsv", summary_fields, region_summary_rows)
    write_tsv(prepared.outdir / "summary.tsv", summary_fields, region_summary_rows)
    recommended_interval = choose_recommended_interval(region_summary_rows)
    write_json(
        prepared.outdir / "recommended_interval.json",
        {
            "recommended_interval": recommended_interval,
            "q50_ne0": prepared.q50_ne0,
        },
    )
    write_scenario_manifest(prepared.outdir / "scenario_manifest.tsv", prepared.scenarios)
    write_final_report(
        prepared.outdir / "final_report.md",
        args=args,
        region=prepared.region,
        sentinels=prepared.sentinels,
        scenarios=prepared.scenarios,
        region_rows=region_summary_rows,
        recommended_interval=recommended_interval,
        q50_ne0=prepared.q50_ne0,
    )


def run_piece_worker(prepared: PreparedBenchmark, args: argparse.Namespace) -> None:
    scenario = find_scenario(prepared.scenarios, str(args.scenario_label))
    if scenario is None:
        raise RuntimeError(f"Unknown scenario label: {args.scenario_label}")
    if args.gamma_purpose == "full" and scenario.smoke_only:
        raise RuntimeError(f"Scenario {scenario.label} is smoke-only and cannot generate full gamma pieces")
    if args.gamma_purpose == "smoke" and not scenario_requires_smoke_prepass(scenario, args):
        raise RuntimeError(f"Scenario {scenario.label} does not require smoke gamma pieces")
    resolved_frequency = resolve_frequency_from_panel(prepared.unique_frequencies, float(args.frequency))
    fragment = generate_piece_fragment(
        prepared=prepared,
        scenario=scenario,
        purpose=str(args.gamma_purpose),
        frequency=resolved_frequency,
        args=args,
        reuse_existing_gamma=not args.no_reuse_existing_gamma,
    )
    print(fragment.piece_path)


def run_chunk_worker(prepared: PreparedBenchmark, args: argparse.Namespace) -> None:
    scenario = find_scenario(prepared.scenarios, str(args.scenario_label))
    if scenario is None:
        raise RuntimeError(f"Unknown scenario label: {args.scenario_label}")
    if args.gamma_purpose == "full" and scenario.smoke_only:
        raise RuntimeError(f"Scenario {scenario.label} is smoke-only and cannot generate full gamma chunks")
    if args.gamma_purpose == "smoke" and not scenario_requires_smoke_prepass(scenario, args):
        raise RuntimeError(f"Scenario {scenario.label} does not require smoke gamma chunks")
    resolved_frequency = resolve_frequency_from_panel(prepared.unique_frequencies, float(args.frequency))
    requested_key = frequency_key(resolved_frequency)
    chunk_root = chunk_root_for(scenario.scenario_dir, str(args.gamma_purpose), resolved_frequency)
    workdir = chunk_workdir_for(chunk_root, int(args.chunk_index))
    run_checked(
        [
            "bash",
            str(RUN_GAMMA_CHUNK_SCRIPT),
            "--pop",
            args.pop,
            "--daf",
            requested_key,
            "--start-rep",
            str(int(args.chunk_start_rep)),
            "--end-rep",
            str(int(args.chunk_end_rep)),
            "--scenario-npz",
            str(scenario.scenario_npz),
            "--present-ne",
            str(int(round(float(scenario.summary["present_diploid_pop_size"])))),
            "--workdir",
            str(workdir),
            "--ms-make-dir",
            str(Path(args.ms_make_dir).resolve()),
            "--ms-binary",
            str(Path(args.ms_binary).resolve()),
            "--backward-script",
            str(Path(args.backward_script).resolve()),
        ],
        cwd=REPO_ROOT,
    )
    print(workdir)


def run_aggregate_worker(prepared: PreparedBenchmark, args: argparse.Namespace) -> None:
    scenario = find_scenario(prepared.scenarios, str(args.scenario_label))
    if scenario is None:
        raise RuntimeError(f"Unknown scenario label: {args.scenario_label}")
    if args.gamma_purpose == "full" and scenario.smoke_only:
        raise RuntimeError(f"Scenario {scenario.label} is smoke-only and cannot aggregate full gamma pieces")
    if args.gamma_purpose == "smoke" and not scenario_requires_smoke_prepass(scenario, args):
        raise RuntimeError(f"Scenario {scenario.label} does not require smoke gamma aggregation")
    resolved_frequency = resolve_frequency_from_panel(prepared.unique_frequencies, float(args.frequency))
    fragment = aggregate_piece_fragment(
        prepared=prepared,
        scenario=scenario,
        purpose=str(args.gamma_purpose),
        frequency=resolved_frequency,
        args=args,
        reuse_existing_gamma=not args.no_reuse_existing_gamma,
    )
    print(fragment.piece_path)


def run_serial_benchmark(prepared: PreparedBenchmark, args: argparse.Namespace) -> None:
    if args.gamma_generation_mode == "chunked":
        raise RuntimeError("Chunked gamma generation must be orchestrated through submit_region_ne0_parallel.py")
    reuse_existing_gamma = not args.no_reuse_existing_gamma
    for scenario in prepared.scenarios:
        try:
            if scenario.smoke_only:
                collect_gamma_fragments(
                    prepared=prepared,
                    scenario=scenario,
                    frequencies=prepared.unique_frequencies,
                    purpose="smoke",
                    args=args,
                    generate_missing=True,
                    reuse_existing_gamma=reuse_existing_gamma,
                    allow_legacy_backfill=True,
                    write_manifest=False,
                )
                continue
            if scenario_requires_smoke_prepass(scenario, args):
                collect_gamma_fragments(
                    prepared=prepared,
                    scenario=scenario,
                    frequencies=prepared.unique_frequencies,
                    purpose="smoke",
                    args=args,
                    generate_missing=True,
                    reuse_existing_gamma=reuse_existing_gamma,
                    allow_legacy_backfill=True,
                    write_manifest=False,
                )
            collect_gamma_fragments(
                prepared=prepared,
                scenario=scenario,
                frequencies=prepared.unique_frequencies,
                purpose="full",
                args=args,
                generate_missing=True,
                reuse_existing_gamma=reuse_existing_gamma,
                allow_legacy_backfill=True,
                write_manifest=False,
            )
        except Exception:
            continue
    finalize_benchmark(prepared, args)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    prepared = prepare_benchmark(args)
    if args.prepare_only:
        print(prepared.outdir)
        return 0
    if args.piece_worker:
        run_piece_worker(prepared, args)
        return 0
    if args.chunk_worker:
        run_chunk_worker(prepared, args)
        return 0
    if args.aggregate_worker:
        run_aggregate_worker(prepared, args)
        return 0
    if args.finalize_only:
        finalize_benchmark(prepared, args)
        print(prepared.outdir)
        return 0
    run_serial_benchmark(prepared, args)
    print(prepared.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
