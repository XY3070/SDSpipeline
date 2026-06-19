#!/data/home/grp-wangyf/intern/miniforge3/envs/sds/bin/python
from __future__ import annotations

import argparse
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence, TextIO


@dataclass(frozen=True)
class PolarizedAlleles:
    ancestral: str
    derived: str
    is_flipped: bool
    derived_carriers: int | None


def open_text(path: Path) -> TextIO:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def split_semicolon_fields(raw_line: str) -> List[str]:
    fields = [field.strip() for field in raw_line.rstrip("\n").split(";")]
    if fields and fields[-1] == "":
        fields = fields[:-1]
    return fields


def split_tab_fields(raw_line: str) -> List[str]:
    fields = raw_line.rstrip("\n").split("\t")
    while fields and fields[-1] == "":
        fields.pop()
    return fields


def parse_sites_line(fields: List[str], sites_path: Path, position: int, expected_states: int) -> List[str]:
    payload = fields[1:]
    if not payload:
        raise RuntimeError(
            f"SampleBranchLengths output incomplete for position {position}: "
            f"{sites_path} has no allele payload on the site row."
        )

    if len(payload) == 1:
        states = list(payload[0].strip().upper())
    else:
        states = [token.strip().upper() for token in payload if token.strip()]

    if len(states) != expected_states:
        raise ValueError(
            f"State count mismatch in {sites_path} for position {position}: "
            f"expected {expected_states} leaves from the NAMES row, found {len(states)} states."
        )
    return states


def read_site_states(sites_path: Path, position: int) -> List[str]:
    with open_text(sites_path) as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]

    if not lines:
        raise RuntimeError(
            f"SampleBranchLengths output incomplete for position {position}: {sites_path} is empty."
        )
    if len(lines) < 2:
        raise RuntimeError(
            f"SampleBranchLengths output incomplete for position {position}: "
            f"{sites_path} must contain at least NAMES and REGION rows."
        )

    name_fields = split_tab_fields(lines[0])
    if not name_fields or name_fields[0] != "NAMES":
        raise ValueError(
            f"Malformed .sites header in {sites_path}: expected first row to start with 'NAMES'."
        )

    region_fields = split_tab_fields(lines[1])
    if not region_fields or region_fields[0] != "REGION":
        raise ValueError(
            f"Malformed .sites header in {sites_path}: expected second row to start with 'REGION'."
        )

    expected_states = len(name_fields) - 1
    if expected_states <= 0:
        raise ValueError(f"Malformed NAMES row in {sites_path}: no leaf identifiers were found.")

    if len(lines) == 2:
        raise RuntimeError(
            f"SampleBranchLengths output incomplete for position {position}: "
            f"{sites_path} only contains NAMES and REGION rows."
        )

    matched_states: List[List[str]] = []
    available_positions: List[int] = []

    for raw_line in lines[2:]:
        fields = split_tab_fields(raw_line)
        if not fields:
            continue
        try:
            row_position = int(fields[0])
        except ValueError as exc:
            raise ValueError(
                f"Malformed site row in {sites_path}: expected an integer position, got {fields[0]!r}."
            ) from exc

        available_positions.append(row_position)
        if row_position != position:
            continue

        matched_states.append(parse_sites_line(fields, sites_path, position, expected_states))

    if not matched_states:
        sample_text = ", ".join(str(pos) for pos in available_positions[:5]) or "none"
        raise FileNotFoundError(
            f"Position {position} not found in {sites_path}; available site rows: {sample_text}."
        )
    if len(matched_states) != 1:
        raise ValueError(
            f"Expected exactly one site row for position {position} in {sites_path}, "
            f"but found {len(matched_states)}."
        )
    return matched_states[0]


def resolve_mut_path(sites_path: Path, mut_arg: str | None) -> Path:
    if mut_arg:
        return Path(mut_arg)
    return sites_path.with_suffix(".mut")


def read_polarized_alleles(mut_path: Path, position: int) -> PolarizedAlleles:
    with open_text(mut_path) as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]

    if not lines:
        raise RuntimeError(
            f"SampleBranchLengths output incomplete for position {position}: {mut_path} is empty."
        )

    header = split_semicolon_fields(lines[0])
    if not header:
        raise ValueError(f"Malformed .mut header in {mut_path}: header row is empty.")

    required_columns = ["pos_of_snp", "is_flipped", "ancestral_allele/alternative_allele"]
    missing_columns = [column for column in required_columns if column not in header]
    if missing_columns:
        raise ValueError(
            f"Malformed .mut header in {mut_path}: missing required column(s) "
            f"{', '.join(missing_columns)}."
        )

    if len(lines) == 1:
        raise RuntimeError(
            f"SampleBranchLengths output incomplete for position {position}: "
            f"{mut_path} contains only the header row."
        )

    idx_position = header.index("pos_of_snp")
    idx_flipped = header.index("is_flipped")
    idx_alleles = header.index("ancestral_allele/alternative_allele")

    matched_rows: List[List[str]] = []
    available_positions: List[int] = []

    for raw_line in lines[1:]:
        fields = split_semicolon_fields(raw_line)
        if len(fields) < len(header):
            raise ValueError(
                f"Malformed .mut row in {mut_path}: expected at least {len(header)} columns "
                f"from the header, found {len(fields)}."
            )

        try:
            row_position = int(fields[idx_position])
        except ValueError as exc:
            raise ValueError(
                f"Malformed pos_of_snp value in {mut_path}: {fields[idx_position]!r} is not an integer."
            ) from exc

        available_positions.append(row_position)
        if row_position == position:
            matched_rows.append(fields)

    if not matched_rows:
        sample_text = ", ".join(str(pos) for pos in available_positions[:5]) or "none"
        raise FileNotFoundError(
            f"Position {position} not found in {mut_path}; available mutation rows: {sample_text}."
        )
    if len(matched_rows) != 1:
        raise ValueError(
            f"Expected exactly one mutation row for position {position} in {mut_path}, "
            f"but found {len(matched_rows)}."
        )

    record = matched_rows[0]
    is_flipped_raw = record[idx_flipped]
    if is_flipped_raw not in {"0", "1"}:
        raise ValueError(
            f"Malformed is_flipped value in {mut_path} for position {position}: {is_flipped_raw!r}."
        )
    is_flipped = is_flipped_raw == "1"

    allele_field = record[idx_alleles].upper()
    allele_parts = [part.strip() for part in allele_field.split("/")]
    if len(allele_parts) != 2 or not allele_parts[0] or not allele_parts[1]:
        raise ValueError(
            f"Malformed ancestral_allele/alternative_allele value in {mut_path} "
            f"for position {position}: {record[idx_alleles]!r}."
        )

    ancestral_candidate, alternative_candidate = allele_parts
    if is_flipped:
        ancestral = alternative_candidate
        derived = ancestral_candidate
    else:
        ancestral = ancestral_candidate
        derived = alternative_candidate

    if ancestral == derived:
        raise ValueError(
            f"Unable to polarize alleles in {mut_path} for position {position}: "
            f"ancestral and derived alleles both resolve to {ancestral!r}."
        )

    derived_carriers: int | None = None
    if "downstream_allele" in header:
        downstream_idx = header.index("downstream_allele")
        population_counts: List[int] = []
        for field in record[downstream_idx + 1 :]:
            if field == "":
                continue
            try:
                population_counts.append(int(field))
            except ValueError as exc:
                raise ValueError(
                    f"Malformed population carrier count in {mut_path} for position {position}: {field!r}."
                ) from exc

        if len(population_counts) == 1:
            derived_carriers = population_counts[0]

    return PolarizedAlleles(
        ancestral=ancestral,
        derived=derived,
        is_flipped=is_flipped,
        derived_carriers=derived_carriers,
    )


def map_states_to_binary(
    raw_states: Sequence[str], polarized: PolarizedAlleles, sites_path: Path, position: int
) -> List[str]:
    binary_states: List[str] = []
    unexpected_states: List[str] = []

    for index, state in enumerate(raw_states):
        normalized = state.upper()
        if normalized == polarized.ancestral:
            binary_states.append("0")
        elif normalized == polarized.derived:
            binary_states.append("1")
        else:
            unexpected_states.append(f"{index}:{normalized}")

    if unexpected_states:
        preview = ", ".join(unexpected_states[:10])
        raise ValueError(
            f"Unexpected allele state(s) in {sites_path} for position {position}: {preview}. "
            f"Expected only ancestral={polarized.ancestral!r} or derived={polarized.derived!r}."
        )

    return binary_states


def read_haps_frequency(haps_path: Path, position: int, is_flipped: bool) -> float:
    with open_text(haps_path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split()
            if len(fields) < 6:
                continue
            try:
                row_position = int(fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"Malformed position field in {haps_path}: {fields[2]!r} is not an integer."
                ) from exc
            if row_position != position:
                continue

            try:
                hap_values = [int(value) for value in fields[5:]]
            except ValueError as exc:
                raise ValueError(
                    f"Malformed haplotype state in {haps_path} for position {position}: "
                    f"expected 0/1 values after column 5."
                ) from exc

            if not hap_values:
                raise ValueError(f"No haplotype states found at position {position} in {haps_path}.")

            alt_frequency = sum(hap_values) / float(len(hap_values))
            return 1.0 - alt_frequency if is_flipped else alt_frequency

    raise FileNotFoundError(f"Position {position} not found in {haps_path}.")


def resolve_pop_frequency(
    polarized: PolarizedAlleles, n_haplotypes: int, haps_path: Path, position: int
) -> float:
    if polarized.derived_carriers is not None:
        if polarized.derived_carriers < 0 or polarized.derived_carriers > n_haplotypes:
            raise ValueError(
                f"Derived carrier count out of range for position {position}: "
                f"{polarized.derived_carriers} carriers in .mut but {n_haplotypes} haplotypes in .sites."
            )
        return polarized.derived_carriers / float(n_haplotypes)
    return read_haps_frequency(haps_path, position, polarized.is_flipped)


def write_lines(path: Path, lines: Iterable[str]) -> None:
    with path.open("wt") as handle:
        for line in lines:
            handle.write(line)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build CLUES2 derived-state and frequency inputs from Relate outputs."
    )
    parser.add_argument("--sites", required=True, help="Relate .sites file from SampleBranchLengths --format n")
    parser.add_argument("--mut", help="Matching Relate .mut file from SampleBranchLengths; defaults to the .sites prefix")
    parser.add_argument("--haps", required=True, help="Prepared Relate .haps(.gz) file with ancestral allele coded as 0")
    parser.add_argument("--position", required=True, type=int, help="Target base position")
    parser.add_argument("--out-prefix", required=True, help="Output prefix")
    args = parser.parse_args()

    sites_path = Path(args.sites)
    mut_path = resolve_mut_path(sites_path, args.mut)
    haps_path = Path(args.haps)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    polarized = read_polarized_alleles(mut_path, args.position)
    raw_states = read_site_states(sites_path, args.position)
    derived_states = map_states_to_binary(raw_states, polarized, sites_path, args.position)
    pop_freq = resolve_pop_frequency(polarized, len(derived_states), haps_path, args.position)

    derived_path = Path(f"{out_prefix}_derived.txt")
    popfreq_path = Path(f"{out_prefix}_popfreq.txt")
    summary_path = Path(f"{out_prefix}_site_summary.tsv")

    n_derived = sum(int(state) for state in derived_states)

    write_lines(derived_path, (f"{state}\n" for state in derived_states))
    popfreq_path.write_text(f"{pop_freq:.12f}\n")
    summary_path.write_text(
        "position\tn_haplotypes\tn_derived\tpopFreq\tis_flipped\tancestral_allele\tderived_allele\tderived_carriers_in_mut\n"
        f"{args.position}\t{len(derived_states)}\t{n_derived}\t{pop_freq:.12f}\t"
        f"{int(polarized.is_flipped)}\t{polarized.ancestral}\t{polarized.derived}\t"
        f"{'' if polarized.derived_carriers is None else polarized.derived_carriers}\n"
    )


if __name__ == "__main__":
    main()
