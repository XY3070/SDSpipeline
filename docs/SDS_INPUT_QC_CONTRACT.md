# SDS Input QC Contract

This document defines the minimum input QC contract for the rebuilt SDS
pipeline. It refines the first-principles workflow in
`SDSlog/logs/sds/2026-06-05_first_principles_rebuild_workflow_compare.md`.

## Scope

The SDS input layer has two distinct universes:

- `test SNP universe`: variants eligible for `t_file`
- `singleton universe`: variants eligible for `s_file`

Both universes must share the same cohort freeze, reference build, chromosome
mask release, and sample order. A production run must record those IDs in a
provenance manifest.

## Variant-Level QC

Default thresholds for the rebuilt pipeline:

| Metric | Default threshold | Implementation note |
| --- | --- | --- |
| Site call rate | `>=0.90` | Equivalent to `F_MISSING<=0.10` after cohort subsetting. |
| VQSR / site filter | `PASS` | If `VQSLOD` is unavailable, use the caller's `PASS` status as the equivalent gate. |
| Read depth | `0.5x-2.0x` of chromosome/site-universe mean | Prefer `INFO/DP_AVG` or equivalent cohort-level DP. |
| QUAL | `>=56` | If the callset uses a different calibrated score, document the equivalent threshold. |
| Hardy-Weinberg / excess heterozygosity | `P>1e-6` | Use `INFO/ExcHet` when available; otherwise compute a cohort-specific HWE exact test. |
| Low-complexity / low-confidence regions | exclude | Use an explicit BED mask. |

This contract does not assume every VCF is a raw GATK VQSR output. The current
Graphtyper/PLINK VCFs expose `PASS`, `QUAL`, `DP_AVG`, `AC/AF`, and `ExcHet`,
so the rebuilt workflow applies equivalent gates from those fields.

## Singleton-Specific QC

Singletons drive SDS directly through singleton-to-test-SNP distance
distributions. They therefore require stricter QC than generic SNPs.

Default singleton gates:

- site is inside the callable analysis genome
- site passes the variant-level QC above
- `AC=1` after cohort subsetting
- `F_MISSING<=0.005`, with an optional absolute missing-sample cap
- site is outside centromere, pericentromeric heterochromatin, subtelomere,
  low-mappability, segmental-duplication, and low-complexity masks
- site is not in a singleton-density outlier window
- carrier sample is not an excess-singleton outlier sample

Default singleton-density filter:

- window size: `20 kb`
- step size: `10 kb`
- remove singleton positions in windows with density `> mean + 4 SD`

Default excess-singleton sample filter:

- compute retained singleton count per sample
- remove samples with count `> mean + 4 SD`
- record removed sample IDs and counts in the run manifest

## Observability

`o_file` must not default to all ones unless the release proves that callable
coverage is effectively uniform.

The minimum rebuilt implementation uses sample-level observability weights:

- `copy`: preserve the upstream `o_file`; valid only for baseline/control arms
- `callrate`: set each sample weight to its retained test-SNP call rate

Future production releases should replace the sample-level proxy with a
position-aware callable-genome observability model when coverage tracks are
available.

## Acrocentric Chromosomes

Do not drop acrocentric chromosomes wholesale. Exclude their short arms,
centromeric/pericentromeric regions, rDNA/satellite arrays, and low-confidence
intervals through the mask layer. Retain q-arm callable regions when they pass
the same QC contract as other autosomes.

## Release Gate

A chr-level or genome-wide SDS release cannot be paper-facing unless it
publishes:

- `site_qc_metrics.tsv`
- `singleton_density_windows.tsv`
- `singleton_sample_counts.tsv`
- `excluded_samples.tsv`
- `observability.tsv`
- `input_qc_summary.tsv`

The release must also include a controlled comparison against the previous
baseline, with one-variable-at-a-time arms before any cumulative final arm.
