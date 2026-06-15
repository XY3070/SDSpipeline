# Repo Layout

## Top level

- `scripts/`
  - Canonical reusable SDS scripts.
  - Kept flat for now because many scripts source helpers from the same directory.
- `benchmark/demography/`
  - Reusable demography benchmark/evaluation helpers used by the SDS pipeline.
- `config/`
  - Path contracts and environment templates.
- `templates/`
  - Templates for run/dataset provenance manifests.
- `manifests/`
  - Human-maintained registries for canonical scripts and migration mapping.
- `docs/`
  - Design, structure, migration rules, and QC contracts.

## Why `scripts/` is still flat

The current script set was copied with relative-path compatibility preserved. A second-stage cleanup can split it into `input/`, `compute/`, `postprocess/`, `demography/`, and `gamma/` once the new repo becomes the true source of edits.

Until then:

- category ownership is recorded in `manifests/canonical_scripts.tsv`
- code edits should happen here, not in the legacy `sds/scripts/`

## QC contracts

`docs/SDS_INPUT_QC_CONTRACT.md` defines the rebuilt SDS input QC contract:
variant-level gates, singleton-specific filters, mask rules, observability
requirements, and controlled comparison requirements.

## Cluster-aware scheduling

Three scripts form the auto-sense scheduling layer:

- `scripts/sense_cluster.sh` — Sourceable bash library that queries LSF
  (`bqueues`, `bhosts`, `bjobs`) and scores queues by
  `fairshare_priority × free_slots / (1 + pending_jobs)`.
  Outputs sourceable `SENSE_*` variables or a human-readable `--report`.
- `scripts/submit_sds_compute_chunked_chr.sh` — Single-chromosome submitter.
  Add `--auto` to enable cluster-aware queue selection and parameter tuning.
  Explicit user flags (`--queue`, `--chunk-rows`, etc.) always override sensed values.
- `scripts/submit_sds_genomewide.sh` — Genome-wide orchestrator that submits
  all (population × chromosome) pairs with dynamic concurrency, re-sensing
  between waves, and automatic failure retry (up to 2× with backoff).
