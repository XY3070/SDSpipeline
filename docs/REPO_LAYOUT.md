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
  - Design, structure, and migration rules.

## Why `scripts/` is still flat

The current script set was copied with relative-path compatibility preserved. A second-stage cleanup can split it into `input/`, `compute/`, `postprocess/`, `demography/`, and `gamma/` once the new repo becomes the true source of edits.

Until then:

- category ownership is recorded in `manifests/canonical_scripts.tsv`
- code edits should happen here, not in the legacy `sds/scripts/`

