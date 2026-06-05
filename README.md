# SDSpipeline

`SDSpipeline` is the canonical git-tracked repository for the reusable SDS pipeline.

It is intentionally separated from:

- `SDSlog/`: report-ready logs, plots, and writing materials
- `SDSworkspace/`: large inputs, outputs, one-off scripts, run roots, and provenance records

## What belongs here

- Stable pipeline entrypoints
- Reusable helper scripts
- Config templates
- Workflow/manifests/contracts
- Small example metadata

## What does not belong here

- Large input/result files
- Scratch directories
- Queue logs from specific runs
- One-off recovery scripts
- Irreproducible ad hoc artifacts

## Current migration status

This repo is seeded from the current workspace so that future edits can move here first.

- Canonical script snapshot: `scripts/`
- Demography support scripts: `benchmark/demography/`
- Structure and migration rules: `docs/`, `config/`, `manifests/`, `templates/`

Edits to reusable pipeline logic should now happen in `SDSpipeline/` first, not under the legacy `sds/` tree.

