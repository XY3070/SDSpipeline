# SDSpipeline

`SDSpipeline` is the canonical git-tracked repository for the reusable SDS pipeline.

It is intentionally separated from:

- `SDSlog/`: report-ready logs, plots, and writing materials
- `SDSworkspace/`: large inputs, outputs, one-off scripts, run roots, and provenance records

## Server portability

This repo should remain portable across servers.

- Reusable logic stays in git.
- Server-local paths stay in `config/paths.env`, which is intentionally untracked.
- Large files, queue logs, run directories, and temporary one-off helpers stay in `SDSworkspace/`.
- If one server has extra data or tools, capture that in `SDSworkspace/provenance/` rather than hardcoding the path into pipeline scripts.

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

Before running the pipeline on a new server:

1. Copy `config/paths.env.example` to `config/paths.env`.
2. Point it at the local `SDSlog/`, `SDSpipeline/`, `SDSworkspace/`, runtime env, and external tool roots.
3. Bootstrap the non-git workspace with `scripts/admin/bootstrap_workspace.sh`.
