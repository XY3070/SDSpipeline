# Workspace Contract

This project now runs across three components:

## 1. `SDSlog/`

Purpose:

- report-facing facts
- paper-writing inputs
- result and plot indexes
- high-level conclusions

Never store here:

- large raw results
- run caches
- queue scratch

## 2. `SDSpipeline/`

Purpose:

- reusable scripts
- stable workflow logic
- templates and schemas
- migration manifests

Never store here:

- big data
- one-off fix scripts
- run-local outputs

## 3. `SDSworkspace/`

Purpose:

- `input/`: raw/frozen input files and derived analysis inputs
- `results/`: finalized and intermediate outputs
- `runs/`: per-run job roots, logs, temporary working directories
- `oneoff/`: throwaway scripts used for a single debug/audit/recovery task
- `provenance/`: run manifests, dataset manifests, and parameter snapshots
- `cache/`, `tmp/`, `external/`: large or transient supporting content

## Required provenance rule

Every production or audit run in `SDSworkspace/` must have a matching manifest under `SDSworkspace/provenance/` that records:

- run ID
- pipeline repo path and commit
- entrypoint script
- input paths
- output paths
- cohort freeze
- masks
- gamma source
- parameters

This is the mechanism that should prevent another `olddefault`-style untraceable artifact.

