# Migration Map

This document defines the first-stage split from the current mixed workspace into `SDSlog/`, `SDSpipeline/`, and `SDSworkspace/`.

## Principles

1. Do not move large historical data before the new boundaries are documented.
2. Copy reusable scripts into `SDSpipeline/` first; retire or relocate legacy copies later.
3. New runs should target `SDSworkspace/` even before all historical data are migrated.

## First-stage mapping

| Current location | New owner | Action |
| --- | --- | --- |
| `SDSlog/` | `SDSlog/` | Keep as-is; remains the reporting repo |
| `sds/scripts/` | `SDSpipeline/scripts/` | Copy reusable/canonical scripts; legacy tree becomes transitional |
| `benchmark/demography/*.py` | `SDSpipeline/benchmark/demography/` | Copy reusable benchmark/evaluation helpers |
| `sds/SDS_SCRIPT_PLAN.md` | `SDSpipeline/docs/LEGACY_SDS_SCRIPT_PLAN.md` | Preserve as design context |
| `sds/data/processed/*` | `SDSworkspace/results/legacy/` | Historical outputs should be treated as workspace artifacts, not repo contents |
| `sds/data/vcf/*`, `raw/*` | `SDSworkspace/input/` | Future canonical inputs belong here |
| ad hoc recovery/debug scripts | `SDSworkspace/oneoff/` | Keep out of git unless promoted into reusable logic |
| queue logs / chunk submit roots / temp dirs | `SDSworkspace/runs/` | Keep out of git |

## Promotion rule

If a one-off script is used more than once or becomes part of a standard workflow, promote it from `SDSworkspace/oneoff/` into `SDSpipeline/` and register it in `manifests/canonical_scripts.tsv`.

