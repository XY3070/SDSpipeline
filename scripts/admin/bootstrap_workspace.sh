#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
WORKSPACE_ROOT="${1:-${SDS_WORKSPACE_ROOT:-$PIPELINE_ROOT/../SDSworkspace}}"

mkdir -p \
    "$WORKSPACE_ROOT/input" \
    "$WORKSPACE_ROOT/input/raw" \
    "$WORKSPACE_ROOT/input/raw/vcf" \
    "$WORKSPACE_ROOT/input/freeze" \
    "$WORKSPACE_ROOT/input/freeze/sample_lists" \
    "$WORKSPACE_ROOT/input/freeze/cohorts" \
    "$WORKSPACE_ROOT/input/derived" \
    "$WORKSPACE_ROOT/input/reference" \
    "$WORKSPACE_ROOT/results" \
    "$WORKSPACE_ROOT/results/production" \
    "$WORKSPACE_ROOT/results/production/sds_input" \
    "$WORKSPACE_ROOT/results/production/sds_output" \
    "$WORKSPACE_ROOT/results/production/demography" \
    "$WORKSPACE_ROOT/results/production/gamma" \
    "$WORKSPACE_ROOT/results/audit" \
    "$WORKSPACE_ROOT/results/audit/raw_header_intersections" \
    "$WORKSPACE_ROOT/results/legacy" \
    "$WORKSPACE_ROOT/results/figures_staging" \
    "$WORKSPACE_ROOT/runs" \
    "$WORKSPACE_ROOT/oneoff" \
    "$WORKSPACE_ROOT/provenance" \
    "$WORKSPACE_ROOT/provenance/runs" \
    "$WORKSPACE_ROOT/provenance/datasets" \
    "$WORKSPACE_ROOT/logs" \
    "$WORKSPACE_ROOT/tmp" \
    "$WORKSPACE_ROOT/cache" \
    "$WORKSPACE_ROOT/external"

cat <<EOF
Workspace bootstrap complete:
  $WORKSPACE_ROOT/input
  $WORKSPACE_ROOT/input/raw
  $WORKSPACE_ROOT/input/raw/vcf
  $WORKSPACE_ROOT/input/freeze
  $WORKSPACE_ROOT/input/freeze/sample_lists
  $WORKSPACE_ROOT/input/freeze/cohorts
  $WORKSPACE_ROOT/input/derived
  $WORKSPACE_ROOT/input/reference
  $WORKSPACE_ROOT/results
  $WORKSPACE_ROOT/results/production
  $WORKSPACE_ROOT/results/production/sds_input
  $WORKSPACE_ROOT/results/production/sds_output
  $WORKSPACE_ROOT/results/production/demography
  $WORKSPACE_ROOT/results/production/gamma
  $WORKSPACE_ROOT/results/audit
  $WORKSPACE_ROOT/results/audit/raw_header_intersections
  $WORKSPACE_ROOT/results/legacy
  $WORKSPACE_ROOT/results/figures_staging
  $WORKSPACE_ROOT/runs
  $WORKSPACE_ROOT/oneoff
  $WORKSPACE_ROOT/provenance
  $WORKSPACE_ROOT/provenance/runs
  $WORKSPACE_ROOT/provenance/datasets
  $WORKSPACE_ROOT/logs
  $WORKSPACE_ROOT/tmp
  $WORKSPACE_ROOT/cache
  $WORKSPACE_ROOT/external
EOF
