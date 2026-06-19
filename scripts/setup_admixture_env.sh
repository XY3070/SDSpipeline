#!/usr/bin/env bash
# ============================================================================
# setup_admixture_env.sh
#
# Bootstrap the uv-managed Python environment for ADAMIXTURE.
# Uses `uv` (not pip/mamba/conda) per project convention.
#
# Idempotent: safe to re-run.
# ============================================================================
set -euo pipefail

if [[ -z "${SDS_PIPELINE_ROOT:-}" ]]; then
  SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  # shellcheck disable=SC1090
  source "${SELF_DIR}/../config/paths.env"
fi

: "${SDS_ADMIXTURE_ENV_PREFIX:?set SDS_ADMIXTURE_ENV_PREFIX in config/paths.env}"

# Pinned Python + ADAMIXTURE versions. Bump deliberately; the lockfile keeps
# runs reproducible until then.
PYTHON_SPEC="3.11"
ADAMIXTURE_SPEC="adamixture>=0.1"

PYPROJECT_DIR="${SDS_ADMIXTURE_ENV_PREFIX%/*}"    # parent of env prefix
PYPROJECT="${PYPROJECT_DIR}/pyproject.toml"
LOCKFILE="${PYPROJECT_DIR}/uv.lock"

mkdir -p "$PYPROJECT_DIR"

if [[ ! -f "$PYPROJECT" ]]; then
  cat > "$PYPROJECT" <<EOF
[project]
name = "admixture-9k-env"
version = "0.1.0"
description = "uv-managed env for ADAMIXTURE (9k Chinese cohort stratification)"
requires-python = ">=${PYTHON_SPEC}"
dependencies = [
    "${ADAMIXTURE_SPEC}",
    "numpy",
    "matplotlib",
    "pandas",
]
EOF
  echo "[setup] wrote $PYPROJECT"
fi

if ! command -v uv >/dev/null 2>&1; then
  echo "[setup] uv not on PATH; installing to user-local bin." >&2
  curl -LsSf https://astral.sh/uv/install.sh | sh
  export PATH="$HOME/.local/bin:$PATH"
fi

cd "$PYPROJECT_DIR"
uv venv --python "$PYTHON_SPEC" "$SDS_ADMIXTURE_ENV_PREFIX"
uv sync --locked 2>/dev/null || uv sync
uv pip install --upgrade "adamixture"

echo ""
echo "[setup] environment ready at: $SDS_ADMIXTURE_ENV_PREFIX"
echo "[setup] entrypoints:"
ls -1 "$SDS_ADMIXTURE_ENV_PREFIX/bin/" | grep -E '^(adamixture|python)' | sed 's/^/    /'
