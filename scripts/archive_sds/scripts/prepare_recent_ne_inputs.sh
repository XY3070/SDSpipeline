#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
# shellcheck source=/data/home/grp-wangyf/xuyuan/sds/scripts/common_env.sh
source "$SCRIPT_DIR/common_env.sh"
activate_relate_runtime

SDS_PYTHON="${SDS_ENV_PREFIX}/bin/python"

exec "$SDS_PYTHON" "$PROJECT_ROOT/benchmark/demography/prepare_recent_ne_inputs.py" "$@"
