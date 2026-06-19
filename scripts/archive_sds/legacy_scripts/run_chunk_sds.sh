#!/bin/bash
set -e

# 参数由 parallel 传入
CHUNK_FILE="$1"
OUT_FILE="$2"
S_FILE="$3"
O_FILE="$4"
B_FILE="$5"
G_FILE="$6"

# 检查 uv 是否可用，如果不可用尝试加载
if ! command -v uv &> /dev/null; then
    export PATH="$HOME/.cargo/bin:$HOME/.local/bin:$PATH"
fi

# 运行计算
# init = 0.00001
# s_file_ncol = 20000
uv run /share/home/grp-wangyf/xuyuan/sdSPY/scripts/compute_SDS.py \
    "$S_FILE" \
    "$CHUNK_FILE" \
    "$O_FILE" \
    "$B_FILE" \
    "$G_FILE" \
    0.00001 \
    20000 \
    --output "$OUT_FILE" > /dev/null 2>&1