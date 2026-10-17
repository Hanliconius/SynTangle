#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
# Recover the two currently runnable cases; retain completed Ficus comparisons.
export SCT_SIMPLE_CASES="1 2"
bash "$repo/validation/real_data/submit_simple_reconstruction.sh" --resume "${1:-$repo/local_results/simple_reconstructed_QI1cVs}"
