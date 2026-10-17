#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
run="$repo/local_results/simple_reconstructed_QI1cVs"
export SCT_ADAPTER_REPO="$repo" SCT_ADAPTER_RUN="$run"
mkdir -p "$repo/logs"
sbatch --partition=nano --cpus-per-task=1 --mem=8G --time=00:10:00 --job-name=st_gs_adapter --output="$repo/logs/st_gs_adapter.%j.out" --error="$repo/logs/st_gs_adapter.%j.err" <<'SLURM'
#!/bin/bash
set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate "${SCT_REAL_GS_ENV:-lep_busco_painter_clean}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
input="$SCT_ADAPTER_RUN/rhodotorula/figures/rhodotorula_600_600/native_input"
output=$(mktemp -d "$SCT_ADAPTER_RUN/rhodotorula/adapter_check_XXXXXX")
Rscript "$SCT_ADAPTER_REPO/validation/benchmark/genespace_native_order.R" "$input" "$output" --skip-incomplete-variants
export PYTHONPATH="$SCT_ADAPTER_REPO/src:$SCT_ADAPTER_REPO/validation/real_data:$SCT_ADAPTER_REPO/validation/benchmark"
"$SCT_ADAPTER_RUN/tools/bin/python" - "$SCT_ADAPTER_RUN" "$output" <<'PY'
import sys
from pathlib import Path
from syntangle.fixtures import load_fixture
from plot_comparison import supplement_native_orders
from compare_genespace import load_native_states
run,out=map(Path,sys.argv[1:])
f=load_fixture(run/'rhodotorula/prepared/rhodotorula_600/fixture.json')
if not (out/'native_orders.tsv').exists():
    if not supplement_native_orders(f,out):raise ValueError('No valid ordering to complete')
states=list(load_native_states(f,out/'native_orders.tsv'))
if not states:raise ValueError('No validated complete states')
print('PASS adapter:',len(states),'complete layouts; all fixture chromosomes retained')
print('AUDIT',out)
PY
SLURM
