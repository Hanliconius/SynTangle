#!/bin/bash
set -euo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$REPO"
export SCT_MAMBA=${MAMBA_EXE:-$(type -P micromamba)}
test -x "$SCT_MAMBA"
export SCT_REPO="$REPO"
export SCT_INPUT_ROOT="${SCT_INPUT_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/SynTangle_pegasus_test/local_results/genespace_stress}"
test -f "$SCT_INPUT_ROOT/benchmark_manifest.tsv"
mkdir -p logs local_results
export SCT_ORDER_ROOT=$(mktemp -d "$REPO/local_results/order_probe_XXXXXX")
git rev-parse HEAD > "$SCT_ORDER_ROOT/commit.txt"
job=$(sbatch --parsable --chdir="$REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_order_probe
#SBATCH --partition=nano
#SBATCH --array=0-8%9
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:05:00
#SBATCH --output=logs/st_order_probe.%A_%a.out
#SBATCH --error=logs/st_order_probe.%A_%a.err
set -euo pipefail
eval "$("$SCT_MAMBA" shell hook --shell bash)"
micromamba activate syntangle_test
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1
export PYTHONPATH="$SCT_REPO/src"
python validation/benchmark/diagnose_order_relaxation.py \
  --history "$SCT_REPO/local_results" --root "$SCT_INPUT_ROOT" --output "$SCT_ORDER_ROOT" --index "$SLURM_ARRAY_TASK_ID"
SLURM
)
printf 'Ordering diagnostic array: %s\nResults: %s\n' "${job%%;*}" "$SCT_ORDER_ROOT"
