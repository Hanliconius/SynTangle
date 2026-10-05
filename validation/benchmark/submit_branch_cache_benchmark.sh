#!/bin/bash
# Submit from an isolated worktree; never modify either running checkout.
set -euo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$REPO"
export SCT_MAMBA=${MAMBA_EXE:-$(type -P micromamba)}
test -x "$SCT_MAMBA"
export SCT_REPO="$REPO"
export SCT_INPUT_ROOT="${SCT_INPUT_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/SynTangle_pegasus_test/local_results/genespace_stress}"
test -f "$SCT_INPUT_ROOT/benchmark_manifest.tsv"
mkdir -p logs local_results
export SCT_PIPE_ROOT=$(mktemp -d "$REPO/local_results/branch_cache_pair_XXXXXX")
mkdir -p "$SCT_PIPE_ROOT/previous"
git archive 5206c63c0700b6dd572a977f9c500e712f5b5af3 src | tar -x -C "$SCT_PIPE_ROOT/previous"
git rev-parse HEAD > "$SCT_PIPE_ROOT/pipeline_commit.txt"
printf '%s\n' 5206c63c0700b6dd572a977f9c500e712f5b5af3 > "$SCT_PIPE_ROOT/previous_commit.txt"
job=$(sbatch --parsable --chdir="$REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_branch_cache
#SBATCH --partition=nano
#SBATCH --array=0-8%9
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:06:00
#SBATCH --output=logs/st_branch_cache.%A_%a.out
#SBATCH --error=logs/st_branch_cache.%A_%a.err
set -euo pipefail
eval "$("$SCT_MAMBA" shell hook --shell bash)"
micromamba activate syntangle_test
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export PYTHONUNBUFFERED=1
python validation/benchmark/benchmark_pipeline.py \
    --root "$SCT_INPUT_ROOT" --baseline "$SCT_PIPE_ROOT/previous" \
    --output "$SCT_PIPE_ROOT" --index "$SLURM_ARRAY_TASK_ID" \
    --budget 120 --nodes 5 --require-identical
SLURM
)
printf 'Pipeline comparison array: %s\nResults: %s\n' "${job%%;*}" "$SCT_PIPE_ROOT"
