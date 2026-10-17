#!/bin/bash
# Run from an isolated checkout: do not update code under the original array.
set -euo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$REPO"
export SCT_MAMBA=${MAMBA_EXE:-$(type -P micromamba)}
test -x "$SCT_MAMBA"
export SCT_REPO="$REPO"
export SCT_INPUT_ROOT="${SCT_INPUT_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/SynTangle_pegasus_test/local_results/genespace_stress}"
export SCT_BENCH_ROOT="$REPO/local_results/scoring_benchmark_$(date +%Y%m%d_%H%M%S)"
mkdir -p logs "$SCT_BENCH_ROOT/baseline"
test -f "$SCT_INPUT_ROOT/benchmark_manifest.tsv"
git archive dbd581b05bb0432fc18a5352a3cb0dee4ec55567 src | tar -x -C "$SCT_BENCH_ROOT/baseline"
git rev-parse HEAD > "$SCT_BENCH_ROOT/cached_commit.txt"
printf '%s\n' dbd581b05bb0432fc18a5352a3cb0dee4ec55567 > "$SCT_BENCH_ROOT/baseline_commit.txt"
job=$(sbatch --parsable --chdir="$REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_score_pair
#SBATCH --partition=nano
#SBATCH --array=0-8%9
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:06:00
#SBATCH --output=logs/st_score_pair.%A_%a.out
#SBATCH --error=logs/st_score_pair.%A_%a.err
set -euo pipefail
eval "$("$SCT_MAMBA" shell hook --shell bash)"
micromamba activate syntangle_test
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export PYTHONUNBUFFERED=1
python validation/benchmark/benchmark_local_scoring.py \
    --root "$SCT_INPUT_ROOT" --baseline "$SCT_BENCH_ROOT/baseline" \
    --output "$SCT_BENCH_ROOT" --index "$SLURM_ARRAY_TASK_ID" \
    --steps 5 --budget 120
SLURM
)
printf 'Paired scoring array: %s\nResults: %s\n' "${job%%;*}" "$SCT_BENCH_ROOT"
