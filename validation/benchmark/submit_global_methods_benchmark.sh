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
export SCT_METHOD_ROOT=$(mktemp -d "$REPO/local_results/global_methods_XXXXXX")
git rev-parse HEAD > "$SCT_METHOD_ROOT/commit.txt"
setup=$(sbatch --parsable --chdir="$REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_methods_setup
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=logs/st_methods_setup.%j.out
#SBATCH --error=logs/st_methods_setup.%j.err
set -euo pipefail
eval "$("$SCT_MAMBA" shell hook --shell bash)"
micromamba activate syntangle_test
python -m venv "$SCT_METHOD_ROOT/venv"
"$SCT_METHOD_ROOT/venv/bin/python" -m pip install -r "$SCT_REPO/validation/benchmark/global_methods_requirements.txt"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export PYTHONPATH="$SCT_REPO/src"
"$SCT_METHOD_ROOT/venv/bin/python" validation/benchmark/benchmark_global_methods.py \
  --prepare --history "$SCT_REPO/local_results" --root "$SCT_INPUT_ROOT" --output "$SCT_METHOD_ROOT"
SLURM
)
setup=${setup%%;*}
job=$(sbatch --parsable --dependency="afterok:$setup" --chdir="$REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_methods
#SBATCH --partition=cpu
#SBATCH --array=0-107
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:12:00
#SBATCH --output=logs/st_methods.%A_%a.out
#SBATCH --error=logs/st_methods.%A_%a.err
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export PYTHONPATH="$SCT_REPO/src"
"$SCT_METHOD_ROOT/venv/bin/python" validation/benchmark/benchmark_global_methods.py \
  --root "$SCT_INPUT_ROOT" --output "$SCT_METHOD_ROOT" --index "$SLURM_ARRAY_TASK_ID"
SLURM
)
job=${job%%;*}
collect=$(sbatch --parsable --dependency="afterany:$job" --chdir="$REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_methods_collect
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=logs/st_methods_collect.%j.out
#SBATCH --error=logs/st_methods_collect.%j.err
set -euo pipefail
"$SCT_METHOD_ROOT/venv/bin/python" validation/benchmark/benchmark_global_methods.py --collect --output "$SCT_METHOD_ROOT"
SLURM
)
printf 'Environment/preparation job: %s\nAll-method array (108 tasks, no throttle): %s\nCollector: %s\nResults: %s\n' "$setup" "$job" "${collect%%;*}" "$SCT_METHOD_ROOT"
