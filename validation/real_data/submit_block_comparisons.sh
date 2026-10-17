#!/bin/bash
set -euo pipefail
SCT_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_REPO"
export SCT_FINAL_ROOT=${1:?Pass completed final validation run directory}
test -x "$SCT_FINAL_ROOT/venv/bin/python"
export SCT_FINAL_CODE="$SCT_FINAL_ROOT/code"
export SCT_MAMBA=${MAMBA_EXE:-$(type -P micromamba)}
export SCT_GENESPACE_ENV=${SCT_GENESPACE_ENV:-lep_busco_painter_clean}
job=$(sbatch --parsable --chdir="$SCT_REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_real_compare
#SBATCH --partition=cpu
#SBATCH --array=0-11
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:25:00
#SBATCH --output=logs/st_real_compare.%A_%a.out
#SBATCH --error=logs/st_real_compare.%A_%a.err
set -euo pipefail
export PYTHONPATH="$SCT_FINAL_CODE/src" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
index=$((SLURM_ARRAY_TASK_ID / 2))
seconds=150
if ((SLURM_ARRAY_TASK_ID % 2)); then seconds=600; fi
"$SCT_FINAL_ROOT/venv/bin/python" "$SCT_FINAL_CODE/validation/real_data/benchmark_blocks.py" --root "$SCT_FINAL_ROOT/real_data" --manifest "$SCT_FINAL_CODE/validation/real_data/datasets.json" --index "$index" --seconds "$seconds"
SLURM
)
job=${job%%;*}
collector=$(sbatch --parsable --dependency="afterany:$job" --chdir="$SCT_REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_real_collect
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=logs/st_real_collect.%j.out
#SBATCH --error=logs/st_real_collect.%j.err
set -euo pipefail
"$SCT_FINAL_ROOT/venv/bin/python" "$SCT_FINAL_CODE/validation/real_data/collect_block_comparisons.py" --root "$SCT_FINAL_ROOT/real_data" --manifest "$SCT_FINAL_CODE/validation/real_data/datasets.json"
SLURM
)
printf 'Collector: %s\n' "${collector%%;*}"
printf 'Block comparison array (12 tasks, no throttle): %s\nResults: %s/real_data\n' "${job%%;*}" "$SCT_FINAL_ROOT"
