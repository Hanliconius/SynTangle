#!/bin/bash
# Retry only the unfinished case on existing biology; preserve nine results.
set -euo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$REPO"
export SCT_MAMBA=${MAMBA_EXE:-$(type -P micromamba)}
test -x "$SCT_MAMBA"
export SCT_REPO="$REPO"
export SCT_ROOT="$REPO/local_results/genespace_pilot"
test -f "$SCT_ROOT/benchmark_manifest.tsv"
mkdir -p logs

retry=$(sbatch --parsable --chdir="$REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=gs_retry9
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/gs_retry9.%j.out
#SBATCH --error=logs/gs_retry9.%j.err
set -euo pipefail
export PYTHONPATH="$SCT_REPO/src${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
case_path="$SCT_ROOT/comparison/gs_8sp_31chr_0ev"
if [ -d "$case_path" ]; then
    cp -a "$case_path" "$case_path.before_retry_$SLURM_JOB_ID"
fi
eval "$("$SCT_MAMBA" shell hook --shell bash)"
micromamba activate syntangle_test
timeout --signal=TERM --kill-after=30s 1500s python \
    validation/benchmark/compare_genespace.py "$SCT_ROOT" \
    --index 9 --micromamba "$SCT_MAMBA"
SLURM
)
retry=${retry%%;*}

collect=$(sbatch --parsable --chdir="$REPO" --dependency="afterany:$retry" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=gs_collect_retry
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:05:00
#SBATCH --output=logs/gs_collect_retry.%j.out
#SBATCH --error=logs/gs_collect_retry.%j.err
set -euo pipefail
export PYTHONPATH="$SCT_REPO/src${PYTHONPATH:+:$PYTHONPATH}"
eval "$("$SCT_MAMBA" shell hook --shell bash)"
micromamba activate syntangle_test
python validation/benchmark/compare_genespace.py "$SCT_ROOT" --collect
SLURM
)
collect=${collect%%;*}
printf 'Case 9 retry: %s\nCollection: %s\n' "$retry" "$collect"
