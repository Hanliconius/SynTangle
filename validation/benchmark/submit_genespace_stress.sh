#!/bin/bash
# Login-node work is limited to paths and Slurm submissions. All R/Python work
# (including generation and collection) runs on compute nodes.
set -euo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$REPO"
export SCT_MAMBA=${MAMBA_EXE:-$(type -P micromamba)}
test -x "$SCT_MAMBA"
export SCT_REPO="$REPO"
export SCT_ROOT="${SCT_STRESS_ROOT:-$REPO/local_results/genespace_stress}"
mkdir -p logs
if [[ -e "$SCT_ROOT/benchmark_manifest.tsv" ]]; then
    printf "Existing stress results at %s; refusing to overwrite. Use SCT_STRESS_ROOT for a new run.\n" "$SCT_ROOT" >&2
    exit 1
fi

prepare=$(sbatch --parsable --chdir="$REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=gs_stress_prepare
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/gs_stress_prepare.%j.out
#SBATCH --error=logs/gs_stress_prepare.%j.err
set -euo pipefail
"$SCT_MAMBA" run -n syntangle_test Rscript \
    validation/benchmark/generate_benchmark_cases.R "$SCT_ROOT" stress
SLURM
)
prepare=${prepare%%;*}

array=$(sbatch --parsable --chdir="$REPO" --dependency="afterok:$prepare" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=gs_stress
#SBATCH --partition=cpu
#SBATCH --array=0-8%3
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/gs_stress.%A_%a.out
#SBATCH --error=logs/gs_stress.%A_%a.err
set -euo pipefail
export PYTHONPATH="$SCT_REPO/src${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
# Keep baselines on disk if the exact solver takes too long. Each case gets
# its own process, outputs and 235-minute external wall-clock limit.
timeout --signal=TERM --kill-after=30s 14100s \
    bash -c 'set -euo pipefail; eval "$("$SCT_MAMBA" shell hook --shell bash)"; micromamba activate syntangle_test; exec python "$@"' gs-python \
    validation/benchmark/compare_genespace.py "$SCT_ROOT" \
    --index "$SLURM_ARRAY_TASK_ID" --micromamba "$SCT_MAMBA" \
    --transition-cap 250000 --branch-node-cap 100000 --local-restarts 4
SLURM
)
array=${array%%;*}

collect=$(sbatch --parsable --chdir="$REPO" --dependency="afterany:$array" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=gs_stress_collect
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:05:00
#SBATCH --output=logs/gs_stress_collect.%j.out
#SBATCH --error=logs/gs_stress_collect.%j.err
set -euo pipefail
export PYTHONPATH="$SCT_REPO/src${PYTHONPATH:+:$PYTHONPATH}"
"$SCT_MAMBA" run -n syntangle_test python \
    validation/benchmark/compare_genespace.py "$SCT_ROOT" --collect
SLURM
)
collect=${collect%%;*}
printf 'Preparation job: %s\nComparison array: %s (9 cases, at most 3 simultaneous)\nCollection job: %s\n' \
    "$prepare" "$array" "$collect"
