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
export SCT_DIAG_ROOT=$(mktemp -d "$REPO/local_results/pipeline_diagnostic_XXXXXX")
git rev-parse HEAD > "$SCT_DIAG_ROOT/diagnostic_commit.txt"
audit=$(sbatch --parsable --chdir="$REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_pipe_audit
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:02:00
#SBATCH --output=logs/st_pipe_audit.%j.out
#SBATCH --error=logs/st_pipe_audit.%j.err
set -euo pipefail
eval "$("$SCT_MAMBA" shell hook --shell bash)"
micromamba activate syntangle_test
python validation/benchmark/diagnose_pipeline.py --audit "$SCT_REPO/local_results"
SLURM
)
audit=${audit%%;*}
array=$(sbatch --parsable --chdir="$REPO" --dependency="afterok:$audit" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_pipe_diag
#SBATCH --partition=nano
#SBATCH --array=3,6,7%3
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:03:00
#SBATCH --output=logs/st_pipe_diag.%A_%a.out
#SBATCH --error=logs/st_pipe_diag.%A_%a.err
set -euo pipefail
eval "$("$SCT_MAMBA" shell hook --shell bash)"
micromamba activate syntangle_test
export PYTHONPATH="$SCT_REPO/src${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
python validation/benchmark/diagnose_pipeline.py \
    --root "$SCT_INPUT_ROOT" --output "$SCT_DIAG_ROOT" \
    --index "$SLURM_ARRAY_TASK_ID" --seconds 90
SLURM
)
printf 'Saved audit job: %s\nDiagnostic array: %s\nResults: %s\n' "$audit" "${array%%;*}" "$SCT_DIAG_ROOT"
