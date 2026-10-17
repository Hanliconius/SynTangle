#!/bin/bash
set -euo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$REPO"
export SCT_MAMBA=${MAMBA_EXE:-$(type -P micromamba)}
test -x "$SCT_MAMBA"
export SCT_INPUT_ROOT="${SCT_INPUT_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/SynTangle_pegasus_test/local_results/genespace_stress}"
export SCT_HISTORY="$REPO/local_results"
export SCT_WARM_RUN="${SCT_WARM_RUN:-}"
test -f "$SCT_INPUT_ROOT/benchmark_manifest.tsv"
mkdir -p logs local_results
export SCT_EFF_ROOT=$(mktemp -d "$REPO/local_results/efficiency_XXXXXX")
export SCT_EFF_CODE="$SCT_EFF_ROOT/code"
mkdir -p "$SCT_EFF_CODE/validation/benchmark" "$SCT_EFF_ROOT/audits"
cp -a src "$SCT_EFF_CODE/"
cp validation/benchmark/*.py validation/benchmark/hybrid_methods_requirements.txt "$SCT_EFF_CODE/validation/benchmark/"
git rev-parse HEAD > "$SCT_EFF_ROOT/commit.txt"
setup=$(sbatch --parsable --chdir="$REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_eff_setup
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=logs/st_eff_setup.%j.out
#SBATCH --error=logs/st_eff_setup.%j.err
set -euo pipefail
eval "$("$SCT_MAMBA" shell hook --shell bash)"
micromamba activate syntangle_test
python -m venv "$SCT_EFF_ROOT/venv"
"$SCT_EFF_ROOT/venv/bin/python" -m pip install -r "$SCT_EFF_CODE/validation/benchmark/hybrid_methods_requirements.txt"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export PYTHONPATH="$SCT_EFF_CODE/src"
args=()
if [[ -n "$SCT_WARM_RUN" ]]; then args+=(--warm-run "$SCT_WARM_RUN"); fi
"$SCT_EFF_ROOT/venv/bin/python" "$SCT_EFF_CODE/validation/benchmark/benchmark_efficiency.py" \
  --prepare --history "$SCT_HISTORY" --root "$SCT_INPUT_ROOT" --output "$SCT_EFF_ROOT" "${args[@]}"
SLURM
)
setup=${setup%%;*}
audit=$(sbatch --parsable --dependency="afterok:$setup" --chdir="$REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_eff_audit
#SBATCH --partition=cpu
#SBATCH --array=0-8
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=logs/st_eff_audit.%A_%a.out
#SBATCH --error=logs/st_eff_audit.%A_%a.err
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export PYTHONPATH="$SCT_EFF_CODE/src"
"$SCT_EFF_ROOT/venv/bin/python" "$SCT_EFF_CODE/validation/benchmark/audit_global_results.py" \
  --history "$SCT_HISTORY" --root "$SCT_INPUT_ROOT" --seconds 120 --case-index "$SLURM_ARRAY_TASK_ID" \
  --output "$SCT_EFF_ROOT/audits/$SLURM_ARRAY_TASK_ID.json"
SLURM
)
audit=${audit%%;*}
job=$(sbatch --parsable --dependency="afterok:$setup" --chdir="$REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_eff
#SBATCH --partition=cpu
#SBATCH --array=0-251
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:12:00
#SBATCH --output=logs/st_eff.%A_%a.out
#SBATCH --error=logs/st_eff.%A_%a.err
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export PYTHONPATH="$SCT_EFF_CODE/src"
"$SCT_EFF_ROOT/venv/bin/python" "$SCT_EFF_CODE/validation/benchmark/benchmark_efficiency.py" \
  --root "$SCT_INPUT_ROOT" --output "$SCT_EFF_ROOT" --index "$SLURM_ARRAY_TASK_ID"
SLURM
)
job=${job%%;*}
collect=$(sbatch --parsable --dependency="afterany:$job:$audit" --chdir="$REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_eff_collect
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=logs/st_eff_collect.%j.out
#SBATCH --error=logs/st_eff_collect.%j.err
set -euo pipefail
"$SCT_EFF_ROOT/venv/bin/python" "$SCT_EFF_CODE/validation/benchmark/benchmark_efficiency.py" --collect --output "$SCT_EFF_ROOT"
SLURM
)
printf 'Setup: %s\nProof audit array (9 tasks): %s\nEfficiency array (252 tasks; no throttle): %s\nCollector: %s\nResults: %s\n' "$setup" "$audit" "$job" "${collect%%;*}" "$SCT_EFF_ROOT"
printf '{"setup":"%s","audit":"%s","array":"%s","collector":"%s"}\n' "$setup" "$audit" "$job" "${collect%%;*}" > "$SCT_EFF_ROOT/jobs.json"
