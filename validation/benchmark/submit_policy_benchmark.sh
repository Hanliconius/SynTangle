#!/bin/bash
set -euo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$REPO"
export SCT_MAMBA=${MAMBA_EXE:-$(type -P micromamba)}
test -x "$SCT_MAMBA"
export SCT_INPUT_ROOT="${SCT_INPUT_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/SynTangle_pegasus_test/local_results/genespace_stress}"
export SCT_HISTORY="$REPO/local_results"
export SCT_REFERENCE_RUN="${SCT_REFERENCE_RUN:-}"
test -f "$SCT_INPUT_ROOT/benchmark_manifest.tsv"
mkdir -p logs local_results
export SCT_POLICY_ROOT=$(mktemp -d "$REPO/local_results/policy_XXXXXX")
export SCT_POLICY_CODE="$SCT_POLICY_ROOT/code"
mkdir -p "$SCT_POLICY_CODE/validation/benchmark" "$SCT_POLICY_ROOT/audits"
cp -a src "$SCT_POLICY_CODE/"
cp validation/benchmark/*.py validation/benchmark/hybrid_methods_requirements.txt "$SCT_POLICY_CODE/validation/benchmark/"
git rev-parse HEAD > "$SCT_POLICY_ROOT/commit.txt"
setup=$(sbatch --parsable --chdir="$REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_policy_setup
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=logs/st_policy_setup.%j.out
#SBATCH --error=logs/st_policy_setup.%j.err
set -euo pipefail
eval "$("$SCT_MAMBA" shell hook --shell bash)"
micromamba activate syntangle_test
python -m venv "$SCT_POLICY_ROOT/venv"
"$SCT_POLICY_ROOT/venv/bin/python" -m pip install -r "$SCT_POLICY_CODE/validation/benchmark/hybrid_methods_requirements.txt"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export PYTHONPATH="$SCT_POLICY_CODE/src"
args=()
if [[ -n "$SCT_REFERENCE_RUN" ]]; then args+=(--reference-run "$SCT_REFERENCE_RUN"); fi
"$SCT_POLICY_ROOT/venv/bin/python" "$SCT_POLICY_CODE/validation/benchmark/benchmark_policy.py" \
  --prepare --history "$SCT_HISTORY" --root "$SCT_INPUT_ROOT" --output "$SCT_POLICY_ROOT" "${args[@]}"
SLURM
)
setup=${setup%%;*}
audit=$(sbatch --parsable --dependency="afterok:$setup" --chdir="$REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_policy_audit
#SBATCH --partition=cpu
#SBATCH --array=0-8
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --output=logs/st_policy_audit.%A_%a.out
#SBATCH --error=logs/st_policy_audit.%A_%a.err
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export PYTHONPATH="$SCT_POLICY_CODE/src"
"$SCT_POLICY_ROOT/venv/bin/python" "$SCT_POLICY_CODE/validation/benchmark/audit_global_results.py" \
  --history "$SCT_HISTORY" --root "$SCT_INPUT_ROOT" --seconds 600 --case-index "$SLURM_ARRAY_TASK_ID" \
  --output "$SCT_POLICY_ROOT/audits/$SLURM_ARRAY_TASK_ID.json"
SLURM
)
audit=${audit%%;*}
job=$(sbatch --parsable --dependency="afterok:$setup" --chdir="$REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_policy
#SBATCH --partition=cpu
#SBATCH --array=0-209
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:12:00
#SBATCH --output=logs/st_policy.%A_%a.out
#SBATCH --error=logs/st_policy.%A_%a.err
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export PYTHONPATH="$SCT_POLICY_CODE/src"
"$SCT_POLICY_ROOT/venv/bin/python" "$SCT_POLICY_CODE/validation/benchmark/benchmark_policy.py" \
  --root "$SCT_INPUT_ROOT" --output "$SCT_POLICY_ROOT" --index "$SLURM_ARRAY_TASK_ID"
SLURM
)
job=${job%%;*}
collect=$(sbatch --parsable --dependency="afterany:$job:$audit" --chdir="$REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_policy_collect
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=logs/st_policy_collect.%j.out
#SBATCH --error=logs/st_policy_collect.%j.err
set -euo pipefail
"$SCT_POLICY_ROOT/venv/bin/python" "$SCT_POLICY_CODE/validation/benchmark/benchmark_policy.py" --collect --output "$SCT_POLICY_ROOT"
SLURM
)
printf 'Setup: %s\nProof audit array (9 tasks): %s\nEfficiency array (210 tasks; no throttle): %s\nCollector: %s\nResults: %s\n' "$setup" "$audit" "$job" "${collect%%;*}" "$SCT_POLICY_ROOT"
printf '{"setup":"%s","audit":"%s","array":"%s","collector":"%s"}\n' "$setup" "$audit" "$job" "${collect%%;*}" > "$SCT_POLICY_ROOT/jobs.json"
