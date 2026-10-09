#!/usr/bin/env bash
set -euo pipefail
SCT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SCT_PIPELINE=/gpfs/automountdir/gpfs/scratch/martinlab/joe/SynTangle_pipeline_test
SCT_PRIOR=$(realpath "${1:-$SCT_PIPELINE/local_results/independent_real_proof_CtcV4i}")
for SCT_TASK in 0 5; do test -f "$SCT_PRIOR/task_$SCT_TASK/COMPLETE"; done
mkdir -p "$SCT_PIPELINE/logs"
SCT_RUN=$(mktemp -d "$SCT_PIPELINE/local_results/independent_order_proof_XXXXXX")
mkdir -p "$SCT_RUN/code"
cp "$SCT_ROOT"/certification/{proof.py,structural.py,order_certificate.py,test_order_certificate.py,test_proof.py} "$SCT_RUN/code/"
SCT_PYTHON=$(type -P python)
export SCT_RUN SCT_PRIOR SCT_PYTHON
printf '%s\n' "$SCT_RUN" > "$SCT_PIPELINE/local_results/latest_independent_order_proof_run.txt"
SCT_SETUP=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=4G --time=00:20:00 --job-name=st_order_setup --output="$SCT_PIPELINE/logs/st_order_setup.%j.out" --error="$SCT_PIPELINE/logs/st_order_setup.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
"$SCT_PYTHON" -m venv "$SCT_RUN/tools"
"$SCT_RUN/tools/bin/python" -m pip install 'scipy==1.17.0' 'numpy==2.3.5' > "$SCT_RUN/dependency_setup.log" 2>&1
cd "$SCT_RUN/code"
"$SCT_RUN/tools/bin/python" -m unittest test_order_certificate test_proof > "$SCT_RUN/checker_tests.log" 2>&1
printf 'PASS checker tests; dependencies ready\n'
SLURM
)
SCT_ARRAY=$(sbatch --parsable --partition=cpu --array=0-1 --dependency="afterok:$SCT_SETUP" --cpus-per-task=1 --mem=8G --time=00:20:00 --job-name=st_order_proof --output="$SCT_PIPELINE/logs/st_order_proof.%A_%a.out" --error="$SCT_PIPELINE/logs/st_order_proof.%A_%a.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
SCT_SOURCE_TASK=$((SLURM_ARRAY_TASK_ID * 5))
SCT_OUT="$SCT_RUN/task_$SLURM_ARRAY_TASK_ID"
mkdir -p "$SCT_OUT"
cp "$SCT_PRIOR/task_$SCT_SOURCE_TASK/"{fixture,layout}.json "$SCT_OUT/"
"$SCT_RUN/tools/bin/python" "$SCT_RUN/code/order_certificate.py" check "$SCT_OUT/fixture.json" "$SCT_OUT/layout.json" "$SCT_OUT/certificate.json" --seconds 600 --max-nodes 101
# Replay in a fresh process with the original Python: no SciPy solver imports.
"$SCT_PYTHON" "$SCT_RUN/code/order_certificate.py" verify "$SCT_OUT/fixture.json" "$SCT_OUT/layout.json" "$SCT_OUT/certificate.json"
printf 'PASS\n' > "$SCT_OUT/COMPLETE"
SLURM
)
printf 'Checker setup: %s\nFull-order audit (two concurrent CPU tasks): %s\nResults: %s\n' "$SCT_SETUP" "$SCT_ARRAY" "$SCT_RUN"
