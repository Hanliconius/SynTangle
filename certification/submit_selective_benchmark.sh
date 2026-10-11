#!/usr/bin/env bash
# One compute node, one CPU core; all proof variants run sequentially.
set -euo pipefail
SCT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SCT_PIPELINE=/gpfs/automountdir/gpfs/scratch/martinlab/joe/SynTangle_pipeline_test
SCT_PRIOR=$(realpath "${1:-$SCT_PIPELINE/local_results/independent_order_proof_oSHLP8}")
test -x "$SCT_PRIOR/tools/bin/python"
for SCT_TASK in 0 1; do test -f "$SCT_PRIOR/task_$SCT_TASK/COMPLETE"; done
SCT_RUN=$(mktemp -d "$SCT_PIPELINE/local_results/proof_selective_XXXXXX")
mkdir -p "$SCT_RUN/code" "$SCT_PIPELINE/logs"
cp "$SCT_ROOT"/certification/{proof.py,structural.py,order_certificate.py,order_strengthened.py,selective_node_benchmark.py,order_selective.py,test_order_selective.py,test_order_strengthened.py,test_order_certificate.py,test_proof.py} "$SCT_RUN/code/"
for SCT_TASK in 0 1; do
    mkdir -p "$SCT_RUN/task_$SCT_TASK"
    cp "$SCT_PRIOR/task_$SCT_TASK/"{fixture,layout}.json "$SCT_RUN/task_$SCT_TASK/"
done
printf '%s\n' "$SCT_RUN" > "$SCT_PIPELINE/local_results/latest_proof_selective_run.txt"
export SCT_RUN SCT_PRIOR
SCT_JOB=$(sbatch --parsable --partition=nano --nodes=1 --cpus-per-task=1 --mem=8G --time=00:25:00 --job-name=st_proof_selective --output="$SCT_PIPELINE/logs/st_proof_selective.%j.out" --error="$SCT_PIPELINE/logs/st_proof_selective.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
cd "$SCT_RUN/code"
"$SCT_PRIOR/tools/bin/python" -m unittest test_order_selective test_order_strengthened test_order_certificate > "$SCT_RUN/tests.log" 2>&1
"$SCT_PRIOR/tools/bin/python" selective_node_benchmark.py --run "$SCT_RUN" --seconds 60
SLURM
)
printf 'Single-core proof benchmark: %s\nResults: %s\nRead: tail -n 16 %s/logs/st_proof_selective.%s.out\n' "$SCT_JOB" "$SCT_RUN" "$SCT_PIPELINE" "${SCT_JOB%%;*}"
