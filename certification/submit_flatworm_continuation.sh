#!/usr/bin/env bash
set -euo pipefail
SCT_FC_CODE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SCT_FC_PIPE=/gpfs/automountdir/gpfs/scratch/martinlab/joe/SynTangle_pipeline_test
SCT_FC_PRIOR=$(realpath "${1:-$SCT_FC_PIPE/local_results/proof_combined_BpUBYJ}")
SCT_FC_PY="$SCT_FC_PIPE/local_results/independent_order_proof_oSHLP8/tools/bin/python"
test -x "$SCT_FC_PY"
test -f "$SCT_FC_PRIOR/COMPLETE"
test -f "$SCT_FC_PRIOR/task_1/selective_frontier_501.json"
SCT_FC_RUN=$(mktemp -d "$SCT_FC_PIPE/local_results/proof_flatworm_continue_XXXXXX")
mkdir -p "$SCT_FC_RUN/code" "$SCT_FC_PIPE/logs"
cp "$SCT_FC_CODE"/{proof.py,structural.py,order_certificate.py,order_strengthened.py,order_selective.py,order_resume.py,test_order_resume.py,test_proof.py} "$SCT_FC_RUN/code/"
cp "$SCT_FC_PRIOR/task_1/"{fixture,layout}.json "$SCT_FC_RUN/"
cp "$SCT_FC_PRIOR/task_1/selective_frontier_501.json" "$SCT_FC_RUN/prior.json"
export SCT_FC_RUN SCT_FC_PY
SCT_FC_JOB=$(sbatch --parsable --partition=nano --nodes=1 --cpus-per-task=1 --mem=8G --time=00:25:00 --job-name=st_flat_continue --output="$SCT_FC_PIPE/logs/st_flat_continue.%j.out" --error="$SCT_FC_PIPE/logs/st_flat_continue.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
cd "$SCT_FC_RUN/code"
"$SCT_FC_PY" -m unittest test_order_resume > "$SCT_FC_RUN/tests.log" 2>&1
"$SCT_FC_PY" order_resume.py "$SCT_FC_RUN/fixture.json" "$SCT_FC_RUN/layout.json" "$SCT_FC_RUN/prior.json" "$SCT_FC_RUN/continued.json" --seconds 900 --max-nodes 10001
printf 'PASS continued certificate replay\n' > "$SCT_FC_RUN/COMPLETE"
SLURM
)
printf 'Single-core flatworm continuation: %s\nResults: %s\n' "$SCT_FC_JOB" "$SCT_FC_RUN"
