#!/usr/bin/env bash
# Rerun the unresolved Fig. 4f obligation with a larger tree budget.
# Previous certificates remain untouched; this is a fresh search, not a resume.
set -euo pipefail
SCT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SCT_PIPELINE=/gpfs/automountdir/gpfs/scratch/martinlab/joe/SynTangle_pipeline_test
SCT_PRIOR=$(realpath "${1:-$SCT_PIPELINE/local_results/independent_order_proof_oSHLP8}")
test -f "$SCT_PRIOR/task_1/COMPLETE"
test -x "$SCT_PRIOR/tools/bin/python"
SCT_RUN=$(mktemp -d "$SCT_PIPELINE/local_results/order_extended_XXXXXX")
mkdir -p "$SCT_RUN/code" "$SCT_RUN/task_1"
cp "$SCT_ROOT"/certification/{proof.py,structural.py,order_certificate.py} "$SCT_RUN/code/"
cp "$SCT_PRIOR/task_1/"{fixture,layout}.json "$SCT_RUN/task_1/"
SCT_PYTHON=$(type -P python)
export SCT_RUN SCT_PRIOR SCT_PYTHON
printf '%s\n' "$SCT_RUN" > "$SCT_PIPELINE/local_results/latest_extended_order_proof_run.txt"
SCT_JOB=$(sbatch --parsable --partition=nano --cpus-per-task=1 --mem=8G --time=00:30:00 --job-name=st_order_more --output="$SCT_PIPELINE/logs/st_order_more.%j.out" --error="$SCT_PIPELINE/logs/st_order_more.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
SCT_OUT="$SCT_RUN/task_1"
"$SCT_PRIOR/tools/bin/python" "$SCT_RUN/code/order_certificate.py" check "$SCT_OUT/fixture.json" "$SCT_OUT/layout.json" "$SCT_OUT/certificate.json" --seconds 600 --max-nodes 1001
"$SCT_PYTHON" "$SCT_RUN/code/order_certificate.py" verify "$SCT_OUT/fixture.json" "$SCT_OUT/layout.json" "$SCT_OUT/certificate.json"
printf 'PASS\n' > "$SCT_OUT/COMPLETE"
SLURM
)
printf 'Extended Fig. 4f proof (nano, 1001 nodes, ten-minute search): %s\nResults: %s\n' "$SCT_JOB" "$SCT_RUN"
