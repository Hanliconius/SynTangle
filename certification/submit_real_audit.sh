#!/usr/bin/env bash
set -euo pipefail
SCT_PROOF_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SCT_PIPELINE=$(realpath "${1:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/SynTangle_pipeline_test}")
SCT_FILTER=$(realpath "${2:-$SCT_PIPELINE/local_results/planarian_filter_9ol27C}")
test -f "$SCT_FILTER/planarians/task_0/report.json"
test -f "$SCT_FILTER/flatworms_fig4f/task_0/report.json"
mkdir -p "$SCT_PIPELINE/logs" "$SCT_PIPELINE/local_results"
SCT_PROOF_RUN=$(mktemp -d "$SCT_PIPELINE/local_results/independent_real_proof_XXXXXX")
mkdir -p "$SCT_PROOF_RUN/code"
cp "$SCT_PROOF_ROOT"/certification/{proof.py,structural.py,audit_real.py} "$SCT_PROOF_RUN/code/"
SCT_PROOF_PYTHON=$(type -P python)
export SCT_PIPELINE SCT_FILTER SCT_PROOF_RUN SCT_PROOF_PYTHON
printf '%s\n' "$SCT_PROOF_RUN" > "$SCT_PIPELINE/local_results/latest_independent_real_proof_run.txt"
SCT_JOBS=$(sbatch --parsable --partition=cpu --array=0-9 --cpus-per-task=1 --mem=4G --time=00:30:00 --job-name=st_real_proof --output="$SCT_PIPELINE/logs/st_real_proof.%A_%a.out" --error="$SCT_PIPELINE/logs/st_real_proof.%A_%a.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
"$SCT_PROOF_PYTHON" "$SCT_PROOF_RUN/code/audit_real.py" --pipeline "$SCT_PIPELINE" --filter-run "$SCT_FILTER" --output "$SCT_PROOF_RUN" --task "$SLURM_ARRAY_TASK_ID"
SLURM
)
printf 'Independent proof audit (ten concurrent tasks): %s\nResults: %s\n' "$SCT_JOBS" "$SCT_PROOF_RUN"
