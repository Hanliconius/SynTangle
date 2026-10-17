#!/usr/bin/env bash
set -euo pipefail
SCT_REPO=$(pwd)
SCT_SOURCE=$(realpath "${1:-local_results/annotated_lep_rjv857}")
test -f "$SCT_SOURCE/figures/annotated_lep_reanalysis/COMPLETE"
test -x "$SCT_SOURCE/tools/bin/python"
mkdir -p logs local_results
SCT_RUN=$(mktemp -d "$SCT_REPO/local_results/block_filter_XXXXXX")
mkdir -p "$SCT_RUN/code"
cp -r src validation "$SCT_RUN/code/"
export SCT_SOURCE SCT_RUN
SCT_JOB=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=12G --time=00:15:00 --array=0-4 --job-name=st_filter --output="$SCT_REPO/logs/st_filter.%A_%a.out" --error="$SCT_REPO/logs/st_filter.%A_%a.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$SCT_RUN/code/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
"$SCT_SOURCE/tools/bin/python" "$SCT_RUN/code/validation/real_data/benchmark_block_filter.py" --source "$SCT_SOURCE" --output "$SCT_RUN" --index "$SLURM_ARRAY_TASK_ID"
SLURM
)
printf 'Filtering array (all 5 concurrent): %s\nResults: %s\n' "$SCT_JOB" "$SCT_RUN"
printf '%s\n' "$SCT_RUN" > local_results/latest_block_filter_run.txt
