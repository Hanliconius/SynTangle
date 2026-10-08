#!/usr/bin/env bash
set -euo pipefail
SCT_REPO=$(pwd)
SCT_SOURCE=$(realpath "${1:-local_results/annotated_lep_rjv857}")
SCT_FILTER=$(realpath "${2:-$(cat local_results/latest_block_filter_run.txt)}")
test -x "$SCT_SOURCE/tools/bin/python"
test -d "$SCT_FILTER/code/src"
mkdir -p "$SCT_FILTER/plot_code/validation/real_data" "$SCT_FILTER/plot_code/validation/benchmark" logs
cp validation/real_data/{plot_block_filter.py,benchmark_block_filter.py,plot_comparison.py} "$SCT_FILTER/plot_code/validation/real_data/"
cp validation/benchmark/{compare_genespace.py,visual_summary.py} "$SCT_FILTER/plot_code/validation/benchmark/"
export SCT_REPO SCT_SOURCE SCT_FILTER
sbatch --partition=nano --cpus-per-task=1 --mem=4G --time=00:10:00 --job-name=st_filter_fig --output="$SCT_REPO/logs/st_filter_fig.%j.out" --error="$SCT_REPO/logs/st_filter_fig.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$SCT_FILTER/code/src"
export MPLBACKEND=Agg
"$SCT_SOURCE/tools/bin/python" "$SCT_FILTER/plot_code/validation/real_data/plot_block_filter.py" --source "$SCT_SOURCE" --run "$SCT_FILTER" --latest-archive "$SCT_REPO/local_results/SynTangle_filter_comparisons.zip"
SLURM
printf 'Archive when complete: %s/local_results/SynTangle_filter_comparisons.zip\n' "$SCT_REPO"
