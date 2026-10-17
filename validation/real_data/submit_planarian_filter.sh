#!/usr/bin/env bash
set -euo pipefail
SCT_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_REPO"
SCT_PLAN=$(realpath "${1:-local_results/planarian_visual_e5OOSr}")
SCT_FLAT=$(realpath "${2:-local_results/flatworms_fig4f_Au9tel}")
SCT_PYTHON="$SCT_PLAN/venv/bin/python"
test -x "$SCT_PYTHON"
test -f "$SCT_PLAN/figures/planarians/COMPLETE"
test -f "$SCT_FLAT/figures/flatworms_fig4f/COMPLETE"
mkdir -p logs local_results
SCT_RUN=$(mktemp -d "$SCT_REPO/local_results/planarian_filter_XXXXXX")
mkdir -p "$SCT_RUN/code" "$SCT_RUN/planarians" "$SCT_RUN/flatworms_fig4f"
cp -a src validation "$SCT_RUN/code/"
export SCT_REPO SCT_PLAN SCT_FLAT SCT_PYTHON SCT_RUN
SCT_ARRAY=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=12G --time=00:15:00 --array=0-9 --job-name=st_plan_filter --output="$SCT_REPO/logs/st_plan_filter.%A_%a.out" --error="$SCT_REPO/logs/st_plan_filter.%A_%a.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$SCT_RUN/code/src" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if ((SLURM_ARRAY_TASK_ID<5)); then
    SCT_CASE=planarians; SCT_SOURCE="$SCT_PLAN"; SCT_INDEX=$SLURM_ARRAY_TASK_ID
else
    SCT_CASE=flatworms_fig4f; SCT_SOURCE="$SCT_FLAT"; SCT_INDEX=$((SLURM_ARRAY_TASK_ID-5))
fi
"$SCT_PYTHON" "$SCT_RUN/code/validation/real_data/benchmark_block_filter.py" --case "$SCT_CASE" --source "$SCT_SOURCE" --output "$SCT_RUN/$SCT_CASE" --index "$SCT_INDEX" --start-panel 2
SLURM
)
SCT_COLLECT=$(sbatch --parsable --dependency="afterany:$SCT_ARRAY" --partition=cpu --cpus-per-task=1 --mem=4G --time=00:15:00 --job-name=st_plan_fplot --output="$SCT_REPO/logs/st_plan_fplot.%j.out" --error="$SCT_REPO/logs/st_plan_fplot.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$SCT_RUN/code/src" MPLBACKEND=Agg
for SCT_CASE in planarians flatworms_fig4f; do
    if [[ "$SCT_CASE" == planarians ]]; then SCT_SOURCE="$SCT_PLAN"; else SCT_SOURCE="$SCT_FLAT"; fi
    "$SCT_PYTHON" "$SCT_RUN/code/validation/real_data/plot_block_filter.py" --case "$SCT_CASE" --source "$SCT_SOURCE" --run "$SCT_RUN/$SCT_CASE" --latest-archive "$SCT_RUN/$SCT_CASE.zip"
done
"$SCT_PYTHON" - "$SCT_RUN" "$SCT_REPO/local_results/SynTangle_planarian_filter_comparisons.zip" <<'PY'
import pathlib,sys,zipfile
run=pathlib.Path(sys.argv[1]);target=pathlib.Path(sys.argv[2]);tmp=target.with_suffix('.partial.zip')
with zipfile.ZipFile(tmp,'w',zipfile.ZIP_DEFLATED) as z:
    for case in ('planarians','flatworms_fig4f'):
        with zipfile.ZipFile(run/(case+'.zip')) as source:
            for item in source.infolist():z.writestr(case+'/'+item.filename,source.read(item))
tmp.replace(target)
print('ARCHIVE:',target)
PY
SLURM
)
printf '%s\n' "$SCT_RUN" > local_results/latest_planarian_filter_run.txt
printf 'Both panels, all 10 tasks concurrently: %s\nPlot/collector: %s\nResults: %s\nArchive when complete: %s/local_results/SynTangle_planarian_filter_comparisons.zip\n' "$SCT_ARRAY" "$SCT_COLLECT" "$SCT_RUN" "$SCT_REPO"
