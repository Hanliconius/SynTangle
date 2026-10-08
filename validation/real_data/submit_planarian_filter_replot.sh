#!/usr/bin/env bash
set -euo pipefail
SCT_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_REPO"
SCT_RUN=$(realpath "${1:-$(cat local_results/latest_planarian_filter_run.txt)}")
SCT_PLAN=$(realpath "${2:-local_results/planarian_visual_e5OOSr}")
SCT_FLAT=$(realpath "${3:-local_results/flatworms_fig4f_Au9tel}")
SCT_PYTHON="$SCT_PLAN/venv/bin/python"
test -x "$SCT_PYTHON"
test -d "$SCT_RUN/code/src"
mkdir -p "$SCT_RUN/replot_code/validation/real_data" "$SCT_RUN/replot_code/validation/benchmark" logs
cp validation/real_data/{plot_block_filter.py,benchmark_block_filter.py,plot_comparison.py} "$SCT_RUN/replot_code/validation/real_data/"
cp validation/benchmark/{visual_summary.py,compare_genespace.py} "$SCT_RUN/replot_code/validation/benchmark/"
export SCT_REPO SCT_RUN SCT_PLAN SCT_FLAT SCT_PYTHON
sbatch --partition=cpu --cpus-per-task=1 --mem=4G --time=00:15:00 --job-name=st_plan_fplot --output="$SCT_REPO/logs/st_plan_fplot.%j.out" --error="$SCT_REPO/logs/st_plan_fplot.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$SCT_RUN/code/src" MPLBACKEND=Agg
for SCT_CASE in planarians flatworms_fig4f; do
    if [[ "$SCT_CASE" == planarians ]]; then SCT_SOURCE="$SCT_PLAN"; else SCT_SOURCE="$SCT_FLAT"; fi
    "$SCT_PYTHON" "$SCT_RUN/replot_code/validation/real_data/plot_block_filter.py" --case "$SCT_CASE" --source "$SCT_SOURCE" --run "$SCT_RUN/$SCT_CASE" --latest-archive "$SCT_RUN/$SCT_CASE.zip"
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

