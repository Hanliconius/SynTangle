#!/usr/bin/env bash
set -euo pipefail
SCT_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_REPO"
SCT_SOURCE=$(realpath "${1:-local_results/flatworms_fig4f_Au9tel}")
SCT_PLAN=$(realpath "${2:-local_results/planarian_visual_e5OOSr}")
SCT_PYTHON="$SCT_PLAN/venv/bin/python"
test -x "$SCT_PYTHON"
test -f "$SCT_SOURCE/figures/flatworms_fig4f/COMPLETE"
SCT_RUN=$(mktemp -d "$SCT_REPO/local_results/fig4f_reconstruction_XXXXXX")
mkdir -p "$SCT_RUN/code/validation/real_data" logs
cp -a src "$SCT_RUN/code/"
cp validation/real_data/{recapitulate_flatworm_fig4f.py,prepare_published_inputs.py} "$SCT_RUN/code/validation/real_data/"
export SCT_REPO SCT_SOURCE SCT_PYTHON SCT_RUN
printf '%s\n' "$SCT_RUN" > local_results/latest_fig4f_reconstruction_run.txt
sbatch --partition=cpu --cpus-per-task=1 --mem=4G --time=00:20:00 --job-name=st_fig4f_rep --output="$SCT_REPO/logs/st_fig4f_rep.%j.out" --error="$SCT_REPO/logs/st_fig4f_rep.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$SCT_RUN/code/src" MPLBACKEND=Agg
"$SCT_PYTHON" "$SCT_RUN/code/validation/real_data/recapitulate_flatworm_fig4f.py" --source "$SCT_SOURCE" --output "$SCT_RUN/figures"
cp "$SCT_RUN/SynTangle_fig4f_reconstruction.zip" "$SCT_REPO/local_results/SynTangle_fig4f_reconstruction.partial.zip"
mv "$SCT_REPO/local_results/SynTangle_fig4f_reconstruction.partial.zip" "$SCT_REPO/local_results/SynTangle_fig4f_reconstruction.zip"
SLURM
printf 'Results: %s\nArchive when complete: %s/local_results/SynTangle_fig4f_reconstruction.zip\n' "$SCT_RUN" "$SCT_REPO"
