#!/usr/bin/env bash
set -euo pipefail
SCT_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_REPO"
mkdir -p logs local_results
SCT_RUN=$(mktemp -d "$SCT_REPO/local_results/rhynchospora_sources_XXXXXX")
mkdir -p "$SCT_RUN/code"
cp validation/real_data/{capture_rhynchospora_fig1a.py,prepare_published_inputs.py} "$SCT_RUN/code/"
SCT_PYTHON=$(type -P python)
export SCT_REPO SCT_RUN SCT_PYTHON
printf '%s\n' "$SCT_RUN" > local_results/latest_rhynchospora_sources_run.txt
SCT_JOB=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=4G --time=00:30:00 --job-name=st_rhync_src --output="$SCT_REPO/logs/st_rhync_source.%j.out" --error="$SCT_REPO/logs/st_rhync_source.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
"$SCT_PYTHON" -m venv "$SCT_RUN/tools"
"$SCT_RUN/tools/bin/pip" install PyMuPDF==1.26.5
"$SCT_RUN/tools/bin/python" "$SCT_RUN/code/capture_rhynchospora_fig1a.py" --run "$SCT_RUN"
cp "$SCT_RUN/SynTangle_rhynchospora_sources.zip" "$SCT_REPO/local_results/SynTangle_rhynchospora_sources.zip"
SLURM
)
printf 'Source capture: %s\nResults: %s\nArchive: %s/local_results/SynTangle_rhynchospora_sources.zip\n' "$SCT_JOB" "$SCT_RUN" "$SCT_REPO"
