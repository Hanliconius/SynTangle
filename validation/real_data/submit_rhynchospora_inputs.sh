#!/usr/bin/env bash
set -euo pipefail
SCT_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_REPO"
mkdir -p logs local_results
SCT_RUN=${1:-$(mktemp -d "$SCT_REPO/local_results/rhynchospora_inputs_XXXXXX")}
mkdir -p "$SCT_RUN/code"
SCT_RUN=$(cd "$SCT_RUN" && pwd)
cp validation/real_data/{prepare_rhynchospora_inputs.py,rhynchospora_input_manifest.json,capture_rhynchospora_fig1a.py,prepare_published_inputs.py} "$SCT_RUN/code/"
SCT_PYTHON=$(type -P python)
"$SCT_PYTHON" -c 'import sys; assert sys.version_info >= (3, 10)'
export SCT_REPO SCT_RUN SCT_PYTHON
printf '%s\n' "$SCT_RUN" > local_results/latest_rhynchospora_inputs_run.txt
SCT_ARRAY=$(sbatch --parsable --partition=cpu --array=0-4 --cpus-per-task=1 --mem=4G --time=04:00:00 --job-name=st_rhync_in --output="$SCT_REPO/logs/st_rhync_input.%A_%a.out" --error="$SCT_REPO/logs/st_rhync_input.%A_%a.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
"$SCT_PYTHON" "$SCT_RUN/code/prepare_rhynchospora_inputs.py" --run "$SCT_RUN" --manifest "$SCT_RUN/code/rhynchospora_input_manifest.json" --task "$SLURM_ARRAY_TASK_ID"
SLURM
)
SCT_ARRAY=${SCT_ARRAY%%;*}
SCT_COLLECT=$(sbatch --parsable --partition=cpu --dependency="afterany:$SCT_ARRAY" --cpus-per-task=1 --mem=2G --time=00:15:00 --job-name=st_rhync_chk --output="$SCT_REPO/logs/st_rhync_check.%j.out" --error="$SCT_REPO/logs/st_rhync_check.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
"$SCT_PYTHON" "$SCT_RUN/code/prepare_rhynchospora_inputs.py" --run "$SCT_RUN" --manifest "$SCT_RUN/code/rhynchospora_input_manifest.json"
SLURM
)
printf 'Downloads (five concurrent tasks): %s\nInput audit: %s\nResults: %s\nRead when finished: cat logs/st_rhync_check.%s.out\n' "$SCT_ARRAY" "$SCT_COLLECT" "$SCT_RUN" "${SCT_COLLECT%%;*}"
