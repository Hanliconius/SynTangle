#!/usr/bin/env bash
set -euo pipefail
source_run=${1:?Usage: submit_pieridae_components.sh EXISTING_TOPOLOGY_VIEW [SECONDS] [--audit-only|--recover|--validate-figures]}
seconds=${2:-600}
mode=${3:---solve}
[[ $mode == --solve || $mode == --audit-only || $mode == --recover || $mode == --validate-figures ]] || { echo 'Invalid mode' >&2; exit 2; }
[[ $seconds =~ ^[1-9][0-9]*$ ]] || { echo 'Seconds must be a positive integer' >&2; exit 2; }
source_run=$(cd "$source_run" && pwd)
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
for name in tables/chromosomes.tsv tables/links.tsv optimization/COMPLETE.json optimization/layout.tsv optimization/input_layout.tsv; do
  [[ -s "$source_run/$name" ]] || { echo "Missing $source_run/$name" >&2; exit 2; }
done
python=${SCT_PYTHON:-$source_run/envs/syntangle/bin/python}
[[ -x $python ]] || { echo "Missing Python: $python; set SCT_PYTHON to the existing SynTangle environment" >&2; exit 2; }
if [[ $mode == --validate-figures ]]; then
  [[ -s "$source_run/code/plotter/analyse.py" ]] || { echo 'Missing frozen plotter' >&2; exit 2; }
  repair=$(mktemp -d "$source_run/validation_repair.XXXXXX")
  cp "$here/probe_pieridae_components.py" "$repair/"
  git -C "$repo" rev-parse HEAD > "$repair/code_commit.txt"
  sha256sum "$repair/probe_pieridae_components.py" > "$repair/code_sha256.txt"
  cat > "$repair/worker.sh" <<'VALIDATE'
#!/usr/bin/env bash
set -euo pipefail
run=$1; repair=$2; python=$3
export PYTHONPATH="$run/code/src"
export PATH="$(dirname "$python"):$PATH"
"$python" "$repair/probe_pieridae_components.py" "$run" "$run" --validate-figures
VALIDATE
  job=$(sbatch --parsable --partition=nano --cpus-per-task=1 --mem=4G --time=15 --job-name=st_pdf_check --output="$repair/check.%j.out" --error="$repair/check.%j.err" "$repair/worker.sh" "$source_run" "$repair" "$python")
  printf 'Validation-only job: %s\nLogs: %s/check.%s.{out,err}\n' "$job" "$repair" "${job%%;*}"
else
out=$(mktemp -d "$source_run/component_test.XXXXXX")
mkdir -p "$out/code/src" "$out/logs"
cp -R "$repo/src/syntangle" "$out/code/src/"
cp "$here/probe_pieridae_components.py" "$out/code/"
render=no
plotter="$source_run/code"
if [[ -s "$source_run/code/plotter/analyse.py" ]]; then plotter="$source_run/code/plotter"; fi
if [[ -s "$plotter/analyse.py" && -s "$plotter/prepare.py" ]]; then
  mkdir -p "$out/code/plotter"
  cp "$plotter/analyse.py" "$plotter/prepare.py" "$out/code/plotter/"
  render=yes
fi
git -C "$repo" rev-parse HEAD > "$out/code_commit.txt"
git -C "$repo" diff -- src validation/real_data/probe_pieridae_components.py validation/real_data/submit_pieridae_components.sh > "$out/source_changes.patch"
find "$out/code" -type f -name '*.py' -exec sha256sum {} + > "$out/code_sha256.txt"
cat > "$out/code/worker.sh" <<'WORKER'
#!/usr/bin/env bash
set -euo pipefail
source_run=$1; out=$2; python=$3; seconds=$4; mode=$5; render=$6
export PYTHONPATH="$out/code/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="$out/matplotlib_cache"
export PATH="$(dirname "$python"):$PATH"
if [[ -d "$source_run/pdf_tools/bin" ]]; then export PATH="$source_run/pdf_tools/bin:$PATH"; fi
args=("$source_run" "$out" --seconds "$seconds")
if [[ $mode == --audit-only ]]; then args+=(--audit-only); fi
if [[ $mode == --recover ]]; then args+=(--recover); fi
if [[ $render == yes && $mode != --audit-only ]]; then args+=(--render); fi
"$python" "$out/code/probe_pieridae_components.py" "${args[@]}"
WORKER
bash -n "$out/code/worker.sh"
wall_minutes=$(( (seconds+59)/60 + 60 ))
job=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=16G --time="$wall_minutes" --job-name=st_components --output="$out/logs/components.%j.out" --error="$out/logs/components.%j.err" "$out/code/worker.sh" "$source_run" "$out" "$python" "$seconds" "$mode" "$render")
job=${job%%;*}
printf 'stage\tjob_id\ncomponent_test\t%s\n' "$job" > "$out/jobs.tsv"
printf 'Component test: %s\nResults: %s\nRead: tail -n 12 %s/logs/components.%s.{out,err}\n' "$job" "$out" "$out" "$job"

fi
