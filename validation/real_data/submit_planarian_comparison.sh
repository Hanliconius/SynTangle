#!/bin/bash
set -euo pipefail
SCT_REAL_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_REAL_REPO"
if ! grep -q -- '--skip-incomplete-variants' validation/benchmark/genespace_native_order.R; then
    printf '%s\n' 'Outdated R helper. Restore validation/benchmark/genespace_native_order.R from the fetched branch before submitting.' >&2
    false
fi
SCT_REAL_FINAL=${1:-$(cat local_results/latest_published_inputs_run.txt)}
SCT_REAL_FINAL=$(cd "$SCT_REAL_FINAL" && pwd)
test -f "$SCT_REAL_FINAL/prepared/planarians/fixture.json"
test -f "$SCT_REAL_FINAL/prepared/planarians/provenance.json"
export SCT_REAL_REPO SCT_REAL_FINAL
export SCT_REAL_SECONDS=${SCT_REAL_SECONDS:-150}
export SCT_REAL_GS_ENV=${SCT_REAL_GS_ENV:-lep_busco_painter_clean}
export SCT_REAL_MAMBA=${MAMBA_EXE:-$(type -P micromamba)}
test -x "$SCT_REAL_MAMBA"
mkdir -p logs local_results
export SCT_REAL_OUTPUT=$(mktemp -d "$SCT_REAL_REPO/local_results/planarian_visual_XXXXXX")
mkdir -p "$SCT_REAL_OUTPUT/code"
cp -a src validation "$SCT_REAL_OUTPUT/code/"
git rev-parse HEAD > "$SCT_REAL_OUTPUT/checkout_head.txt"
mkdir -p "$SCT_REAL_OUTPUT/prepared"
cp -a "$SCT_REAL_FINAL/prepared/planarians" "$SCT_REAL_OUTPUT/prepared/planarians"
printf '%s\n' '{"datasets":[{"id":"planarians"}]}' > "$SCT_REAL_OUTPUT/comparison_manifest.json"
printf '%s\n' "$SCT_REAL_OUTPUT" > local_results/latest_planarian_visual_run.txt
SCT_REAL_SETUP=$(sbatch --parsable --chdir="$SCT_REAL_REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_planarian_fig_setup
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=logs/st_planarian_fig_setup.%j.out
#SBATCH --error=logs/st_planarian_fig_setup.%j.err
set -euo pipefail
eval "$("$SCT_REAL_MAMBA" shell hook --shell bash)"
micromamba activate syntangle_test
python -m venv "$SCT_REAL_OUTPUT/venv"
"$SCT_REAL_OUTPUT/venv/bin/python" -m pip install -r "$SCT_REAL_OUTPUT/code/validation/benchmark/hybrid_methods_requirements.txt" matplotlib==3.10.8 pypdf==6.1.1
SLURM
)
SCT_REAL_ARRAY=$(sbatch --parsable --dependency="afterok:${SCT_REAL_SETUP%%;*}" --chdir="$SCT_REAL_REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_planarian_fig
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=00:20:00
#SBATCH --output=logs/st_planarian_fig.%j.out
#SBATCH --error=logs/st_planarian_fig.%j.err
set -euo pipefail
eval "$("$SCT_REAL_MAMBA" shell hook --shell bash)"
micromamba activate "$SCT_REAL_GS_ENV"
export PYTHONPATH="$SCT_REAL_OUTPUT/code/src" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
SCT_REAL_CASE=planarians
"$SCT_REAL_OUTPUT/venv/bin/python" "$SCT_REAL_OUTPUT/code/validation/real_data/plot_comparison.py" \
  --prepared "$SCT_REAL_OUTPUT/prepared/$SCT_REAL_CASE" \
  --output "$SCT_REAL_OUTPUT/figures/$SCT_REAL_CASE" --seconds "$SCT_REAL_SECONDS"
SLURM
)
SCT_REAL_COLLECT=$(sbatch --parsable --dependency="afterany:${SCT_REAL_ARRAY%%;*}" --chdir="$SCT_REAL_REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_planarian_fig_collect
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/st_planarian_fig_collect.%j.out
#SBATCH --error=logs/st_planarian_fig_collect.%j.err
set -euo pipefail
"$SCT_REAL_OUTPUT/venv/bin/python" "$SCT_REAL_OUTPUT/code/validation/real_data/collect_comparisons.py" \
  --run "$SCT_REAL_OUTPUT" --manifest "$SCT_REAL_OUTPUT/comparison_manifest.json" \
  --latest-archive "$SCT_REAL_REPO/local_results/SynTangle_planarian_comparisons.zip"
SLURM
)
printf 'Setup: %s\nPlanarian comparison: %s\nCollector: %s\nPDF directory: %s/figures\nDownload archive: %s/local_results/SynTangle_planarian_comparisons.zip\n' "${SCT_REAL_SETUP%%;*}" "${SCT_REAL_ARRAY%%;*}" "${SCT_REAL_COLLECT%%;*}" "$SCT_REAL_OUTPUT" "$SCT_REAL_REPO"
