#!/usr/bin/env bash
set -euo pipefail
SCT_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_REPO"
SCT_PRIOR=$(realpath "${1:-local_results/rhynchospora_followup_jY8c5Y}")
test -f "$SCT_PRIOR/annotation/LIFTOFF_COMPLETE"
test -x "$SCT_PRIOR/tools/bin/liftoff"
if ! grep -q -- '--skip-incomplete-variants' validation/benchmark/genespace_native_order.R; then
    printf 'Restore the current genespace_native_order.R before submitting.\n' >&2; exit 1
fi
SCT_MAMBA=${MAMBA_EXE:-$(type -P micromamba)}
SCT_PYTHON=$(type -P python)
SCT_GS_ENV=${SCT_REAL_GS_ENV:-lep_busco_painter_clean}
mkdir -p logs local_results
SCT_RUN=$(mktemp -d "$SCT_REPO/local_results/rhynchospora_visual_XXXXXX")
mkdir -p "$SCT_RUN/code"
cp -a src validation "$SCT_RUN/code/"
printf '%s\n' "$SCT_RUN" > local_results/latest_rhynchospora_visual_run.txt
printf '%s\n' '{"datasets":[{"id":"rhynchospora_transfer_reanalysis"}]}' > "$SCT_RUN/comparison_manifest.json"
git rev-parse HEAD > "$SCT_RUN/checkout_head.txt"
export SCT_REPO SCT_PRIOR SCT_RUN SCT_MAMBA SCT_PYTHON SCT_GS_ENV
SCT_SETUP=$(sbatch --parsable --partition=nano --cpus-per-task=1 --mem=4G --time=00:30:00 --job-name=st_rh_fig_set --output="$SCT_REPO/logs/st_rh_fig_setup.%j.out" --error="$SCT_REPO/logs/st_rh_fig_setup.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
cd "$SCT_RUN/code/validation/real_data"
"$SCT_PYTHON" -m unittest test_rhynchospora_comparison > "$SCT_RUN/importer_tests.log" 2>&1
"$SCT_PYTHON" -m venv "$SCT_RUN/venv"
"$SCT_RUN/venv/bin/python" -m pip install -r "$SCT_RUN/code/validation/benchmark/hybrid_methods_requirements.txt" matplotlib==3.10.8 pypdf==6.1.1 > "$SCT_RUN/dependency_setup.log" 2>&1
"$SCT_PRIOR/tools/bin/python" "$SCT_RUN/code/validation/real_data/rhynchospora_comparison.py" --stage prepare --run "$SCT_RUN" --prior "$SCT_PRIOR"
SLURM
)
SCT_LIFT=$(sbatch --parsable --partition=nano --dependency="afterok:${SCT_SETUP%%;*}" --cpus-per-task=8 --mem=32G --time=00:30:00 --job-name=st_rh_fig_lift --output="$SCT_REPO/logs/st_rh_fig_lift.%j.out" --error="$SCT_REPO/logs/st_rh_fig_lift.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PATH="$SCT_PRIOR/tools/bin:$PATH" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python "$SCT_RUN/code/validation/real_data/rhynchospora_comparison.py" --stage lift --run "$SCT_RUN" --prior "$SCT_PRIOR" --threads "$SLURM_CPUS_PER_TASK"
SLURM
)
SCT_COMPARE=$(sbatch --parsable --partition=nano --dependency="afterok:${SCT_LIFT%%;*}" --cpus-per-task=1 --mem=12G --time=00:30:00 --job-name=st_rh_fig --output="$SCT_REPO/logs/st_rh_fig.%j.out" --error="$SCT_REPO/logs/st_rh_fig.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$SCT_RUN/code/src" MPLBACKEND=Agg PYTHONUNBUFFERED=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$SCT_PRIOR/tools/bin/python" "$SCT_RUN/code/validation/real_data/rhynchospora_comparison.py" --stage fixture --run "$SCT_RUN" --prior "$SCT_PRIOR"
eval "$("$SCT_MAMBA" shell hook --shell bash)"
micromamba activate "$SCT_GS_ENV"
"$SCT_RUN/venv/bin/python" "$SCT_RUN/code/validation/real_data/plot_comparison.py" --prepared "$SCT_RUN/prepared/rhynchospora_transfer_reanalysis" --output "$SCT_RUN/figures/rhynchospora_transfer_reanalysis" --seconds 300
SLURM
)
SCT_COLLECT=$(sbatch --parsable --partition=nano --dependency="afterany:${SCT_COMPARE%%;*}" --cpus-per-task=1 --mem=4G --time=00:10:00 --job-name=st_rh_fig_zip --output="$SCT_REPO/logs/st_rh_fig_collect.%j.out" --error="$SCT_REPO/logs/st_rh_fig_collect.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
"$SCT_RUN/venv/bin/python" "$SCT_RUN/code/validation/real_data/collect_comparisons.py" --run "$SCT_RUN" --manifest "$SCT_RUN/comparison_manifest.json" --latest-archive "$SCT_REPO/local_results/SynTangle_rhynchospora_comparisons.zip"
SLURM
)
printf 'Setup: %s\nAustro transfer: %s\nComparison: %s\nCollector: %s\nResults: %s\nArchive: %s/local_results/SynTangle_rhynchospora_comparisons.zip\n' "$SCT_SETUP" "$SCT_LIFT" "$SCT_COMPARE" "$SCT_COLLECT" "$SCT_RUN" "$SCT_REPO"
