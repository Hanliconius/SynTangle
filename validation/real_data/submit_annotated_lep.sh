#!/bin/bash
set -euo pipefail
SCT_LEP_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_LEP_REPO"
if ! grep -q -- '--skip-incomplete-variants' validation/benchmark/genespace_native_order.R; then
    printf '%s\n' 'Restore the current validation/benchmark/genespace_native_order.R before submitting.' >&2
    false
fi
export SCT_LEP_REPO SCT_LEP_MAMBA=${MAMBA_EXE:-$(type -P micromamba)}
export SCT_REAL_GS_ENV=${SCT_REAL_GS_ENV:-lep_busco_painter_clean}
export SCT_REAL_SECONDS=${SCT_REAL_SECONDS:-150}
mkdir -p logs local_results
export SCT_LEP_RUN=$(mktemp -d "$SCT_LEP_REPO/local_results/annotated_lep_XXXXXX")
mkdir -p "$SCT_LEP_RUN/code"
cp -a src validation "$SCT_LEP_RUN/code/"
git rev-parse HEAD > "$SCT_LEP_RUN/checkout_head.txt"
printf '%s\n' '{"datasets":[{"id":"annotated_lep_reanalysis"}]}' > "$SCT_LEP_RUN/comparison_manifest.json"
printf '%s\n' "$SCT_LEP_RUN" > local_results/latest_annotated_lep_run.txt
SCT_LEP_SETUP=$(sbatch --parsable --chdir="$SCT_LEP_REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_lep_setup
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/st_lep_setup.%j.out
#SBATCH --error=logs/st_lep_setup.%j.err
set -euo pipefail
"$SCT_LEP_MAMBA" create -y -p "$SCT_LEP_RUN/tools" -c conda-forge -c bioconda python=3.11 pip blast mcscanx
"$SCT_LEP_MAMBA" list -p "$SCT_LEP_RUN/tools" --json > "$SCT_LEP_RUN/tool_versions.json"
"$SCT_LEP_RUN/tools/bin/python" -m pip install -r "$SCT_LEP_RUN/code/validation/benchmark/hybrid_methods_requirements.txt" matplotlib==3.10.8 pypdf==6.1.1
export PYTHONPATH="$SCT_LEP_RUN/code/src" PYTHONUNBUFFERED=1
"$SCT_LEP_RUN/tools/bin/python" "$SCT_LEP_RUN/code/validation/real_data/annotated_lep_pipeline.py" --stage prepare --run "$SCT_LEP_RUN" --manifest "$SCT_LEP_RUN/code/validation/real_data/annotated_lep.json"
SLURM
)
SCT_LEP_PAIRS=$(sbatch --parsable --dependency="afterok:${SCT_LEP_SETUP%%;*}" --chdir="$SCT_LEP_REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_lep_blocks
#SBATCH --partition=cpu
#SBATCH --array=0-5
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/st_lep_blocks.%A_%a.out
#SBATCH --error=logs/st_lep_blocks.%A_%a.err
set -euo pipefail
export PATH="$SCT_LEP_RUN/tools/bin:$PATH" PYTHONPATH="$SCT_LEP_RUN/code/src" PYTHONUNBUFFERED=1
python "$SCT_LEP_RUN/code/validation/real_data/annotated_lep_pipeline.py" --stage pair --index "$SLURM_ARRAY_TASK_ID" --run "$SCT_LEP_RUN" --manifest "$SCT_LEP_RUN/code/validation/real_data/annotated_lep.json"
SLURM
)
SCT_LEP_COMPARE=$(sbatch --parsable --dependency="afterok:${SCT_LEP_PAIRS%%;*}" --chdir="$SCT_LEP_REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_lep_compare
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=logs/st_lep_compare.%j.out
#SBATCH --error=logs/st_lep_compare.%j.err
set -euo pipefail
export PYTHONPATH="$SCT_LEP_RUN/code/src" PYTHONUNBUFFERED=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
"$SCT_LEP_RUN/tools/bin/python" "$SCT_LEP_RUN/code/validation/real_data/annotated_lep_pipeline.py" --stage collect --run "$SCT_LEP_RUN" --manifest "$SCT_LEP_RUN/code/validation/real_data/annotated_lep.json"
eval "$("$SCT_LEP_MAMBA" shell hook --shell bash)"
micromamba activate "$SCT_REAL_GS_ENV"
"$SCT_LEP_RUN/tools/bin/python" "$SCT_LEP_RUN/code/validation/real_data/plot_comparison.py" --prepared "$SCT_LEP_RUN/prepared/annotated_lep_reanalysis" --output "$SCT_LEP_RUN/figures/annotated_lep_reanalysis" --seconds "$SCT_REAL_SECONDS"
SLURM
)
SCT_LEP_COLLECT=$(sbatch --parsable --dependency="afterany:${SCT_LEP_COMPARE%%;*}" --chdir="$SCT_LEP_REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_lep_collect
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/st_lep_collect.%j.out
#SBATCH --error=logs/st_lep_collect.%j.err
set -euo pipefail
"$SCT_LEP_RUN/tools/bin/python" "$SCT_LEP_RUN/code/validation/real_data/collect_comparisons.py" --run "$SCT_LEP_RUN" --manifest "$SCT_LEP_RUN/comparison_manifest.json" --latest-archive "$SCT_LEP_REPO/local_results/SynTangle_annotated_lep_comparisons.zip"
SLURM
)
printf 'Annotation setup: %s\nAll six pairwise block jobs concurrently: %s\nComparison: %s\nCollector: %s\nResults: %s\n' "${SCT_LEP_SETUP%%;*}" "${SCT_LEP_PAIRS%%;*}" "${SCT_LEP_COMPARE%%;*}" "${SCT_LEP_COLLECT%%;*}" "$SCT_LEP_RUN"
