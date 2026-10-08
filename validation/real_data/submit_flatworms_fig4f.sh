#!/bin/bash
set -euo pipefail
SCT_F_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_F_REPO"
export SCT_F_PREVIOUS=${1:-$(cat local_results/latest_planarian_visual_run.txt)}
SCT_F_PREVIOUS=$(cd "$SCT_F_PREVIOUS" && pwd)
export SCT_F_PREVIOUS
test -x "$SCT_F_PREVIOUS/venv/bin/python"
export SCT_F_MAMBA=${MAMBA_EXE:-$(type -P micromamba)}
export SCT_REAL_GS_ENV=${SCT_REAL_GS_ENV:-lep_busco_painter_clean}
export SCT_REAL_SECONDS=${SCT_REAL_SECONDS:-150}
mkdir -p logs local_results
export SCT_F_RUN=$(mktemp -d "$SCT_F_REPO/local_results/flatworms_fig4f_XXXXXX")
mkdir -p "$SCT_F_RUN/code"
cp -a src validation "$SCT_F_RUN/code/"
git rev-parse HEAD > "$SCT_F_RUN/checkout_head.txt"
printf '%s\n' "$SCT_F_RUN" > local_results/latest_flatworms_fig4f_run.txt
printf '%s\n' '{"datasets":[{"id":"flatworms_fig4f"}]}' > "$SCT_F_RUN/comparison_manifest.json"
sbatch --chdir="$SCT_F_REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_flatworm_f
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=logs/st_flatworm_f.%j.out
#SBATCH --error=logs/st_flatworm_f.%j.err
set -euo pipefail
export PYTHONPATH="$SCT_F_RUN/code/src" PYTHONUNBUFFERED=1 MPLBACKEND=Agg OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$SCT_F_PREVIOUS/venv/bin/python" "$SCT_F_RUN/code/validation/real_data/prepare_planarian_fig4f.py" --output "$SCT_F_RUN/prepared/flatworms_fig4f"
eval "$("$SCT_F_MAMBA" shell hook --shell bash)"
micromamba activate "$SCT_REAL_GS_ENV"
"$SCT_F_PREVIOUS/venv/bin/python" "$SCT_F_RUN/code/validation/real_data/plot_comparison.py" --prepared "$SCT_F_RUN/prepared/flatworms_fig4f" --output "$SCT_F_RUN/figures/flatworms_fig4f" --seconds "$SCT_REAL_SECONDS"
"$SCT_F_PREVIOUS/venv/bin/python" "$SCT_F_RUN/code/validation/real_data/collect_comparisons.py" --run "$SCT_F_RUN" --manifest "$SCT_F_RUN/comparison_manifest.json" --latest-archive "$SCT_F_REPO/local_results/SynTangle_flatworms_fig4f.zip"
SLURM
printf 'Results: %s\nArchive when complete: %s/local_results/SynTangle_flatworms_fig4f.zip\n' "$SCT_F_RUN" "$SCT_F_REPO"
