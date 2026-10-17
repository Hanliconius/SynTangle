#!/bin/bash
set -euo pipefail
SCT_PAPER_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_PAPER_REPO"
export SCT_PAPER_RUN=${1:-$(cat local_results/latest_planarian_visual_run.txt)}
SCT_PAPER_RUN=$(cd "$SCT_PAPER_RUN" && pwd)
export SCT_PAPER_RUN
test -f "$SCT_PAPER_RUN/figures/planarians/COMPLETE"
test -x "$SCT_PAPER_RUN/venv/bin/python"
mkdir -p "$SCT_PAPER_RUN/paper_order_code/validation/real_data" logs
cp validation/real_data/plot_planarian_paper_order.py validation/real_data/plot_comparison.py "$SCT_PAPER_RUN/paper_order_code/validation/real_data/"
cp -a validation/benchmark "$SCT_PAPER_RUN/paper_order_code/validation/"
sbatch --chdir="$SCT_PAPER_REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_paper_order
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/st_paper_order.%j.out
#SBATCH --error=logs/st_paper_order.%j.err
set -euo pipefail
export PYTHONPATH="$SCT_PAPER_RUN/code/src" MPLBACKEND=Agg PYTHONUNBUFFERED=1
"$SCT_PAPER_RUN/venv/bin/python" "$SCT_PAPER_RUN/paper_order_code/validation/real_data/plot_planarian_paper_order.py" --run "$SCT_PAPER_RUN"
SLURM
