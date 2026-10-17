#!/bin/bash
set -euo pipefail
SCT_PAPER_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_PAPER_REPO"
mkdir -p logs local_results
export SCT_PAPER_RUN=$(mktemp -d "$SCT_PAPER_REPO/local_results/published_inputs_XXXXXX")
mkdir -p "$SCT_PAPER_RUN/code"
cp -a src validation "$SCT_PAPER_RUN/code/"
git rev-parse HEAD > "$SCT_PAPER_RUN/checkout_head.txt"
printf '%s\n' "$SCT_PAPER_RUN" > local_results/latest_published_inputs_run.txt
SCT_PAPER_SETUP=$(sbatch --parsable --chdir="$SCT_PAPER_REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_paper_setup
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=logs/st_paper_setup.%j.out
#SBATCH --error=logs/st_paper_setup.%j.err
set -euo pipefail
python -m venv "$SCT_PAPER_RUN/venv"
"$SCT_PAPER_RUN/venv/bin/python" -m pip install openpyxl==3.1.5
SLURM
)
SCT_PAPER_ARRAY=$(sbatch --parsable --dependency="afterok:${SCT_PAPER_SETUP%%;*}" --chdir="$SCT_PAPER_REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_paper_prep
#SBATCH --partition=cpu
#SBATCH --array=0-1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/st_paper_prep.%A_%a.out
#SBATCH --error=logs/st_paper_prep.%A_%a.err
set -euo pipefail
export PYTHONPATH="$SCT_PAPER_RUN/code/src" PYTHONUNBUFFERED=1
SCT_PAPER_DATASETS=(planarians leptidea)
SCT_PAPER_CASE=${SCT_PAPER_DATASETS[$SLURM_ARRAY_TASK_ID]}
"$SCT_PAPER_RUN/venv/bin/python" "$SCT_PAPER_RUN/code/validation/real_data/prepare_published_inputs.py" \
  --dataset "$SCT_PAPER_CASE" --output "$SCT_PAPER_RUN/prepared/$SCT_PAPER_CASE"
SLURM
)
printf 'Setup: %s\nPreparation array (both concurrently): %s\nResults: %s\n' "${SCT_PAPER_SETUP%%;*}" "${SCT_PAPER_ARRAY%%;*}" "$SCT_PAPER_RUN"
