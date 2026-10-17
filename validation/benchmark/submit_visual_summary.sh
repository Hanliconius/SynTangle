#!/bin/bash
set -euo pipefail
SCT_VIS_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_VIS_REPO"
export SCT_VIS_REPO
export SCT_VIS_ZIP=${SCT_VIS_ZIP:-$SCT_VIS_REPO/validation/benchmark/visual_examples/SynTangle_visual_examples.zip}
test -f "$SCT_VIS_ZIP"
mkdir -p logs local_results
export SCT_VIS_OUTPUT=$(mktemp -d "$SCT_VIS_REPO/local_results/visual_XXXXXX")
export SCT_VIS_MAMBA=${MAMBA_EXE:-$(type -P micromamba)}
test -x "$SCT_VIS_MAMBA"
mkdir -p "$SCT_VIS_OUTPUT/code"
cp -a src "$SCT_VIS_OUTPUT/code/"
cp validation/benchmark/visual_summary.py "$SCT_VIS_OUTPUT/code/"
cp "$SCT_VIS_ZIP" "$SCT_VIS_OUTPUT/input.zip"
git rev-parse HEAD > "$SCT_VIS_OUTPUT/commit.txt"
SCT_VIS_JOB=$(sbatch --parsable --chdir="$SCT_VIS_REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_visual
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=logs/st_visual.%j.out
#SBATCH --error=logs/st_visual.%j.err
set -euo pipefail
eval "$("$SCT_VIS_MAMBA" shell hook --shell bash)"
micromamba activate syntangle_test
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
export PYTHONPATH="$SCT_VIS_OUTPUT/code/src"
python -m venv --system-site-packages "$SCT_VIS_OUTPUT/venv"
"$SCT_VIS_OUTPUT/venv/bin/python" -m pip install matplotlib==3.10.8
"$SCT_VIS_OUTPUT/venv/bin/python" - <<'PY'
import os
from pathlib import Path
from zipfile import ZipFile
root=Path(os.environ['SCT_VIS_OUTPUT'])/'data';root.mkdir()
with ZipFile(root.parent/'input.zip') as z:
    for name in z.namelist():
        if not (root/name).resolve().is_relative_to(root.resolve()):raise ValueError('Unsafe archive path')
    z.extractall(root)
PY
"$SCT_VIS_OUTPUT/venv/bin/python" "$SCT_VIS_OUTPUT/code/visual_summary.py" \
  --data-root "$SCT_VIS_OUTPUT/data" --output "$SCT_VIS_OUTPUT/figures"
SLURM
)
printf 'Plotting job: %s\nFigures: %s/figures\n' "${SCT_VIS_JOB%%;*}" "$SCT_VIS_OUTPUT"
