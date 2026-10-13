#!/usr/bin/env bash
set -euo pipefail
SCT_RECON_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_RECON_REPO"
SCT_RECON_CACHE=$(realpath "${1:-$(cat local_results/latest_simple_papers_run.txt)}")
SCT_RECON_RUN=$(mktemp -d "$SCT_RECON_REPO/local_results/simple_reconstructed_XXXXXX")
mkdir -p "$SCT_RECON_RUN/code" logs
cp -a src validation "$SCT_RECON_RUN/code/"
printf '%s\n' "$SCT_RECON_RUN" > local_results/latest_simple_reconstruction_run.txt
SCT_RECON_TOOLS="$SCT_RECON_RUN/tools"
export SCT_RECON_RUN SCT_RECON_CACHE SCT_RECON_TOOLS
SCT_RECON_SETUP=$(sbatch --parsable --partition=cpu --cpus-per-task=2 --mem=8G --time=01:00:00 --job-name=st_simple_tools --output=logs/st_simple_tools.%j.out --error=logs/st_simple_tools.%j.err <<'SLURM'
#!/bin/bash
set -euo pipefail
micromamba create -y -p "$SCT_RECON_TOOLS" -c conda-forge -c bioconda --strict-channel-priority python=3.11 pip blast mcscanx minimap2 gffread
"$SCT_RECON_TOOLS/bin/python" -m pip install -r "$SCT_RECON_RUN/code/validation/benchmark/hybrid_methods_requirements.txt" matplotlib==3.10.8 pypdf==6.1.1
printf 'PASS reconstruction tools ready\n'
SLURM
)
SCT_RECON_TERMINALS=()
for SCT_RECON_CASE in 0 1 2 3; do
  read -r SCT_RECON_ID SCT_RECON_NSP SCT_RECON_NPAIR < <(python - "$SCT_RECON_CASE" "$SCT_RECON_RUN" <<'PY'
import json,sys
from pathlib import Path
c=json.loads((Path(sys.argv[2])/'code/validation/real_data/simple_reconstruction_manifest.json').read_text())['cases'][int(sys.argv[1])]
print(c['id'],len(c['species']),len(c['pairs']))
PY
)
  export SCT_RECON_CASE SCT_RECON_ID
  SCT_RECON_PREP=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=16G --time=04:00:00 --array="0-$((SCT_RECON_NSP-1))" --dependency="afterok:${SCT_RECON_SETUP%%;*}" --job-name="st_${SCT_RECON_ID}_prep" --output="logs/st_simple_${SCT_RECON_ID}_prep.%A_%a.out" --error="logs/st_simple_${SCT_RECON_ID}_prep.%A_%a.err" <<'SLURM'
#!/bin/bash
set -euo pipefail
export PATH="$SCT_RECON_TOOLS/bin:$PATH" PYTHONPATH="$SCT_RECON_RUN/code/src:$SCT_RECON_RUN/code/validation/real_data" OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1
python "$SCT_RECON_RUN/code/validation/real_data/simple_reconstruction.py" prepare --run "$SCT_RECON_RUN" --case "$SCT_RECON_CASE" --index "$SLURM_ARRAY_TASK_ID" --cache "$SCT_RECON_CACHE/$SCT_RECON_ID/sources"
SLURM
)
  SCT_RECON_PAIR=$(sbatch --parsable --partition=cpu --cpus-per-task=4 --mem=16G --time=04:00:00 --array="0-$((SCT_RECON_NPAIR-1))" --dependency="afterany:${SCT_RECON_PREP%%;*}" --job-name="st_${SCT_RECON_ID}_pairs" --output="logs/st_simple_${SCT_RECON_ID}_pair.%A_%a.out" --error="logs/st_simple_${SCT_RECON_ID}_pair.%A_%a.err" <<'SLURM'
#!/bin/bash
set -euo pipefail
export PATH="$SCT_RECON_TOOLS/bin:$PATH" PYTHONPATH="$SCT_RECON_RUN/code/src:$SCT_RECON_RUN/code/validation/real_data" OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1
python "$SCT_RECON_RUN/code/validation/real_data/simple_reconstruction.py" pair --run "$SCT_RECON_RUN" --case "$SCT_RECON_CASE" --index "$SLURM_ARRAY_TASK_ID"
SLURM
)
  SCT_RECON_NFIG=3
  if [[ "$SCT_RECON_ID" == ficus ]]; then SCT_RECON_NFIG=12; fi
  SCT_RECON_FIG=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=16G --time=00:30:00 --array="0-$((SCT_RECON_NFIG-1))" --dependency="afterany:${SCT_RECON_PAIR%%;*}" --job-name="st_${SCT_RECON_ID}_fig" --output="logs/st_simple_${SCT_RECON_ID}_fig.%A_%a.out" --error="logs/st_simple_${SCT_RECON_ID}_fig.%A_%a.err" <<'SLURM'
#!/bin/bash
set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate "${SCT_REAL_GS_ENV:-lep_busco_painter_clean}"
export PATH="$SCT_RECON_TOOLS/bin:$PATH" PYTHONPATH="$SCT_RECON_RUN/code/src:$SCT_RECON_RUN/code/validation/real_data" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
"$SCT_RECON_TOOLS/bin/python" "$SCT_RECON_RUN/code/validation/real_data/simple_reconstruction.py" compare --run "$SCT_RECON_RUN" --case "$SCT_RECON_CASE" --index "$SLURM_ARRAY_TASK_ID" --python "$SCT_RECON_TOOLS/bin/python"
SLURM
)
  SCT_RECON_TERMINALS+=("${SCT_RECON_FIG%%;*}")
  printf '%s: inputs=%s pairs=%s comparisons=%s\n' "$SCT_RECON_ID" "$SCT_RECON_PREP" "$SCT_RECON_PAIR" "$SCT_RECON_FIG"
done
SCT_RECON_DEP=$(IFS=:;printf '%s' "${SCT_RECON_TERMINALS[*]}")
SCT_RECON_COLLECT=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=4G --time=00:20:00 --dependency="afterany:$SCT_RECON_DEP" --job-name=st_simple_summary --output=logs/st_simple_summary.%j.out --error=logs/st_simple_summary.%j.err <<'SLURM'
#!/bin/bash
set -euo pipefail
export PYTHONPATH="$SCT_RECON_RUN/code/src:$SCT_RECON_RUN/code/validation/real_data"
"$SCT_RECON_TOOLS/bin/python" "$SCT_RECON_RUN/code/validation/real_data/simple_reconstruction.py" collect --run "$SCT_RECON_RUN"
SLURM
)
printf 'Setup: %s\nCollector: %s\nResults: %s\n' "$SCT_RECON_SETUP" "$SCT_RECON_COLLECT" "$SCT_RECON_RUN"
