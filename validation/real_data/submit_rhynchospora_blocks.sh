#!/usr/bin/env bash
set -euo pipefail
SCT_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_REPO"
SCT_SOURCE=$(realpath "${1:-local_results/rhynchospora_visual_PTiYgr}")
test -f "$SCT_SOURCE/prepared/rhynchospora_transfer_reanalysis/fixture.json"
test -x "$SCT_SOURCE/venv/bin/python"
SCT_MAMBA=${MAMBA_EXE:-$(type -P micromamba)}
SCT_GS_ENV=${SCT_REAL_GS_ENV:-lep_busco_painter_clean}
SCT_RUN=$(mktemp -d "$SCT_REPO/local_results/rhynchospora_blocks_XXXXXX")
mkdir -p "$SCT_RUN/code" logs
cp -a src validation "$SCT_RUN/code/"
printf '%s\n' "$SCT_RUN" > local_results/latest_rhynchospora_blocks_run.txt
printf '%s\n' '{"datasets":[{"id":"rhynchospora_collinear_runs"}]}' > "$SCT_RUN/comparison_manifest.json"
export SCT_REPO SCT_SOURCE SCT_RUN SCT_MAMBA SCT_GS_ENV
SCT_JOB=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=16G --time=02:00:00 --job-name=st_rh_blocks --output="$SCT_REPO/logs/st_rh_blocks.%j.out" --error="$SCT_REPO/logs/st_rh_blocks.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$SCT_RUN/code/src" MPLBACKEND=Agg PYTHONUNBUFFERED=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
SCT_PY="$SCT_SOURCE/venv/bin/python"
"$SCT_PY" "$SCT_RUN/code/validation/real_data/rhynchospora_blocks.py" --source "$SCT_SOURCE/prepared/rhynchospora_transfer_reanalysis" --output "$SCT_RUN/prepared/rhynchospora_collinear_runs"
eval "$("$SCT_MAMBA" shell hook --shell bash)"
micromamba activate "$SCT_GS_ENV"
"$SCT_PY" "$SCT_RUN/code/validation/real_data/plot_rhynchospora_blocks.py" --prepared "$SCT_RUN/prepared/rhynchospora_collinear_runs" --genes "$SCT_SOURCE/prepared/rhynchospora_transfer_reanalysis" --output "$SCT_RUN/figures/rhynchospora_collinear_runs" --seconds 300
SLURM
)
SCT_ZIP=$(sbatch --parsable --partition=nano --cpus-per-task=1 --mem=4G --time=00:10:00 --dependency="afterany:${SCT_JOB%%;*}" --job-name=st_rh_block_zip --output="$SCT_REPO/logs/st_rh_blocks_zip.%j.out" --error="$SCT_REPO/logs/st_rh_blocks_zip.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
"$SCT_SOURCE/venv/bin/python" "$SCT_RUN/code/validation/real_data/collect_comparisons.py" --run "$SCT_RUN" --manifest "$SCT_RUN/comparison_manifest.json" --latest-archive "$SCT_REPO/local_results/SynTangle_rhynchospora_blocks.zip"
SLURM
)
printf 'Block comparison: %s\nCollector: %s\nResults: %s\n' "$SCT_JOB" "$SCT_ZIP" "$SCT_RUN"
