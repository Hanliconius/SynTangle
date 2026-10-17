#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
run=$(mktemp -d "$repo/local_results/ribbon_objective_XXXXXX")
mkdir -p "$run/code" "$repo/logs"
cp -a "$repo/src" "$repo/validation" "$run/code/"
export SCT_RIBBON_RUN="$run" SCT_RIBBON_REPO="$repo"
array=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=16G --time=00:20:00 --array=0-3 --job-name=st_ribbon --output="$repo/logs/st_ribbon.%A_%a.out" --error="$repo/logs/st_ribbon.%A_%a.err" <<'SLURM'
#!/bin/bash
set -euo pipefail
export PYTHONPATH="$SCT_RIBBON_RUN/code/src:$SCT_RIBBON_RUN/code/validation/real_data:$SCT_RIBBON_RUN/code/validation/benchmark" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg
"$SCT_RIBBON_REPO/local_results/simple_reconstructed_QI1cVs/tools/bin/python" "$SCT_RIBBON_RUN/code/validation/real_data/ribbon_experiment.py" run --run "$SCT_RIBBON_RUN" --index "$SLURM_ARRAY_TASK_ID"
SLURM
)
collector=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=4G --time=00:10:00 --dependency="afterok:${array%%;*}" --job-name=st_rib_zip --output="$repo/logs/st_ribbon_zip.%j.out" --error="$repo/logs/st_ribbon_zip.%j.err" <<'SLURM'
#!/bin/bash
set -euo pipefail
export PYTHONPATH="$SCT_RIBBON_RUN/code/src"
"$SCT_RIBBON_REPO/local_results/simple_reconstructed_QI1cVs/tools/bin/python" "$SCT_RIBBON_RUN/code/validation/real_data/ribbon_experiment.py" collect --run "$SCT_RIBBON_RUN"
SLURM
)
printf 'Ribbon experiments: %s\nCollector: %s\nResults: %s\n' "$array" "$collector" "$run"
