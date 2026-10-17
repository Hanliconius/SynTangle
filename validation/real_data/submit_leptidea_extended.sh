#!/usr/bin/env bash
set -euo pipefail
SCT_LX_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_LX_REPO"
SCT_LX_SOURCE=$(realpath "${1:-local_results/leptidea_paper_SaryLn}")
SCT_LX_PY=$(realpath local_results/annotated_lep_rjv857/tools/bin/python)
SCT_LX_DEP=$(cat local_results/latest_preservation_job.txt)
[[ "$SCT_LX_DEP" =~ ^[0-9]+$ ]]
SCT_LX_RUN=$(mktemp -d "$SCT_LX_REPO/local_results/leptidea_extended_XXXXXX")
mkdir -p "$SCT_LX_RUN/code" logs
cp -a src validation "$SCT_LX_RUN/code/"
export SCT_LX_REPO SCT_LX_SOURCE SCT_LX_PY SCT_LX_RUN
SCT_LX_JOB=$(sbatch --parsable --partition=cpu --array=0-2 --cpus-per-task=1 --mem=16G --time=00:50:00 --dependency="afterok:$SCT_LX_DEP" --job-name=st_lep_extend --output=logs/st_lep_extend.%A_%a.out --error=logs/st_lep_extend.%A_%a.err <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$SCT_LX_RUN/code/src:$SCT_LX_RUN/code/validation/real_data" MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$SCT_LX_PY" "$SCT_LX_RUN/code/validation/real_data/leptidea_extended.py" run --source "$SCT_LX_SOURCE" --run "$SCT_LX_RUN" --index "$SLURM_ARRAY_TASK_ID"
SLURM
)
SCT_LX_COLLECT=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=4G --time=00:10:00 --dependency="afterany:${SCT_LX_JOB%%;*}" --job-name=st_lep_ext_zip --output=logs/st_lep_ext_zip.%j.out --error=logs/st_lep_ext_zip.%j.err <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$SCT_LX_RUN/code/src:$SCT_LX_RUN/code/validation/real_data" MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$SCT_LX_PY" "$SCT_LX_RUN/code/validation/real_data/leptidea_extended.py" collect --source "$SCT_LX_SOURCE" --run "$SCT_LX_RUN"
SLURM
)
printf 'Extended Leptidea, all three concurrent: %s\nCollector: %s\nResults: %s\n' "$SCT_LX_JOB" "$SCT_LX_COLLECT" "$SCT_LX_RUN"
