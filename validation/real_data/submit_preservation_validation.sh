#!/usr/bin/env bash
set -euo pipefail
SCT_PV_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_PV_REPO"
SCT_PV_RUN=$(mktemp -d "$SCT_PV_REPO/local_results/preservation_XXXXXX")
SCT_PV_PY=$(realpath local_results/annotated_lep_rjv857/tools/bin/python)
SCT_PV_LEP=$(realpath local_results/leptidea_paper_SaryLn)
SCT_PV_RH=$(realpath local_results/rhynchospora_repaired_e8yRYz)
mkdir -p "$SCT_PV_RUN/code" logs
cp -a src validation "$SCT_PV_RUN/code/"
export SCT_PV_REPO SCT_PV_RUN SCT_PV_PY SCT_PV_LEP SCT_PV_RH
SCT_PV_JOB=$(sbatch --parsable --partition=cpu --array=0-2 --cpus-per-task=1 --mem=4G --time=00:20:00 --job-name=st_preserve --output=logs/st_preserve.%A_%a.out --error=logs/st_preserve.%A_%a.err <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$SCT_PV_RUN/code/src:$SCT_PV_RUN/code/validation/real_data" MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
case "$SLURM_ARRAY_TASK_ID" in
0) "$SCT_PV_PY" -m unittest test_preservation_audit test_ribbon_whole_flip test_leptidea_paper_comparison test_rhynchospora_native_regions > "$SCT_PV_RUN/tests.log" 2>&1 ;;
1) "$SCT_PV_PY" "$SCT_PV_RUN/code/validation/real_data/replot_leptidea_endpoints.py" --run "$SCT_PV_LEP" ;;
2) "$SCT_PV_PY" "$SCT_PV_RUN/code/validation/real_data/replot_rhynchospora_endpoints.py" --run "$SCT_PV_RH" ;;
esac
printf 'PASS preservation task %s\n' "$SLURM_ARRAY_TASK_ID"
SLURM
)
printf '%s\n' "${SCT_PV_JOB%%;*}" > local_results/latest_preservation_job.txt
printf 'Preservation validation, all three concurrent: %s\nResults: %s\n' "$SCT_PV_JOB" "$SCT_PV_RUN"
