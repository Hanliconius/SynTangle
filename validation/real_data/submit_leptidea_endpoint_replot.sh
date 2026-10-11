#!/usr/bin/env bash
set -euo pipefail
SCT_LE_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_LE_REPO"
SCT_LE_RUN=$(realpath "${1:-local_results/leptidea_paper_SaryLn}")
SCT_LE_PY=$(realpath local_results/annotated_lep_rjv857/tools/bin/python)
test -f "$SCT_LE_RUN/figures/leptidea_paper_assemblies/comparison_audit.json"
SCT_LE_CODE=$(mktemp -d "$SCT_LE_RUN/endpoint_code_XXXXXX")
cp -a src validation "$SCT_LE_CODE/"
export SCT_LE_RUN SCT_LE_PY SCT_LE_CODE
sbatch --partition=nano --cpus-per-task=1 --mem=4G --time=00:10:00 --job-name=st_lep_endpoint --output=logs/st_lep_endpoint.%j.out --error=logs/st_lep_endpoint.%j.err <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$SCT_LE_CODE/src:$SCT_LE_CODE/validation/real_data" MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$SCT_LE_PY" -m unittest test_ribbon_whole_flip
"$SCT_LE_PY" "$SCT_LE_CODE/validation/real_data/replot_leptidea_endpoints.py" --run "$SCT_LE_RUN"
SLURM
