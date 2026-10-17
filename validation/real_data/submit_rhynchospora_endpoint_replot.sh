#!/usr/bin/env bash
set -euo pipefail
SCT_EP_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_EP_REPO"
SCT_EP_RUN=$(realpath "${1:-local_results/rhynchospora_repaired_e8yRYz}")
SCT_EP_PY=$(realpath local_results/annotated_lep_rjv857/tools/bin/python)
test -f "$SCT_EP_RUN/figures/rhynchospora_native_regions/comparison_audit.json"
SCT_EP_CODE=$(mktemp -d "$SCT_EP_RUN/endpoint_code_XXXXXX")
cp -a src validation "$SCT_EP_CODE/"
export SCT_EP_REPO SCT_EP_RUN SCT_EP_PY SCT_EP_CODE
sbatch --partition=nano --cpus-per-task=1 --mem=4G --time=00:10:00 --job-name=st_rh_endpoint --output=logs/st_rh_endpoint.%j.out --error=logs/st_rh_endpoint.%j.err <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$SCT_EP_CODE/src:$SCT_EP_CODE/validation/real_data" MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$SCT_EP_PY" -m unittest test_ribbon_whole_flip
"$SCT_EP_PY" "$SCT_EP_CODE/validation/real_data/replot_rhynchospora_endpoints.py" --run "$SCT_EP_RUN"
SLURM
