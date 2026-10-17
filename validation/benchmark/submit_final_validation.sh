#!/bin/bash
set -euo pipefail
SCT_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_REPO"
mkdir -p logs local_results
export SCT_FINAL_ROOT=$(mktemp -d "$SCT_REPO/local_results/final_XXXXXX")
if [[ -n "${SCT_REUSE_PREP:-}" ]]; then
    test -d "$SCT_REUSE_PREP"
    cp -a "$SCT_REUSE_PREP" "$SCT_FINAL_ROOT/real_data"
fi
export SCT_FINAL_CODE="$SCT_FINAL_ROOT/code"
export SCT_MAMBA=${MAMBA_EXE:-$(type -P micromamba)}
test -x "$SCT_MAMBA"
mkdir -p "$SCT_FINAL_CODE"
cp -a src tests examples validation "$SCT_FINAL_CODE/"
git rev-parse HEAD > "$SCT_FINAL_ROOT/checkout_head.txt"
python - "$SCT_FINAL_CODE" "$SCT_FINAL_ROOT" <<'PY'
import hashlib,json,sys
from pathlib import Path
root=Path(sys.argv[1]);files={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in root.rglob('*') if p.is_file() and p.suffix in ('.py','.R','.json','.sh','.txt','.md') and '__pycache__' not in p.parts}
Path(sys.argv[2],'code_sha256.json').write_text(json.dumps(files,indent=2)+'\n')
PY
setup=$(sbatch --parsable --chdir="$SCT_REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_final_setup
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=logs/st_final_setup.%j.out
#SBATCH --error=logs/st_final_setup.%j.err
set -euo pipefail
eval "$("$SCT_MAMBA" shell hook --shell bash)"
micromamba activate syntangle_test
python -m venv "$SCT_FINAL_ROOT/venv"
"$SCT_FINAL_ROOT/venv/bin/python" -m pip install -r "$SCT_FINAL_CODE/validation/benchmark/hybrid_methods_requirements.txt" cvxpy==1.7.5
SLURM
)
setup=${setup%%;*}
checks=$(sbatch --parsable --dependency="afterok:$setup" --chdir="$SCT_REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_final_checks
#SBATCH --partition=cpu
#SBATCH --array=0-2
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/st_final_checks.%A_%a.out
#SBATCH --error=logs/st_final_checks.%A_%a.err
set -euo pipefail
export PYTHONPATH="$SCT_FINAL_CODE/src:$SCT_FINAL_CODE/tests" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
cd "$SCT_FINAL_CODE"
case "$SLURM_ARRAY_TASK_ID" in
 0) "$SCT_FINAL_ROOT/venv/bin/python" -m unittest discover -s tests -v ;;
 1) "$SCT_FINAL_ROOT/venv/bin/python" -m unittest test_final_exhaustive -v ;;
 2) "$SCT_FINAL_ROOT/venv/bin/python" -m unittest test_real_data_import test_constraint_components test_search_pipeline test_search_deadline test_tangledness_refinement -v ;;
esac
printf 'PASS final validation group %s\n' "$SLURM_ARRAY_TASK_ID"
SLURM
)
prep=$(sbatch --parsable --dependency="afterok:$setup" --chdir="$SCT_REPO" <<'SLURM'
#!/bin/bash
#SBATCH --job-name=st_real_prep
#SBATCH --partition=nano
#SBATCH --array=0-5
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=logs/st_real_prep.%A_%a.out
#SBATCH --error=logs/st_real_prep.%A_%a.err
set -euo pipefail
export PYTHONPATH="$SCT_FINAL_CODE/src" PYTHONUNBUFFERED=1
"$SCT_FINAL_ROOT/venv/bin/python" "$SCT_FINAL_CODE/validation/real_data/prepare_genespace_blocks.py" --manifest "$SCT_FINAL_CODE/validation/real_data/datasets.json" --index "$SLURM_ARRAY_TASK_ID" --output "$SCT_FINAL_ROOT/real_data"
SLURM
)
printf 'Setup: %s\nFinal checks (all 3 concurrently): %s\nReal-data preparation (all 6 concurrently): %s\nResults: %s\n' "$setup" "${checks%%;*}" "${prep%%;*}" "$SCT_FINAL_ROOT"
printf '%s\n' "$SCT_FINAL_ROOT" > local_results/latest_final_run.txt
