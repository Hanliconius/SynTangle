#!/usr/bin/env bash
set -euo pipefail
SCT_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_REPO"
SCT_PRIOR=$(realpath "${1:-local_results/rhynchospora_native_mNeWPw}")
SCT_SOURCE=$(realpath local_results/rhynchospora_visual_PTiYgr)
test -f "$SCT_PRIOR/input_audit.json"
test -x "$SCT_SOURCE/venv/bin/python"
SCT_MAMBA=${MAMBA_EXE:-$(type -P micromamba)}
SCT_GS_ENV=${SCT_REAL_GS_ENV:-lep_busco_painter_clean}
SCT_RUN=$(mktemp -d "$SCT_REPO/local_results/rhynchospora_repaired_XXXXXX")
mkdir -p "$SCT_RUN/code" logs
cp -a src validation "$SCT_RUN/code/"
export SCT_REPO SCT_PRIOR SCT_SOURCE SCT_RUN SCT_MAMBA SCT_GS_ENV
SCT_SETUP=$(sbatch --parsable --partition=nano --cpus-per-task=1 --mem=8G --time=00:25:00 --job-name=st_rh_fix_setup --output="$SCT_REPO/logs/st_rh_fix_setup.%j.out" --error="$SCT_REPO/logs/st_rh_fix_setup.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
"$SCT_MAMBA" create -y -p "$SCT_RUN/tools" -c conda-forge -c bioconda python=3.11 numpy=1.26.4 scipy=1.13.1 orthofinder=2.5.5 diamond mcscanx gffread > "$SCT_RUN/environment_setup.log" 2>&1
"$SCT_RUN/tools/bin/python" "$SCT_RUN/code/validation/real_data/repair_rhynchospora_orthofinder.py" prepare --prior "$SCT_PRIOR" --run "$SCT_RUN"
SLURM
)
SCT_JOB=$(sbatch --parsable --partition=cpu --cpus-per-task=4 --mem=32G --time=04:00:00 --dependency="afterok:${SCT_SETUP%%;*}" --job-name=st_rh_fix --output="$SCT_REPO/logs/st_rh_fix.%j.out" --error="$SCT_REPO/logs/st_rh_fix.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$SCT_RUN/code/src" MPLBACKEND=Agg PYTHONUNBUFFERED=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export PATH="$SCT_RUN/tools/bin:$PATH"
cd "$SCT_RUN/code/validation/real_data"
"$SCT_RUN/tools/bin/python" -m unittest test_rhynchospora_orthofinder_repair > "$SCT_RUN/repair_tests.log" 2>&1
"$SCT_RUN/tools/bin/python" repair_rhynchospora_orthofinder.py resume --run "$SCT_RUN" --threads "$SLURM_CPUS_PER_TASK"
eval "$("$SCT_MAMBA" shell hook --shell bash)"
micromamba activate "$SCT_GS_ENV"
export PATH="$SCT_RUN/tools/bin:$PATH"
Rscript "$SCT_RUN/code/validation/real_data/rhynchospora_repaired_regions.R" "$SCT_RUN" "$SCT_RUN/tools/bin"
"$SCT_SOURCE/venv/bin/python" "$SCT_RUN/code/validation/real_data/rhynchospora_native_regions.py" compare --run "$SCT_RUN" --seconds 300
"$SCT_SOURCE/venv/bin/python" "$SCT_RUN/code/validation/real_data/collect_rhynchospora_repaired.py"
SLURM
)
printf 'Repair setup: %s\nAudited discovery and comparison: %s\nResults: %s\nArchive when complete: %s/local_results/SynTangle_rhynchospora_native_repaired.zip\n' "$SCT_SETUP" "$SCT_JOB" "$SCT_RUN" "$SCT_REPO"
