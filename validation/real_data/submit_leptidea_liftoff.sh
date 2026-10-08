#!/usr/bin/env bash
set -euo pipefail
SCT_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_REPO"
SCT_ASSEMBLIES=$(realpath "${1:-$(cat local_results/latest_published_inputs_run.txt)/prepared/leptidea}")
test -f "$SCT_ASSEMBLIES/provenance.json"
SCT_MAMBA=${MAMBA_EXE:-$(type -P micromamba)}
test -x "$SCT_MAMBA"
mkdir -p logs local_results
SCT_RUN=$(mktemp -d "$SCT_REPO/local_results/leptidea_liftoff_XXXXXX")
mkdir -p "$SCT_RUN/code/validation/real_data"
cp validation/real_data/{leptidea_liftoff.py,prepare_published_inputs.py} "$SCT_RUN/code/validation/real_data/"
export SCT_REPO SCT_ASSEMBLIES SCT_MAMBA SCT_RUN
printf '%s\n' "$SCT_RUN" > local_results/latest_leptidea_liftoff_run.txt
SCT_SETUP=$(sbatch --parsable --partition=cpu --cpus-per-task=2 --mem=8G --time=02:00:00 --job-name=st_lift_setup --output="$SCT_REPO/logs/st_lift_setup.%j.out" --error="$SCT_REPO/logs/st_lift_setup.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
"$SCT_MAMBA" create -y -p "$SCT_RUN/tools" -c conda-forge -c bioconda python=3.10 liftoff=1.6.3 minimap2
"$SCT_MAMBA" list -p "$SCT_RUN/tools" --json > "$SCT_RUN/tool_versions.json"
"$SCT_RUN/tools/bin/python" "$SCT_RUN/code/validation/real_data/leptidea_liftoff.py" --stage prepare --run "$SCT_RUN" --assemblies "$SCT_ASSEMBLIES"
SLURM
)
SCT_ARRAY=$(sbatch --parsable --dependency="afterok:${SCT_SETUP%%;*}" --partition=cpu --cpus-per-task=8 --mem=32G --time=12:00:00 --array=0-3 --job-name=st_liftoff --output="$SCT_REPO/logs/st_liftoff.%A_%a.out" --error="$SCT_REPO/logs/st_liftoff.%A_%a.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PATH="$SCT_RUN/tools/bin:$PATH" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python "$SCT_RUN/code/validation/real_data/leptidea_liftoff.py" --stage lift --run "$SCT_RUN" --index "$SLURM_ARRAY_TASK_ID" --threads "$SLURM_CPUS_PER_TASK"
SLURM
)
SCT_COLLECT=$(sbatch --parsable --dependency="afterok:${SCT_ARRAY%%;*}" --partition=cpu --cpus-per-task=1 --mem=4G --time=00:15:00 --job-name=st_lift_audit --output="$SCT_REPO/logs/st_lift_audit.%j.out" --error="$SCT_REPO/logs/st_lift_audit.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
"$SCT_RUN/tools/bin/python" "$SCT_RUN/code/validation/real_data/leptidea_liftoff.py" --stage collect --run "$SCT_RUN"
SLURM
)
printf 'Setup: %s\nLiftoff (all four concurrently): %s\nAnnotation audit: %s\nResults: %s\n' "$SCT_SETUP" "$SCT_ARRAY" "$SCT_COLLECT" "$SCT_RUN"
