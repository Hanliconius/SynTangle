#!/usr/bin/env bash
set -euo pipefail
SCT_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_REPO"
SCT_INPUTS=$(realpath "${1:-local_results/rhynchospora_inputs_UkAKU9}")
SCT_PUBLISHED=$(realpath "${2:-local_results/rhynchospora_sources_ACs9FQ}")
test -f "$SCT_INPUTS/input_readiness.json"
test -f "$SCT_PUBLISHED/published_chromosome_geometry.tsv"
SCT_MAMBA=${MAMBA_EXE:-$(type -P micromamba)}
mkdir -p logs local_results
SCT_RUN=$(mktemp -d "$SCT_REPO/local_results/rhynchospora_followup_XXXXXX")
mkdir -p "$SCT_RUN/code"
cp validation/real_data/{rhynchospora_followup.py,prepare_rhynchospora_inputs.py,capture_rhynchospora_fig1a.py,prepare_published_inputs.py} "$SCT_RUN/code/"
export SCT_REPO SCT_INPUTS SCT_PUBLISHED SCT_MAMBA SCT_RUN
printf '%s\n' "$SCT_RUN" > local_results/latest_rhynchospora_followup_run.txt
SCT_SETUP=$(sbatch --parsable --partition=cpu --cpus-per-task=2 --mem=8G --time=02:00:00 --job-name=st_rhync_set --output="$SCT_REPO/logs/st_rhync_setup.%j.out" --error="$SCT_REPO/logs/st_rhync_setup.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
"$SCT_MAMBA" create -y -p "$SCT_RUN/tools" -c conda-forge -c bioconda python=3.10 liftoff=1.6.3 minimap2 > "$SCT_RUN/environment_setup.log" 2>&1
"$SCT_MAMBA" list -p "$SCT_RUN/tools" --json > "$SCT_RUN/tool_versions.json"
"$SCT_RUN/tools/bin/python" "$SCT_RUN/code/rhynchospora_followup.py" --stage prepare --run "$SCT_RUN" --inputs "$SCT_INPUTS" --published "$SCT_PUBLISHED"
SLURM
)
SCT_LIFT=$(sbatch --parsable --dependency="afterok:${SCT_SETUP%%;*}" --partition=cpu --cpus-per-task=8 --mem=32G --time=12:00:00 --job-name=st_rhync_lift --output="$SCT_REPO/logs/st_rhync_lift.%j.out" --error="$SCT_REPO/logs/st_rhync_lift.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PATH="$SCT_RUN/tools/bin:$PATH" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python "$SCT_RUN/code/rhynchospora_followup.py" --stage lift --run "$SCT_RUN" --threads "$SLURM_CPUS_PER_TASK"
SLURM
)
SCT_AUDIT=$(sbatch --parsable --dependency="afterany:${SCT_LIFT%%;*}" --partition=cpu --cpus-per-task=1 --mem=4G --time=00:30:00 --job-name=st_rhync_audit --output="$SCT_REPO/logs/st_rhync_annotation.%j.out" --error="$SCT_REPO/logs/st_rhync_annotation.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
"$SCT_RUN/tools/bin/python" "$SCT_RUN/code/rhynchospora_followup.py" --stage collect --run "$SCT_RUN"
SLURM
)
printf 'Setup/mapping candidates: %s\nAnnotation transfer: %s\nAudit: %s\nResults: %s\nRead when finished: cat logs/st_rhync_annotation.%s.out\n' "$SCT_SETUP" "$SCT_LIFT" "$SCT_AUDIT" "$SCT_RUN" "${SCT_AUDIT%%;*}"
