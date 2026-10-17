#!/usr/bin/env bash
set -euo pipefail
SCT_PC_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_PC_REPO"
SCT_PC_TOOLS=$(realpath local_results/annotated_lep_rjv857/tools)
SCT_PC_GFFREAD=$(realpath local_results/rhynchospora_native_mNeWPw/tools/bin/gffread)
SCT_PC_MAMBA=${MAMBA_EXE:-$(type -P micromamba)}
SCT_PC_GS_ENV=${SCT_REAL_GS_ENV:-lep_busco_painter_clean}
for exe in python blastp makeblastdb MCScanX; do test -x "$SCT_PC_TOOLS/bin/$exe"; done
test -x "$SCT_PC_GFFREAD"
SCT_PC_ROOT=$(mktemp -d "$SCT_PC_REPO/local_results/codfish_carex_XXXXXX")
mkdir -p "$SCT_PC_ROOT/code" logs
cp -a src validation "$SCT_PC_ROOT/code/"
SCT_PC_MANIFEST="$SCT_PC_ROOT/code/validation/real_data/codfish_carex_manifest.json"
printf '%s\n' "$SCT_PC_ROOT" > local_results/latest_codfish_carex_run.txt
export SCT_PC_REPO SCT_PC_TOOLS SCT_PC_GFFREAD SCT_PC_MAMBA SCT_PC_GS_ENV SCT_PC_ROOT SCT_PC_MANIFEST
for SCT_PC_CASE in 0 1; do
  if [[ "$SCT_PC_CASE" == 0 ]]; then SCT_PC_SPECIES=5; SCT_PC_PAIRS=14; SCT_PC_NAME=codfish; else SCT_PC_SPECIES=3; SCT_PC_PAIRS=5; SCT_PC_NAME=carex; fi
  SCT_PC_RUN="$SCT_PC_ROOT/$SCT_PC_NAME"
  mkdir -p "$SCT_PC_RUN"
  export SCT_PC_CASE SCT_PC_RUN SCT_PC_NAME
  SCT_PC_PREP=$(sbatch --parsable --partition=cpu --array="0-$SCT_PC_SPECIES" --cpus-per-task=1 --mem=8G --time=01:00:00 --job-name="st_${SCT_PC_NAME}_prep" --output="logs/st_${SCT_PC_NAME}_prep.%A_%a.out" --error="logs/st_${SCT_PC_NAME}_prep.%A_%a.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$SCT_PC_ROOT/code/src:$SCT_PC_ROOT/code/validation/real_data" PYTHONUNBUFFERED=1 OPENBLAS_NUM_THREADS=1
"$SCT_PC_TOOLS/bin/python" "$SCT_PC_ROOT/code/validation/real_data/paper_genome_reanalysis.py" prepare --run "$SCT_PC_RUN" --manifest "$SCT_PC_MANIFEST" --case "$SCT_PC_CASE" --index "$SLURM_ARRAY_TASK_ID" --gffread "$SCT_PC_GFFREAD"
SLURM
)
  SCT_PC_PAIR=$(sbatch --parsable --partition=cpu --array="0-$SCT_PC_PAIRS" --cpus-per-task=4 --mem=8G --time=03:00:00 --dependency="afterok:${SCT_PC_PREP%%;*}" --job-name="st_${SCT_PC_NAME}_pair" --output="logs/st_${SCT_PC_NAME}_pair.%A_%a.out" --error="logs/st_${SCT_PC_NAME}_pair.%A_%a.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PATH="$SCT_PC_TOOLS/bin:$PATH" PYTHONPATH="$SCT_PC_ROOT/code/src:$SCT_PC_ROOT/code/validation/real_data" PYTHONUNBUFFERED=1 OPENBLAS_NUM_THREADS=1
python "$SCT_PC_ROOT/code/validation/real_data/paper_genome_reanalysis.py" pair --run "$SCT_PC_RUN" --manifest "$SCT_PC_MANIFEST" --case "$SCT_PC_CASE" --index "$SLURM_ARRAY_TASK_ID"
SLURM
)
  SCT_PC_FIG=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=16G --time=01:00:00 --dependency="afterok:${SCT_PC_PAIR%%;*}" --job-name="st_${SCT_PC_NAME}_fig" --output="logs/st_${SCT_PC_NAME}_fig.%j.out" --error="logs/st_${SCT_PC_NAME}_fig.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$SCT_PC_ROOT/code/src:$SCT_PC_ROOT/code/validation/real_data" PYTHONUNBUFFERED=1 MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$SCT_PC_TOOLS/bin/python" "$SCT_PC_ROOT/code/validation/real_data/paper_genome_reanalysis.py" collect --run "$SCT_PC_RUN" --manifest "$SCT_PC_MANIFEST" --case "$SCT_PC_CASE"
SCT_PC_ID=$("$SCT_PC_TOOLS/bin/python" -c 'import json,sys; print(json.load(open(sys.argv[1]))["cases"][int(sys.argv[2])]["id"])' "$SCT_PC_MANIFEST" "$SCT_PC_CASE")
eval "$("$SCT_PC_MAMBA" shell hook --shell bash)"
micromamba activate "$SCT_PC_GS_ENV"
"$SCT_PC_TOOLS/bin/python" "$SCT_PC_ROOT/code/validation/real_data/plot_comparison.py" --prepared "$SCT_PC_RUN/prepared/$SCT_PC_ID" --output "$SCT_PC_RUN/figures/$SCT_PC_ID" --seconds 900
printf '{"datasets":[{"id":"%s"}]}\n' "$SCT_PC_ID" > "$SCT_PC_RUN/comparison_manifest.json"
"$SCT_PC_TOOLS/bin/python" "$SCT_PC_ROOT/code/validation/real_data/collect_comparisons.py" --run "$SCT_PC_RUN" --manifest "$SCT_PC_RUN/comparison_manifest.json" --latest-archive "$SCT_PC_REPO/local_results/SynTangle_${SCT_PC_NAME}_comparisons.zip"
SLURM
)
  printf '%s preparation: %s; all pairs concurrently: %s; comparison: %s\n' "$SCT_PC_NAME" "$SCT_PC_PREP" "$SCT_PC_PAIR" "$SCT_PC_FIG"
done
SCT_PC_INVENTORY=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=2G --time=00:30:00 --job-name=st_carex_sources --output=logs/st_carex_sources.%j.out --error=logs/st_carex_sources.%j.err <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$SCT_PC_ROOT/code/src:$SCT_PC_ROOT/code/validation/real_data" PYTHONUNBUFFERED=1 OPENBLAS_NUM_THREADS=1
"$SCT_PC_TOOLS/bin/python" "$SCT_PC_ROOT/code/validation/real_data/paper_genome_reanalysis.py" inventory --run "$SCT_PC_ROOT" --manifest "$SCT_PC_MANIFEST"
SLURM
)
printf 'Carex author-source inventory: %s\nResults: %s\n' "$SCT_PC_INVENTORY" "$SCT_PC_ROOT"
