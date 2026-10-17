#!/usr/bin/env bash
set -euo pipefail
SCT_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_REPO"
SCT_RUN=$(realpath "${1:-local_results/codfish_carex_kdt7Zz/codfish}")
SCT_TOOLS=$(realpath local_results/annotated_lep_rjv857/tools)
SCT_WORK="$SCT_RUN/chunked_blast_v1"
mkdir -p "$SCT_WORK/code" logs
# Refuse concurrent writes to pair directories from the obsolete array.
SCT_OLD=$(scontrol show job -o 138436841 2>/dev/null || true)
if [[ "$SCT_OLD" == *JobName=st_codfish_pair* ]]; then scancel 138436841; fi
cp -a validation/real_data/. "$SCT_WORK/code/"
export SCT_RUN SCT_TOOLS SCT_WORK
SCT_SETUP=$(sbatch --parsable -p cpu -c 1 --mem=4G -t 00:30:00 -J st_cod_chunk_setup -o logs/st_cod_chunk_setup.%j.out -e logs/st_cod_chunk_setup.%j.err <<'SLURM'
#!/bin/bash
set -euo pipefail
export PATH="$SCT_TOOLS/bin:$PATH" PYTHONPATH="$SCT_WORK/code"
python "$SCT_WORK/code/codfish_chunked.py" setup --run "$SCT_RUN" --work "$SCT_WORK" --manifest "$SCT_WORK/code/codfish_carex_manifest.json"
SLURM
)
SCT_SEARCH=$(sbatch --parsable -p cpu -c 2 --mem=4G -t 04:00:00 --array=0-119 --dependency="afterok:${SCT_SETUP%%;*}" -J st_cod_chunk -o logs/st_cod_chunk.%A_%a.out -e logs/st_cod_chunk.%A_%a.err <<'SLURM'
#!/bin/bash
set -euo pipefail
export PATH="$SCT_TOOLS/bin:$PATH" PYTHONPATH="$SCT_WORK/code" OPENBLAS_NUM_THREADS=1
python "$SCT_WORK/code/codfish_chunked.py" search --run "$SCT_RUN" --work "$SCT_WORK" --index "$SLURM_ARRAY_TASK_ID"
SLURM
)
SCT_BLOCK=$(sbatch --parsable -p cpu -c 1 --mem=8G -t 01:00:00 --array=0-14 --dependency="afterok:${SCT_SEARCH%%;*}" -J st_cod_blocks -o logs/st_cod_blocks.%A_%a.out -e logs/st_cod_blocks.%A_%a.err <<'SLURM'
#!/bin/bash
set -euo pipefail
export PATH="$SCT_TOOLS/bin:$PATH" PYTHONPATH="$SCT_WORK/code" OPENBLAS_NUM_THREADS=1
python "$SCT_WORK/code/codfish_chunked.py" blocks --run "$SCT_RUN" --work "$SCT_WORK" --index "$SLURM_ARRAY_TASK_ID"
SLURM
)
scontrol update JobId=138436842 Dependency="afterok:${SCT_BLOCK%%;*}"
printf 'Setup: %s\n120 concurrent BLAST chunks: %s\n15 block tasks: %s\nExisting comparison: 138436842\n' "$SCT_SETUP" "$SCT_SEARCH" "$SCT_BLOCK"
