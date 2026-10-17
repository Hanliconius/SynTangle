#!/usr/bin/env bash
set -euo pipefail
SCT_SIMPLE_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_SIMPLE_REPO"
SCT_SIMPLE_MODE=${1:-intake}
SCT_SIMPLE_TOOLS=$(realpath local_results/annotated_lep_rjv857/tools)
if [[ "$SCT_SIMPLE_MODE" == intake ]]; then
  SCT_SIMPLE_RUN=$(mktemp -d "$SCT_SIMPLE_REPO/local_results/simple_papers_XXXXXX")
  mkdir -p "$SCT_SIMPLE_RUN/code" logs
  cp -a src validation "$SCT_SIMPLE_RUN/code/"
  printf '%s\n' "$SCT_SIMPLE_RUN" > local_results/latest_simple_papers_run.txt
elif [[ "$SCT_SIMPLE_MODE" == compare ]]; then
  SCT_SIMPLE_RUN=$(realpath "${2:-$(cat local_results/latest_simple_papers_run.txt)}")
  export PYTHONPATH="$SCT_SIMPLE_RUN/code/src:$SCT_SIMPLE_RUN/code/validation/real_data"
  "$SCT_SIMPLE_TOOLS/bin/python" - "$SCT_SIMPLE_RUN" <<'PY'
import sys,json
from pathlib import Path
from simple_papers import gate
r=Path(sys.argv[1]);m=json.loads((r/'code/validation/real_data/simple_papers_manifest.json').read_text())
for c in m['cases']:gate(r/c['id']/'prepared')
print('PASS four source/layout gates')
PY
else
  printf 'Usage: bash validation/real_data/submit_simple_papers.sh intake | compare [RUN]\n' >&2;exit 2
fi
export SCT_SIMPLE_REPO SCT_SIMPLE_RUN SCT_SIMPLE_TOOLS SCT_SIMPLE_MODE
SCT_SIMPLE_ARRAY=0-3
SCT_SIMPLE_LIMIT=01:00:00
if [[ "$SCT_SIMPLE_MODE" == compare ]]; then SCT_SIMPLE_ARRAY=0-11; SCT_SIMPLE_LIMIT=00:30:00; fi
SCT_SIMPLE_JOB=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=8G --time="$SCT_SIMPLE_LIMIT" --array="$SCT_SIMPLE_ARRAY" --job-name="st_simple_$SCT_SIMPLE_MODE" --output="logs/st_simple_$SCT_SIMPLE_MODE.%A_%a.out" --error="logs/st_simple_$SCT_SIMPLE_MODE.%A_%a.err" <<'SLURM'
#!/bin/bash
set -euo pipefail
export PATH="$SCT_SIMPLE_TOOLS/bin:$PATH" PYTHONPATH="$SCT_SIMPLE_RUN/code/src:$SCT_SIMPLE_RUN/code/validation/real_data" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg PYTHONUNBUFFERED=1
if [[ "$SCT_SIMPLE_MODE" == compare ]]; then
  eval "$(micromamba shell hook --shell bash)"
  micromamba activate "${SCT_REAL_GS_ENV:-lep_busco_painter_clean}"
  export PATH="$SCT_SIMPLE_TOOLS/bin:$PATH"
fi
/usr/bin/time -o "$SCT_SIMPLE_RUN/${SCT_SIMPLE_MODE}_${SLURM_ARRAY_TASK_ID}.resources.tsv" -f 'wall_seconds=%e\nmax_rss_kb=%M' "$SCT_SIMPLE_TOOLS/bin/python" "$SCT_SIMPLE_RUN/code/validation/real_data/simple_papers.py" "$SCT_SIMPLE_MODE" --run "$SCT_SIMPLE_RUN" --index "$SLURM_ARRAY_TASK_ID"
SLURM
)
SCT_SIMPLE_COLLECT=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=2G --time=00:15:00 --dependency="afterany:${SCT_SIMPLE_JOB%%;*}" --job-name=st_simple_collect --output=logs/st_simple_collect.%j.out --error=logs/st_simple_collect.%j.err <<'SLURM'
#!/bin/bash
set -euo pipefail
export PYTHONPATH="$SCT_SIMPLE_RUN/code/src:$SCT_SIMPLE_RUN/code/validation/real_data" OPENBLAS_NUM_THREADS=1
"$SCT_SIMPLE_TOOLS/bin/python" "$SCT_SIMPLE_RUN/code/validation/real_data/simple_papers.py" collect --run "$SCT_SIMPLE_RUN"
SLURM
)
printf 'Concurrent %s array: %s\nCollector: %s\nResults: %s\n' "$SCT_SIMPLE_MODE" "$SCT_SIMPLE_JOB" "$SCT_SIMPLE_COLLECT" "$SCT_SIMPLE_RUN"
