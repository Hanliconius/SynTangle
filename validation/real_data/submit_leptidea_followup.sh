#!/usr/bin/env bash
set -euo pipefail
SCT_LF_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_LF_REPO"
SCT_LF_SOURCE=$(realpath "${1:-local_results/leptidea_paper_SaryLn}")
SCT_LF_PY=$(realpath local_results/annotated_lep_rjv857/tools/bin/python)
SCT_LF_PRIOR=${2:-$(cat local_results/latest_leptidea_extended_run.txt 2>/dev/null || true)}
if [[ -z "$SCT_LF_PRIOR" ]]; then
  SCT_LF_PRIOR=$("$SCT_LF_PY" - "$SCT_LF_REPO/local_results" <<'PYFIND'
import json,sys
from pathlib import Path
found=[]
for root in Path(sys.argv[1]).glob('leptidea_extended_*'):
    summary=root/'figures/leptidea_paper_assemblies/summary.json'
    if summary.exists():
        s=json.loads(summary.read_text())
        if s['upper_bound']<=3468 and s['lower_bound']>=359:found.append((summary.stat().st_mtime,root))
if not found:raise SystemExit('No completed extended run with validated U<=3468 and L>=359 found')
print(max(found,key=lambda x:x[0])[1])
PYFIND
)
fi
test -n "$SCT_LF_PRIOR"
SCT_LF_PRIOR=$(realpath "$SCT_LF_PRIOR")
"$SCT_LF_PY" -c 'import json,sys; from pathlib import Path; p=Path(sys.argv[1]); s=json.loads((p/"figures/leptidea_paper_assemblies/summary.json").read_text()); print("Starting from saved U",s["upper_bound"],"L",s["lower_bound"])' "$SCT_LF_PRIOR"
SCT_LF_RUN=$(mktemp -d "$SCT_LF_REPO/local_results/leptidea_followup_XXXXXX")
mkdir -p "$SCT_LF_RUN/code" logs
cp -a src validation "$SCT_LF_RUN/code/"
export SCT_LF_REPO SCT_LF_SOURCE SCT_LF_PY SCT_LF_RUN SCT_LF_PRIOR
SCT_LF_JOB=$(sbatch --parsable --partition=cpu --array=0-1 --cpus-per-task=1 --mem=16G --time=01:15:00 --job-name=st_lep_follow --output=logs/st_lep_follow.%A_%a.out --error=logs/st_lep_follow.%A_%a.err <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$SCT_LF_RUN/code/src:$SCT_LF_RUN/code/validation/real_data" MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$SCT_LF_PY" "$SCT_LF_RUN/code/validation/real_data/leptidea_followup.py" run --source "$SCT_LF_SOURCE" --prior "$SCT_LF_PRIOR" --run "$SCT_LF_RUN" --index "$SLURM_ARRAY_TASK_ID"
SLURM
)
SCT_LF_COLLECT=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=4G --time=00:10:00 --dependency="afterany:${SCT_LF_JOB%%;*}" --job-name=st_lep_follow_zip --output=logs/st_lep_follow_zip.%j.out --error=logs/st_lep_follow_zip.%j.err <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$SCT_LF_RUN/code/src:$SCT_LF_RUN/code/validation/real_data" MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$SCT_LF_PY" "$SCT_LF_RUN/code/validation/real_data/leptidea_followup.py" collect --source "$SCT_LF_SOURCE" --prior "$SCT_LF_PRIOR" --run "$SCT_LF_RUN"
SLURM
)
printf 'Adaptive Leptidea follow-up, both concurrent: %s\nCollector: %s\nResults: %s\n' "$SCT_LF_JOB" "$SCT_LF_COLLECT" "$SCT_LF_RUN"
