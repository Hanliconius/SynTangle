#!/usr/bin/env bash
set -euo pipefail
SCT_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_REPO"
SCT_SOURCE=$(realpath "${1:-local_results/rhynchospora_visual_PTiYgr}")
SCT_TE_ID=${2:-}
test -f "$SCT_SOURCE/prepared/rhynchospora_transfer_reanalysis/fixture.json"
test -x "$SCT_SOURCE/venv/bin/python"
SCT_MAMBA=${MAMBA_EXE:-$(type -P micromamba)}
SCT_GS_ENV=${SCT_REAL_GS_ENV:-lep_busco_painter_clean}
SCT_RUN=$(mktemp -d "$SCT_REPO/local_results/rhynchospora_TE_XXXXXX")
mkdir -p "$SCT_RUN/code" logs
cp -a src validation "$SCT_RUN/code/"
printf '%s\n' "$SCT_RUN" > local_results/latest_rhynchospora_TE_run.txt
printf '%s\n' '{"datasets":[{"id":"all_anchors"},{"id":"CDS_TE_screened"}]}' > "$SCT_RUN/comparison_manifest.json"
export SCT_REPO SCT_SOURCE SCT_RUN SCT_MAMBA SCT_GS_ENV SCT_TE_ID
SCT_AUDIT=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=8G --time=01:00:00 --job-name=st_rh_TEaudit --output="$SCT_REPO/logs/st_rh_TEaudit.%j.out" --error="$SCT_REPO/logs/st_rh_TEaudit.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$SCT_RUN/code/src" PYTHONUNBUFFERED=1
cd "$SCT_RUN/code/validation/real_data"
"$SCT_SOURCE/venv/bin/python" -m unittest test_rhynchospora_te_audit > "$SCT_RUN/tests.log" 2>&1
SCT_ARGS=()
if [[ -n "$SCT_TE_ID" ]]; then SCT_ARGS=(--repeat-file-id "$SCT_TE_ID"); fi
"$SCT_SOURCE/venv/bin/python" rhynchospora_te_audit.py --source "$SCT_SOURCE" --run "$SCT_RUN" "${SCT_ARGS[@]}"
SLURM
)
SCT_COMPARE=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=16G --time=02:00:00 --dependency="afterok:${SCT_AUDIT%%;*}" --job-name=st_rh_TEpair --output="$SCT_REPO/logs/st_rh_TEpair.%j.out" --error="$SCT_REPO/logs/st_rh_TEpair.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
if [[ ! -f "$SCT_RUN/AUDIT_COMPLETE" ]]; then
    printf 'No matched repeat source; inspect audit log. No filtered comparison claimed.\n'; exit 1
fi
export PYTHONPATH="$SCT_RUN/code/src" MPLBACKEND=Agg PYTHONUNBUFFERED=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
eval "$("$SCT_MAMBA" shell hook --shell bash)"
micromamba activate "$SCT_GS_ENV"
SCT_PY="$SCT_SOURCE/venv/bin/python"
for SCT_CASE in all_anchors CDS_TE_screened; do
    SCT_GENES="$SCT_SOURCE/prepared/rhynchospora_transfer_reanalysis"
    if [[ "$SCT_CASE" == CDS_TE_screened ]]; then SCT_GENES="$SCT_RUN/filtered_genes"; fi
    "$SCT_PY" "$SCT_RUN/code/validation/real_data/rhynchospora_blocks.py" --source "$SCT_GENES" --output "$SCT_RUN/prepared/$SCT_CASE"
    "$SCT_PY" "$SCT_RUN/code/validation/real_data/plot_rhynchospora_blocks.py" --prepared "$SCT_RUN/prepared/$SCT_CASE" --genes "$SCT_SOURCE/prepared/rhynchospora_transfer_reanalysis" --output "$SCT_RUN/figures/$SCT_CASE" --seconds 300
done
SLURM
)
SCT_ZIP=$(sbatch --parsable --partition=nano --cpus-per-task=1 --mem=4G --time=00:10:00 --dependency="afterany:${SCT_COMPARE%%;*}" --job-name=st_rh_TEzip --output="$SCT_REPO/logs/st_rh_TEzip.%j.out" --error="$SCT_REPO/logs/st_rh_TEzip.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
"$SCT_SOURCE/venv/bin/python" "$SCT_RUN/code/validation/real_data/collect_comparisons.py" --run "$SCT_RUN" --manifest "$SCT_RUN/comparison_manifest.json" --latest-archive "$SCT_REPO/local_results/SynTangle_rhynchospora_TE_comparisons.zip"
"$SCT_SOURCE/venv/bin/python" - <<'PY'
import os,shutil,zipfile
from pathlib import Path
run=Path(os.environ['SCT_RUN']);archive=run/'SynTangle_real_comparisons.zip'
with zipfile.ZipFile(archive,'a',compression=zipfile.ZIP_DEFLATED) as z:
    for name in ('te_anchor_audit.json','repeat_candidates.json','AUDIT_COMPLETE','AUDIT_NEEDS_SOURCE','sources/IXRT5Y.inventory.json'):
        p=run/name
        if p.is_file():z.write(p,name)
latest=Path(os.environ['SCT_REPO'])/'local_results/SynTangle_rhynchospora_TE_comparisons.zip'
tmp=latest.with_suffix('.partial');shutil.copyfile(archive,tmp);tmp.replace(latest)
PY
SLURM
)
printf 'TE audit: %s\nComparison: %s\nCollector: %s\nResults: %s\n' "$SCT_AUDIT" "$SCT_COMPARE" "$SCT_ZIP" "$SCT_RUN"
