#!/usr/bin/env bash
set -euo pipefail
SCT_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$SCT_REPO"
SCT_SOURCE=$(realpath "${1:-local_results/rhynchospora_visual_PTiYgr}")
test -f "$SCT_SOURCE/source_manifest.json"
test -x "$SCT_SOURCE/venv/bin/python"
SCT_MAMBA=${MAMBA_EXE:-$(type -P micromamba)}
SCT_GS_ENV=${SCT_REAL_GS_ENV:-lep_busco_painter_clean}
SCT_RUN=$(mktemp -d "$SCT_REPO/local_results/rhynchospora_native_XXXXXX")
mkdir -p "$SCT_RUN/code" logs
cp -a src validation "$SCT_RUN/code/"
printf '%s\n' "$SCT_RUN" > local_results/latest_rhynchospora_native_run.txt
export SCT_REPO SCT_SOURCE SCT_RUN SCT_MAMBA SCT_GS_ENV
SCT_SETUP=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=8G --time=01:00:00 --job-name=st_rh_native_setup --output="$SCT_REPO/logs/st_rh_native_setup.%j.out" --error="$SCT_REPO/logs/st_rh_native_setup.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
"$SCT_MAMBA" create -y -p "$SCT_RUN/tools" -c conda-forge -c bioconda python=3.11 orthofinder=2.5.5 diamond mcscanx gffread
export PATH="$SCT_RUN/tools/bin:$PATH"
test -x "$SCT_RUN/tools/bin/MCScanX_h"
"$SCT_RUN/tools/bin/python" "$SCT_RUN/code/validation/real_data/rhynchospora_native_regions.py" prepare --run "$SCT_RUN" --source "$SCT_SOURCE"
SLURM
)
SCT_JOB=$(sbatch --parsable --partition=cpu --cpus-per-task=4 --mem=32G --time=08:00:00 --dependency="afterok:${SCT_SETUP%%;*}" --job-name=st_rh_native --output="$SCT_REPO/logs/st_rh_native.%j.out" --error="$SCT_REPO/logs/st_rh_native.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="$SCT_RUN/code/src" MPLBACKEND=Agg PYTHONUNBUFFERED=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
eval "$("$SCT_MAMBA" shell hook --shell bash)"
micromamba activate "$SCT_GS_ENV"
export PATH="$SCT_RUN/tools/bin:$PATH"
Rscript "$SCT_RUN/code/validation/real_data/rhynchospora_native_regions.R" "$SCT_RUN" "$SCT_RUN/tools/bin"
"$SCT_SOURCE/venv/bin/python" "$SCT_RUN/code/validation/real_data/rhynchospora_native_regions.py" compare --run "$SCT_RUN" --seconds 300
SLURM
)
SCT_COLLECT=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=4G --time=00:15:00 --dependency="afterany:${SCT_JOB%%;*}" --job-name=st_rh_native_zip --output="$SCT_REPO/logs/st_rh_native_zip.%j.out" --error="$SCT_REPO/logs/st_rh_native_zip.%j.err" <<'SLURM'
#!/usr/bin/env bash
set -euo pipefail
"$SCT_SOURCE/venv/bin/python" - <<'PY'
import json,os,shutil
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from pypdf import PdfWriter
run=Path(os.environ['SCT_RUN']); repo=Path(os.environ['SCT_REPO'])
out=run/'figures/rhynchospora_native_regions'
status=dict(native_complete=(run/'NATIVE_COMPLETE').exists(),comparison_complete=(out/'COMPLETE').exists(),scope='Three-genome protein-based GENESPACE reanalysis; transferred breviuscula annotation; not exact original twenty-genome discovery')
(run/'status.json').write_text(json.dumps(status,indent=2)+'\n')
with PdfWriter() as writer:
    for path in [run/'native_regions.pdf',run/'native_blocks.pdf',out/'comparison.pdf']:
        if path.exists():writer.append(str(path))
    if writer.pages:writer.write(str(run/'SynTangle_native_regions.pdf'))
archive=run/'SynTangle_rhynchospora_native_regions.zip'
with ZipFile(archive,'w',ZIP_DEFLATED) as z:
    paths=list(out.rglob('*'))+[run/name for name in ('native_regions.pdf','native_blocks.pdf','native_regions.rds','native_blocks.rds','genespace_parameters.rds','native_regions.tsv','native_blocks.tsv','native_regions_chromosomes.tsv','native_blocks_chromosomes.tsv','input_audit.json','status.json','SynTangle_native_regions.pdf')]
    for path in paths:
        if path.is_file():z.write(path,str(path.relative_to(run)))
latest=repo/'local_results/SynTangle_rhynchospora_native_regions.zip'
tmp=latest.with_suffix('.partial');shutil.copyfile(archive,tmp);tmp.replace(latest)
print(json.dumps(status,indent=2));print('ARCHIVE:',latest)
PY
SLURM
)
printf 'Setup/translation: %s\nNative discovery/comparison: %s\nCollector: %s\nResults: %s\n' "$SCT_SETUP" "$SCT_JOB" "$SCT_COLLECT" "$SCT_RUN"
