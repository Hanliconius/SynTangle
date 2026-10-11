import json,os,shutil
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from pypdf import PdfWriter
from preservation_audit import require_audit
run=Path(os.environ['SCT_RUN']); repo=Path(os.environ['SCT_REPO'])
out=run/'figures/rhynchospora_native_regions'
require_audit(out)
status=dict(native_complete=(run/'NATIVE_COMPLETE').exists(),comparison_complete=(out/'COMPLETE').exists(),orthofinder_complete=(run/'ORTHOFINDER_COMPLETE').exists(),scope='Three-genome protein-based GENESPACE reanalysis; transferred breviuscula annotation; not exact original twenty-genome discovery')
(run/'status.json').write_text(json.dumps(status,indent=2)+'\n')
with PdfWriter() as writer:
    for path in [run/'native_regions.pdf',run/'native_blocks.pdf',out/'comparison.pdf']:
        if path.exists():writer.append(str(path))
    if writer.pages:writer.write(str(run/'SynTangle_native_regions.pdf'))
archive=run/'SynTangle_rhynchospora_native_repaired.zip'
with ZipFile(archive,'w',ZIP_DEFLATED) as z:
    paths=list(out.rglob('*'))+[run/name for name in ('native_regions.pdf','native_blocks.pdf','native_regions.rds','native_blocks.rds','genespace_parameters.rds','native_regions.tsv','native_blocks.tsv','native_regions_chromosomes.tsv','native_blocks_chromosomes.tsv','input_audit.json','orthofinder_repair_audit.json','repair_inputs.json','repair_versions.json','orthofinder_resume.log','orthogroup_mode.txt','status.json','SynTangle_native_regions.pdf')]
    for path in paths:
        if path.is_file():z.write(path,str(path.relative_to(run)))
latest=repo/'local_results/SynTangle_rhynchospora_native_repaired.zip'
tmp=latest.with_suffix('.partial');shutil.copyfile(archive,tmp);tmp.replace(latest)
print(json.dumps(status,indent=2));print('ARCHIVE:',latest)
