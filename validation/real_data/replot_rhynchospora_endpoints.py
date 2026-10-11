"""Replot saved states with endpoint-faithful whole-chromosome reflections."""
import argparse,json,shutil
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from syntangle.fixtures import load_fixture
from syntangle.saved_layout import decode_saved_layout,validate_saved_layout
from syntangle.orientation_space import orientation_basis
from syntangle.layout import score_crossings
from plot_comparison import render,atomic_json
import plot_comparison,visual_summary

def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args()
    source=a.run/'figures/rhynchospora_native_regions';out=a.run/'endpoint_checked';out.mkdir(exist_ok=True)
    fixture=load_fixture(source/'fixture.json');audit=json.loads((source/'comparison_audit.json').read_text())
    basis=orientation_basis(fixture,fixture.chromosome_refs);panels=[];states=[]
    for item in audit['panels']:
        state=decode_saved_layout(fixture,item['state']);validate_saved_layout(fixture,state,basis)
        score=score_crossings(fixture,state).crossings
        if score!=item['crossings']:raise ValueError('Saved score mismatch')
        panels.append((state,item['title'],'Directed endpoints reflected with every whole-chromosome flip'))
        states.append(state)
    if len(states)!=4:raise ValueError('Expected four saved panels')
    native,final=states[1],states[3]
    flips=[dict(species=r.species_id,chromosome=r.chromosome_id,before=native.chromosome_orientation[r],after=final.chromosome_orientation[r]) for r in fixture.chromosome_refs if native.chromosome_orientation[r]!=final.chromosome_orientation[r]]
    evidence=[dict(species=c.ref.species_id,chromosome=c.ref.chromosome_id,length=c.length,blocks=[dict(id=b.occurrence_id,homology=b.homology_id,start=b.start,end=b.end,strand=b.strand) for b in c.blocks]) for c in fixture.chromosomes]
    plot_comparison.riparian=lambda *args,**kw:visual_summary.riparian(*args,**kw,respect_block_strand=True)
    render(fixture,panels,out,audit['provenance'])
    atomic_json(out/'endpoint_audit.json',dict(whole_chromosome_flips_native_to_final=flips,evidence=evidence,internal_moves=0,rule='All panels use one fixed fixture. Endpoints x transform as x -> L-x for negative chromosome orientation.',panel_scores=[x['crossings'] for x in audit['panels']]))
    shutil.copy2(source/'fixture.json',out/'fixture.json')
    shutil.copy2(source/'comparison_audit.json',out/'saved_layouts.json')
    archive=a.run.parent/'SynTangle_rhynchospora_endpoint_checked.zip'
    temporary=archive.with_suffix('.partial')
    with ZipFile(temporary,'w',ZIP_DEFLATED) as z:
        for f in sorted(out.iterdir()):
            if f.is_file():z.write(f,f.name)
    temporary.replace(archive)
    print('WHOLE CHROMOSOME FLIPS',json.dumps(flips),flush=True)
    print('Saved midpoint scores unchanged:',[x['crossings'] for x in audit['panels']],flush=True)
    print('PDF:',out/'comparison.pdf',flush=True);print('ARCHIVE:',archive,flush=True)
if __name__=='__main__':main()
