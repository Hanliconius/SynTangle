"""Check saved chromosome-only states and replot native block directions."""
import argparse,json,shutil
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from syntangle.fixtures import load_fixture
from syntangle.saved_layout import decode_saved_layout,validate_saved_layout
from syntangle.orientation_space import orientation_basis
from syntangle.layout import score_crossings
from syntangle.visualize import _unambiguous_links
from plot_comparison import render,atomic_json

def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args()
    case='leptidea_paper_assemblies';source=a.run/'figures'/case;prepared=a.run/'prepared'/case
    out=a.run/'endpoint_checked';out.mkdir(exist_ok=True)
    fixture=load_fixture(prepared/'fixture.json');audit=json.loads((source/'comparison_audit.json').read_text());basis=orientation_basis(fixture,fixture.chromosome_refs)
    panels=[];states=[]
    for item in audit['panels']:
        state=decode_saved_layout(fixture,item['state']);validate_saved_layout(fixture,state,basis)
        score=score_crossings(fixture,state).crossings
        if score!=item['crossings']:raise ValueError('Saved score differs')
        # Independent direct O(n^2) count on fixed block midpoints.
        from itertools import combinations
        direct=0
        for a1,a2 in zip(fixture.species_ids,fixture.species_ids[1:]):
            points=[]
            for h,(c1,b1),(c2,b2) in _unambiguous_links(fixture,a1,a2):
                coords=[]
                for c,b in ((c1,b1),(c2,b2)):
                    x=(b.start+b.end)/(2*c.length)
                    if state.chromosome_orientation[c.ref]<0:x=1-x
                    coords.append((state.chromosome_order[c.ref.species_id].index(c.ref),x))
                points.append(coords)
            direct+=sum((u[0]<v[0] and u[1]>v[1]) or (u[0]>v[0] and u[1]<v[1]) for u,v in combinations(points,2))
        if direct!=score:raise ValueError('Independent midpoint count differs')
        states.append(state);panels.append((state,item['title'],'Fixed blocks; directed endpoints preserved under whole-chromosome reflection'))
    if len(states)!=4:raise ValueError('Expected four saved panels')
    provenance=audit['provenance'];provenance['respect_block_strand']=True
    render(fixture,panels,out,provenance)
    changes=[]
    for i,j in ((1,2),(2,3)):
        flips=[r.label for r in fixture.chromosome_refs if states[i].chromosome_orientation[r]!=states[j].chromosome_orientation[r]]
        order=[sp for sp in fixture.species_ids if states[i].chromosome_order[sp]!=states[j].chromosome_order[sp]]
        changes.append(dict(from_panel=i,to_panel=j,whole_chromosome_flips=flips,changed_row_orders=order))
        print(f'PANELS {i}->{j}: whole_chromosome_flips={len(flips)} rows_reordered={len(order)}',flush=True)
    atomic_json(out/'endpoint_audit.json',dict(changes=changes,panel_scores=[p['crossings'] for p in audit['panels']],rule='One immutable fixture; only chromosome permutations and whole-chromosome reflections. Source block strands are now drawn.',scoring='Independent pairwise midpoint counts match; does not certify global optimality'))
    shutil.copy2(prepared/'fixture.json',out/'fixture.json');shutil.copy2(source/'comparison_audit.json',out/'saved_layouts.json');shutil.copy2(source/'result.json',out/'result.json')
    archive=a.run.parent/'SynTangle_Leptidea_endpoint_checked.zip';temporary=archive.with_suffix('.partial')
    with ZipFile(temporary,'w',ZIP_DEFLATED) as z:
        for f in sorted(out.iterdir()):
            if f.is_file():z.write(f,f.name)
    temporary.replace(archive)
    print('PDF:',out/'comparison.pdf',flush=True);print('ARCHIVE:',archive,flush=True)
if __name__=='__main__':main()
