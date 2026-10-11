"""Fail-closed checks for fixed evidence and directed rendered ribbon endpoints."""
import hashlib,json,math
from pathlib import Path
from syntangle.orientation_space import orientation_basis
from syntangle.saved_layout import validate_saved_layout
from syntangle.visualize import _unambiguous_links

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def evidence_signature(fixture):
    rows=[(c.ref.species_id,c.ref.chromosome_id,c.length,[(b.occurrence_id,b.homology_id,b.start,b.end,b.strand) for b in c.blocks]) for c in fixture.chromosomes]
    return hashlib.sha256(json.dumps(rows,separators=(',',':')).encode()).hexdigest()
def check_panel(fixture,state,axes,geometry):
    validate_saved_layout(fixture,state,orientation_basis(fixture,fixture.chromosome_refs))
    links=[link for a,b in zip(fixture.species_ids,fixture.species_ids[1:]) for link in _unambiguous_links(fixture,a,b)]
    from matplotlib.patches import PathPatch
    paths=[p for p in axes.patches if isinstance(p,PathPatch)]
    if len(paths)!=len(links):raise ValueError('Ribbon identity/count differs from fixed evidence')
    records=[]
    for patch,(hid,(c1,b1),(c2,b2)) in zip(paths,links):
        endpoints=[]
        for c,b in ((c1,b1),(c2,b2)):
            if not 0<=b.start<b.end<=c.length:raise ValueError('Invalid immutable block coordinates')
            x,width,y=geometry[c.ref]
            # Preserve biological endpoint identity, then reflect both endpoints.
            start,end=(b.end,b.start) if b.strand=='-' else (b.start,b.end)
            coords=[x+width*(v/c.length if state.chromosome_orientation[c.ref]==1 else 1-v/c.length) for v in (start,end)]
            endpoints.append((*coords,y))
        l1,r1,y1=endpoints[0];l2,r2,y2=endpoints[1];mid=(y1+y2)/2
        expected=[(l1,y1+.065),(l1,mid),(l2,mid),(l2,y2-.065),(r2,y2-.065),(r2,mid),(r1,mid),(r1,y1+.065),(l1,y1+.065)]
        actual=patch.get_path().vertices
        if len(actual)!=len(expected) or any(not math.isfinite(float(v)) or abs(float(v)-e)>1e-10 for row,exp in zip(actual,expected) for v,e in zip(row,exp)):
            raise ValueError('Illegal directed ribbon endpoint transform: '+hid)
        records.append(dict(homology_id=hid,source_occurrence=b1.occurrence_id,target_occurrence=b2.occurrence_id,source_chromosome=c1.ref.label,target_chromosome=c2.ref.label,directed_endpoints=endpoints,source_strand=b1.strand,target_strand=b2.strand))
    return dict(status='PASS',links=len(records),links_sha256=hashlib.sha256(json.dumps(records,sort_keys=True).encode()).hexdigest(),state=state.to_dict(),directed_links=records)

def require_audit(output,min_panels=4):
    output=Path(output);data=json.loads((output/'preservation_audit.json').read_text())
    if data.get('status')!='PASS' or len(data['panels'])<min_panels or any(p['status']!='PASS' for p in data['panels']):raise ValueError('Missing or failed preservation audit')
    if data['pdf_sha256']!=digest(output/'comparison.pdf'):raise ValueError('PDF changed after preservation audit')
    return data
