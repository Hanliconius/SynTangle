"""Import archived directed GENESPACE block links without rerunning orthology.

Each row is one observed pairwise ribbon, NOT a merged orthogroup. Adjacent
species pairs are retained; the manifest fixes the row chain before solving.
"""
import argparse
from collections import Counter, defaultdict
import csv, gzip, hashlib, json, math, re, urllib.request
from pathlib import Path
from syntangle.fixtures import fixture_from_dict


def natural(s):
    return tuple((1,int(x)) if x.isdigit() else (0,x) for x in re.split(r'(\d+)',s))


def signature(row, reverse=False):
    sides=[]
    for i in (2,1) if reverse else (1,2):
        a,b=float(row[f'startBp{i}']),float(row[f'endBp{i}'])
        if not math.isfinite(a+b) or min(a,b)<0 or a==b:
            raise ValueError('Invalid/nonpositive block span')
        sides.append((row[f'chr{i}'],min(a,b),max(a,b)))
    orient=row['orient']
    if orient not in ('+','-'): raise ValueError('Unknown block orientation; inspect source')
    return tuple(sides)+(orient,)


def convert(rows, species, case_id):
    if len(set(species))!=len(species) or len(species)<2: raise ValueError('Distinct ordered species required')
    chromosomes={s:defaultdict(list) for s in species}; edges=[]; audit=[]
    for a,b in zip(species,species[1:]):
        direct=[(i,x) for i,x in enumerate(rows,2) if (x['gen1'],x['gen2'])==(a,b)]
        reverse=[(i,x) for i,x in enumerate(rows,2) if (x['gen1'],x['gen2'])==(b,a)]
        if not direct or not reverse: raise ValueError(f'Missing directional pair {a}/{b}')
        if Counter(signature(x)[:-1] for _,x in direct)!=Counter(signature(x,True)[:-1] for _,x in reverse):
            raise ValueError(f'Reciprocal evidence differs for {a}/{b}; manual inspection required, no silent deduplication')
        reciprocal=defaultdict(Counter)
        for _,x in reverse:reciprocal[signature(x,True)[:-1]][x['orient']]+=1
        direct_annotations=defaultdict(Counter)
        for _,x in direct:direct_annotations[signature(x)[:-1]][x['orient']]+=1
        conflicts=[dict(geometry=key,direct_orientations=dict(value),reciprocal_orientations=dict(reciprocal[key])) for key,value in direct_annotations.items() if value!=reciprocal[key]]
        audit.append(dict(pair=[a,b],retained=len(direct),verified_reciprocal_rows=len(reverse),orientation_conflicts=conflicts))
        for i,x in direct:
            sig=signature(x); hid=f'row_{i}_{a}_{b}'
            for sp,side in zip((a,b),sig[:2]):
                chrom,start,end=side
                chromosomes[sp][chrom].append(dict(occurrence_id=hid+'_'+sp,homology_id=hid,start=start,end=end,strand='?'))
            edges.append(dict(homology_id=hid,source_line=i,source_block_id=x['blkID'],source_orientation=x['orient'],reciprocal_orientation_counts=dict(reciprocal[sig[:-1]]),orientation_conflict=direct_annotations[sig[:-1]]!=reciprocal[sig[:-1]]))
    data=dict(fixture_version=1,id=case_id,title=case_id,purpose='Archived observed block-link chain; chromosome movements only',species=[],orientation_constraints=[])
    for sp in species:
        data['species'].append(dict(id=sp,chromosomes=[dict(id=ch,length=max(x['end'] for x in chromosomes[sp][ch]),display_rank=j+1,blocks=chromosomes[sp][ch]) for j,ch in enumerate(sorted(chromosomes[sp],key=natural))]))
    fixture_from_dict(data)  # schema + coordinates checked before publishing input
    provenance=dict(pair_audit=audit,edges=edges,source_rows=len(rows),retained_links=len(edges),
        source_rows_outside_selected_directed_chain=len(rows)-len(edges),
        chromosome_extent='Maximum observed block endpoint; NOT assembly chromosome lengths',
        display_order='Natural chromosome-name order; NOT recovered published drawing',
        objective='Unweighted crossings between block midpoints of adjacent species; NOT gene-anchor crossings',
        ambiguity='Overlapping, repeated and duplicated-coordinate block links retained as separate observed edges; no orthogroup equivalence inferred',
        orientation='Both directional orientation annotations retained in edge provenance; conflicts explicit; neither imposed as a hard whole-chromosome parity constraint',
        scope='Linked chromosomes only; self-comparisons, reciprocal copies and nonadjacent species comparisons explicitly outside this benchmark')
    return data,provenance


def main():
    p=argparse.ArgumentParser();p.add_argument('--manifest',type=Path,required=True);p.add_argument('--index',type=int,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    m=json.loads(a.manifest.read_text());d=m['datasets'][a.index];out=a.output/d['id'];out.mkdir(parents=True,exist_ok=True)
    target=out/'syntenicBlocks.txt.gz'
    if not target.exists():
        url=f"https://raw.githubusercontent.com/{m['repository']}/{m['commit']}/{d['path']}"
        tmp=target.with_suffix('.partial')
        try:
            with urllib.request.urlopen(url,timeout=120) as response,tmp.open('wb') as f:
                while chunk:=response.read(1024*1024): f.write(chunk)
            tmp.replace(target)
        finally: tmp.unlink(missing_ok=True)
    raw=target.read_bytes();blob=hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()
    if blob!=d['git_blob_sha']: raise ValueError('Downloaded source does not match pinned Git blob')
    with gzip.open(target,'rt') as f: rows=list(csv.DictReader(f,delimiter='\t'))
    data,provenance=convert(rows,d['species'],d['id'])
    provenance.update(source=d,repository=m['repository'],commit=m['commit'],paper=m['paper'],sha256=hashlib.sha256(raw).hexdigest())
    (out/'fixture.json').write_text(json.dumps(data,indent=2)+'\n');(out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    conflicts=sum(len(x['orientation_conflicts']) for x in provenance['pair_audit'])
    print(f"PREPARED {d['id']}: {len(provenance['edges'])} observed links; {sum(len(s['chromosomes']) for s in data['species'])} linked chromosomes; orientation_conflicts={conflicts}",flush=True)

if __name__=='__main__':main()
