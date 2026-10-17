"""Auditable monotone-run compression; not the paper's original block caller.

Every adjacent gene link belongs to exactly one run, including singletons.
Block-midpoint optimization is a different objective from gene crossings.
"""
import argparse
from collections import defaultdict
import copy
import hashlib
import json
from pathlib import Path


def build(data, max_rank_gap=20, max_bp_gap=2000000):
    output=copy.deepcopy(data)
    output['id']='rhynchospora_collinear_runs'
    tables={s['id']:{b['homology_id']:(c,b) for c in s['chromosomes'] for b in c['blocks']} for s in data['species']}
    ranks={}
    for s in data['species']:
        for c in s['chromosomes']:
            for rank,b in enumerate(sorted(c['blocks'],key=lambda b:(b['start']+b['end'],b['homology_id']))):
                ranks[s['id'],b['homology_id']]=rank
    targets={(s['id'],c['id']):c for s in output['species'] for c in s['chromosomes']}
    for c in targets.values():c['blocks']=[]
    audit=[]
    for left,right in zip(data['species'],data['species'][1:]):
        a,z=left['id'],right['id']; groups=defaultdict(list)
        shared=set(tables[a])&set(tables[z])
        for hid in sorted(shared):
            ca,ba=tables[a][hid];cz,bz=tables[z][hid]
            groups[ca['id'],cz['id']].append(hid)
        for (ac,zc),ids in sorted(groups.items()):
            ids.sort(key=lambda h:(ranks[a,h],h))
            runs=[];run=[];direction=0
            for hid in ids:
                if run:
                    prev=run[-1];delta=ranks[z,hid]-ranks[z,prev]
                    sign=1 if delta>0 else -1
                    close=(0<ranks[a,hid]-ranks[a,prev]<=max_rank_gap and 0<abs(delta)<=max_rank_gap)
                    for sp in (a,z):
                        b=tables[sp][hid][1];p=tables[sp][prev][1]
                        close=close and abs((b['start']+b['end'])-(p['start']+p['end']))<=2*max_bp_gap
                    if not close or (direction and sign!=direction):
                        runs.append(run);run=[];direction=0
                    else:direction=sign
                run.append(hid)
            if run:runs.append(run)
            for members in runs:
                bid=f'run_{len(audit):06d}'
                record=dict(block_id=bid,pair=[a,z],chromosomes=[ac,zc],members=members)
                audit.append(record)
                for sp,chrom in ((a,ac),(z,zc)):
                    bs=[tables[sp][h][1] for h in members]
                    targets[sp,chrom]['blocks'].append(dict(occurrence_id=sp+'_'+bid,homology_id=bid,
                        start=min(b['start'] for b in bs),end=max(b['end'] for b in bs),strand='+'))
        observed=[h for r in audit if r['pair']==[a,z] for h in r['members']]
        if len(observed)!=len(shared) or set(observed)!=shared:raise AssertionError('Anchor partition failed')
    return output,audit


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--rank-gap',type=int,default=20);p.add_argument('--bp-gap',type=int,default=2000000)
    a=p.parse_args()
    if a.rank_gap<1 or a.bp_gap<1:p.error('Gap limits must be positive')
    source=a.source/'fixture.json';data=json.loads(source.read_text())
    blocks,membership=build(data,a.rank_gap,a.bp_gap)
    from syntangle.fixtures import fixture_from_dict
    fixture_from_dict(blocks)
    provenance=json.loads((a.source/'provenance.json').read_text())
    provenance.update(source_kind='Reconstructed monotone gene runs; NOT published synteny blocks',
        retained_links=len(membership),crossing_unit='run',anchor_description='reconstructed monotone runs',
        input_panel_title='Candidate paper order: reconstructed runs',
        objective='Unweighted run-midpoint crossings; differs from full gene-anchor objective',
        scope='All adjacent gene anchors represented once in run membership, including singletons; compression is not proven lossless',
        gene_fixture_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        pair_audit=[dict(pair=[x['id'],y['id']],retained=sum(r['pair']==[x['id'],y['id']] for r in membership)) for x,y in zip(data['species'],data['species'][1:])])
    # A-B runs receive the donor chromosome only for unanimous donor membership.
    donor={b['homology_id']:c['id'] for c in data['species'][0]['chromosomes'] for b in c['blocks']}
    provenance['homology_colour_reference']={r['block_id']:next(iter(values)) for r in membership if len(values:={donor[h] for h in r['members']})==1}
    a.output.mkdir(parents=True,exist_ok=True)
    for name,value in [('fixture.json',blocks),('provenance.json',provenance),('block_membership.json',dict(rank_gap=a.rank_gap,bp_gap=a.bp_gap,blocks=membership,total_gene_links=sum(len(r['members']) for r in membership),singletons=sum(len(r['members'])==1 for r in membership)))]:
        (a.output/name).write_text(json.dumps(value,indent=2)+'\n')
    print(f'BLOCKS={len(membership)} GENE_LINKS={sum(len(r["members"]) for r in membership)} SINGLETONS={sum(len(r["members"])==1 for r in membership)}',flush=True)


if __name__=='__main__':main()
