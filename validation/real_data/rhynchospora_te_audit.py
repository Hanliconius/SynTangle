"""Test donor coding-sequence overlap with explicitly classified TE intervals.

Repeat overlap is a screening flag, not proof that a protein is TE-derived.
Never infer that ordinary genes are TEs from intronic repeat overlap alone.
"""
import argparse
from bisect import bisect_left
from collections import defaultdict
import copy
import gzip
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import unquote
from prepare_published_inputs import download
from prepare_rhynchospora_inputs import verify


def opener(path):
    return gzip.open(path,'rt') if path.suffix=='.gz' else path.open()


def merge(intervals):
    out=[]
    for start,end in sorted(intervals):
        if start<0 or end<=start:raise ValueError('Invalid interval')
        if out and start<=out[-1][1]:out[-1]=(out[-1][0],max(out[-1][1],end))
        else:out.append((start,end))
    return out


def overlap(cds,repeats):
    return sum(max(0,min(b,d)-max(a,c)) for a,b in merge(cds) for c,d in merge(repeats))


def attrs(text):
    return {k:unquote(v) for k,v in (x.split('=',1) for x in text.split(';') if '=' in x)}


def coding_genes(path):
    parents={}; genes={}; coding=defaultdict(list)
    with opener(path) as f:
        for line in f:
            if line.startswith('##FASTA'):break
            if not line.strip() or line.startswith('#'):continue
            r=line.rstrip().split('\t')
            if len(r)!=9:raise ValueError('Malformed donor GFF')
            a=attrs(r[8]); ps=a.get('Parent','').split(',') if a.get('Parent') else []
            if 'ID' in a:
                if a['ID'] in parents and parents[a['ID']]!=ps:raise ValueError('Conflicting GFF parents')
                parents[a['ID']]=ps
            if r[2]=='gene':genes[a['ID']]=(r[0],int(r[3])-1,int(r[4]))
            if r[2]=='CDS':
                if not ps:raise ValueError('CDS without parent')
                for parent in ps:coding[parent].append((r[0],int(r[3])-1,int(r[4])))
    def gene_for(node,seen):
        if node in seen:raise ValueError('GFF parent cycle')
        if node in genes:return {node}
        return set().union(*(gene_for(p,seen|{node}) for p in parents.get(node,[])))
    candidates=defaultdict(list)
    for transcript,rows in coding.items():
        gs=gene_for(transcript,set())
        if len(gs)!=1:raise ValueError('CDS has missing or ambiguous gene ancestor')
        gene=next(iter(gs));chrom,lo,hi=genes[gene]
        if any(c!=chrom or a<lo or b>hi for c,a,b in rows):raise ValueError('CDS outside gene')
        spans=merge([(a,b) for c,a,b in rows])
        candidates[gene].append((sum(b-a for a,b in spans),transcript,spans))
    return {g:dict(chromosome=genes[g][0],transcript=sorted(v,key=lambda x:(-x[0],x[1]))[0][1],
                   cds=sorted(v,key=lambda x:(-x[0],x[1]))[0][2]) for g,v in candidates.items()}


def read_te(path,lengths):
    # Positive classification only: simple repeats/satellites are never TE flags.
    positive=re.compile(r'transposable[_ -]?element|retrotransposon|\btransposon\b|(?:LTR|LINE|SINE|DNA|RC|MITE)[/_:]|\bhelitron\b',re.I)
    negative=re.compile(r'simple[_ -]?repeat|low[_ -]?complexity|satellite|\btyba\b',re.I)
    intervals=defaultdict(list);counts=defaultdict(int)
    isbed=path.name.lower().removesuffix('.gz').endswith('.bed')
    with opener(path) as f:
        for line in f:
            if line.startswith('##FASTA'):break
            if not line.strip() or line.startswith(('#','track ','browser ')):continue
            r=line.rstrip().split('\t')
            if isbed:
                if len(r)<4:raise ValueError('TE BED needs a classification column')
                chrom,start,end=r[0],int(r[1]),int(r[2]);label=' '.join(r[3:])
            else:
                if len(r)!=9:raise ValueError('Expected GFF3 repeat annotation')
                chrom,start,end=r[0],int(r[3])-1,int(r[4]);label=r[2]+' '+unquote(r[8])
            if chrom not in lengths or not 0<=start<end<=lengths[chrom]:raise ValueError('Repeat coordinates do not match exact donor chromosomes')
            counts['rows']+=1
            if negative.search(label):counts['non_TE_repeat_rows']+=1;continue
            if not positive.search(label):counts['unclassified_rows']+=1;continue
            intervals[chrom].append((start,end));counts['TE_classified_rows']+=1
    if not intervals:raise ValueError('No explicitly TE-classified intervals; no empty negative screen accepted')
    return {c:merge(v) for c,v in intervals.items()},dict(counts)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--run',type=Path,required=True)
    p.add_argument('--repeat-file-id',type=int);p.add_argument('--fraction',type=float,default=.5)
    a=p.parse_args()
    if not 0<a.fraction<=1:p.error('Fraction must be in (0,1]')
    a.run.mkdir(parents=True,exist_ok=True);src=a.run/'sources';src.mkdir(exist_ok=True)
    inventory=src/'IXRT5Y.inventory.json'
    download('https://edmond.mpg.de/api/datasets/:persistentId/?persistentId=doi:10.17617/3.IXRT5Y',inventory)
    raw=json.loads(inventory.read_text())
    if raw.get('status')!='OK':raise ValueError('Deposit inventory failed')
    records=[]
    for entry in raw['data']['latestVersion']['files']:
        d=entry['dataFile'];records.append(dict(id=d['id'],filename=d['filename'],bytes=d['filesize'],restricted=entry.get('restricted',False),checksum=d.get('checksum')))
    candidates=sorted([r for r in records if re.search(r'repeat|dante|edta|transpos|(?:^|[._-])TE(?:[._-]|$)',r['filename'],re.I)],key=lambda r:(not r['filename'].startswith('Rhync_tenuis_REC.hap1.chr.'),r['filename']))
    (a.run/'repeat_candidates.json').write_text(json.dumps(candidates,indent=2)+'\n')
    for r in candidates[:25]:print(f'REPEAT_CANDIDATE id={r["id"]} {r["filename"]}',flush=True)
    if len(candidates)>25:print(f'{len(candidates)-25} additional candidates saved in repeat_candidates.json',flush=True)
    exact=[r for r in candidates if r['filename'].startswith('Rhync_tenuis_REC.hap1.chr.') and re.search(r'\.(gff3?|bed)(\.gz)?$',r['filename'],re.I)]
    selected=[r for r in records if r['id']==a.repeat_file_id] if a.repeat_file_id else exact
    if len(selected)!=1:
        (a.run/'AUDIT_NEEDS_SOURCE').write_text('Select the exact donor TE interval annotation after inspecting repeat_candidates.json.\n')
        print('AUDIT_NEEDS_SOURCE: repeat track ambiguous or absent; no filtered plot produced.',flush=True);return
    record=selected[0]
    if not record['filename'].startswith('Rhync_tenuis_REC.hap1.chr.') or record['restricted']:raise ValueError('Repeat track does not identify exact donor representation')
    target=src/record['filename'];download(f'https://edmond.mpg.de/api/access/datafile/{record["id"]}',target);verify(target,record)
    source=a.source/'prepared/rhynchospora_transfer_reanalysis'
    data=json.loads((source/'fixture.json').read_text());provenance=json.loads((source/'provenance.json').read_text())
    prior=Path(provenance['prior_run']);donor_gff=prior/'sources/donor.gff'
    donor_sha=hashlib.sha256(donor_gff.read_bytes()).hexdigest()
    if donor_sha!=provenance['donor_gff_sha256']:raise ValueError('Donor annotation identity changed')
    donor=data['species'][0]
    if donor['id']!='R_tenuis':raise ValueError('Unexpected donor row')
    lengths={c['id']:c['length'] for c in donor['chromosomes']}
    te,track_counts=read_te(target,lengths);coding=coding_genes(donor_gff)
    ends={c:[b for a,b in spans] for c,spans in te.items()}
    rows=[];flagged=set();without_cds=0
    for c in donor['chromosomes']:
        for block in c['blocks']:
            hid=block['homology_id']
            if hid not in coding:without_cds+=1;rows.append(dict(gene=hid,status='NO_CDS_UNASSESSED'));continue
            g=coding[hid];spans=g['cds'];hits=[];r=te.get(g['chromosome'],[])
            for lo,hi in spans:
                i=bisect_left(ends.get(g['chromosome'],[]),lo+1)
                while i<len(r) and r[i][0]<hi:hits.append((max(lo,r[i][0]),min(hi,r[i][1])));i+=1
            covered=sum(hi-lo for lo,hi in merge(hits));total=sum(hi-lo for lo,hi in spans)
            fraction=covered/total
            if fraction>=a.fraction:flagged.add(hid)
            rows.append(dict(gene=hid,transcript=g['transcript'],cds_bp=total,TE_overlap_bp=covered,fraction=fraction,flagged=hid in flagged))
    audit=dict(repeat_source=record,repeat_track_counts=track_counts,threshold=a.fraction,flagged_genes=len(flagged),without_cds=without_cds,
        scope='Donor CDS TE-overlap screen only; not protein-domain classification or independent screening of target genes',
        warning='TE overlap does not establish TE-derived protein function. Genes without CDS remain unassessed and retained.',genes=rows)
    (a.run/'te_anchor_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    filtered=copy.deepcopy(data)
    for s in filtered['species']:
        for c in s['chromosomes']:c['blocks']=[b for b in c['blocks'] if b['homology_id'] not in flagged]
    out=a.run/'filtered_genes';out.mkdir(exist_ok=True)
    provenance.update(TE_screen=audit|{'genes':'See te_anchor_audit.json'},source_kind='Transfer anchors excluding donor CDS TE-overlap flags',
        scope=audit['scope'],workflow_note='TE-overlap sensitivity test; not exact published evidence. '+audit['warning'])
    pairs=[]
    for x,y in zip(filtered['species'],filtered['species'][1:]):
        ids=lambda s:{b['homology_id'] for c in s['chromosomes'] for b in c['blocks']}
        pairs.append(dict(pair=[x['id'],y['id']],retained=len(ids(x)&ids(y))))
    provenance.update(pair_audit=pairs,retained_links=sum(r['retained'] for r in pairs))
    from syntangle.fixtures import fixture_from_dict
    from syntangle.layout import initial_layout_state,score_crossings
    original_fixture=fixture_from_dict(data);filtered_fixture=fixture_from_dict(filtered)
    original_c=score_crossings(original_fixture,initial_layout_state(original_fixture)).crossings
    filtered_c=score_crossings(filtered_fixture,initial_layout_state(filtered_fixture)).crossings
    audit.update(original_gene_links=sum(r['retained'] for r in json.loads((source/'provenance.json').read_text())['pair_audit']),
        retained_gene_links=provenance['retained_links'],initial_full_C=original_c,initial_screened_C=filtered_c,
        initial_crossing_pairs_touching_flagged_anchor=original_c-filtered_c)
    (a.run/'te_anchor_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    provenance['TE_screen']=audit|{'genes':'See te_anchor_audit.json'}
    (out/'fixture.json').write_text(json.dumps(filtered,indent=2)+'\n')
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    (a.run/'AUDIT_COMPLETE').write_text('PASS\n')
    print(f'TE_AUDIT flagged={len(flagged)} no_CDS_unassessed={without_cds} retained_links={provenance["retained_links"]}',flush=True)
    print(f'INITIAL gene_C={original_c} screened_C={filtered_c} crossing_pairs_touching_flags={original_c-filtered_c}',flush=True)


if __name__=='__main__':main()
