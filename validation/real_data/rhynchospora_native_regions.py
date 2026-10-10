"""Protein-based GENESPACE reanalysis and exact native-ribbon import, Pegasus only."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
from urllib.parse import unquote

IDS = ['R_tenuis', 'R_austrobrasiliensis', 'R_breviuscula']
REF_ORDER = ['Chr2_h1','Chr5_h1','Chr1_h1','Chr4_h1','Chr3_h1']
PALETTE = ['#F9AC60','#307BB5','#D61F27','#ADD8E7','#FCF7BF']


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''): h.update(chunk)
    return h.hexdigest()


def annotations(path, allowed):
    records={}; rows=[]
    with path.open() as f:
        for line in f:
            if line.startswith('##FASTA'):break
            if line.startswith('#') or not line.strip():continue
            r=line.rstrip().split('\t')
            if len(r)!=9:raise ValueError('Malformed GFF')
            if r[0] not in allowed:continue
            attrs={k:unquote(v) for k,v in (a.split('=',1) for a in r[8].split(';') if '=' in a)}
            if 'ID' in attrs and r[2] in ('gene','mRNA','transcript'):
                key=attrs['ID']
                if key in records:raise ValueError('Duplicate selected GFF ID: '+key)
                records[key]=(r,attrs)
            rows.append(line)
    def gene(identifier, trail=()):
        if identifier in trail:raise ValueError('Cyclic GFF parents')
        if identifier not in records:raise ValueError('Unknown GFF parent '+identifier)
        r,a=records[identifier]
        if r[2]=='gene':return identifier
        parents=a.get('Parent','').split(',')
        found={gene(x,trail+(identifier,)) for x in parents if x}
        if len(found)!=1:raise ValueError('Ambiguous GFF gene parent '+identifier)
        return found.pop()
    transcripts={k:gene(k) for k,(r,a) in records.items() if r[2] in ('mRNA','transcript')}
    return records,transcripts,rows


def proteins(path):
    name=None;parts=[]
    with path.open() as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:yield name,''.join(parts)
                name=line[1:].split()[0];parts=[]
            else:parts.append(line.strip())
    if name is not None:yield name,''.join(parts)


def prepare(run, source):
    manifest=json.loads((source/'source_manifest.json').read_text())
    prior=Path(manifest['prior_run']); inputs=Path(manifest['inputs'])
    original=json.loads((source/'prepared/rhynchospora_transfer_reanalysis/fixture.json').read_text())
    paths=[(prior/'sources/donor.fasta',prior/'sources/donor.gff'),
           (source/'sources/selected_austro.fasta',inputs/'sources/R_austrobrasiliensis_paper_deposit/Rhync_austrobrasiliensis_6344D.hic.hap1.chr.helixer.gff3'),
           (prior/'sources/target.fasta',prior/'annotation/lifted.gff')]
    gs=run/'genespace';(gs/'bed').mkdir(parents=True);(gs/'peptide').mkdir()
    raw=run/'translation';raw.mkdir(); audit=[]
    for sp,(fasta,gff),spec in zip(IDS,paths,original['species']):
        if sp!=spec['id']:raise ValueError('Species mismatch')
        lengths={c['id']:c['length'] for c in spec['chromosomes']}
        records,transcripts,rows=annotations(gff,set(lengths))
        selected=raw/(sp+'.gff3');selected.write_text('##gff-version 3\n'+''.join(rows))
        private=raw/(sp+'.fasta');private.symlink_to(fasta.resolve())
        pep=raw/(sp+'.fa')
        print('TRANSLATING '+sp,flush=True)
        subprocess.run(['gffread',str(selected),'-g',str(private),'-y',str(pep)],check=True)
        best={}; excluded=0;seen=set()
        for tid,seq in proteins(pep):
            if tid in seen:raise ValueError('Duplicate translated transcript')
            seen.add(tid)
            if tid not in transcripts:raise ValueError('Unmapped translated transcript '+tid)
            seq=seq.rstrip('*').upper()
            if not seq or '*' in seq or not set(seq)<=set('ACDEFGHIKLMNPQRSTVWYBXZJUO'):excluded+=1;continue
            gid=transcripts[tid]
            if gid not in best or (-len(seq),tid)<(-len(best[gid][1]),best[gid][0]):best[gid]=(tid,seq)
        with (gs/'bed'/(sp+'.bed')).open('w') as bed,(gs/'peptide'/(sp+'.fa')).open('w') as fa:
            for n,(gid,(tid,seq)) in enumerate(sorted(best.items())):
                r,_=records[gid];start,end=int(r[3]),int(r[4])
                if not 1<=start<=end<=lengths[r[0]]:raise ValueError('Gene outside chromosome')
                identifier=f'{sp}_g{n:07d}'
                # GENESPACE's parser preserves GFF's 1-based inclusive endpoints.
                bed.write(f'{r[0]}\t{start}\t{end}\t{identifier}\n');fa.write(f'>{identifier}\n{seq}\n')
        if {records[g][0][0] for g in best}!=set(lengths):raise ValueError('A selected chromosome lacks usable peptides')
        audit.append(dict(species=sp,fasta=str(fasta),gff=str(gff),fasta_sha256=sha(fasta),gff_sha256=sha(gff),genes=len(best),rejected_translations=excluded,annotation='Liftoff transfer from Tenuis' if sp==IDS[-1] else 'Deposited Helixer annotation'))
    (run/'input_audit.json').write_text(json.dumps(dict(inputs=audit,homology='OrthoFinder peptide sequence inference; donor ID equality is not used',scope='Three-genome reanalysis; not original twenty-genome run',chromosomes=[dict(id=sp['id'],chromosomes=[dict(id=c['id'],length=c['length'],blocks=[]) for c in sp['chromosomes']]) for sp in original['species']]),indent=2)+'\n')
    print('PREPARED native GENESPACE peptide and BED inputs',flush=True)


def import_regions(rows, species):
    import copy
    data=dict(fixture_version=1,id='rhynchospora_native_regions',title='Rhynchospora native GENESPACE regions',purpose='Fixed native region-midpoint comparison',species=copy.deepcopy(species))
    targets={(s['id'],c['id']):c for s in data['species'] for c in s['chromosomes']}
    for c in targets.values():c['blocks']=[]
    colours={};refmap={};audit=[]
    pairs=set(zip(IDS,IDS[1:]))
    for i,r in enumerate(rows):
        a,b=r['genome1'],r['genome2']
        if (a,b) not in pairs:raise ValueError('Native plotted ribbon is not an adjacent forward pair')
        hid=f'native_region_{i:07d}'
        for side,sp in ((1,a),(2,b)):
            key=(sp,r[f'chr{side}'])
            if key not in targets:raise ValueError('Unknown native chromosome')
            x,y=float(r[f'startBp{side}']),float(r[f'endBp{side}'])
            lo,hi=min(x,y)-1,max(x,y)
            if not 0<=lo<hi<=targets[key]['length']:raise ValueError('Native region out of bounds')
            targets[key]['blocks'].append(dict(occurrence_id=sp+'_'+hid,homology_id=hid,start=lo,end=hi,strand='-' if x>y else '+'))
        if r['refChr'] not in REF_ORDER:raise ValueError('Unknown native reference chromosome')
        colours[hid]=r['color'];refmap[hid]=r['refChr'];audit.append(dict(homology_id=hid,source_row=i,native=r))
    if not audit:raise ValueError('No native plotted regions')
    return data,colours,refmap,audit


def compare(run, seconds):
    from syntangle.fixtures import fixture_from_dict
    from syntangle.layout import initial_layout_state,score_crossings,LayoutState
    from syntangle.hybrid_experiment import run_hybrid
    from syntangle.saved_layout import decode_saved_layout
    from plot_comparison import render,atomic_json,improve_flips
    import plot_comparison,visual_summary
    plot_comparison.riparian=lambda *args,**kw: visual_summary.riparian(*args,**kw,respect_block_strand=True)
    with (run/'native_regions.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    audit=json.loads((run/'input_audit.json').read_text())
    data,colours,refmap,membership=import_regions(rows,audit['chromosomes'])
    fixture=fixture_from_dict(data);initial=initial_layout_state(fixture)
    with (run/'native_regions_chromosomes.tsv').open() as f:chroms=list(csv.DictReader(f,delimiter='\t'))
    refs={(r.species_id,r.chromosome_id):r for r in fixture.chromosome_refs}
    actual={(r['genome'],r['chr']) for r in chroms}
    if len(chroms)!=len(actual) or actual!=set(refs):raise ValueError('Native plot omitted/duplicated chromosomes')
    order={sp:tuple(refs[sp,r['chr']] for r in sorted((r for r in chroms if r['genome']==sp),key=lambda r:float(r['plotOrd']))) for sp in IDS}
    orientation={r:(-1 if r.species_id==IDS[-1] and r.chromosome_id=='Chr3_h1' else 1) for r in refs.values()}
    native=LayoutState(order,orientation)
    out=run/'figures/rhynchospora_native_regions';out.mkdir(parents=True)
    provenance=dict(pair_audit=[dict(pair=[a,b],retained=sum(r['genome1']==a and r['genome2']==b for r in rows)) for a,b in zip(IDS,IDS[1:])],orientation_note='Native ribbon endpoint directions retained; whole-chromosome flips only.',source_kind='Native GENESPACE useRegions=TRUE plotted ribbons; three-genome reanalysis',colour_reference=IDS[-1],reference_palette=dict(zip(REF_ORDER,PALETTE)),homology_colour_reference=refmap,retained_links=len(rows),scope='Same fixed native region evidence in every panel; full selected FASTA lengths; whole-chromosome moves only',objective='Unweighted native region-midpoint crossings; not gene-level or polygon-overlap optimality',input_panel_title='Selected chromosome-name order',chromosome_extent='Full selected FASTA lengths',crossing_unit='native region',input_layout_note='Selected sequence order; exact Fig.1a order unverified.',anchor_description='sequence-inferred orthology and native regions',gs_workflow_note='GS baseline comes from the full native discovery and region plotting workflow, with the author reference order and inversion.',workflow_note='Deposited T/A annotations; transferred B annotation; three-genome discovery rather than the original twenty genomes.')
    (out/'fixture.json').write_text(json.dumps(data,indent=2)+'\n')
    atomic_json(out/'native_region_import.json',dict(regions=membership,source_sha256=sha(run/'native_regions.tsv')))
    score=lambda state:score_crossings(fixture,state).crossings
    panels=[(initial,'Native regions: initial display','Full chromosome lengths'),(native,'Native GENESPACE region display','Author reference order / Chr3 inversion')]
    render(fixture,panels,out,provenance)
    print('FLIP ASSISTANCE START',flush=True)
    assisted,_=improve_flips(fixture,native,0)
    panels.append((assisted,'Native GENESPACE + our flips','Same native region evidence'))
    start=min((initial,native,assisted),key=score);best=start
    def checkpoint(state):
        nonlocal best
        if score(state)<=score(best):best=state;atomic_json(out/'incumbent.json',dict(crossings=score(state),state=state.to_dict()))
    checkpoint(start);render(fixture,panels,out,provenance)
    print('SYNTANGLE START',flush=True)
    result=run_hybrid(fixture,start,'hybrid_adaptive',seconds,progress=checkpoint,mirror_symmetry=True,neighborhood_min_decisions=512,neighborhood_fraction=.15,neighborhood_max_seconds=30,deduplicate_neighborhoods=True)
    final=decode_saved_layout(fixture,result['optimized_state'])
    if score(final)!=result['upper_bound'] or score(final)>score(start):raise AssertionError('Final score mismatch')
    panels.append((final,'SynTangle: same native regions',f"L={result['lower_bound']:,}; U={result['upper_bound']:,}; {result['optimality_status']}"))
    render(fixture,panels,out,provenance);atomic_json(out/'result.json',result)
    atomic_json(out/'comparison_audit.json',dict(provenance=provenance,panels=[dict(title=t,crossings=score(s),state=s.to_dict()) for s,t,n in panels]))
    (out/'COMPLETE').write_text('PASS\n')
    print(f"NATIVE REGIONS={len(rows)} GS={score(native)} GS+flips={score(assisted)} ST={score(final)} L={result['lower_bound']}",flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','compare']);p.add_argument('--run',type=Path,required=True);p.add_argument('--source',type=Path);p.add_argument('--seconds',type=float,default=300);a=p.parse_args()
    if a.stage=='prepare':prepare(a.run,a.source)
    else:compare(a.run,a.seconds)
