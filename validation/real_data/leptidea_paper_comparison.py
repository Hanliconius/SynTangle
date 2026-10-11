"""Four original Leptidea assemblies: translated Liftoff proteins -> MCScanX blocks."""
import argparse,csv,json,subprocess,sys
from collections import defaultdict
from pathlib import Path
from annotated_lep_pipeline import pair,write_json
from rhynchospora_native_regions import annotations,proteins,sha
from leptidea_liftoff import NAMES,LEPTIDEA

CASE='leptidea_paper_assemblies'

def best_proteins(peptides,transcripts,allowed):
    best={};seen=set();rejected=0
    for tid,seq in peptides:
        if tid in seen:raise ValueError('Duplicate transcript '+tid)
        seen.add(tid)
        if tid not in transcripts:raise ValueError('Unknown translated transcript '+tid)
        gid=transcripts[tid];seq=seq.rstrip('*').upper()
        if gid not in allowed:continue
        if not seq or '*' in seq or not set(seq)<=set('ACDEFGHIKLMNPQRSTVWYBXZJUO'):
            rejected+=1;continue
        if gid not in best or (-len(seq),tid)<(-len(best[gid][1]),best[gid][0]):best[gid]=(tid,seq)
    return best,rejected

def prepare(run,source,gffread):
    if not (source/'homology/COMPLETE').exists():raise ValueError('Recovery incomplete')
    mapping=json.loads((source/'chromosome_map.json').read_text())
    accepted=defaultdict(dict)
    with (source/'homology/candidate_gene_anchors.tsv').open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['donor_gene_id'] in accepted[r['assembly']]:raise ValueError('Duplicate accepted donor gene')
            accepted[r['assembly']][r['donor_gene_id']]=r
    audit=[]
    for name in NAMES:
        m=mapping['assemblies'][name]
        if m['accession']!=LEPTIDEA[name]:raise ValueError('Wrong paper assembly')
        report=source/'targets'/name/'sequence_report.jsonl'
        if sha(report)!=m['sequence_report_sha256']:raise ValueError('Changed sequence report')
        chroms={r['accession']:dict(id=r['sequence_name'],length=int(r['length'])) for r in m['chromosomes']}
        gff=source/'annotations'/name/'lifted.gff';genome=source/'targets'/name/'genome.fna'
        records,transcripts,rows=annotations(gff,set(chroms))
        out=run/'sources'/name;out.mkdir(parents=True)
        selected=out/'translation.gff';selected.write_text('##gff-version 3\n'+''.join(rows))
        private=out/'genome.fna';private.symlink_to(genome.resolve())
        subprocess.run([str(gffread),str(selected),'-g',str(private),'-y',str(out/'translated.faa')],check=True)
        best,rejected=best_proteins(proteins(out/'translated.faa'),transcripts,accepted[name])
        genes={}
        with (out/'proteins.faa').open('w') as f:
            for gid,(tid,seq) in sorted(best.items()):
                r=accepted[name][gid];identifier=name+'__'+gid
                start,end=int(r['start0']),int(r['end'])
                if not 0<=start<end<=int(r['chromosome_length']):raise ValueError('Invalid gene bounds')
                genes[identifier]=dict(species=name,chromosome=r['chromosome'],start=start+1,end=end)
                f.write('>'+identifier+'\n'+seq+'\n')
        if not genes:raise ValueError('No usable translated proteins: '+name)
        write_json(out/'genes.json',genes);write_json(out/'chromosomes.json',chroms)
        audit.append(dict(assembly=name,accession=LEPTIDEA[name],accepted=len(accepted[name]),proteins=len(genes),rejected_translations=rejected,gff_sha256=sha(gff)))
        print('PROTEINS READY',name,len(genes),flush=True)
    write_json(run/'manifest.json',dict(id=CASE,species=[dict(id=n) for n in NAMES],threads=4))
    write_json(run/'input_audit.json',dict(assemblies=audit,recovery=str(source),chromosome_map=mapping,anchor_sha256=sha(source/'homology/candidate_gene_anchors.tsv'),scope='Original paper assemblies; new Liftoff-derived proteins and sequence-based BLAST/MCScanX discovery. Not original published annotation or ribbons.'))

def collect(run):
    from syntangle.fixtures import fixture_from_dict
    rows=[]
    for index in range(6):
        root=run/'pairs'/str(index)
        if not (root/'COMPLETE').exists():raise ValueError('Missing pair '+str(index))
        rows.extend(json.loads((root/'blocks.json').read_text()))
    blocks={n:defaultdict(list) for n in NAMES};pairs=[]
    for a,b in zip(NAMES,NAMES[1:]):
        selected=[r for r in rows if {r['genome1'],r['genome2']}=={a,b}]
        if not selected:raise ValueError('No adjacent blocks')
        pairs.append(dict(pair=[a,b],retained=len(selected)))
        for r in selected:
            for sp in (a,b):
                i=1 if r['genome1']==sp else 2
                blocks[sp][r[f'chr{i}']].append(dict(occurrence_id=r['blkID']+'_'+sp,homology_id=r['blkID'],start=int(r[f'startBp{i}'])-1,end=int(r[f'endBp{i}']),strand='+' if i==1 or r['orient']=='+' else '-'))
    mapping=json.loads((run/'input_audit.json').read_text())['chromosome_map']['assemblies']
    data=dict(fixture_version=1,id=CASE,title='Four paper Leptidea assemblies',purpose='New protein-based block reanalysis',species=[],orientation_constraints=[])
    unlinked={}
    for sp in NAMES:
        chroms=mapping[sp]['chromosomes'];unlinked[sp]=[c['sequence_name'] for c in chroms if c['sequence_name'] not in blocks[sp]]
        data['species'].append(dict(id=sp,chromosomes=[dict(id=c['sequence_name'],length=c['length'],display_rank=k+1,blocks=blocks[sp][c['sequence_name']]) for k,c in enumerate(chroms) if c['sequence_name'] in blocks[sp]]))
    fixture_from_dict(data)
    out=run/'prepared'/CASE;out.mkdir(parents=True)
    write_json(out/'fixture.json',data)
    # The existing comparison importer expects 0-based half-open spans here.
    projected=[{**r,**{f'startBp{i}':int(r[f'startBp{i}'])-1 for i in (1,2)}} for r in rows]
    with (out/'published_block_coordinates.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(projected[0]),delimiter='\t');w.writeheader();w.writerows(projected)
    write_json(out/'provenance.json',dict(scope='Four original male paper assemblies; new Liftoff/protein reanalysis',source_kind='New reciprocal BLASTP/MCScanX blocks',retained_links=sum(x['retained'] for x in pairs),pair_audit=pairs,chromosome_extent='Full recovered assembly chromosome lengths',colour_reference=NAMES[0],row_order=NAMES,unlinked_chromosomes=unlinked,exact_published_display=False,input_layout_note='Original assembly chromosome-number order; published figure order is not yet verified.',orientation_note='MCScanX block directions retained; no subchromosomal rearrangement permitted.',workflow_note='New protein homology searches; shared donor IDs are not used to define matches.',inference_parameters=dict(blastp='Reciprocal default searches; five distinct top subjects',mcscanx='Default except maximum gene gap -m 10'),input_audit=json.loads((run/'input_audit.json').read_text())))
    from plot_comparison import direct_reference_colours
    provenance=json.loads((out/'provenance.json').read_text())
    direct_reference_colours(fixture_from_dict(data),out,provenance)
    provenance['colour_note']='Colours track dominant direct overlap with Spanish L. sinapis chromosomes; ambiguous links grey.'
    write_json(out/'provenance.json',provenance)
    print('PREPARED',CASE,sum(x['retained'] for x in pairs),'scored blocks',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','pair','collect']);p.add_argument('--run',type=Path,required=True);p.add_argument('--source',type=Path);p.add_argument('--gffread',type=Path);p.add_argument('--index',type=int);a=p.parse_args()
    if a.stage=='prepare':prepare(a.run,a.source,a.gffread)
    elif a.stage=='pair':pair(a.run,json.loads((a.run/'manifest.json').read_text()),a.index)
    else:collect(a.run)
