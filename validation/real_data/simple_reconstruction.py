"""Reconstruct fixed evidence from versioned genomes; never claim author evidence."""
import argparse,csv,gzip,hashlib,json,re,shutil,subprocess,sys,tempfile
from collections import defaultdict
from pathlib import Path
from urllib.parse import urljoin,unquote
from zipfile import ZipFile
from simple_papers import fetch,sha

def dump(p,d):
    p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix(p.suffix+'.partial');t.write_text(json.dumps(d,indent=2)+'\n');t.replace(p)
def textopen(p):return gzip.open(p,'rt') if p.suffix=='.gz' else p.open()
def fasta(p):
    name=None;seq=[]
    with textopen(p) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:yield name,''.join(seq)
                name=line[1:].split()[0];seq=[]
            else:seq.append(line.strip())
    if name is not None:yield name,''.join(seq)
def lengths(p):
    d={};name=None
    with textopen(p) as f:
        for line in f:
            if line.startswith('>'):
                name=line[1:].split()[0]
                if name in d:raise ValueError('Duplicate sequence '+name)
                d[name]=0
            elif line.strip():
                if name is None:raise ValueError('FASTA sequence before header')
                d[name]+=len(line.strip())
    if not d or any(n<=0 for n in d.values()):raise ValueError('Empty FASTA')
    return d

def asset(record,out,cache):
    if 'directory' in record:
        page=out/('directory_'+hashlib.sha256(record['directory'].encode()).hexdigest()[:16]+'.html');fetch(record['directory'],page)
        choices=sorted({urljoin(record['directory'],unquote(x)) for x in re.findall(r'href="([^"]+)"',page.read_text()) if x.endswith(record['filename_suffix'])})
        if len(choices)!=1:raise ValueError('Pinned directory lacks unique '+record['filename_suffix']+' '+repr(choices))
        record={**record,'url':choices[0]}
    name=record.get('cache_name') or record['url'].split('?',1)[0].rsplit('/',1)[-1]
    p=out/name
    if cache and (cache/name).exists() and not p.exists():
        # Validate reuse against the published checksum; no unverified cache copy.
        expected=record.get('md5')
        if expected:
            h=hashlib.md5()
            with (cache/name).open('rb') as f:
                for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
            if h.hexdigest()==expected:shutil.copy2(cache/name,p)
    audit=fetch(record['url'],p,record.get('md5'))
    return p,audit

def choose(lengths,spec,report=None):
    if report is not None:
        names={};wanted=spec.get('chromosome_names')
        for r in report:
            if r.get('role')!='assembled-molecule' or r.get('assignedMoleculeLocationType')!='Chromosome':continue
            label=r.get('chrName')
            if label in ('MT','Un') or wanted and label not in wanted:continue
            matches={r.get(k) for k in ('genbankAccession','refseqAccession','sequenceName')}&lengths.keys()
            if len(matches)!=1:raise ValueError('Sequence report identity not unique')
            seq=matches.pop()
            if lengths[seq]!=int(r['length']):raise ValueError('Report length mismatch')
            names[seq]=label
        if wanted and set(names.values())!=set(wanted):raise ValueError('Missing named paper chromosome arms')
        selected={k:lengths[k] for k in names}
    else:
        selected={k:n for k,n in lengths.items() if re.fullmatch(spec['chromosome_regex'],k)};names={k:k for k in selected}
    if len(selected)!=spec['expected_chromosomes']:
        raise ValueError(f'{spec["id"]}: chromosome selector matched {len(selected)}, expected {spec["expected_chromosomes"]}; inspect headers, never choose longest N. Headers: '+repr(list(lengths)[:25]))
    return selected,names

def gene_records(gff,selected):
    records={};parents={};aliases={};excluded=defaultdict(int);gene_instances=defaultdict(list)
    with textopen(gff) as f:
        for line in f:
            if line.startswith('##FASTA'):break
            if line.startswith('#') or not line.strip():continue
            fields=line.rstrip().split('\t')
            if len(fields)!=9:raise ValueError('Malformed GFF')
            if fields[0] not in selected:excluded[fields[0]]+=1;continue
            if fields[2] not in ('gene','mRNA','transcript','CDS'):continue
            start,end=int(fields[3]),int(fields[4])
            if not 1<=start<=end<=selected[fields[0]]:raise ValueError('GFF outside chromosome')
            a={k:unquote(v) for k,v in (x.split('=',1) for x in fields[8].split(';') if '=' in x)}
            ident=a.get('ID')
            if ident:
                if ident in records and fields[2]!='CDS':
                    if records[ident]==(fields,a):continue
                    if fields[2]=='gene':
                        ident=ident+'@'+fields[0]+':'+fields[3]+'-'+fields[4]+':'+fields[6]
                    else:
                        raise ValueError('Distinct placements share feature ID '+ident+'; manual biological review required')
                records[ident]=(fields,a)
                if fields[2]=='gene':gene_instances[a['ID']].append(ident)
                if a.get('Parent'):parents[ident]=a['Parent'].split(',')
            for value in [a.get('protein_id'),a.get('transcript_id')]:
                if value:
                    owner=ident or a.get('Parent','').split(',')[0]
                    aliases.setdefault(value,set()).add(owner)
    def ancestor(x,seen=None):
        seen=set() if seen is None else seen
        if x in seen:raise ValueError('GFF parent cycle')
        if x not in records:return set()
        row,_=records[x]
        if row[2]=='gene':return {x}
        result=set()
        for p in parents.get(x,[]):
            candidates=gene_instances.get(p)
            if candidates:
                candidates=[g for g in candidates if records[g][0][0]==row[0] and records[g][0][6]==row[6] and int(records[g][0][3])<=int(row[3])<=int(row[4])<=int(records[g][0][4])]
                if len(candidates)!=1:raise ValueError('Ambiguous gene placement for '+x)
                p=candidates[0]
            result|=ancestor(p,seen|{x})
        return result
    mapping={}
    for x in records:
        gs=ancestor(x)
        if len(gs)==1:mapping[x]=next(iter(gs))
    for x,owners in aliases.items():
        gs=set().union(*(ancestor(p) for p in owners))
        if len(gs)==1:mapping[x]=next(iter(gs))
    return records,mapping,dict(excluded)

def prepare(run,case,index,cache):
    spec=case['species'][index];out=run/'sources'/spec['id'];out.mkdir(parents=True,exist_ok=True)
    if (out/'COMPLETE').exists() and (out/'input_audit.json').exists():
        prior=json.loads((out/'input_audit.json').read_text())
        if (out/'COMPLETE').read_text().strip()==sha(out/'input_audit.json') and prior['source']==spec and (out/'chromosomes.fna').exists() and sha(out/'chromosomes.fna')==prior['chromosome_fasta_sha256']:
            if case['method']=='minimap2' or all((out/n).exists() for n in ('genes.json','proteins.faa')):
                print('REUSED INPUT',case['id'],spec['id'],flush=True);return
    (out/'COMPLETE').unlink(missing_ok=True)
    if spec.get('blocked'):
        dump(out/'status.json',dict(status='SOURCE_BLOCKED',reason=spec['blocked']));raise ValueError(spec['blocked'])
    audits=[];report=None
    if spec.get('accession'):
        acc=spec['accession'];u=f'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=SEQUENCE_REPORT'
        if case['method']!='minimap2':u+='&include_annotation_type=GENOME_GFF'
        package=out/(acc+'.zip');audits.append(fetch(u,package))
        with ZipFile(package) as z:
            if z.testzip():raise ValueError('Corrupt NCBI ZIP')
            def extract(suffix,target):
                members=[n for n in z.namelist() if n.endswith(suffix)]
                if len(members)!=1:raise ValueError('NCBI package needs one '+suffix)
                with z.open(members[0]) as src,target.open('wb') as dst:shutil.copyfileobj(src,dst)
            genome=out/'genome.fna';extract('.fna',genome)
            rp=out/'sequence_report.jsonl';extract('sequence_report.jsonl',rp)
            report=[json.loads(l) for l in rp.read_text().splitlines() if l.strip()]
            if any(r.get('assemblyAccession')!=acc for r in report):raise ValueError('Wrong NCBI accession')
            if case['method']!='minimap2':gff=out/'genomic.gff';extract('genomic.gff',gff)
    else:
        genome,audit=asset(spec['genome'],out,cache);audits.append(audit)
        if case['method']!='minimap2':gff,audit=asset(spec['gff'],out,cache);audits.append(audit)
    all_lengths=lengths(genome)
    if spec.get('gwh_header_metadata'):
        report=[]
        with textopen(genome) as f:
            for line in f:
                if not line.startswith('>'):continue
                seq=line[1:].split()[0];label=re.search(r'OriSeqID=([^\s]+)',line)
                declared=re.search(r'Len=(\d+)',line)
                if not label or not declared:raise ValueError('GWH FASTA lacks sequence identity metadata')
                if re.fullmatch(r'(?i)Chr0*[1-9][0-9]*',label.group(1)):
                    report.append(dict(role='assembled-molecule',assignedMoleculeLocationType='Chromosome',chrName=label.group(1),genbankAccession=seq,length=int(declared.group(1))))
    selected,names=choose(all_lengths,spec,report)
    # Keep only the documented chromosome scope in alignment input; full sizes retained.
    canonical=out/'chromosomes.fna';temp=canonical.with_suffix('.partial')
    with textopen(genome) as src,temp.open('w') as dst:
        keep=False
        for line in src:
            if line.startswith('>'):
                name=line[1:].split()[0];keep=name in selected
                if keep:dst.write('>'+name+'\n')
            elif keep:dst.write(line)
    temp.replace(canonical)
    if case['method']!='minimap2':
        records,mapping,excluded=gene_records(gff,selected)
        if spec.get('protein'):pep,audit=asset(spec['protein'],out,cache);audits.append(audit)
        else:
            # Translate supplied CDS, never transfer/re-predict annotations silently.
            expanded=out/'annotation.gff'
            if gff!=expanded:
                with textopen(gff) as src,expanded.open('w') as dst:shutil.copyfileobj(src,dst)
            pep=out/'translated.faa';subprocess.run(['gffread',str(expanded),'-g',str(canonical),'-y',str(pep)],check=True)
        best={};unmatched=0
        for pid,seq in fasta(pep):
            gid=mapping.get(pid)
            if gid is None:unmatched+=1;continue
            if gid not in best or (-len(seq),pid)<(-len(best[gid][1]),best[gid][0]):best[gid]=(pid,seq)
        if not best:raise ValueError('No supplied proteins match GFF parent graph; no guessed identifiers')
        genes={}
        with (out/'proteins.faa').open('w') as f:
            for gid,(pid,seq) in sorted(best.items()):
                fields,_=records[gid];key=spec['id']+'__'+gid
                genes[key]=dict(species=spec['id'],chromosome=fields[0],start=int(fields[3]),end=int(fields[4]),strand=fields[6]);f.write('>'+key+'\n'+seq+'\n')
        dump(out/'genes.json',genes)
        audits.append(dict(protein_unmatched=unmatched,selected_genes=len(genes),excluded_annotation_sequences=excluded))
    dump(out/'chromosomes.json',{k:dict(id=k,length=v,display_name=names[k]) for k,v in selected.items()})
    dump(out/'input_audit.json',dict(source=spec,assets=audits,selected=selected,excluded_sequences={k:v for k,v in all_lengths.items() if k not in selected},chromosome_fasta_sha256=sha(canonical),method=case['method']))
    (out/'COMPLETE').write_text(sha(out/'input_audit.json')+'\n');print('PREPARED',case['id'],spec['id'],len(selected),flush=True)

def paf_rows(path,a,b,ca,cb,index,minbp,mapq):
    rows=[]
    with path.open() as f:
        for line in f:
            c=line.rstrip().split('\t')
            if len(c)<12:raise ValueError('Malformed PAF')
            tags={x[:2]:x[5:] for x in c[12:] if len(x)>=5}
            if tags.get('tp')!='P' or int(c[10])<minbp or int(c[11])<mapq:continue
            q,t=c[0],c[5]
            if q not in ca or t not in cb:raise ValueError('PAF chromosome missing')
            if int(c[1])!=ca[q]['length'] or int(c[6])!=cb[t]['length']:raise ValueError('PAF chromosome length changed')
            qs,qe,ts,te=map(int,[c[2],c[3],c[7],c[8]])
            if c[4] not in ('+','-') or not 0<=qs<qe<=int(c[1]) or not 0<=ts<te<=int(c[6]):raise ValueError('Invalid signed endpoints')
            rows.append(dict(blkID=f'pair{index}_alignment{len(rows)}',orient=c[4],genome1=a,chr1=q,startBp1=qs+1,endBp1=qe,genome2=b,chr2=t,startBp2=ts+1,endBp2=te))
    return rows

def pair(run,case,index):
    a,b=case['pairs'][index];out=run/'pairs'/str(index);out.mkdir(parents=True,exist_ok=True)
    if (out/'COMPLETE').exists() and (out/'blocks.json').exists() and (out/'pair_audit.json').exists():
        prior=json.loads((out/'pair_audit.json').read_text())
        if {k:v for k,v in prior.get('parameters',{}).items() if k!='pairs'}=={k:v for k,v in case.items() if k!='pairs'} and prior.get('pair')==[a,b] and (out/'COMPLETE').read_text().strip()==sha(out/'blocks.json'):
            print('REUSED PAIR',a,b,flush=True);return
    (out/'COMPLETE').unlink(missing_ok=True)
    for sp in (a,b):
        p=run/'sources'/sp
        if not (p/'COMPLETE').exists() or (p/'COMPLETE').read_text().strip()!=sha(p/'input_audit.json'):raise ValueError('Unvalidated input '+sp)
    if case['method']=='minimap2':
        path=out/'alignment.paf';temp=out/'alignment.paf.partial'
        try:
            with temp.open('w') as f:subprocess.run(['minimap2','-x','asm20','--secondary=no','-t',str(case['threads']),str(run/'sources'/b/'chromosomes.fna'),str(run/'sources'/a/'chromosomes.fna')],stdout=f,check=True)
            temp.replace(path)
        finally:temp.unlink(missing_ok=True)
        ca=json.loads((run/'sources'/a/'chromosomes.json').read_text());cb=json.loads((run/'sources'/b/'chromosomes.json').read_text())
        rows=paf_rows(path,a,b,ca,cb,index,case['min_alignment_bp'],case['min_mapq'])
    else:
        genes={};hits=[]
        for sp in (a,b):genes.update(json.loads((run/'sources'/sp/'genes.json').read_text()))
        for q,t in ((a,b),(b,a)):
            db=out/('db_'+t);subprocess.run(['makeblastdb','-in',str(run/'sources'/t/'proteins.faa'),'-dbtype','prot','-out',str(db)],stdout=subprocess.DEVNULL,check=True)
            raw=out/(q+'.tsv');temp=raw.with_suffix('.partial')
            try:
                subprocess.run(['blastp','-query',str(run/'sources'/q/'proteins.faa'),'-db',str(db),'-evalue','1e-5','-outfmt','6','-num_threads',str(case['threads']),'-out',str(temp)],check=True);temp.replace(raw)
            finally:temp.unlink(missing_ok=True)
            grouped=defaultdict(list)
            with raw.open() as f:
                for l in f:cols=l.rstrip().split('\t');grouped[cols[0]].append(cols)
            for values in grouped.values():
                seen=set()
                for c in sorted(values,key=lambda c:(-float(c[11]),float(c[10]),c[1])):
                    if c[1] in seen:continue
                    if len(seen)>=case['blast_best_subjects']:break
                    seen.add(c[1]);hits.append('\t'.join(c))
        prefix=out/'synteny';prefix.with_suffix('.blast').write_text('\n'.join(hits)+'\n')
        with prefix.with_suffix('.gff').open('w') as f:
            for gid,g in sorted(genes.items()):f.write(f'{g["species"]}__{g["chromosome"]}\t{gid}\t{g["start"]}\t{g["end"]}\n')
        subprocess.run(['MCScanX',str(prefix),'-s',str(case['mcscanx_min_anchors'])],check=True)
        rows=[];block=[];sign=None;bid=None
        def finish():
            if not block:return
            r=dict(blkID=f'pair{index}_block{bid}',orient=sign)
            for i in (1,2):
                gs=[genes[v[i-1]] for v in block]
                if len({(g['species'],g['chromosome']) for g in gs})!=1:raise ValueError('Block crosses chromosome')
                r.update({f'genome{i}':gs[0]['species'],f'chr{i}':gs[0]['chromosome'],f'startBp{i}':min(g['start'] for g in gs),f'endBp{i}':max(g['end'] for g in gs)})
            if {r['genome1'],r['genome2']}=={a,b}:rows.append(r)
        for l in prefix.with_suffix('.collinearity').read_text().splitlines():
            if l.startswith('## Alignment'):
                finish();block=[];bid=l.split(':')[0].split()[-1];orient=l.split()[-1]
                if orient not in ('plus','minus'):raise ValueError('Unknown orientation')
                sign='+' if orient=='plus' else '-'
            elif l.strip() and not l.startswith('#'):block.append(tuple(l.split(':',1)[1].split()[:2]))
        finish()
    if not rows:raise ValueError('No retained blocks; inspect alignment, no empty comparison')
    dump(out/'blocks.json',rows);dump(out/'pair_audit.json',dict(pair=[a,b],method=case['method'],parameters=case,blocks_sha256=sha(out/'blocks.json'),fixed_links=len(rows)))
    (out/'COMPLETE').write_text(sha(out/'blocks.json')+'\n');print('PAIR',a,b,len(rows),flush=True)

def fixture(run,case,pair_index=None,tag=""):
    from syntangle.fixtures import fixture_from_dict
    from prepare_genespace_blocks import natural
    from plot_comparison import direct_reference_colours
    indices=range(len(case['pairs'])) if pair_index is None else [pair_index];rows=[]
    for i in indices:
        p=run/'pairs'/str(i)
        if not (p/'COMPLETE').exists() or (p/'COMPLETE').read_text().strip()!=sha(p/'blocks.json'):raise ValueError('Missing/tampered pair '+str(i))
        rows+=json.loads((p/'blocks.json').read_text())
    all_rows=rows
    ids=[s['id'] for s in case['species']] if pair_index is None else case['pairs'][pair_index]
    adjacent={frozenset((a,b)) for a,b in zip(ids,ids[1:])}
    rows=[r for r in rows if frozenset((r['genome1'],r['genome2'])) in adjacent]
    blocks={sp:defaultdict(list) for sp in ids}
    for r in rows:
        for i in (1,2):
            sp=r[f'genome{i}'];blocks[sp][r[f'chr{i}']].append(dict(occurrence_id=r['blkID']+'_'+sp,homology_id=r['blkID'],start=r[f'startBp{i}']-1,end=r[f'endBp{i}'],strand='+' if i==1 else r['orient']))
    data=dict(fixture_version=1,id=case['id']+('_pair'+str(pair_index) if pair_index is not None else ''),title=case['id'],purpose='Reconstructed homology, fixed endpoints',species=[],orientation_constraints=[])
    for sp in ids:
        chroms=json.loads((run/'sources'/sp/'chromosomes.json').read_text());ordered=sorted(chroms,key=natural)
        data['species'].append(dict(id=sp,chromosomes=[dict(id=c,length=chroms[c]['length'],display_rank=i+1,blocks=blocks[sp][c]) for i,c in enumerate(ordered)]))
    obj=fixture_from_dict(data);out=run/'prepared'/(data['id']+'_'+tag);dump(out/'fixture.json',data)
    converted=[{**r,**{f'startBp{i}':r[f'startBp{i}']-1 for i in (1,2)}} for r in all_rows]
    with (out/'published_block_coordinates.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(converted[0]),delimiter='\t');w.writeheader();w.writerows(converted)
    provenance=dict(pair_audit=[dict(pair=case['pairs'][i],retained=len(json.loads((run/'pairs'/str(i)/'blocks.json').read_text()))) for i in indices],paper=case['paper'],figure=case['figure'],exact_published_display=False,source_kind='Reconstructed '+case['method']+' evidence',scope=case['reconstruction_note'],retained_links=len(rows),colour_reference=case['colour_reference'],input_panel_title='Reconstructed input',input_layout_note='Natural chromosome labels; published orientation not yet verified',chromosome_extent='Full fixed assembly chromosome lengths',row_order=ids)
    provenance['chromosome_labels']={sp:{k:v['display_name'] for k,v in json.loads((run/'sources'/sp/'chromosomes.json').read_text()).items()} for sp in ids}
    direct_reference_colours(obj,out,provenance);dump(out/'provenance.json',provenance)
    return out

def compare(run,case,index,python):
    # Ficus is four independent edges, never a fabricated five-species chain.
    pair_index=index//3 if case['id']=='ficus' else None;budget=[30,120,600][index%3]
    prepared=fixture(run,case,pair_index,str(budget));out=run/'figures'/(prepared.name+'_'+str(budget));out.mkdir(parents=True,exist_ok=True)
    subprocess.run([python,str(Path(__file__).with_name('plot_comparison.py')),'--prepared',str(prepared),'--output',str(out),'--seconds',str(budget)],check=True)
    from preservation_audit import require_audit
    require_audit(out);(out/'COMPLETE').write_text('PASS endpoint preservation\n');print('COMPARISON',case['id'],budget,out,flush=True)

def collect(run):
    reports=[]
    for c in json.loads((run/'code/validation/real_data/simple_reconstruction_manifest.json').read_text())['cases']:
        root=run/c['id'];reports.append(dict(case=c['id'],prepared_species=[p.parent.name for p in (root/'sources').glob('*/COMPLETE')],completed_pairs=[p.parent.name for p in (root/'pairs').glob('*/COMPLETE')],completed_comparisons=[str(p.parent.relative_to(run)) for p in (root/'figures').glob('*/COMPLETE')],blocked=[json.loads(p.read_text()) for p in (root/'sources').glob('*/status.json')]))
    dump(run/'reconstruction_summary.json',reports)
    # Archive figures/audits/fixtures only; multi-GB genomes stay outside download bundle.
    archive=run/'SynTangle_simple_reconstructed.zip'
    with ZipFile(archive.with_suffix('.partial'),'w') as z:
        for p in run.rglob('*'):
            if p.is_file() and ('figures' in p.parts or 'prepared' in p.parts or p.name.endswith('audit.json') or p.name in ('reconstruction_summary.json','simple_reconstruction_manifest.json')):z.write(p,p.relative_to(run))
    archive.with_suffix('.partial').replace(archive);print(json.dumps(reports,indent=2));print('ARCHIVE',archive)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','pair','compare','collect']);p.add_argument('--run',type=Path,required=True);p.add_argument('--case',type=int,default=0);p.add_argument('--index',type=int,default=0);p.add_argument('--cache',type=Path);p.add_argument('--python',default=sys.executable);a=p.parse_args()
    if a.stage=='collect':collect(a.run)
    else:
        case=json.loads(Path(__file__).with_name('simple_reconstruction_manifest.json').read_text())['cases'][a.case];root=a.run/case['id']
        if a.stage=='prepare':prepare(root,case,a.index,a.cache)
        elif a.stage=='pair':pair(root,case,a.index)
        else:compare(root,case,a.index,a.python)
