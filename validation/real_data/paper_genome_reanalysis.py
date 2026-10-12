"""Versioned codfish/Carex preparation; new protein homology, never paper-link claims."""
import argparse,csv,gzip,hashlib,itertools,json,re,shutil,subprocess,xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path
from urllib.parse import urlencode
from zipfile import ZipFile
from prepare_published_inputs import download,digest
from annotated_lep_pipeline import pair,write_json
from rhynchospora_native_regions import annotations,proteins
from leptidea_paper_comparison import best_proteins

def verified_download(record,path):
    download(record['url'],path)
    expected=record.get('md5')
    if expected:
        h=hashlib.md5()
        with path.open('rb') as f:
            for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
        if h.hexdigest()!=expected:
            path.unlink();raise ValueError('Deposited checksum mismatch: '+str(path))
    return dict(**record,sha256=digest(path))

def unzip_member(z,suffix,out):
    members=[n for n in z.namelist() if n.endswith(suffix)]
    if len(members)!=1:raise ValueError('NCBI package needs exactly one '+suffix+'; found '+repr(members))
    with z.open(members[0]) as src,out.open('wb') as dst:shutil.copyfileobj(src,dst)
    return members[0]

def expand(path,out):
    with gzip.open(path,'rb') as src,out.open('wb') as dst:shutil.copyfileobj(src,dst)

def resolve_accession(nucleotide,out):
    base='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    search=out/'nucleotide_search.xml'
    download(base+'esearch.fcgi?'+urlencode(dict(db='nuccore',term=nucleotide+'[Accession]',retmode='xml')),search)
    ids=[r.text for r in ET.parse(search).findall('./IdList/Id')]
    if len(ids)!=1:raise ValueError('Nucleotide does not resolve uniquely: '+nucleotide)
    links=out/'assembly_link.xml';download(base+'elink.fcgi?'+urlencode(dict(dbfrom='nuccore',db='assembly',id=ids[0],retmode='xml')),links)
    assembly_ids=sorted({r.text for r in ET.parse(links).findall('.//LinkSetDb/Link/Id')})
    if len(assembly_ids)!=1:raise ValueError('Published nucleotide has ambiguous/missing assembly link')
    summary=out/'assembly_summary.json';download(base+'esummary.fcgi?'+urlencode(dict(db='assembly',id=assembly_ids[0],retmode='json')),summary)
    accession=json.loads(summary.read_text())['result'][assembly_ids[0]]['assemblyaccession']
    if not re.fullmatch(r'GC[AF]_\d+\.\d+',accession):raise ValueError('Unversioned assembly')
    return accession

def fasta_lengths(path):
    lengths={};name=None
    with path.open() as f:
        for line in f:
            if line.startswith('>'):
                name=line[1:].split()[0]
                if name in lengths:raise ValueError('Duplicate FASTA sequence '+name)
                lengths[name]=0
            elif line.strip():
                if name is None:raise ValueError('Sequence before FASTA header')
                lengths[name]+=len(line.strip())
    if not lengths or any(n<=0 for n in lengths.values()):raise ValueError('Empty FASTA sequence')
    return lengths

def select_chromosomes(lengths,spec,report=None):
    if report is None:
        minimum=spec['min_chromosome_length'];selected={k:n for k,n in lengths.items() if n>minimum}
        rule=f'Deposit sequences longer than {minimum} bp; count verified against paper'
    else:
        if spec.get('resolve_from_nucleotide'):
            # The author explicitly names CM022079 through CM022107.
            chosen=[r for r in report if r.get('genbankAccession','').split('.')[0] in {f'CM{i:06d}' for i in range(22079,22108)}]
            rule='Published CM022079-CM022107 accession set'
        else:
            chosen=[r for r in report if r.get('role')=='assembled-molecule' and r.get('assignedMoleculeLocationType')=='Chromosome' and r.get('chrName') not in ('MT','Un')]
            rule='NCBI assembled nuclear chromosome molecules'
        selected={}
        for row in chosen:
            keys=[row.get('refseqAccession'),row.get('genbankAccession'),row.get('sequenceName')]
            found=[k for k in keys if k in lengths]
            if len(set(found))!=1:raise ValueError('Sequence report/FASTA identity ambiguous')
            k=found[0]
            if lengths[k]!=int(row['length']):raise ValueError('Sequence report length differs')
            selected[k]=lengths[k]
    expected=spec.get('expected_chromosomes')
    if not selected or expected is not None and len(selected)!=expected:
        raise ValueError(f'{spec["id"]}: selected {len(selected)} sequences; expected {expected}; inspect source, no longest-N fallback')
    return selected,rule

def prepare_species(run,case,index,gffread):
    spec=case['species'][index];out=run/'sources'/spec['id'];out.mkdir(parents=True,exist_ok=True)
    files={};report=None
    genome=out/'genome.fna';gff=out/'annotation.gff'
    if 'genome' in spec:
        for key,record,target in [('genome',spec['genome'],genome),('gff',spec['gff'],gff)]:
            compressed=out/(key+'.gz');files[key]=verified_download(record,compressed);expand(compressed,target)
    else:
        accession=spec.get('accession') or resolve_accession(spec['resolve_from_nucleotide'],out)
        package=out/(accession+'.zip');types='GENOME_FASTA&include_annotation_type=SEQUENCE_REPORT'
        if 'gff' not in spec:types+='&include_annotation_type=GENOME_GFF'
        url=f'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download?include_annotation_type={types}'
        files['ncbi_package']=verified_download(dict(url=url),package)
        with ZipFile(package) as z:
            if z.testzip() is not None:raise ValueError('Corrupt NCBI package')
            unzip_member(z,'.fna',genome);unzip_member(z,'sequence_report.jsonl',out/'sequence_report.jsonl')
            if 'gff' not in spec:unzip_member(z,'genomic.gff',gff)
        report=[json.loads(l) for l in (out/'sequence_report.jsonl').read_text().splitlines() if l.strip()]
        if not report or any(r.get('assemblyAccession')!=accession for r in report):raise ValueError('NCBI returned different accession')
        if 'gff' in spec:
            compressed=out/'gff.gz';files['gff']=verified_download(spec['gff'],compressed);expand(compressed,gff)
        files['versioned_accession']=accession
    lengths=fasta_lengths(genome);chosen,rule=select_chromosomes(lengths,spec,report)
    records,transcripts,rows=annotations(gff,set(chosen))
    selected=out/'translation.gff';selected.write_text('##gff-version 3\n'+''.join(rows))
    pep=out/'translated.faa';subprocess.run([str(gffread),str(selected),'-g',str(genome),'-y',str(pep)],check=True)
    allowed={gid for gid,(r,attrs) in records.items() if r[2]=='gene'}
    best,rejected=best_proteins(proteins(pep),transcripts,allowed)
    if not best:raise ValueError('No matching translated genes: '+spec['id'])
    genes={}
    with (out/'proteins.faa').open('w') as f:
        for gid,(tid,seq) in sorted(best.items()):
            r,_=records[gid];start,end=int(r[3]),int(r[4]);identifier=spec['id']+'__'+gid
            if not 1<=start<=end<=chosen[r[0]]:raise ValueError('Gene outside fixed chromosome')
            genes[identifier]=dict(species=spec['id'],chromosome=r[0],start=start,end=end)
            f.write('>'+identifier+'\n'+seq+'\n')
    unannotated=sorted(set(chosen)-{g['chromosome'] for g in genes.values()})
    if unannotated:raise ValueError('Selected chromosomes lack usable proteins: '+repr(unannotated))
    write_json(out/'genes.json',genes);write_json(out/'chromosomes.json',{k:dict(id=k,length=n) for k,n in chosen.items()})
    write_json(out/'input_audit.json',dict(species=spec,files=files,chromosome_selection=rule,selected=chosen,excluded_sequences={k:n for k,n in lengths.items() if k not in chosen},gene_count=len(genes),rejected_translations=rejected,genome_sha256=digest(genome),gff_sha256=digest(gff)))
    (out/'COMPLETE').write_text('PASS versioned inputs, chromosome count and annotation coordinates\n')
    print('PREPARED',spec['id'],len(chosen),'chromosomes',len(genes),'proteins',flush=True)

def collect(run,case):
    from syntangle.fixtures import fixture_from_dict
    from prepare_genespace_blocks import natural
    from plot_comparison import direct_reference_colours
    ids=[d['id'] for d in case['species']];blocks={s:defaultdict(list) for s in ids};rows=[]
    for i in range(len(list(itertools.combinations(ids,2)))):
        root=run/'pairs'/str(i)
        if not (root/'COMPLETE').exists():raise ValueError('Missing completed pair '+str(i))
        rows.extend(json.loads((root/'blocks.json').read_text()))
    pairs=[];adjacent=[]
    for a,b in zip(ids,ids[1:]):
        selected=[r for r in rows if {r['genome1'],r['genome2']}=={a,b}]
        if not selected:raise ValueError('Missing adjacent homology pair')
        pairs.append(dict(pair=[a,b],retained=len(selected)));adjacent.extend(selected)
        for r in selected:
            for sp in (a,b):
                i=1 if r['genome1']==sp else 2
                blocks[sp][r[f'chr{i}']].append(dict(occurrence_id=r['blkID']+'_'+sp,homology_id=r['blkID'],start=int(r[f'startBp{i}'])-1,end=int(r[f'endBp{i}']),strand='+' if i==1 or r['orient']=='+' else '-'))
    data=dict(fixture_version=1,id=case['id'],title=case['id'],purpose='New protein-based homology reanalysis',species=[],orientation_constraints=[]);unlinked={}
    for sp in ids:
        chroms=json.loads((run/'sources'/sp/'chromosomes.json').read_text())
        ordered=sorted(chroms,key=(lambda c:(-chroms[c]['length'],c)) if case['size_order'] and sp!=case['colour_reference'] else natural)
        unlinked[sp]=[c for c in ordered if c not in blocks[sp]]
        # Every selected full chromosome remains present, including unlinked ones.
        data['species'].append(dict(id=sp,chromosomes=[dict(id=c,length=chroms[c]['length'],display_rank=i+1,blocks=blocks[sp][c]) for i,c in enumerate(ordered)]))
    fixture=fixture_from_dict(data);out=run/'prepared'/case['id'];out.mkdir(parents=True,exist_ok=True)
    write_json(out/'fixture.json',data)
    projected=[{**r,**{f'startBp{i}':int(r[f'startBp{i}'])-1 for i in (1,2)}} for r in rows]
    with (out/'published_block_coordinates.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(projected[0]),delimiter='\t');w.writeheader();w.writerows(projected)
    provenance=dict(paper=case['paper'],figure=case['figure'],exact_published_display=False,scope='Original study species; deposited assembly/annotation or accession-matched reconstruction; newly inferred fixed blocks',source_kind='New reciprocal BLASTP/MCScanX block inference',retained_links=len(adjacent),pair_audit=pairs,colour_reference=case['colour_reference'],chromosome_extent='Full selected assembly sequence lengths',row_order=ids,input_layout_note=case['input_layout_note'],orientation_note='Signed MCScanX block endpoints; only whole chromosome moves allowed',workflow_note='Reanalysis, not original paper evidence. GENESPACE ordering function evaluated on block proxy anchors; not a gene-level GENESPACE discovery rerun.',unlinked_chromosomes=unlinked,source_audits=[json.loads((run/'sources'/sp/'input_audit.json').read_text()) for sp in ids])
    direct_reference_colours(fixture,out,provenance)
    provenance['colour_note']='Colours follow dominant direct overlap with '+case['colour_reference']+' chromosomes; ambiguous links remain grey.'
    write_json(out/'provenance.json',provenance);print('FIXTURE',case['id'],len(adjacent),'fixed signed blocks',flush=True)

def inventory(run,manifest):
    out=run/'carex_source_inventory';out.mkdir(parents=True,exist_ok=True);records=[]
    for record in manifest['carex_author_records']+[manifest['carex_annotation_record']]:
        path=out/(str(record)+'.json');download(f'https://zenodo.org/api/records/{record}',path);data=json.loads(path.read_text())
        entries=[dict(name=f['key'],size=f['size'],checksum=f.get('checksum'),url=f['links']['self']) for f in data['files']]
        records.append(dict(record=record,files=entries,metadata_sha256=digest(path)))
        for f in entries:print('SOURCE',record,f['name'],f['size'],flush=True)
    path=out/'littledalei_figshare.json';download(f'https://api.figshare.com/v2/articles/{manifest["littledalei_annotation_record"]}/versions/1',path)
    data=json.loads(path.read_text());records.append(dict(figshare=manifest['littledalei_annotation_record'],files=data['files'],metadata_sha256=digest(path)))
    write_json(out/'inventory.json',records)
    print('INVENTORY',out/'inventory.json',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','pair','collect','inventory']);p.add_argument('--run',type=Path,required=True);p.add_argument('--manifest',type=Path,required=True);p.add_argument('--case',type=int,choices=[0,1],default=0);p.add_argument('--index',type=int);p.add_argument('--gffread',type=Path);a=p.parse_args();manifest=json.loads(a.manifest.read_text());case=manifest['cases'][a.case]
    if a.stage=='prepare':prepare_species(a.run,case,a.index,a.gffread)
    elif a.stage=='pair':
        for d in case['species']:
            if not (a.run/'sources'/d['id']/'COMPLETE').exists():raise ValueError('Unvalidated species input')
        pair(a.run,case,a.index)
    elif a.stage=='collect':collect(a.run,case)
    else:inventory(a.run,manifest)
