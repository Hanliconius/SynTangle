"""Transfer one annotated L. sinapis donor onto the four paper assemblies."""
import argparse,csv,hashlib,json,shutil,subprocess,zipfile
from collections import Counter,defaultdict
from pathlib import Path
from urllib.parse import unquote
from prepare_published_inputs import download, LEPTIDEA

DONOR='GCF_905404315.1'
NAMES=list(LEPTIDEA)

def save(path,data):
    tmp=path.with_suffix('.partial.json');tmp.write_text(json.dumps(data,indent=2)+'\n');tmp.replace(path)

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for part in iter(lambda:f.read(1024*1024),b''):h.update(part)
    return h.hexdigest()

def unpack(archive,output,required):
    output.mkdir(parents=True,exist_ok=True);files={}
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None:raise ValueError('Corrupt NCBI archive')
        for key,suffix in required.items():
            matches=[n for n in z.namelist() if n.endswith(suffix)]
            if len(matches)!=1:raise ValueError(f'Expected one {suffix} in {archive}: {matches}')
            dest=output/key;tmp=dest.with_suffix('.partial')
            with z.open(matches[0]) as src,tmp.open('wb') as dst:shutil.copyfileobj(src,dst)
            tmp.replace(dest);files[key]={'archive_member':matches[0],'sha256':sha(dest)}
    return files

def genes(path):
    result=[]
    for line in path.read_text().splitlines():
        if not line or line.startswith('#'):continue
        f=line.split('\t')
        if len(f)!=9:raise ValueError('Malformed GFF')
        if f[2]!='gene':continue
        attrs={k:unquote(v) for k,v in (a.split('=',1) for a in f[8].split(';') if '=' in a)}
        result.append(dict(id=attrs['ID'],sequence=f[0],start=int(f[3])-1,end=int(f[4]),strand=f[6],attrs=attrs))
    return result

def prepare(a):
    provenance=json.loads((a.assemblies/'provenance.json').read_text())
    records={r['name']:r for r in provenance['assemblies']}
    if {n:r['accession'] for n,r in records.items()}!=LEPTIDEA:raise ValueError('Paper assembly manifest differs')
    manifest={'donor':DONOR,'paper':'https://doi.org/10.1007/s10577-023-09713-z',
              'scope':'New Liftoff annotation transfer; not the original MAKER annotation or published ribbon evidence',
              'parameters':{'coverage':.9,'identity':.9,'extra_copy_identity':.95},'targets':[]}
    donor=a.run/'donor';donor.mkdir(exist_ok=True)
    url=f'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{DONOR}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=SEQUENCE_REPORT'
    archive=donor/'dataset.zip';download(url,archive)
    manifest['donor_files']=unpack(archive,donor,{'genome.fna':'.fna','annotation.gff':'genomic.gff','sequence_report.jsonl':'sequence_report.jsonl'})
    manifest['donor_archive']={'url':url,'sha256':sha(archive)}
    if not genes(donor/'annotation.gff'):raise ValueError('Donor annotation contains no genes')
    # Prebuild shared donor FASTA index before concurrent Liftoff workers.
    from pyfaidx import Fasta
    indexed=Fasta(str(donor/'genome.fna'))
    missing={g['sequence'] for g in genes(donor/'annotation.gff')}-set(indexed.keys())
    indexed.close()
    if missing:raise ValueError('Donor GFF and FASTA sequence IDs differ')
    for name in NAMES:
        record=records[name];archive=a.assemblies/name/(record['accession']+'.zip')
        if sha(archive)!=record['sha256']:raise ValueError('Cached paper assembly checksum differs')
        files=unpack(archive,a.run/'targets'/name,{'genome.fna':'.fna','sequence_report.jsonl':'sequence_report.jsonl'})
        manifest['targets'].append({'name':name,'accession':record['accession'],'files':files,'archive_sha256':record['sha256']})
    save(a.run/'manifest.json',manifest)
    print('PREPARED original four Leptidea assemblies and annotated donor',flush=True)

def lift(a):
    name=NAMES[a.index];out=a.run/'annotations'/name;out.mkdir(parents=True,exist_ok=True)
    # Each worker gets its own GFF/database; shared reference sequence is preindexed.
    shutil.copy2(a.run/'donor/annotation.gff',out/'donor.gff')
    cmd=['liftoff',str(a.run/'targets'/name/'genome.fna'),str(a.run/'donor/genome.fna'),
         '-g',str(out/'donor.gff'),'-o',str(out/'lifted.gff'),'-u',str(out/'unmapped.txt'),
         '-dir',str(out/'intermediate'),'-p',str(a.threads),'-a','0.9','-s','0.9',
         '-exclude_partial','-copies','-sc','0.95']
    save(out/'command.json',cmd);print('START',name,flush=True);subprocess.run(cmd,check=True,cwd=out)
    if not genes(out/'lifted.gff'):raise ValueError('No genes transferred')
    (out/'COMPLETE').write_text('PASS\n');print('DONE',name,flush=True)

def collect(a):
    donor_ids={g['id'] for g in genes(a.run/'donor/annotation.gff')};accepted={};audits=[]
    output=a.run/'homology';output.mkdir(exist_ok=True)
    for name in NAMES:
        root=a.run/'annotations'/name
        if not (root/'COMPLETE').exists():raise ValueError('Incomplete transfer: '+name)
        sequences={}
        for line in (a.run/'targets'/name/'sequence_report.jsonl').read_text().splitlines():
            r=json.loads(line)
            if r.get('assemblyAccession')!=LEPTIDEA[name]:raise ValueError('Target sequence report accession differs')
            if r.get('role')=='assembled-molecule' and r.get('assignedMoleculeLocationType')=='Chromosome' and r.get('chrName') not in ('MT','M'):
                for k in ('genbankAccession','refseqAccession'):
                    if r.get(k):sequences[r[k]]={'chromosome':r['chrName'],'length':int(r['length'])}
        if not sequences:raise ValueError('No nuclear chromosome report: '+name)
        groups=defaultdict(list);excluded=Counter()
        for g in genes(root/'lifted.gff'):
            identifier=g['id']
            if identifier not in donor_ids:
                base=identifier.rsplit('_',1)[0]
                if base in donor_ids and g['attrs'].get('extra_copy_number','0')!='0':identifier=base
                else:excluded['unknown_donor_ID']+=1;continue
            groups[identifier].append(g)
        keep={}
        for identifier,rows in groups.items():
            if len(rows)!=1 or rows[0]['attrs'].get('extra_copy_number','0')!='0':
                excluded['detected_extra_or_duplicate_copy']+=1;continue
            g=rows[0];attrs=g['attrs']
            if float(attrs.get('coverage','0'))<.9 or float(attrs.get('sequence_ID','0'))<.9:
                excluded['low_or_missing_mapping_quality']+=1;continue
            if g['sequence'] not in sequences:excluded['outside_nuclear_chromosomes']+=1;continue
            if not 0<=g['start']<g['end']<=sequences[g['sequence']]['length']:raise ValueError('Invalid target coordinates')
            keep[identifier]={**g,**sequences[g['sequence']]}
        accepted[name]=keep
        audit={'assembly':name,'donor_genes':len(donor_ids),'mapped_donor_IDs':len(groups),'accepted_unique_chromosomal':len(keep),
               'excluded':dict(excluded),'gff_sha256':sha(root/'lifted.gff'),
               'unmapped_lines':len((root/'unmapped.txt').read_text().splitlines()),
               'scope':'Unique among detected transfers; donor ID is candidate homology, not independent orthology proof'}
        audits.append(audit);print('ANNOTATION',name,'accepted',len(keep),'of',len(donor_ids),flush=True)
    shared=set.intersection(*(set(accepted[n]) for n in NAMES))
    fields=['donor_gene_id','assembly','sequence','chromosome','chromosome_length','start0','end','strand','coverage','identity','shared_all_four']
    with (output/'candidate_gene_anchors.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(fields)
        for name in NAMES:
            for identifier,g in sorted(accepted[name].items()):
                w.writerow([identifier,name,g['sequence'],g['chromosome'],g['length'],g['start'],g['end'],g['strand'],g['attrs']['coverage'],g['attrs']['sequence_ID'],int(identifier in shared)])
    (output/'shared_all_four_gene_ids.txt').write_text(''.join(h+'\n' for h in sorted(shared)))
    save(output/'annotation_audit.json',{'assemblies':audits,'shared_all_four':len(shared),
          'thresholds':{'coverage':.9,'identity':.9,'extra_copy_identity':.95},
          'next':'Check mapping recovery before collinear block construction; report shared-set exclusions explicitly'})
    if not shared:raise ValueError('No accepted gene shared across four targets')
    (output/'COMPLETE').write_text('PASS\n');print('SHARED_ALL_FOUR',len(shared),'ANCHORS',output/'candidate_gene_anchors.tsv',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=['prepare','lift','collect'],required=True)
    p.add_argument('--run',type=Path,required=True);p.add_argument('--assemblies',type=Path)
    p.add_argument('--index',type=int,choices=range(4));p.add_argument('--threads',type=int,default=8)
    a=p.parse_args();globals()[a.stage](a)
