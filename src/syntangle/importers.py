"""Experimental, provenance-preserving importers for external synteny formats."""
import csv
import hashlib
import json
from pathlib import Path
from .table_input import CHROM_COLUMNS, LINK_COLUMNS, load_tables


def records(path):
    with Path(path).open() as handle:
        for number,line in enumerate(handle,1):
            if line.strip() and not line.startswith('#'):
                yield number,line.strip().split()


def sizes(path, species):
    rows=[];seen=set()
    for number,fields in records(path):
        if len(fields)<2: raise ValueError(f'{path}:{number}: expected chromosome and length')
        name=fields[0];length=int(fields[1])
        if name in seen or length<=0: raise ValueError(f'{path}:{number}: duplicate chromosome or invalid length')
        seen.add(name);rows.append(dict(species=species,chromosome=name,length=length,rank=len(rows)+1,orientation=1))
    if not rows:raise ValueError(f'{path}: no chromosome sizes')
    return rows


def bed(path, lengths):
    genes={}
    for number,fields in records(path):
        if len(fields)<4:raise ValueError(f'{path}:{number}: BED needs chromosome, start, end, gene ID')
        chrom,start,end,gene=fields[:4];start=int(start);end=int(end)
        if gene in genes:raise ValueError(f'{path}:{number}: duplicate gene ID {gene}')
        if chrom not in lengths or not 0<=start<end<=lengths[chrom]:
            raise ValueError(f'{path}:{number}: invalid BED interval or unknown chromosome')
        genes[gene]=(chrom,start,end)
    return genes


def import_paf(args, qlengths, tlengths):
    links=[];dropped=[];total=0
    if args.min_mapq<0 or args.min_mapq>254 or args.min_alignment_length<0:
        raise ValueError('Invalid PAF filter threshold')
    for number,f in records(args.paf):
        total+=1
        if len(f)<12:raise ValueError(f'PAF line {number}: expected at least 12 columns')
        q,ql,qs,qe,strand,t,tl,ts,te,matches,span,mapq=f[:12]
        ql,qs,qe,tl,ts,te,matches,span,mapq=map(int,(ql,qs,qe,tl,ts,te,matches,span,mapq))
        if q not in qlengths or t not in tlengths or qlengths[q]!=ql or tlengths[t]!=tl:
            raise ValueError(f'PAF line {number}: chromosome identity/length differs from supplied sizes')
        if strand not in ('+','-') or not 0<=qs<qe<=ql or not 0<=ts<te<=tl or not 0<=matches<=span or span<=0 or not 0<=mapq<=255:
            raise ValueError(f'PAF line {number}: invalid coordinates, strand or alignment statistics')
        reason=[]
        if args.min_mapq and (mapq==255 or mapq<args.min_mapq):reason.append('mapq_below_threshold_or_unknown')
        if span<args.min_alignment_length:reason.append('alignment_length')
        if reason:
            dropped.append(dict(source_line=number,reason=';'.join(reason)));continue
        # Query direction is positive. A negative relative strand reverses target endpoints.
        a,b=(ts,te) if strand=='+' else (te,ts)
        links.append(dict(zip(LINK_COLUMNS,[f'paf_line_{number}',args.query_species,q,qs,qe,args.target_species,t,a,b])))
    return links,dropped,total


def import_jcvi(args,qlengths,tlengths):
    qgenes=bed(args.query_bed,qlengths);tgenes=bed(args.target_bed,tlengths)
    links=[];total=0
    for number,f in records(args.simple):
        if f[0]=='StartGeneA':continue
        # JCVI allows a display-color prefix such as red*geneA.
        f[0]=f[0].split('*',1)[-1]
        if len(f)!=6:raise ValueError(f'JCVI line {number}: expected six-column .simple, not --rich or --coords output')
        a,b,c,d,score,orientation=f
        if orientation not in ('+','-'):raise ValueError(f'JCVI line {number}: missing explicit + or - orientation')
        if int(score)<0:raise ValueError(f'JCVI line {number}: negative block score')
        try:qa,qb=qgenes[a],qgenes[b];ta,tb=tgenes[c],tgenes[d]
        except KeyError as error:raise ValueError(f'JCVI line {number}: missing BED gene {error.args[0]}') from error
        if qa[0]!=qb[0] or ta[0]!=tb[0]:raise ValueError(f'JCVI line {number}: block endpoints span chromosomes')
        if qa[1]>qb[1] or ta[1]>tb[1]:raise ValueError(f'JCVI line {number}: expected ascending endpoint genes before applying orientation')
        qs,qe=min(qa[1],qb[1]),max(qa[2],qb[2]);ts,te=min(ta[1],tb[1]),max(ta[2],tb[2])
        if orientation=='-':ts,te=te,ts
        total+=1
        links.append(dict(zip(LINK_COLUMNS,[f'jcvi_line_{number}',args.query_species,qa[0],qs,qe,args.target_species,ta[0],ts,te])))
    return links,[],total


def run_import(args):
    out=Path(args.output)
    if out.exists() and (not out.is_dir() or any(out.iterdir())):raise ValueError('Output must be new or empty')
    if not args.query_species or not args.target_species or args.query_species==args.target_species:
        raise ValueError('Supply two distinct species identifiers')
    chromosomes=sizes(args.query_sizes,args.query_species)+sizes(args.target_sizes,args.target_species)
    qlengths={r['chromosome']:r['length'] for r in chromosomes if r['species']==args.query_species}
    tlengths={r['chromosome']:r['length'] for r in chromosomes if r['species']==args.target_species}
    if args.command=='import-paf':
        links,dropped,total=import_paf(args,qlengths,tlengths)
        sources=[args.paf,args.query_sizes,args.target_sizes]
        policy=dict(min_mapq=args.min_mapq,min_alignment_length=args.min_alignment_length,
                    unknown_mapq='retained with min_mapq=0; excluded otherwise',secondary_alignments='retained',representation='one ribbon per retained alignment span; CIGAR not interpreted')
    else:
        links,dropped,total=import_jcvi(args,qlengths,tlengths)
        sources=[args.simple,args.query_bed,args.target_bed,args.query_sizes,args.target_sizes]
        policy=dict(representation='one ribbon per simple block; envelope of endpoint BED genes; explicit block orientation',filtering='none',colors='JCVI highlight colors not imported')
    if not links:raise ValueError('No retained links; no dataset written')
    # Validate through the same reader before publishing any output directory.
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        stage=Path(tmp)
        for name,columns,rows in [('chromosomes.tsv',CHROM_COLUMNS,chromosomes),('links.tsv',LINK_COLUMNS,links),('excluded.tsv',['source_line','reason'],dropped)]:
            with (stage/name).open('w') as handle:
                writer=csv.DictWriter(handle,columns,delimiter='\t');writer.writeheader();writer.writerows(rows)
        load_tables(stage/'chromosomes.tsv',stage/'links.tsv')
        audit=dict(schema_version=1,format=args.command,validation_status='experimental_not_validated_on_real_datasets',
                   source_records=total,retained_links=len(links),excluded_records=len(dropped),
                   chromosome_count=len(chromosomes),coordinate_convention='zero-based interval boundaries',
                   policy=policy,input_sha256={str(p):hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in sources})
        (stage/'import_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
        out.mkdir(parents=True,exist_ok=True)
        for path in stage.iterdir():(out/path.name).write_bytes(path.read_bytes())
    print(f'EXPERIMENTAL {args.command}: {len(links)} links; {len(dropped)} excluded; {out.resolve()}')
    return 0
