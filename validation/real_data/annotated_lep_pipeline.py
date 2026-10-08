"""Pegasus stages for an explicitly new annotated Lepidoptera comparison."""
import argparse, csv, gzip, hashlib, itertools, json, re, subprocess
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote
from prepare_published_inputs import download, digest


def write_json(p, data):
    p.write_text(json.dumps(data, indent=2)+'\n')


def fasta(path):
    name, pieces = None, []
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if name is not None: yield name, ''.join(pieces)
                name, pieces = line[1:].split()[0], []
            else: pieces.append(line.strip())
    if name is not None: yield name, ''.join(pieces)


def prepare(run, manifest):
    for d in manifest['species']:
        out = run/'sources'/d['id']; out.mkdir(parents=True, exist_ok=True)
        acc = d['accession']; digits = acc.split('_')[1].split('.')[0]
        prefix = acc+'_'+d['assembly_name']
        base = f'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/{digits[:3]}/{digits[3:6]}/{digits[6:]}/{prefix}/'
        files = {}
        for suffix in ['_genomic.gff.gz', '_protein.faa.gz', '_assembly_report.txt']:
            p = out/(prefix+suffix); download(base+p.name, p)
            files[suffix] = dict(url=base+p.name, sha256=digest(p))
        # RefSeq accession and chromosome lengths come from the assembly report.
        chroms = {}
        for line in (out/(prefix+'_assembly_report.txt')).read_text().splitlines():
            if line.startswith('#'): continue
            cols = line.split('\t')
            if len(cols) >= 9 and cols[1] == 'assembled-molecule' and cols[3] == 'Chromosome' and cols[2] != 'MT':
                chroms[cols[6]] = dict(id=cols[2], length=int(cols[8]))
        if not chroms: raise ValueError('No chromosome sequences in '+acc)
        genes, protein_gene, excluded = {}, {}, defaultdict(int)
        with gzip.open(out/(prefix+'_genomic.gff.gz'), 'rt') as f:
            for line in f:
                if line.startswith('#'): continue
                fields = line.rstrip().split('\t')
                if len(fields) != 9: raise ValueError('Malformed GFF')
                seq, kind = fields[0], fields[2]
                if seq not in chroms:
                    excluded['features_outside_nuclear_chromosomes'] += 1; continue
                attrs = dict(x.split('=', 1) for x in fields[8].split(';') if '=' in x)
                match = re.search(r'(?:^|,)GeneID:(\d+)(?:,|$)', unquote(attrs.get('Dbxref','')))
                if not match: continue
                gene = d['id']+'_GeneID_'+match[1]
                start, end = int(fields[3]), int(fields[4])
                if start < 1 or end < start or end > chroms[seq]['length']: raise ValueError('Out-of-bounds GFF')
                if kind == 'gene':
                    if gene in genes: raise ValueError('Duplicated GeneID placement')
                    genes[gene] = dict(species=d['id'], chromosome=chroms[seq]['id'], start=start, end=end)
                if kind == 'CDS' and attrs.get('protein_id'):
                    pid = unquote(attrs['protein_id'])
                    if pid in protein_gene and protein_gene[pid] != gene: raise ValueError('Ambiguous protein/GeneID mapping')
                    protein_gene[pid] = gene
        selected = {}
        for pid, seq in fasta(out/(prefix+'_protein.faa.gz')):
            gene = protein_gene.get(pid)
            if gene not in genes:
                excluded['proteins_without_retained_gene'] += 1; continue
            if gene not in selected or (-len(seq), pid) < (-len(selected[gene][1]), selected[gene][0]):
                selected[gene] = (pid, seq)
        if not selected: raise ValueError('No proteins matched retained gene coordinates')
        with (out/'proteins.faa').open('w') as f:
            for gene in sorted(selected): f.write('>'+gene+'\n'+selected[gene][1]+'\n')
        write_json(out/'genes.json', {gene:genes[gene] for gene in selected})
        write_json(out/'chromosomes.json', chroms)
        write_json(out/'provenance.json', dict(source=d, files=files, selected_genes=len(selected),
            protein_selection='Longest supplied protein per GeneID; protein ID breaks length ties',
            original_protein_ids={g:p for g,(p,s) in selected.items()}, excluded=dict(excluded),
            chromosome_scope='Nuclear assembled chromosomes only; unplaced sequences and MT excluded'))
        print('ANNOTATIONS READY', d['id'], len(selected), 'genes', len(chroms), 'chromosomes', flush=True)


def pair(run, manifest, index):
    a, b = list(itertools.combinations(manifest['species'], 2))[index]
    out = run/'pairs'/str(index); out.mkdir(parents=True, exist_ok=True)
    blast = []
    for query, target in [(a,b),(b,a)]:
        q = run/'sources'/query['id']; t = run/'sources'/target['id']
        db = out/('db_'+target['id'])
        subprocess.run(['makeblastdb','-in',str(t/'proteins.faa'),'-dbtype','prot','-out',str(db)],check=True)
        raw = out/(query['id']+'.blast.tsv')
        subprocess.run(['blastp','-query',str(q/'proteins.faa'),'-db',str(db),'-outfmt','6',
                        '-num_threads','8','-out',str(raw)],check=True)
        hits = defaultdict(list)
        with raw.open() as f:
            for line in f:
                cols = line.rstrip().split('\t'); hits[cols[0]].append(cols)
        for gene in sorted(hits):
            seen = set()
            for cols in sorted(hits[gene],key=lambda c:(-float(c[11]),float(c[10]),c[1],int(c[6]))):
                if cols[1] in seen: continue
                if len(seen) == 5: break
                seen.add(cols[1]); blast.append('\t'.join(cols))
    prefix = out/'synteny'
    with prefix.with_suffix('.blast').open('w') as f: f.write('\n'.join(blast)+'\n')
    all_genes = {}
    for d in [a,b]: all_genes.update(json.loads((run/'sources'/d['id']/'genes.json').read_text()))
    with prefix.with_suffix('.gff').open('w') as f:
        for gene, v in sorted(all_genes.items(),key=lambda x:(x[1]['species'],x[1]['chromosome'],x[1]['start'],x[0])):
            f.write(f'{v["species"]}__{v["chromosome"]}\t{gene}\t{v["start"]}\t{v["end"]}\n')
    subprocess.run(['MCScanX',str(prefix),'-m','10'],check=True)
    rows, current, orientation, blk = [], [], None, None
    def finish():
        if not current: return
        sides = [all_genes[x] for x,y in current], [all_genes[y] for x,y in current]
        row = dict(blkID=f'pair{index}_block{blk}',orient='+' if orientation=='plus' else '-')
        for i, side in enumerate(sides,1):
            if len({(g['species'],g['chromosome']) for g in side}) != 1: raise ValueError('Block crosses chromosome boundary')
            row.update({f'genome{i}':side[0]['species'],f'chr{i}':side[0]['chromosome'],
                        f'startBp{i}':min(g['start'] for g in side),f'endBp{i}':max(g['end'] for g in side)})
        rows.append(row)
    with prefix.with_suffix('.collinearity').open() as f:
        for line in f:
            if line.startswith('## Alignment'):
                finish(); current=[]; blk=line.split(':')[0].split()[-1]; orientation=line.split()[-1]
                if orientation not in ('plus','minus'): raise ValueError('Unknown block orientation')
            elif line.strip() and not line.startswith('#'):
                vals=line.split(':',1)[1].split(); x,y=vals[:2]
                if x not in all_genes or y not in all_genes: raise ValueError('Unknown collinearity gene')
                current.append((x,y))
        finish()
    if not rows: raise ValueError('No syntenic blocks; inspect pair outputs')
    write_json(out/'blocks.json',rows)
    (out/'COMPLETE').write_text('PASS\n')
    print('BLOCKS READY',a['id'],b['id'],len(rows),flush=True)


def collect(run, manifest):
    from syntangle.fixtures import fixture_from_dict
    from prepare_genespace_blocks import natural
    rows=[]
    for i in range(len(list(itertools.combinations(manifest['species'],2)))):
        p=run/'pairs'/str(i)
        if not (p/'COMPLETE').exists(): raise ValueError('Missing completed pair '+str(i))
        rows.extend(json.loads((p/'blocks.json').read_text()))
    out=run/'prepared'/manifest['id'];out.mkdir(parents=True,exist_ok=True)
    with (out/'published_block_coordinates.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    ids=[d['id'] for d in manifest['species']]; blocks={s:defaultdict(list) for s in ids}; edges=[]; audit=[]
    for a,b in zip(ids,ids[1:]):
        selected=[r for r in rows if {r['genome1'],r['genome2']}=={a,b}]
        if not selected: raise ValueError('Missing adjacent block pair')
        audit.append(dict(pair=[a,b],retained=len(selected)))
        for r in selected:
            hid=r['blkID']
            for sp in [a,b]:
                i=1 if r['genome1']==sp else 2
                blocks[sp][r[f'chr{i}']].append(dict(occurrence_id=hid+'_'+sp,homology_id=hid,
                    start=r[f'startBp{i}'],end=r[f'endBp{i}'],strand='?'))
            edges.append(dict(homology_id=hid,source_orientation=r['orient']))
    data=dict(fixture_version=1,id=manifest['id'],title=manifest['id'],purpose=manifest['scope'],orientation_constraints=[],species=[])
    lengths={d['id']:{x['id']:x['length'] for x in json.loads((run/'sources'/d['id']/'chromosomes.json').read_text()).values()} for d in manifest['species']}
    for sp in ids:
        data['species'].append(dict(id=sp,chromosomes=[dict(id=c,length=lengths[sp][c],display_rank=i+1,blocks=blocks[sp][c]) for i,c in enumerate(sorted(blocks[sp],key=natural))]))
    fixture_from_dict(data)
    write_json(out/'fixture.json',data)
    write_json(out/'provenance.json',dict(scope=manifest['scope'],retained_links=len(edges),edges=edges,pair_audit=audit,
        inference_parameters={'blastp':'Default settings; reciprocal protein searches; eight threads',
            'hit_selection':'Five distinct best subjects per query by bitscore, evalue and deterministic ID tie break',
            'mcscanx':'Default settings except -m 10; separately inferred for all six species pairs'},
        tool_versions=json.loads((run/'tool_versions.json').read_text()),
        chromosome_extent='Full nuclear chromosome lengths from NCBI assembly reports',exact_published_display=False,
        source_kind='New BLASTP/MCScanX block inference; not deposited published block evidence',
        row_order=ids, objective='Unweighted adjacent block-midpoint crossings',
        unlinked_chromosomes={s:sorted(set(lengths[s])-set(blocks[s])) for s in ids},
        sources=[json.loads((run/'sources'/d['id']/'provenance.json').read_text()) for d in manifest['species']]))
    print('PREPARED',manifest['id'],len(edges),'scored links',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=['prepare','pair','collect'],required=True)
    p.add_argument('--run',type=Path,required=True);p.add_argument('--manifest',type=Path,required=True);p.add_argument('--index',type=int)
    args=p.parse_args();m=json.loads(args.manifest.read_text())
    if args.stage=='pair':pair(args.run,m,args.index)
    else:globals()[args.stage](args.run,m)
