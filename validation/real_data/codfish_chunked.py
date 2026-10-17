"""Checkpointed query chunks against full protein databases; BLAST defaults unchanged."""
import argparse, hashlib, itertools, json, os, subprocess
from pathlib import Path
import annotated_lep_pipeline as legacy

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()

def records(path):
    record=[]
    for line in path.read_text().splitlines(True):
        if line.startswith('>') and record:
            yield ''.join(record); record=[]
        record.append(line)
    if record: yield ''.join(record)

def valid(raw,marker,key):
    return raw.exists() and marker.exists() and json.loads(marker.read_text())==dict(key=key,sha256=digest(raw))

def setup(run,work,case,chunks):
    work.mkdir(parents=True,exist_ok=True); species=case['species']; tasks=[]
    for s in species:
        source=run/'sources'/s['id']/'proteins.faa'
        if not (source.parent/'COMPLETE').exists(): raise ValueError('Unvalidated input')
        db=work/('db_'+s['id'])
        subprocess.run(['makeblastdb','-in',str(source),'-dbtype','prot','-out',str(db)],check=True)
        seqs=list(records(source)); parts=[[] for _ in range(chunks)]
        for i,r in enumerate(seqs): parts[i%chunks].append(r)
        for i,part in enumerate(parts):
            if not part: raise ValueError('Empty query chunk')
            (work/(s['id']+f'.{i}.faa')).write_text(''.join(part))
    for pair,(a,b) in enumerate(itertools.combinations(species,2)):
        for q,t in [(a,b),(b,a)]:
            for c in range(chunks):
                tasks.append(dict(pair=pair,query=q['id'],target=t['id'],chunk=c,query_sha=digest(work/(q['id']+f'.{c}.faa')),target_sha=digest(run/'sources'/t['id']/'proteins.faa')))
    (work/'tasks.json').write_text(json.dumps(tasks)); (work/'case.json').write_text(json.dumps(case))
    print('READY',len(tasks),'chunks',flush=True)

def search(work,index):
    task=json.loads((work/'tasks.json').read_text())[index]; raw=work/f'{index}.tsv'; marker=work/f'{index}.complete.json'
    key=dict(task=task,blast_version=subprocess.check_output(['blastp','-version'],text=True),options=['-outfmt','6'])
    if valid(raw,marker,key): print('REUSED',index,flush=True); return
    partial=work/f'{index}.partial'; partial.unlink(missing_ok=True)
    subprocess.run(['blastp','-query',str(work/(task['query']+f".{task['chunk']}.faa")),'-db',str(work/('db_'+task['target'])),'-outfmt','6','-num_threads',os.environ.get('SLURM_CPUS_PER_TASK','2'),'-out',str(partial)],check=True)
    partial.replace(raw); tmp=marker.with_suffix('.partial'); tmp.write_text(json.dumps(dict(key=key,sha256=digest(raw)))); tmp.replace(marker)
    print('COMPLETE',index,flush=True)

def blocks(run,work,index):
    case=json.loads((work/'case.json').read_text()); tasks=json.loads((work/'tasks.json').read_text()); selected=[(i,t) for i,t in enumerate(tasks) if t['pair']==index]
    out=run/'pairs'/str(index); out.mkdir(parents=True,exist_ok=True)
    (out/'COMPLETE').unlink(missing_ok=True)
    # Legacy files without success markers are never trusted, even if nonempty.
    expected={}
    for q in sorted({t['query'] for _,t in selected}):
        dest=work/f'pair{index}_{q}.merged.tsv'; temp=dest.with_suffix('.partial')
        with temp.open('wb') as f:
            for i,t in selected:
                if t['query']!=q: continue
                raw=work/f'{i}.tsv'; marker=work/f'{i}.complete.json'
                if not marker.exists(): raise ValueError('Incomplete chunk')
                saved=json.loads(marker.read_text())
                if saved['key']['task']!=t or saved['sha256']!=digest(raw): raise ValueError('Changed chunk')
                with raw.open('rb') as source:
                    import shutil
                    shutil.copyfileobj(source,f)
        temp.replace(dest); expected[q]=dest
    original=subprocess.run
    def reuse(cmd,*args,**kwargs):
        if cmd[0]=='makeblastdb': return subprocess.CompletedProcess(cmd,0)
        if cmd[0]=='blastp':
            query=Path(cmd[cmd.index('-query')+1]).parent.name
            import shutil
            shutil.copyfile(expected[query],cmd[cmd.index('-out')+1])
            return subprocess.CompletedProcess(cmd,0)
        return original(cmd,*args,**kwargs)
    subprocess.run=reuse
    try: legacy.pair(run,case,index)
    finally: subprocess.run=original

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('stage',choices=['setup','search','blocks']); p.add_argument('--run',type=Path,required=True); p.add_argument('--work',type=Path,required=True); p.add_argument('--manifest',type=Path); p.add_argument('--index',type=int); a=p.parse_args()
    if a.stage=='setup': setup(a.run,a.work,json.loads(a.manifest.read_text())['cases'][0],4)
    elif a.stage=='search': search(a.work,a.index)
    else: blocks(a.run,a.work,a.index)
