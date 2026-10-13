"""Acquire exact paper deposits; gate comparison on reviewed fixed evidence and layout."""
import argparse, hashlib, json, shutil, subprocess, urllib.request, time
from pathlib import Path

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def fetch(url,path,md5=None):
    path.parent.mkdir(parents=True,exist_ok=True)
    def good():
        if not path.exists():return False
        if md5:
            h=hashlib.md5()
            with path.open('rb') as f:
                for b in iter(lambda:f.read(1048576),b''):h.update(b)
            return h.hexdigest()==md5
        return (path.with_suffix(path.suffix+'.sha256')).exists() and path.with_suffix(path.suffix+'.sha256').read_text().strip()==sha(path)
    if not good():
        temp=path.with_suffix(path.suffix+'.partial');temp.unlink(missing_ok=True)
        try:
            command=['curl','--http1.1','--user-agent','Mozilla/5.0','--fail','--location','--silent','--show-error','--retry','3','--connect-timeout','30','--max-time','1800','--output',str(temp),url]
            for attempt in range(3):
                try:
                    subprocess.run(command,check=True);break
                except subprocess.CalledProcessError as exc:
                    temp.unlink(missing_ok=True)
                    if exc.returncode not in (18,28,35,52,56) or attempt==2:raise
                    time.sleep(2*(attempt+1))
            if md5:
                h=hashlib.md5()
                with temp.open('rb') as f:
                    for b in iter(lambda:f.read(1048576),b''):h.update(b)
                if h.hexdigest()!=md5:raise ValueError('Deposit checksum mismatch')
            temp.replace(path);path.with_suffix(path.suffix+'.sha256').write_text(sha(path)+'\n')
        finally:temp.unlink(missing_ok=True)
    return {'url':url,'path':str(path),'sha256':sha(path),'bytes':path.stat().st_size}

def intake(run,case):
    out=run/case['id'];out.mkdir(parents=True,exist_ok=True); files=[]; errors=[]
    assets=list(case['assets'])
    urls=[('paper.html',case['source_page'])]
    if case.get('inventory_api'):urls.append(('deposit_inventory.json',case['inventory_api']))
    if case.get('code_inventory'):urls.append(('code_inventory.json',case['code_inventory']))
    urls += [(f'supplement_{i}'+Path(u).suffix,u) for i,u in enumerate(case['supplementary_urls'])]
    urls += [(f'figure_{i}'+Path(u).suffix,u) for i,u in enumerate(case['figure_urls'])]
    for name,url in urls:
        try:files.append(fetch(url,out/'sources'/name))
        except Exception as e:errors.append({'url':url,'error':str(e)})
    for asset in assets:
        try:files.append(fetch(asset['url'],out/'sources'/asset['name'],asset['md5']))
        except Exception as e:errors.append({'url':asset['url'],'error':str(e)})
    for accession in case.get('assembly_candidates',[]):
        url=f'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/dataset_report'
        # Metadata only: do not substitute a strain for an unverified figure row.
        try:files.append(fetch(url,out/'sources'/(accession+'.json')))
        except Exception as e:errors.append({'url':url,'error':str(e)})
    report=dict(case=case,files=files,download_errors=errors,status='NEEDS_EVIDENCE_AND_LAYOUT_REVIEW',exact_replication=False,
                remaining=['Identify exact comparator assembly versions','Locate published anchors/blocks or document reconstruction','Verify chromosome order, orientation and displayed exclusions','Create prepared fixture and reviewed layout'])
    (out/'source_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(case['id'],'files=',len(files),'download_errors=',len(errors),'status=',report['status'],flush=True)
    if errors:raise RuntimeError('Source downloads incomplete; cached successes retained')

def gate(prepared):
    from syntangle.fixtures import load_fixture
    # Use the same validated decoder as all real-data comparisons.
    from plot_comparison import decode_saved_layout
    f=load_fixture(prepared/'fixture.json'); review=json.loads((prepared/'source_review.json').read_text())
    for name in ['fixture.json','provenance.json','input_layout.json']:
        if review['sha256'].get(name)!=sha(prepared/name):raise ValueError('Changed reviewed input: '+name)
    if review.get('status')!='PASS' or not review.get('assembly_versions_verified'):raise ValueError('Source review incomplete')
    if review.get('source_kind') not in ['published_evidence','reconstructed_homology']:raise ValueError('Unknown evidence origin')
    if review.get('published_layout_verified') is not True:raise ValueError('Published layout has not been recovered')
    decode_saved_layout(f,json.loads((prepared/'input_layout.json').read_text()))
    return review

def compare(run,case,budget):
    prepared=run/case['id']/'prepared';gate(prepared)
    out=run/case['id']/'comparisons'/str(budget)
    subprocess.run(['python',str(Path(__file__).with_name('plot_comparison.py')),'--prepared',str(prepared),'--output',str(out),'--seconds',str(budget)],check=True)
    from preservation_audit import require_audit
    require_audit(out)
    (out/'COMPLETE').write_text('PASS reviewed evidence and endpoint preservation\n')

def collect(run,cases):
    rows=[]
    for case in cases:
        for seconds in [30,120,600]:
            out=run/case['id']/'comparisons'/str(seconds)
            row=dict(case=case['id'],seconds=seconds,status='NOT_RUN')
            if (out/'COMPLETE').exists():
                from preservation_audit import require_audit
                require_audit(out)
                result=json.loads((out/'result.json').read_text()); row.update(status='COMPLETE',upper=result['upper_bound'],lower=result['lower_bound'],optimality=result['optimality_status'])
            rows.append(row)
    (run/'summary.json').write_text(json.dumps(rows,indent=2)+'\n')
    import tempfile
    with tempfile.TemporaryDirectory(dir=run.parent) as tmp:
        archive=shutil.make_archive(str(Path(tmp)/'bundle'),'zip',run)
        shutil.move(archive,run/'SynTangle_simple_papers.zip')
    print('ARCHIVE:',run/'SynTangle_simple_papers.zip',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['intake','compare','collect']);p.add_argument('--run',type=Path,required=True);p.add_argument('--manifest',type=Path,default=Path(__file__).with_name('simple_papers_manifest.json'));p.add_argument('--index',type=int,default=0);a=p.parse_args();cases=json.loads(a.manifest.read_text())['cases']
    if a.stage=='intake':intake(a.run,cases[a.index])
    elif a.stage=='compare':compare(a.run,cases[a.index//3],[30,120,600][a.index%3])
    else:collect(a.run,cases)
