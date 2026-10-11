"""Resume cached OrthoFinder searches in a compatible isolated environment."""
import argparse,csv,gzip,hashlib,json,shutil,subprocess
from collections import Counter
from pathlib import Path
from rhynchospora_native_regions import IDS,proteins,sha


def read_ids(path):
    rows={}
    for line in path.read_text().splitlines():
        if not line.strip():continue
        key,value=line.split(':',1);value=value.strip()
        if key.strip() in rows:raise ValueError('Duplicate species ID')
        rows[key.strip()]=Path(value).stem
    if set(rows.values())!=set(IDS) or len(rows)!=3:raise ValueError('Cached searches use different species')
    return rows


def assignment(table,expected,seen):
    groups=[]
    with table.open() as f:
        reader=csv.DictReader(f,delimiter='\t')
        if set(reader.fieldnames or [])!={'Orthogroup',*IDS}:raise ValueError('Unexpected orthogroup columns')
        for row in reader:
            group=[]
            for species in IDS:
                identifiers=[x.strip() for x in row[species].split(',') if x.strip()]
                if any(x not in expected[species] for x in identifiers):raise ValueError('Unknown or wrong-species peptide in orthogroup')
                for identifier in identifiers:
                    if identifier in seen:raise ValueError('Peptide assigned more than once')
                    seen.add(identifier);group.append(identifier)
            if group:groups.append(tuple(sorted(group)))
    return sorted(groups)


def validate_results(result,peptides):
    expected={sp:{tid for tid,seq in proteins(peptides/(sp+'.fa'))} for sp in IDS}
    if any(not value for value in expected.values()):raise ValueError('Empty peptide input')
    if sum(map(len,expected.values()))!=len(set.union(*expected.values())):raise ValueError('Peptide IDs overlap across species')
    seen=set();groups=assignment(result/'Orthogroups/Orthogroups.tsv',expected,seen)
    unassigned=assignment(result/'Orthogroups/Orthogroups_UnassignedGenes.tsv',expected,seen)
    if seen!=set.union(*expected.values()):raise ValueError('Incomplete peptide assignment')
    for name in ('Species_Tree/SpeciesTree_rooted.txt','Phylogenetic_Hierarchical_Orthogroups/N0.tsv'):
        path=result/name
        if not path.is_file() or path.stat().st_size==0:raise ValueError('Missing completed tree/hierarchical result: '+name)
    with (result/'Phylogenetic_Hierarchical_Orthogroups/N0.tsv').open() as f:
        reader=csv.DictReader(f,delimiter='\t')
        if not set(IDS)<=set(reader.fieldnames or []):raise ValueError('Hierarchical group table omits species')
        hog_rows=0
        for row in reader:
            for species in IDS:
                if any(g.strip() not in expected[species] for g in row[species].split(',') if g.strip()):raise ValueError('Unknown peptide in hierarchical group table')
            hog_rows+=1
        if not hog_rows:raise ValueError('Empty hierarchical group table')
    tree=(result/'Species_Tree/SpeciesTree_rooted.txt').read_text()
    if any(sp not in tree for sp in IDS):raise ValueError('Species tree omits input species')
    return dict(genes={sp:len(x) for sp,x in expected.items()},groups=len(groups),unassigned_groups=len(unassigned),
                membership_sha256=hashlib.sha256(json.dumps(groups).encode()).hexdigest())


def prepare(prior,run):
    for name in ('bed','peptide'):
        shutil.copytree(prior/'genespace'/name,run/'genespace'/name)
    shutil.copy2(prior/'input_audit.json',run/'input_audit.json')
    candidates=list((prior/'genespace/orthofinder').rglob('WorkingDirectory/SpeciesIDs.txt'))
    if len(candidates)!=1:raise ValueError(f'Expected one cached search directory, found {len(candidates)}')
    cached=candidates[0].parent;ids=read_ids(cached/'SpeciesIDs.txt')
    for i in ids:
        if not (cached/f'Species{i}.fa').is_file():raise ValueError('Missing cached FASTA')
        for j in ids:
            paths=[cached/f'Blast{i}_{j}.txt',cached/f'Blast{i}_{j}.txt.gz']
            if sum(p.is_file() for p in paths)!=1:raise ValueError('Missing or ambiguous cached search pair')
    # Copy rather than link: OrthoFinder may write into its continuation directory.
    seed=run/'seed/WorkingDirectory';shutil.copytree(cached,seed)
    files=[p for p in seed.iterdir() if p.is_file() and (p.name.startswith(('Blast','Species','SequenceIDs')))]
    audit=dict(prior=str(prior),cache=str(cached),searches_reused=True,cache_files=[dict(name=p.name,sha256=sha(p),bytes=p.stat().st_size) for p in sorted(files)],
      peptide_files={p.name:sha(p) for p in sorted((run/'genespace/peptide').glob('*.fa'))},
      cause='OrthoFinder 2.5.5 DendroBLAST array-to-scalar failure; NumPy compatibility repair, not an increased RAM request')
    (run/'repair_inputs.json').write_text(json.dumps(audit,indent=2)+'\n')
    print('REUSING nine cached DIAMOND searches; prior run preserved',flush=True)


def resume(run,threads):
    import numpy,scipy,sys
    if numpy.__version__!='1.26.4':raise ValueError('Repair requires NumPy 1.26.4')
    command=['orthofinder','-b',str(run/'seed/WorkingDirectory'),'-t',str(threads),'-a','1','-X']
    versions=dict(python=sys.version,numpy=numpy.__version__,scipy=scipy.__version__,command=command)
    (run/'repair_versions.json').write_text(json.dumps(versions,indent=2)+'\n')
    with (run/'orthofinder_resume.log').open('w') as log:
        subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    results=list((run/'seed').rglob('Orthogroups/Orthogroups.tsv'))
    if len(results)!=1:raise ValueError('Ambiguous completed result directory')
    result=results[0].parent.parent
    audit=validate_results(result,run/'genespace/peptide')
    # GENESPACE expects search files directly inside the new WorkingDirectory,
    # whereas OrthoFinder continuation may refer back to the copied seed.
    wd=result/'WorkingDirectory'
    for source in (run/'seed/WorkingDirectory').iterdir():
        if source.is_file() and source.name.startswith(('Blast','Species','SequenceIDs')):
            target=wd/source.name
            if target.exists():
                if sha(target)!=sha(source):raise ValueError('Continuation search/ID cache changed')
            else:shutil.copy2(source,target)
    ofdir=run/'genespace/orthofinder';ofdir.mkdir()
    (ofdir/result.name).symlink_to(result.resolve(),target_is_directory=True)
    audit.update(versions=versions,result=str(result),exit_code=0)
    old=Path(json.loads((run/'repair_inputs.json').read_text())['prior'])
    oldtables=list((old/'genespace/orthofinder').glob('Results*/Orthogroups/Orthogroups.tsv'))
    if len(oldtables)==1:
        expected={sp:{tid for tid,seq in proteins(run/'genespace/peptide'/(sp+'.fa'))} for sp in IDS}
        oldgroups=assignment(oldtables[0],expected,set());audit['previous_ordinary_groups_identical']=hashlib.sha256(json.dumps(oldgroups).encode()).hexdigest()==audit['membership_sha256']
    (run/'orthofinder_repair_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    (run/'ORTHOFINDER_COMPLETE').write_text('PASS exit status, complete input assignments, rooted tree and hierarchical results\n')
    print('ORTHOFINDER COMPLETE; ordinary groups identical to prior:',audit.get('previous_ordinary_groups_identical','unknown'),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','resume']);p.add_argument('--prior',type=Path);p.add_argument('--run',type=Path,required=True);p.add_argument('--threads',type=int,default=4);a=p.parse_args()
    if a.stage=='prepare':prepare(a.prior,a.run)
    else:resume(a.run,a.threads)
