"""Fresh matched-budget searches from a validated saved full-evidence incumbent."""
import argparse,json,shutil,time
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from syntangle.fixtures import load_fixture
from syntangle.saved_layout import decode_saved_layout,validate_saved_layout
from syntangle.orientation_space import orientation_basis
from syntangle.hybrid_experiment import run_hybrid
from syntangle.layout import score_crossings
from preservation_audit import require_audit,evidence_signature,digest
from plot_comparison import render,atomic_json
CASE='leptidea_paper_assemblies'
POLICIES=[('hybrid_adaptive',900),('milp_hint',900),('hybrid_adaptive',1800)]
def run_task(source,run,index):
    require_audit(source/'endpoint_checked')
    fixture=load_fixture(source/'prepared'/CASE/'fixture.json');fixed=evidence_signature(fixture);basis=orientation_basis(fixture,fixture.chromosome_refs)
    audit=json.loads((source/'figures'/CASE/'comparison_audit.json').read_text())
    if audit.get('fixture_sha256')!=digest(source/'prepared'/CASE/'fixture.json'):raise ValueError('Saved comparison belongs to different evidence')
    initial=min(audit['panels'],key=lambda p:p['crossings'])
    state=decode_saved_layout(fixture,initial['state']);validate_saved_layout(fixture,state,basis)
    if score_crossings(fixture,state).crossings!=initial['crossings']:raise ValueError('Invalid saved incumbent')
    method,seconds=POLICIES[index];out=run/f'task_{index}';out.mkdir();best=state
    def checkpoint(candidate):
        nonlocal best
        validate_saved_layout(fixture,candidate,basis)
        if evidence_signature(fixture)!=fixed:raise ValueError('Solver changed fixed evidence')
        score=score_crossings(fixture,candidate).crossings
        if score<=score_crossings(fixture,best).crossings:
            best=candidate;atomic_json(out/'incumbent.json',dict(state=candidate.to_dict(),crossings=score,evidence_sha256=fixed))
    checkpoint(state);print('START',index,method,seconds,'C',initial['crossings'],flush=True)
    result=run_hybrid(fixture,state,method,seconds,progress=checkpoint,mirror_symmetry=True,neighborhood_min_decisions=512,neighborhood_fraction=.15,neighborhood_max_seconds=60,deduplicate_neighborhoods=True)
    final=decode_saved_layout(fixture,result['optimized_state']);checkpoint(final)
    if score_crossings(fixture,final).crossings!=result['upper_bound']:raise ValueError('Result score differs')
    atomic_json(out/'result.json',dict(**result,fixture_sha256=digest(source/'prepared'/CASE/'fixture.json'),evidence_sha256=fixed,requested_seconds=seconds))
    (out/'COMPLETE').write_text('PASS full evidence and legal chromosome-only state\n')
    print('DONE',index,'C',result['upper_bound'],'L',result['lower_bound'],'seconds',result['seconds'],flush=True)
def collect(source,run):
    fixture=load_fixture(source/'prepared'/CASE/'fixture.json');sha=digest(source/'prepared'/CASE/'fixture.json')
    old=json.loads((source/'figures'/CASE/'comparison_audit.json').read_text());previous=json.loads((source/'figures'/CASE/'result.json').read_text())
    if old.get('fixture_sha256')!=sha:raise ValueError('Prior bound belongs to different evidence')
    previous_state=decode_saved_layout(fixture,previous['optimized_state'])
    if score_crossings(fixture,previous_state).crossings!=previous['upper_bound']:raise ValueError('Prior bound/candidate mismatch')
    results=[]
    for i in range(len(POLICIES)):
        out=run/f'task_{i}'
        if not (out/'COMPLETE').exists():raise ValueError('Incomplete extended task '+str(i)+'; no complete archive promoted')
        r=json.loads((out/'result.json').read_text())
        if r['fixture_sha256']!=sha:raise ValueError('Different evidence; bounds cannot be combined')
        results.append(r)
    # Bound transfer requires the same unchanged fixture/objective; conditional
    # neighborhood and filtered bounds are never used here.
    lower=max([previous['lower_bound']]+[r['lower_bound'] for r in results])
    candidates=[(p['crossings'],p['state']) for p in old['panels']]+[(r['upper_bound'],r['optimized_state']) for r in results]
    upper,saved=min(candidates,key=lambda x:x[0]);best=decode_saved_layout(fixture,saved)
    if score_crossings(fixture,best).crossings!=upper or lower>upper:raise ValueError('Combined scores or bounds inconsistent')
    out=run/'figures'/CASE;out.mkdir(parents=True)
    panels=[(decode_saved_layout(fixture,p['state']),p['title'],'Validated fixed block evidence') for p in old['panels'][:3]]
    panels.append((best,'SynTangle: extended searches',f'L={lower}; U={upper}; '+('solver-reported optimum' if lower==upper else 'bounded best known')))
    render(fixture,panels,out,old['provenance']);require_audit(out)
    atomic_json(out/'summary.json',dict(lower_bound=lower,upper_bound=upper,gap=upper-lower,fixture_sha256=sha,prior_lower_bound=previous['lower_bound'],prior_upper_bound=previous['upper_bound'],trials=results,scope='Fresh searches; strongest same-evidence global bound and incumbent retained. No previous search frontier is claimed to resume.'))
    shutil.copy2(source/'prepared'/CASE/'fixture.json',out/'fixture.json')
    archive=run.parent/'SynTangle_Leptidea_extended.zip';temp=archive.with_suffix('.partial')
    with ZipFile(temp,'w',ZIP_DEFLATED) as z:
        for directory in (run/'figures',):
            for f in sorted(directory.rglob('*')):
                if f.is_file():z.write(f,str(f.relative_to(run)))
        for i in range(3):z.write(run/f'task_{i}/result.json',f'task_{i}/result.json')
    temp.replace(archive);print('COMBINED L',lower,'U',upper,'gap',upper-lower,flush=True);print('ARCHIVE:',archive,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['run','collect']);p.add_argument('--source',type=Path,required=True);p.add_argument('--run',type=Path,required=True);p.add_argument('--index',type=int,choices=range(3));a=p.parse_args()
    if a.stage=='run':run_task(a.source,a.run,a.index)
    else:collect(a.source,a.run)
