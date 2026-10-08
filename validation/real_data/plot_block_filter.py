"""Audit and render saved block-filter results; never invoke a solver."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import zipfile
from syntangle.fixtures import load_fixture
from syntangle.saved_layout import decode_saved_layout, validate_saved_layout
from syntangle.orientation_space import orientation_basis
from benchmark_block_filter import metrics, atomic
from plot_comparison import render, direct_reference_colours


def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--run',type=Path,required=True);p.add_argument('--latest-archive',type=Path,required=True)
    p.add_argument('--case',default='annotated_lep_reanalysis')
    a=p.parse_args();prepared=a.source/'prepared'/a.case
    fixture=load_fixture(prepared/'fixture.json');basis=orientation_basis(fixture,fixture.chromosome_refs)
    sha=hashlib.sha256((prepared/'fixture.json').read_bytes()).hexdigest()
    old=json.loads((a.source/'figures'/a.case/'comparison_audit.json').read_text())
    if old['fixture_sha256']!=sha:raise ValueError('Original comparison fixture differs')
    output=a.run/'visual_summary';output.mkdir(exist_ok=True)
    reports=[];incomplete=[]
    for i in range(5):
        task=a.run/f'task_{i}'
        if not (task/'COMPLETE').exists():
            incomplete.append(i);continue
        d=json.loads((task/'report.json').read_text())
        if d['fixture_sha256']!=sha:raise ValueError('Cannot combine bounds from different evidence')
        for statekey,metrickey in [('best_state','final'),('visual_best_state','visual_best')]:
            state=decode_saved_layout(fixture,d[statekey]);validate_saved_layout(fixture,state,basis)
            actual=metrics(fixture,state)
            if actual['crossings']!=d[metrickey]['crossings']:raise ValueError('Saved score mismatch')
            if abs(actual['span_weighted_crossings']-d[metrickey]['span_weighted_crossings'])>1e-7:
                raise ValueError('Saved weighted metric mismatch')
        reports.append((i,d))
    if not reports:raise ValueError('No completed filtering tasks')
    best_i,best=min(reports,key=lambda item:(item[1]['final']['crossings'],item[0]))
    lower=max(d['full_lower_bound'] for i,d in reports)
    upper=best['final']['crossings']
    if not 0<=lower<=upper:raise ValueError('Inconsistent combined bounds')
    summary={'fixture_sha256':sha,'complete_tasks':[i for i,d in reports],
             'incomplete_tasks':incomplete,'best_task':best_i,'lower_bound':lower,'upper_bound':upper,
             'gap':upper-lower,'optimality_status':'proven optimum' if lower==upper else 'bounded best known',
             'colour_scope':'Bombyx membership only; grey indicates unknown reference assignment, not solver status',
             'tasks':[dict(index=i,filter_fraction=d['filter_fraction'],initial=d['initial'],final=d['final'],
                           visual_best=d['visual_best'],wall_seconds=d['wall_seconds'],
                           component_sizes=[s['component_sizes'] for s in d['stages']],events=d['events']) for i,d in reports]}
    panels=[]
    for original in old['panels'][:3]:
        state=decode_saved_layout(fixture,original['state']);validate_saved_layout(fixture,state,basis)
        m=metrics(fixture,state)
        if m['crossings']!=original['crossings']:raise ValueError('Original panel score mismatch')
        panels.append((state,original['title'],f"all links; span-weighted score={m['span_weighted_crossings']:.2f}"))
    panels.append((decode_saved_layout(fixture,best['best_state']),'SynTangle: staged filtering',
                   f"all links restored; L={lower:,}; U={upper:,}; {summary['optimality_status']}"))
    provenance=dict(old['provenance'])
    colour_audit=(direct_reference_colours(fixture,prepared,provenance)
                  if a.case=='annotated_lep_reanalysis' else
                  {'rule':'Saved source reference membership', 'reference':provenance.get('colour_reference')})
    provenance['independent_row_scaling']=True
    summary['display_scaling']='Independent row scales; proportional chromosome lengths within species; equal total row width'
    atomic(output/'colour_audit.json',colour_audit)
    summary['colour_scope']=provenance.get('colour_note','Saved source reference membership; unassigned links grey')
    summary['original_syntangle_crossings']=min(p['crossings'] for p in old['panels'])
    summary['recovered_saved_best']=upper<=summary['original_syntangle_crossings']
    original_output=output/'saved_comparison';original_output.mkdir(exist_ok=True)
    render(fixture,[(decode_saved_layout(fixture,p['state']),p['title'],'saved full-evidence layout')
                    for p in old['panels']],original_output,provenance)
    render(fixture,panels,output,provenance)
    atomic(output/'summary.json',summary)
    lines=['# Full-evidence filtering comparison','',f'Full-data bound: {lower} <= C* <= {upper}; {summary["optimality_status"]}.','',
           '| Removed initially | Full C | Weighted score | Best weighted score | Wall seconds |',
           '|---|---:|---:|---:|---:|']
    for i,d in reports:
        lines.append(f"| {d['filter_fraction']:.0%} | {d['final']['crossings']} | {d['final']['span_weighted_crossings']:.3f} | {d['visual_best']['span_weighted_crossings']:.3f} | {d['wall_seconds']:.2f} |")
        print('FILTER',d['filter_fraction'],'C',d['final']['crossings'],'WEIGHTED',round(d['final']['span_weighted_crossings'],3),
              'COMPONENT_SIZES',[s['component_sizes'] for s in d['stages']],flush=True)
    lines += ['', 'The weighted score is a diagnostic of block span, not a certified weighted objective.',
              'Filtered certificates never enter the combined full-data bound.',f'Incomplete tasks: {incomplete}.']
    (output/'summary.md').write_text('\n'.join(lines)+'\n')
    shutil.copy2(prepared/'fixture.json',output/'fixture.json')
    shutil.copy2(prepared/'provenance.json',output/'provenance.json')
    archive=a.run/'SynTangle_filter_comparisons.zip';temp=archive.with_suffix('.partial.zip')
    with zipfile.ZipFile(temp,'w',zipfile.ZIP_DEFLATED) as z:
        for f in sorted(output.rglob('*')):
            if f.is_file():z.write(f,'visual_summary/'+str(f.relative_to(output)))
        for i,d in reports:z.write(a.run/f'task_{i}/report.json',f'task_{i}/report.json')
        z.write(Path(__file__),'code/plot_block_filter.py')
        z.write(Path(__file__).with_name('benchmark_block_filter.py'),'code/benchmark_block_filter.py')
    temp.replace(archive)
    a.latest_archive.parent.mkdir(parents=True,exist_ok=True)
    latesttemp=a.latest_archive.with_suffix('.partial.zip');shutil.copy2(archive,latesttemp);latesttemp.replace(a.latest_archive)
    print('PDF:',output/'comparison.pdf');print('ARCHIVE:',a.latest_archive)

if __name__=='__main__':main()
