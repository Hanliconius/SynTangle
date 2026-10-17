"""Controlled scheduling/backend/hint/neighborhood ablations on frozen starts."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import sys
import time
from statistics import median
from benchmark_pipeline import atomic

WORKER_ENTRY=Path(__file__).resolve()

VARIANTS = {
    'equal_nohint': ('milp_reclaim', 'equal', 'highs'),
    'weighted_nohint': ('milp_reclaim', 'weighted', 'highs'),
    'weighted_hint': ('milp_hint', 'weighted', 'highs'),
    'weighted_three': ('hybrid_3', 'weighted', 'highs'),
    'weighted_adaptive': ('hybrid_adaptive', 'weighted', 'highs'),
    'weighted_scipy': ('milp_reclaim', 'weighted', 'scipy'),
}


def design(entries):
    tasks = []
    for budget, repeats, indices in ((150, 2, range(9)), (600, 1, range(6, 9))):
        for repeat in range(repeats):
            for index in indices:
                for start in ('public', 'pre_hybrid'):
                    for variant in VARIANTS:
                        tasks.append(dict(case_index=index, case_id=entries[index]['case_id'],
                                          budget=budget, repeat=repeat, start=start, variant=variant))
    return tasks


def fingerprint(directory):
    return {str(p.relative_to(directory)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(directory.rglob('*')) if p.is_file()}


def prepare(args):
    from syntangle import load_validation_bundle, score_crossings
    from syntangle.layout import initial_layout_state, LayoutState, canonicalize_component_order
    from syntangle.orientation_space import orientation_basis
    from syntangle.saved_layout import decode_saved_layout, validate_saved_layout
    root = Path(args.root); output = Path(args.output)
    with (root/'benchmark_manifest.tsv').open() as handle:
        entries = list(csv.DictReader(handle, delimiter='\t'))
    if len(entries) != 9: raise ValueError('Expected the unchanged nine-case stress manifest')
    reference=Path(args.reference_run) if args.reference_run else None
    if reference:
        warm_run=reference
        reference_hashes={r['case_id']:r['input_files'] for r in json.loads((reference/'inputs.json').read_text())}
    elif args.warm_run:
        warm_run = Path(args.warm_run)
    else:
        runs = [p for p in Path(args.history).glob('hybrid_methods_*')
                if all((p/(e['case_id']+'.warm.json')).exists() for e in entries)]
        if not runs: raise ValueError('No prior hybrid frozen-start directory found')
        warm_run = min(runs, key=lambda p: (p.stat().st_mtime_ns, str(p)))
    provenance = []
    for entry in entries:
        fixture_dir = root/entry['case_dir']; fixture = load_validation_bundle(fixture_dir)
        basis = orientation_basis(fixture, fixture.chromosome_refs)
        raw = initial_layout_state(fixture)
        public = canonicalize_component_order(fixture, LayoutState(raw.chromosome_order, basis.base_assignment))
        if reference:
            if fingerprint(fixture_dir)!=reference_hashes[entry['case_id']]:
                raise AssertionError('Input fingerprint changed since reference run')
            public_frozen=json.loads((reference/(entry['case_id']+'.public.json')).read_text())
            public=decode_saved_layout(fixture,public_frozen['state'])
            if score_crossings(fixture,public).crossings!=public_frozen['score']:
                raise AssertionError('Reference public score mismatch')
        warm_path = warm_run/(entry['case_id']+('.pre_hybrid.json' if reference else '.warm.json'))
        frozen = json.loads(warm_path.read_text())
        warm = decode_saved_layout(fixture, frozen['state'])
        for name, state in (('public', public), ('pre_hybrid', warm)):
            validate_saved_layout(fixture, state, basis)
            score = score_crossings(fixture, state).crossings
            if name == 'pre_hybrid' and score != frozen['score']:
                raise AssertionError('Historical frozen start score mismatch')
            atomic(output/(entry['case_id']+'.'+name+'.json'), dict(
                score=score, state=state.to_dict(), source=str(warm_path) if name=='pre_hybrid' else
                (str(reference/(entry['case_id']+'.public.json')) if reference else
                 'public chromosome order + propagated hard orientation base assignment')))
            print('FROZEN', entry['case_id'], name, 'C='+str(score), flush=True)
        provenance.append(dict(case_id=entry['case_id'], input_files=fingerprint(fixture_dir),
                               frozen_warm_source=str(warm_path), historical_source=frozen.get('source')))
    atomic(output/'inputs.json', provenance)
    atomic(output/'tasks.json', design(entries))
    print('Prepared',len(design(entries)),'tasks; frozen start source:',warm_run,flush=True)


def worker(args):
    from syntangle import load_validation_bundle, score_crossings
    from syntangle.saved_layout import decode_saved_layout
    from syntangle.hybrid_experiment import run_hybrid
    import scipy, highspy, numpy
    fixture = load_validation_bundle(args.case)
    input_record=Path(args.result).parents[2]/'inputs.json'
    if input_record.exists():
        records=json.loads(input_record.read_text())
        expected=next(r['input_files'] for r in records if r['case_id']==Path(args.result).parents[1].name)
        if fingerprint(Path(args.case))!=expected:raise AssertionError('Input evidence changed after setup')
    frozen = json.loads(Path(args.warm).read_text()); state = decode_saved_layout(fixture, frozen['state'])
    started = time.perf_counter(); trajectory = []
    def checkpoint(candidate):
        score = score_crossings(fixture, candidate).crossings
        if trajectory and score >= trajectory[-1]['crossings']: return
        trajectory.append(dict(seconds=time.perf_counter()-started, crossings=score))
        atomic(Path(args.result).with_suffix('.incumbent.json'), dict(
            crossings=score, state=candidate.to_dict(), trajectory=trajectory))
    checkpoint(state)
    specification=VARIANTS[args.variant]
    method,scheduling,backend=specification[:3]
    options=specification[3] if len(specification)>3 else {}
    if method=='budget_selective':method='hybrid_adaptive' if args.seconds<300 else 'hybrid_3'
    result = run_hybrid(fixture, state, method, args.seconds, args.seed, checkpoint,
                        scheduling=scheduling, backend=backend,**options)
    checkpoint(decode_saved_layout(fixture, result['optimized_state']))
    diagnostics = result['component_diagnostics']
    result['profile'] = dict(
        preparation_seconds=result['preparation_seconds'],
        graph_model_seconds=sum(d.get('graph_preparation_seconds', 0) for d in diagnostics),
        linear_matrix_seconds=sum(d.get('linear_model_build_seconds', 0) for d in diagnostics),
        neighborhood_seconds=sum(d.get('neighborhood_seconds', 0) for d in diagnostics),
        global_seconds=sum(d.get('global_seconds', 0) for d in diagnostics),
        global_solver_seconds=sum(d.get('solver_seconds', 0) for d in diagnostics),
        global_nodes=sum(d.get('nodes', 0) for d in diagnostics),
        mirror_reduced_components=sum(d.get('mirror_reduction',{}).get('applied',False) for d in diagnostics),
        neighborhood_subsolves=sum(d.get('neighborhood_search', {}).get('subsolves', 0) for d in diagnostics),
        neighborhood_solver_seconds=sum(n.get('solver_seconds',0) for d in diagnostics
            for n in d.get('neighborhood_search',{}).get('neighborhoods',[])),
        peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        numpy_version=numpy.__version__, scipy_version=scipy.__version__, highs_version=highspy.Highs().version(),
        incumbent_trajectory=trajectory,
        timing_scope='Inclusive stage measurements overlap; matrix build may occur inside neighborhood/global work. Do not sum all columns.',
        scipy_callback_scope='SciPy variant saves returned incumbents only; direct HiGHS also exports improving callbacks.')
    atomic(Path(args.result), result)


def task(args):
    output = Path(args.output); row = json.loads((output/'tasks.json').read_text())[args.index]
    with (Path(args.root)/'benchmark_manifest.tsv').open() as handle:
        entries = list(csv.DictReader(handle, delimiter='\t'))
    entry = entries[row['case_index']]
    destination = output/row['case_id']/f"{row['variant']}_{row['start']}_{row['budget']}_r{row['repeat']}"
    destination.mkdir(parents=True, exist_ok=False)
    warm = output/(row['case_id']+'.'+row['start']+'.json'); result_path = destination/'result.json'
    report = dict(row, starting_crossings=json.loads(warm.read_text())['score'], status='started')
    atomic(destination/'task.json', report)
    command = [sys.executable, str(WORKER_ENTRY), '--worker', '--case',
               str(Path(args.root)/entry['case_dir']), '--warm', str(warm), '--result', str(result_path),
               '--variant', row['variant'], '--seconds', str(row['budget']), '--seed', entry['seed']]
    started = time.perf_counter()
    with (destination/'worker.log').open('w') as log:
        try:
            run = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=row['budget']+30)
            report['status'] = 'complete' if run.returncode == 0 else f'failed:{run.returncode}'
        except subprocess.TimeoutExpired:
            report['status'] = 'wall cutoff; final bound unavailable'
    report['wall_seconds'] = time.perf_counter()-started
    if report['status'] == 'complete':
        result = report['result'] = json.loads(result_path.read_text())
        if result['upper_bound'] > report['starting_crossings']: raise AssertionError('Frozen incumbent lost')
    else:
        saved = result_path.with_suffix('.incumbent.json')
        if saved.exists(): report['incumbent'] = json.loads(saved.read_text())
    atomic(destination/'task.json', report)
    print(row['case_id'], row['variant'], row['start'], 'budget='+str(row['budget']),
          report['status'], flush=True)
    if report['status'].startswith('failed:'): raise RuntimeError('Worker failed; inspect '+str(destination/'worker.log'))


def collect(args):
    output = Path(args.output); records = []
    lines = ['# Efficiency ablation', '',
             '| Case | Start | Variant | Budget | Repeat | C | Lower | Proven | Prep s | LNS s | Global s | Nodes | RSS MiB | Mirror comps | Status |',
             '|---|---|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---|']
    for p in sorted(output.glob('*/*/task.json')):
        row = json.loads(p.read_text()); records.append(row)
        r = row.get('result', {}); profile = r.get('profile', {})
        upper = r.get('upper_bound', row.get('incumbent', {}).get('crossings', 'missing'))
        lower = r.get('lower_bound', 'unavailable')
        proven = bool(r) and upper == lower
        lines.append(f"| {row['case_id']} | {row['start']} | {row['variant']} | {row['budget']} | {row['repeat']} | "
                     f"{upper} | {lower} | {proven} | {profile.get('preparation_seconds',0):.2f} | "
                     f"{profile.get('neighborhood_seconds',0):.2f} | {profile.get('global_seconds',0):.2f} | "
                     f"{profile.get('global_nodes',0)} | {profile.get('peak_rss_kib',0)/1024:.1f} | "
                     f"{profile.get('mirror_reduced_components',0)} | {row['status']} |")
    audits = []
    for path in sorted((output/'audits').glob('*.json')): audits.extend(json.loads(path.read_text()))
    lines += ['', f'Expected {len(json.loads((output/"tasks.json").read_text()))} tasks; found {len(records)}. Audit case records: {len(audits)}/9.',
              'Completion is not proof; only matching global bounds count as reported optima.',
              'Equal vs weighted allocation both retain proven-zero components and visit small unresolved models first.',
              'The pre_hybrid start is frozen historical evidence, not the latest known optimum.',
              'Repeat counts and fixture subsets are recorded in tasks.json. Timing differences need repeated confirmation.',
              'Stage timings overlap. Incumbent trajectories are improvements, not global bound trajectories.']
    for audit in audits: lines.append(f"AUDIT {audit['case_id']}: {audit['status']}")
    expected_groups={}
    for task in json.loads((output/'tasks.json').read_text()):
        key=(task['budget'],task['start'],task['variant'])
        expected_groups[key]=expected_groups.get(key,0)+1
    groups={}
    for row in records:
        groups.setdefault((row['budget'],row['start'],row['variant']),[]).append(row)
    overview=['# Efficiency overview','',
        '| Budget | Start | Variant | Complete/expected | Proven | Median wall s | Median LNS s | Median global s |',
        '|---:|---|---|---:|---:|---:|---:|---:|']
    for (budget,start,variant),rows in sorted(groups.items()):
        good=[r for r in rows if r['status']=='complete' and 'result' in r]
        expected=expected_groups[(budget,start,variant)]
        proven=sum(r['result']['upper_bound']==r['result']['lower_bound'] for r in good)
        if not good:continue
        overview.append(f"| {budget} | {start} | {variant} | {len(good)}/{expected} | {proven} | "
            f"{median(r['wall_seconds'] for r in good):.2f} | "
            f"{median(r['result']['profile']['neighborhood_seconds'] for r in good):.2f} | "
            f"{median(r['result']['profile']['global_seconds'] for r in good):.2f} |")
    overview += ['', 'Medians describe this fixed mixture of fixtures; compare matched individual cases before ranking methods.', '']
    lines=overview+lines
    (output/'summary.md').write_text('\n'.join(lines)+'\n')
    atomic(output/'summary.json', records); atomic(output/'audit_summary.json', audits)
    print('\n'.join(lines))


def main():
    p = argparse.ArgumentParser()
    for name in ('prepare','worker','collect'): p.add_argument('--'+name, action='store_true')
    for name in ('root','output','history','warm-run','reference-run','case','warm','result','variant'): p.add_argument('--'+name)
    p.add_argument('--index', type=int); p.add_argument('--seconds', type=float)
    p.add_argument('--seed', type=int, default=1)
    args = p.parse_args()
    if args.prepare: prepare(args)
    elif args.worker: worker(args)
    elif args.collect: collect(args)
    else: task(args)

if __name__ == '__main__': main()
