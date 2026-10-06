"""Independent method/budget tasks with identical frozen saved incumbents."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys
import time
from benchmark_pipeline import atomic

METHODS=('pipeline','milp','lns','sdp')
BUDGETS=(30,150,600)


def prepare(args):
    from benchmark_bound_pair import choose_saved
    from syntangle import load_validation_bundle
    root=Path(args.root);output=Path(args.output)
    entries=list(csv.DictReader((root/'benchmark_manifest.tsv').open(),delimiter='\t'))
    for entry in entries:
        fixture=load_validation_bundle(root/entry['case_dir'])
        score,source,state=choose_saved(fixture,entry['case_id'],Path(args.history))
        atomic(output/(entry['case_id']+'.warm.json'),dict(score=score,source=source,state=state))
        print(entry['case_id'],'FROZEN_C='+str(score),flush=True)


def worker(args):
    from syntangle import load_validation_bundle, score_crossings
    from syntangle.saved_layout import decode_saved_layout
    from syntangle.global_experiment import run_experiment
    fixture=load_validation_bundle(args.case)
    warm=json.loads(Path(args.warm).read_text())
    state=decode_saved_layout(fixture,warm['state'])
    started=time.perf_counter()
    def checkpoint(state):
        score=score_crossings(fixture,state).crossings
        atomic(Path(args.result).with_suffix('.incumbent.json'),dict(
            crossings=score,state=state.to_dict(),seconds=time.perf_counter()-started))
    checkpoint(state)
    result=run_experiment(fixture,state,args.method,args.seconds,args.seed,checkpoint)
    atomic(Path(args.result),result)


def task(args):
    root=Path(args.root);output=Path(args.output)
    entries=list(csv.DictReader((root/'benchmark_manifest.tsv').open(),delimiter='\t'))
    entry=entries[args.index%9];method=METHODS[(args.index//9)%4];budget=BUDGETS[args.index//36]
    destination=output/entry['case_id']/f'{method}_{budget}'
    destination.mkdir(parents=True,exist_ok=False)
    warm_path=output/(entry['case_id']+'.warm.json')
    warm=json.loads(warm_path.read_text());result_path=destination/'result.json'
    report=dict(case_id=entry['case_id'],method=method,budget=budget,starting_crossings=warm['score'],
                starting_source=warm['source'],status='started')
    atomic(destination/'task.json',report)
    if method=='pipeline':
        state_path=destination/'state.json';atomic(state_path,warm['state'])
        command=[sys.executable,str(Path(__file__).with_name('benchmark_pipeline.py')),
                 '--worker','--case',str(root/entry['case_dir']),'--output',str(result_path),
                 '--warm-layout',str(state_path),'--seed',entry['seed'],'--nodes','100000000',
                 '--solve-seconds',str(budget)]
    else:
        command=[sys.executable,str(Path(__file__).resolve()),'--worker','--case',str(root/entry['case_dir']),
                 '--result',str(result_path),'--warm',str(warm_path),'--method',method,
                 '--seconds',str(budget),'--seed',entry['seed']]
    started=time.perf_counter()
    with (destination/'worker.log').open('w') as log:
        try:
            run=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,timeout=budget+20)
            report['status']='complete' if run.returncode==0 else f'failed:{run.returncode}'
        except subprocess.TimeoutExpired:
            report['status']='wall cutoff; final bounds unavailable'
    report['wall_seconds']=time.perf_counter()-started
    if report['status']=='complete':
        result=report['result']=json.loads(result_path.read_text())
        if result['upper_bound']>warm['score']:raise AssertionError('Frozen incumbent lost')
        print(f"{entry['case_id']} {method} budget={budget} START_C={warm['score']} "
              f"C={result['upper_bound']} lower={result['lower_bound']} "
              f"gap={result['optimality_gap']} wall={report['wall_seconds']:.3f}s",flush=True)
        if method=='sdp':
            print('SDP_NUMERICAL_COMPONENT_ESTIMATES='+str([
                d.get('numerical_relaxation_objective') for d in result['component_diagnostics']]),flush=True)
    else:
        checkpoint=result_path.with_suffix('.incumbent.json')
        if checkpoint.exists():report['incumbent']=json.loads(checkpoint.read_text())
        print(entry['case_id'],method,'budget='+str(budget),report['status'],
              'saved_C='+str(report.get('incumbent',{}).get('crossings')),flush=True)
    baseline=root/'comparison'/entry['case_id']/'results.tsv'
    if baseline.exists():
        report['genespace_baselines']=[r for r in csv.DictReader(baseline.open(),delimiter='\t')
                                     if r['method'] in ('GENESPACE','GENESPACE_plus_flips')]
    atomic(destination/'task.json',report)
    if report['status'].startswith('failed:'):raise RuntimeError('Worker failed; inspect '+str(destination/'worker.log'))


def collect(args):
    output=Path(args.output)
    lines=['# Global method comparison','',
           '| Case | Seconds allowed | Method | Starting C | Final/saved C | Lower | Status | Wall seconds |',
           '|---|---:|---|---:|---:|---:|---|---:|']
    records=[]
    for path in sorted(output.glob('*/*/task.json')):
        row=json.loads(path.read_text());records.append(row)
        result=row.get('result',{});incumbent=row.get('incumbent',{})
        lines.append(f"| {row['case_id']} | {row['budget']} | {row['method']} | {row['starting_crossings']} | "
                     f"{result.get('upper_bound',incumbent.get('crossings','missing'))} | "
                     f"{result.get('lower_bound','unavailable')} | {row['status']} | {row.get('wall_seconds',0):.2f} |")
    lines+=['','SDP estimates are numerical diagnostics, not certified bounds. Neighborhood bounds are not global bounds.',
            f'Expected 108 tasks; found {len(records)} records. Missing/failed/cutoff tasks are not wins.',
            'All methods start independently from the same frozen legal layout; prior branch exclusions are not resumed.']
    (output/'summary.md').write_text('\n'.join(lines)+'\n')
    atomic(output/'summary.json',records)
    print('\n'.join(lines))


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--prepare',action='store_true');p.add_argument('--worker',action='store_true');p.add_argument('--collect',action='store_true')
    for name in ('root','history','output','case','warm','result','method'):p.add_argument('--'+name)
    p.add_argument('--index',type=int);p.add_argument('--seconds',type=float);p.add_argument('--seed',type=int,default=1)
    args=p.parse_args()
    if args.prepare:prepare(args)
    elif args.worker:worker(args)
    elif args.collect:collect(args)
    else:task(args)

if __name__=='__main__':main()
