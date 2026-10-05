"""Short old/new solver comparison, with prior GENESPACE baselines on same inputs."""
import argparse
import csv
import inspect
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def atomic(path, value):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2))
    temporary.replace(path)


def worker(args):
    from syntangle import load_validation_bundle, optimize_auto, score_crossings
    fixture = load_validation_bundle(args.case)
    started = time.perf_counter()
    def checkpoint(state, crossings):
        if score_crossings(fixture, state).crossings != crossings:
            raise AssertionError('Exported incumbent fails canonical scoring')
        atomic(Path(args.output).with_suffix('.incumbent.json'), dict(
            crossings=crossings, seconds=time.perf_counter()-started,
            state=state.to_dict(), status='feasible incumbent; proof pending'))
        print(f'INCUMBENT C={crossings} seconds={time.perf_counter()-started:.4f}', flush=True)
    kwargs = dict(transition_cap_per_component=250000,
                  branch_node_cap_per_component=args.nodes, local_restarts=1,
                  local_max_improving_steps=5, component_workers=1, seed=args.seed)
    if 'progress_callback' in inspect.signature(optimize_auto).parameters:
        kwargs['progress_callback'] = checkpoint
    result = optimize_auto(fixture, **kwargs)
    if score_crossings(fixture, result.layout.optimized_state).crossings != result.layout.optimized_score.crossings:
        raise AssertionError('Final layout fails canonical scoring')
    payload = result.to_dict()
    payload['seconds'] = time.perf_counter()-started
    atomic(Path(args.output), payload)
    print(f'DONE C={result.layout.optimized_score.crossings} status={result.layout.optimality_status} '
          f'gap={result.details.get("optimality_gap", 0)} seconds={payload["seconds"]:.4f}', flush=True)


def paired(args):
    root = Path(args.root).resolve()
    with (root/'benchmark_manifest.tsv').open() as handle:
        entry = list(csv.DictReader(handle, delimiter='\t'))[args.index]
    output = Path(args.output).resolve()/entry['case_id']
    output.mkdir(parents=True, exist_ok=False)
    saved = root/'comparison'/entry['case_id']/'results.tsv'
    baselines = []
    if saved.exists():
        with saved.open() as handle:
            baselines = [r for r in csv.DictReader(handle, delimiter='\t')
                         if r['method'] in ('GENESPACE', 'GENESPACE_plus_flips')]
    methods = []
    for name, source in [('previous', Path(args.baseline).resolve()),
                         ('pipeline', Path(__file__).resolve().parents[2])]:
        started = time.perf_counter()
        env = dict(os.environ, PYTHONPATH=str(source/'src'), PYTHONUNBUFFERED='1')
        command = [sys.executable, str(Path(__file__).resolve()), '--worker',
                   '--case', str(root/entry['case_dir']), '--output', str(output/f'{name}.json'),
                   '--seed', entry['seed'], '--nodes', str(args.nodes)]
        with (output/f'{name}.log').open('w') as log:
            try:
                result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT,
                                        timeout=args.budget)
                status = 'complete' if result.returncode == 0 else f'failed:{result.returncode}'
            except subprocess.TimeoutExpired:
                status = 'timeout'
        row = dict(method=name, status=status, process_seconds=time.perf_counter()-started)
        if status == 'complete':
            row['result'] = json.loads((output/f'{name}.json').read_text())
        checkpoint = output/f'{name}.incumbent.json'
        if checkpoint.exists():
            row['incumbent'] = json.loads(checkpoint.read_text())
        methods.append(row)
        atomic(output/'paired.json', dict(case_id=entry['case_id'], budget_seconds=args.budget,
               branch_nodes=args.nodes, methods=methods, prior_genespace=baselines,
               note=('Exact cache comparison; complete results must match apart from elapsed time.'
                     if args.require_identical else
                     'Heuristic strategy changed; compare score, bounds and elapsed time, not identical search.')))
        final = row.get('result')
        if final:
            print(entry['case_id'], name, status, f'C={final["optimized_score"]["crossings"]}',
                  final['optimality_status'], f'{final["seconds"]:.3f}s', flush=True)
        else:
            print(entry['case_id'], name, status,
                  f'incumbent={row.get("incumbent", {}).get("crossings", "missing")}', flush=True)
    if args.require_identical:
        complete = [row for row in methods if row['status'] == 'complete']
        if len(complete) == 2:
            records = [{key: value for key, value in row['result'].items()
                        if key != 'seconds'} for row in complete]
            identical = records[0] == records[1]
            speedup = complete[0]['result']['seconds'] / complete[1]['result']['seconds']
            print(f'IDENTICAL_RESULT_AND_SEARCH={identical} SPEEDUP={speedup:.3f}x', flush=True)
            if not identical:
                raise AssertionError('Cache changed layout, bounds, reductions or search diagnostics')
        else:
            print('IDENTITY_UNVERIFIED: at least one solve did not complete', flush=True)
    for baseline in baselines:
        print(entry['case_id'], baseline['method'], f'PRIOR_C={baseline["crossings"]}', flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--require-identical', action='store_true')
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--root')
    parser.add_argument('--case')
    parser.add_argument('--baseline')
    parser.add_argument('--output', required=True)
    parser.add_argument('--index', type=int)
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--nodes', type=int, default=1000)
    parser.add_argument('--budget', type=int, default=120)
    args = parser.parse_args()
    worker(args) if args.worker else paired(args)


if __name__ == '__main__':
    main()
