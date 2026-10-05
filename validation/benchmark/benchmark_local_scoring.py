"""Paired, time-limited local-search benchmark on unchanged stress inputs."""
import argparse
import csv
import hashlib
import inspect
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def atomic_json(path, payload):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(payload, indent=2))
    temporary.replace(path)


def worker(args):
    from syntangle import load_validation_bundle, score_crossings
    from syntangle.heuristic import optimize_local_search
    fixture = load_validation_bundle(args.case)
    started = time.perf_counter()
    def checkpoint(state, score):
        # Independently validate every exported incumbent.
        assert score_crossings(fixture, state).crossings == score
        atomic_json(Path(args.output).with_suffix('.incumbent.json'), dict(
            crossings=score, seconds=time.perf_counter()-started,
            state=state.to_dict(), status='heuristic'))
        print(f'INCUMBENT C={score} seconds={time.perf_counter()-started:.4f}', flush=True)
    kwargs = dict(restarts=1, max_improving_steps=args.steps, seed=args.seed)
    if 'progress_callback' in inspect.signature(optimize_local_search).parameters:
        kwargs['progress_callback'] = checkpoint
    result = optimize_local_search(fixture, **kwargs)
    seconds = time.perf_counter()-started
    payload = result.to_dict()
    state = result.layout.optimized_state.to_dict()
    payload.update(seconds=seconds, state_sha256=hashlib.sha256(
        json.dumps(state, sort_keys=True).encode()).hexdigest())
    atomic_json(Path(args.output), payload)
    print(f'DONE C={result.layout.optimized_score.crossings} seconds={seconds:.4f}', flush=True)


def paired(args):
    with (Path(args.root) / 'benchmark_manifest.tsv').open() as handle:
        entry = list(csv.DictReader(handle, delimiter='\t'))[args.index]
    case = Path(args.root).resolve() / entry['case_dir']
    output = Path(args.output).resolve() / entry['case_id']
    output.mkdir(parents=True, exist_ok=False)
    rows = []
    for method, source in [('baseline', Path(args.baseline).resolve()),
                           ('cached', Path(__file__).resolve().parents[2])]:
        env = dict(os.environ, PYTHONPATH=str(source / 'src'), PYTHONUNBUFFERED='1')
        command = [sys.executable, str(Path(__file__).resolve()), '--worker',
                   '--case', str(case), '--output', str(output / f'{method}.json'),
                   '--steps', str(args.steps), '--seed', entry['seed']]
        started = time.perf_counter()
        with (output / f'{method}.log').open('w') as log:
            try:
                completed = subprocess.run(command, env=env, stdout=log,
                    stderr=subprocess.STDOUT, timeout=args.budget)
                status = 'complete' if completed.returncode == 0 else f'failed:{completed.returncode}'
            except subprocess.TimeoutExpired:
                status = 'timeout'
        row = dict(method=method, status=status, process_seconds=time.perf_counter()-started)
        if status == 'complete':
            row['result'] = json.loads((output / f'{method}.json').read_text())
        rows.append(row)
        atomic_json(output / 'paired.json', dict(case_id=entry['case_id'], steps=args.steps,
            restarts=1, budget_seconds=args.budget, methods=rows))
        print(entry['case_id'], method, status, flush=True)
    if all(r['status'] == 'complete' for r in rows):
        a, b = (r['result'] for r in rows)
        identical = (a['state_sha256'] == b['state_sha256'] and a['diagnostics'] == b['diagnostics'])
        speedup = a['seconds'] / b['seconds']
        print(f'IDENTICAL_SEARCH={identical} SPEEDUP={speedup:.3f}x '
              f'BASELINE={a["seconds"]:.4f}s CACHED={b["seconds"]:.4f}s', flush=True)
        if not identical:
            raise AssertionError('Scoring optimization changed search outcome or evaluation counts')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--case')
    parser.add_argument('--root')
    parser.add_argument('--baseline')
    parser.add_argument('--output', required=True)
    parser.add_argument('--index', type=int)
    parser.add_argument('--steps', type=int, default=5)
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--budget', type=int, default=120)
    args = parser.parse_args()
    worker(args) if args.worker else paired(args)

if __name__ == '__main__':
    main()
