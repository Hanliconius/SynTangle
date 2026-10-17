"""Same saved incumbent and time/node budgets; compare independent vs coupled bounds."""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import time


from benchmark_bound_pair import choose_saved


def main():
    from syntangle import load_validation_bundle
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', required=True)
    parser.add_argument('--history', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--index', type=int, required=True)
    args = parser.parse_args()
    root = Path(args.root).resolve()
    with (root/'benchmark_manifest.tsv').open() as handle:
        entry = list(csv.DictReader(handle, delimiter='\t'))[args.index]
    fixture = load_validation_bundle(root/entry['case_dir'])
    output = Path(args.output).resolve()/entry['case_id']
    output.mkdir(parents=True, exist_ok=False)
    score, source, state = choose_saved(fixture, entry['case_id'], Path(args.history))
    (output/'warm.json').write_text(json.dumps(state))
    report = dict(case_id=entry['case_id'], starting_crossings=score, starting_source=source,
                  node_cap_per_component=5000, cooperative_deadline_seconds=150, methods=[])
    for cluster in (0, 6):
        result_path = output/f'coupled_{cluster}.json'
        command = [sys.executable, str(Path(__file__).with_name('benchmark_pipeline.py')),
            '--worker', '--case', str(root/entry['case_dir']), '--output', str(result_path),
            '--seed', entry['seed'], '--nodes', '5000', '--solve-seconds', '150',
            '--coupled-bound-size', str(cluster), '--warm-layout', str(output/'warm.json')]
        started = time.perf_counter()
        with (output/f'coupled_{cluster}.log').open('w') as log:
            try:
                run = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=210)
                status = 'complete' if run.returncode == 0 else f'failed:{run.returncode}'
            except subprocess.TimeoutExpired:
                status = 'emergency timeout; bounds unavailable'
        row = dict(coupled_size=cluster, status=status, seconds=time.perf_counter()-started)
        if status == 'complete':
            result = row['result'] = json.loads(result_path.read_text())
            if result['upper_bound'] > score:
                raise AssertionError('Saved feasible incumbent was lost')
            print(f"{entry['case_id']} coupled={cluster} START_C={score} "
                  f"C={result['upper_bound']} lower={result['lower_bound']} "
                  f"gap={result['optimality_gap']} branch_nodes={result['branch_search_nodes_evaluated']} "
                  f"table_entries={result['factor_table_entries_evaluated']} "
                  f"continuation_nodes={result['factor_continuation_nodes_evaluated']} "
                  f"deadline={result['deadline_reached']} {row['seconds']:.3f}s", flush=True)
            roots=[(d['component_id'],d.get('root_independent_lower_bound'),d.get('root_lower_bound'))
                   for d in result['component_diagnostics'] if d['method']=='implicit-branch-and-bound']
            print('ROOT_BOUNDS component,independent,coupled='+str(roots),flush=True)
        else:
            checkpoint = result_path.with_suffix('.incumbent.json')
            if checkpoint.exists():
                row['incumbent'] = json.loads(checkpoint.read_text())
            print(entry['case_id'], 'coupled='+str(cluster), status, flush=True)
        report['methods'].append(row)
        (output/'paired.json').write_text(json.dumps(report, indent=2))
        if status.startswith('failed:'):
            raise RuntimeError('Worker failed; inspect '+str(output))


if __name__ == '__main__':
    main()
