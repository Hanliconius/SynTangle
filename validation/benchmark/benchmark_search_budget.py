"""Cached solver budget ladder on unchanged stress inputs; preserve timeout incumbents."""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--index', required=True, type=int)
    args = parser.parse_args()
    if not 0 <= args.index < 18:
        parser.error('index must be 0..17')
    case_index, level = divmod(args.index, 2)
    nodes = (50, 500)[level]
    root = Path(args.root).resolve()
    with (root/'benchmark_manifest.tsv').open() as handle:
        entry = list(csv.DictReader(handle, delimiter='\t'))[case_index]
    output = Path(args.output).resolve()/entry['case_id']/f'nodes_{nodes}'
    output.mkdir(parents=True, exist_ok=False)
    result_path = output/'result.json'
    command = [sys.executable, str(Path(__file__).with_name('benchmark_pipeline.py')),
               '--worker', '--case', str(root/entry['case_dir']),
               '--output', str(result_path), '--seed', entry['seed'], '--nodes', str(nodes)]
    started = time.perf_counter()
    with (output/'solver.log').open('w') as log:
        try:
            run = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=180)
            status = 'complete' if run.returncode == 0 else f'failed:{run.returncode}'
        except subprocess.TimeoutExpired:
            status = 'timeout; solve unresolved'
    report = dict(case_id=entry['case_id'], nodes_per_component=nodes,
                  cutoff_seconds=180, status=status, seconds=time.perf_counter()-started)
    if status == 'complete':
        report['result'] = json.loads(result_path.read_text())
        result = report['result']
        details = result
        print(f"{entry['case_id']} nodes={nodes} COMPLETE C={result['optimized_score']['crossings']} "
              f"lower={details.get('lower_bound')} upper={details.get('upper_bound')} "
              f"gap={details.get('optimality_gap')} {result['optimality_status']} "
              f"{report['seconds']:.3f}s", flush=True)
    else:
        checkpoint = result_path.with_suffix('.incumbent.json')
        if checkpoint.exists():
            report['incumbent'] = json.loads(checkpoint.read_text())
        print(f"{entry['case_id']} nodes={nodes} {status} "
              f"last_saved_incumbent={report.get('incumbent', {}).get('crossings', 'missing')} "
              "final_bound=unavailable; no optimality claim", flush=True)
    saved = root/'comparison'/entry['case_id']/'results.tsv'
    if saved.exists():
        with saved.open() as handle:
            report['prior_genespace'] = [r for r in csv.DictReader(handle, delimiter='\t')
                if r['method'] in ('GENESPACE', 'GENESPACE_plus_flips')]
        for row in report['prior_genespace']:
            print(f"{entry['case_id']} {row['method']} PRIOR_C={row['crossings']}", flush=True)
    (output/'summary.json').write_text(json.dumps(report, indent=2))
    if status.startswith('failed:'):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
