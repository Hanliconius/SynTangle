"""Bind saved real-data layouts to original fixtures and check independent bounds."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

from proof import Evidence, digest, read_json
from structural import report

CASES = ('planarians', 'flatworms_fig4f')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--pipeline', type=Path, required=True)
    p.add_argument('--filter-run', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--task', type=int, required=True)
    a = p.parse_args()
    if not 0 <= a.task < 10:
        raise ValueError('Task must be 0..9')
    case, index = CASES[a.task // 5], a.task % 5
    saved = a.filter_run / case / f'task_{index}' / 'report.json'
    if not (saved.parent / 'COMPLETE').exists():
        raise ValueError('Saved filtering run is incomplete')
    source = json.loads(saved.read_text())
    matches = sorted(path for path in (a.pipeline / 'local_results').glob(f'*/prepared/{case}/fixture.json')
                     if digest(path) == source['fixture_sha256'])
    if not matches:
        raise ValueError(f'Cannot find original full-evidence fixture matching {case} saved hash')
    out = a.output / f'task_{a.task}'
    out.mkdir(parents=True, exist_ok=False)
    fixture, layout = out / 'fixture.json', out / 'layout.json'
    shutil.copy2(matches[0], fixture)
    layout.write_text(json.dumps(source['best_state'], indent=2) + '\n')
    evidence = Evidence(read_json(fixture))
    measured = report(evidence, read_json(layout))
    if measured['upper'] != source['final']['crossings']:
        raise ValueError('Independent rational score differs from saved production score')
    certificate = out / 'certificate.json'
    checker = Path(__file__).with_name('structural.py')
    subprocess.run([sys.executable, str(checker), 'check', str(fixture), str(layout), str(certificate)],
                   check=True, stdout=subprocess.DEVNULL)
    replayed = False
    if measured['status'] == 'EXACT_STRUCTURAL':
        subprocess.run([sys.executable, str(checker), 'verify', str(fixture), str(layout), str(certificate)],
                       check=True, stdout=subprocess.DEVNULL)
        replayed = True
    result = dict(case=case, filter_index=index, independent_upper=measured['upper'],
                  independent_lower=measured['lower'], independent_gap=measured['gap'],
                  status=measured['status'], replay_verified=replayed,
                  production_lower_not_assumed=source['full_lower_bound'],
                  fixture_sha256=digest(fixture), layout_sha256=digest(layout),
                  source_report_sha256=digest(saved), original_fixture=str(matches[0]),
                  scope='Original full imported evidence; whole-chromosome permutations/flips; unweighted strict midpoint crossings')
    (out / 'audit.json').write_text(json.dumps(result, indent=2) + '\n')
    (out / 'COMPLETE').write_text('PASS\n')
    print(f"{case} filter_index={index} U={measured['upper']} L={measured['lower']} gap={measured['gap']} {measured['status']} replay={replayed}")


if __name__ == '__main__':
    main()
