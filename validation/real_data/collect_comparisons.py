"""Bundle Pegasus comparison output and identify incomplete cases explicitly."""
import argparse
import json
import shutil
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
from pypdf import PdfWriter


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--run', type=Path, required=True)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--latest-archive', type=Path, required=True)
    a = p.parse_args()
    rows = []
    writer = PdfWriter()
    for entry in json.loads(a.manifest.read_text())['datasets']:
        directory = a.run / 'figures' / entry['id']
        complete = (directory / 'COMPLETE').exists()
        pdf = directory / 'comparison.pdf'
        rows.append({'case': entry['id'], 'complete_comparison': complete,
                     'pdf_exists': pdf.exists(),
                     'note': 'Four-panel comparison' if complete else 'Incomplete: any PDF may contain input/baselines only'})
        if complete:
            if not pdf.exists():
                raise ValueError(f'Completed case has no PDF: {entry["id"]}')
            writer.append(str(pdf))
    if any(r['complete_comparison'] for r in rows):
        writer.write(str(a.run / 'SynTangle_real_comparisons.pdf'))
    writer.close()
    (a.run / 'comparison_status.json').write_text(json.dumps(rows, indent=2) + '\n')
    (a.run / 'DOWNLOAD_README.txt').write_text(
        'See comparison_status.json before interpreting any figure.\n'
        'Combined PDF contains completed comparisons only; individual files may be partial.\n'
        'These are real block projections, not reconstructed published layouts or gene-level GENESPACE runs.\n'
        'Completion means the comparison returned and was rescored; consult result.json for optimality bounds.\n')
    archive = a.run / 'SynTangle_real_comparisons.zip'
    with ZipFile(archive, 'w', ZIP_DEFLATED) as z:
        for directory in [a.run / 'figures', a.run / 'prepared']:
            for path in sorted(directory.rglob('*')):
                if path.is_file():
                    z.write(path, str(path.relative_to(a.run)))
        for name in ['SynTangle_real_comparisons.pdf', 'comparison_status.json',
                     'DOWNLOAD_README.txt', 'checkout_head.txt']:
            path = a.run / name
            if path.exists():
                z.write(path, name)
    temporary = a.latest_archive.with_suffix('.partial')
    shutil.copyfile(archive, temporary)
    temporary.replace(a.latest_archive)
    print(json.dumps(rows, indent=2))
    print(f'Completed comparisons: {sum(r["complete_comparison"] for r in rows)}/{len(rows)}')
    print(f'Archive: {archive}')
    print(f'Latest archive: {a.latest_archive}')


if __name__ == '__main__':
    main()
