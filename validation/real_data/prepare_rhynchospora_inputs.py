"""Verify deposited Fig. 1a candidates before selecting chromosomes or homology.

This stage produces no optimized PDF. Run its five downloads concurrently on
Pegasus; detailed inventories stay in files rather than filling Slurm logs.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from capture_rhynchospora_fig1a import ASSETS, inventory_sequences
from prepare_published_inputs import download


def verify(path, record):
    if path.stat().st_size != record['bytes']:
        raise ValueError(f"Size mismatch: {path.name}; move/remove the invalid cache before retrying")
    checksum = record['checksum']
    if checksum['type'].upper() != 'MD5':
        raise ValueError('Unexpected deposit checksum algorithm')
    h = hashlib.md5()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    if h.hexdigest() != checksum['value']:
        raise ValueError(f"Checksum mismatch: {path.name}; move/remove the invalid cache before retrying")


def asset_path(run, record):
    group, filename = ASSETS[record['id']]
    if filename != record['filename'] or record['restricted']:
        raise ValueError('Unexpected or restricted source asset')
    return run / 'sources' / group / filename


def prepare_asset(run, record):
    path = asset_path(run, record)
    download(f"https://edmond.mpg.de/api/access/datafile/{record['id']}", path)
    verify(path, record)
    receipt = {'record': record, 'path': str(path), 'verified': True}
    if path.suffix == '.fasta':
        receipt['sequences'] = inventory_sequences(path)
        names = [r['sequence'] for r in receipt['sequences']]
        if len(set(names)) != len(names):
            raise ValueError('Duplicate FASTA sequence identifiers')
    target = run / 'receipts' / f"{record['id']}.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(receipt, indent=2) + '\n')
    print(f"VERIFIED {path.name}", flush=True)


def annotation_inventory(path, lengths):
    counts, sequences = Counter(), Counter()
    with path.open() as f:
        for number, line in enumerate(f, 1):
            if line.startswith('##FASTA'):
                break
            if not line.strip() or line.startswith('#'):
                continue
            row = line.rstrip('\n').split('\t')
            if len(row) != 9:
                raise ValueError(f'Malformed GFF row {number}')
            name, start, end = row[0], int(row[3]), int(row[4])
            if name not in lengths or not 1 <= start <= end <= lengths[name]:
                raise ValueError(f'GFF/FASTA coordinate mismatch at row {number}: {name}')
            counts[row[2]] += 1
            sequences[name] += 1
    if not counts['CDS']:
        raise ValueError('Annotation has no CDS features')
    return {'features': dict(counts), 'features_by_sequence': dict(sequences)}


def collect(run, manifest):
    groups, missing = {}, []
    for record in manifest['assets']:
        receipt = run / 'receipts' / f"{record['id']}.json"
        if not receipt.exists():
            missing.append(record['filename'])
            continue
        data = json.loads(receipt.read_text())
        if data['record'] != record or data.get('verified') is not True:
            raise ValueError('Receipt differs from pinned manifest')
        group = ASSETS[record['id']][0]
        groups.setdefault(group, {})['fasta' if 'sequences' in data else 'gff'] = data
    report = {'paper': manifest['paper'], 'groups': {}, 'missing_downloads': missing,
              'comparison_complete': False, 'original_block_table_available': False,
              'next_requirements': ['Validate deposited sequence identities and their mapping to published chromosome bars.',
                                    'Obtain or transfer an annotation onto the exact breviuscula hap1 FASTA.',
                                    'Recompute homology if the original GENESPACE block table remains unavailable; label this as reconstructed evidence.'],
              'no_scaffolds_joined_or_excluded': True}
    for group, data in sorted(groups.items()):
        if 'fasta' not in data:
            continue
        rows = data['fasta']['sequences']
        lengths = {r['sequence']: r['length'] for r in rows}
        item = {'sequences': rows, 'expected_published_bars': manifest['expected_published_bars'][group],
                'sequence_count_matches_bars': len(rows) == manifest['expected_published_bars'][group],
                'annotation_available': 'gff' in data,
                'published_identity_confirmed': False}
        if 'gff' in data:
            item['annotation'] = annotation_inventory(Path(data['gff']['path']), lengths)
        report['groups'][group] = item
        names = ', '.join(r['sequence'] for r in rows[:6])
        if len(rows) > 6:
            names += f", ... (+{len(rows)-6})"
        print(f"{group}: sequences={len(rows)} paper_bars={item['expected_published_bars']} annotation={item['annotation_available']} names=[{names}]", flush=True)
    (run / 'input_readiness.json').write_text(json.dumps(report, indent=2) + '\n')
    print('AUDIT:', run / 'input_readiness.json', flush=True)
    print('No optimized PDF yet: chromosome mapping and breviuscula annotation remain to be established.', flush=True)
    if missing:
        raise ValueError(f'{len(missing)} downloads failed; inspect the corresponding input task logs')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--run', required=True, type=Path)
    p.add_argument('--manifest', required=True, type=Path)
    p.add_argument('--task', type=int)
    args = p.parse_args()
    manifest = json.loads(args.manifest.read_text())
    run = args.run.resolve()
    run.mkdir(parents=True, exist_ok=True)
    if args.task is None:
        collect(run, manifest)
    else:
        if not 0 <= args.task < len(manifest['assets']):
            raise ValueError('Task index out of range')
        prepare_asset(run, manifest['assets'][args.task])


if __name__ == '__main__':
    main()
