"""Annotation transfer and chromosome-identity diagnostics, not a figure replica."""
import argparse
from collections import Counter
import csv
import itertools
import json
import math
from pathlib import Path
import shutil
import subprocess
from urllib.parse import unquote

from prepare_rhynchospora_inputs import verify, asset_path


def genes(path):
    with path.open() as f:
        for line in f:
            if line.startswith('##FASTA'):
                break
            if line.startswith('#') or not line.strip():
                continue
            row = line.rstrip('\n').split('\t')
            if len(row) != 9:
                raise ValueError('Malformed GFF row')
            if row[2] != 'gene':
                continue
            attrs = dict((k, unquote(v)) for k, v in
                         (a.split('=', 1) for a in row[8].split(';') if '=' in a))
            yield {'id': attrs['ID'], 'sequence': row[0], 'start': int(row[3]),
                   'end': int(row[4]), 'strand': row[6], 'attributes': attrs}


def prepare(a):
    manifest = json.loads((a.inputs / 'code/rhynchospora_input_manifest.json').read_text())
    records = {r['id']: r for r in manifest['assets']}
    audit = json.loads((a.inputs / 'input_readiness.json').read_text())
    if audit['missing_downloads']:
        raise ValueError('Input downloads are incomplete')
    source = a.run / 'sources'
    source.mkdir(exist_ok=True)
    for key, identifier in [('donor.fasta', 333736), ('donor.gff', 333735), ('target.fasta', 336260)]:
        original = asset_path(a.inputs, records[identifier])
        verify(original, records[identifier])
        shutil.copy2(original, source / key)
    names = {r['sequence'] for r in audit['groups']['R_tenuis_REC_h1']['sequences']}
    donor_genes = list(genes(source / 'donor.gff'))
    if not donor_genes or {g['sequence'] for g in donor_genes} - names:
        raise ValueError('Donor genes do not match donor FASTA')
    if len({g['id'] for g in donor_genes}) != len(donor_genes):
        raise ValueError('Duplicate donor gene IDs')
    # Index privately copied FASTAs, avoiding changes to prior runs.
    from pyfaidx import Fasta
    for name in ('donor.fasta', 'target.fasta'):
        fasta = Fasta(str(source / name)); fasta.close()
    geometry = a.published / 'published_chromosome_geometry.tsv'
    with geometry.open() as f:
        bars = [r for r in csv.DictReader(f, delimiter='\t') if r['species'] == 'R_austrobrasiliensis']
    bars.sort(key=lambda r: int(r['published_chromosome_label']))
    if len(bars) != 3:
        raise ValueError('Expected the three audited published bars')
    widths = [float(r['x1_pdf_points']) - float(r['x0_pdf_points']) for r in bars]
    if min(widths) <= 0:
        raise ValueError('Invalid published geometry')
    observed = [w / sum(widths) for w in widths]
    sequences = audit['groups']['R_austrobrasiliensis_paper_deposit']['sequences']
    ranked = []
    for combination in itertools.permutations(sequences, 3):
        lengths = [r['length'] for r in combination]
        fractions = [x / sum(lengths) for x in lengths]
        ranked.append({'names': [r['sequence'] for r in combination],
                       'lengths': lengths, 'fraction_difference_L1': sum(abs(x-y) for x, y in zip(fractions, observed))})
    ranked.sort(key=lambda r: (r['fraction_difference_L1'], r['names']))
    mapping = dict(all_sequences=sequences, published_width_fractions=observed,
                   top_candidates=ranked[:10], confirmed_mapping=None,
                   warning='Length similarity is not chromosome identity. No sequences selected, joined, excluded or renamed.')
    (a.run / 'chromosome_mapping_candidates.json').write_text(json.dumps(mapping, indent=2) + '\n')
    provenance = dict(inputs=str(a.inputs), published_source=str(a.published), source_records=records,
                      donor_gene_count=len(donor_genes), target='Exact companion breviuscula hap1 assembly',
                      annotation_scope='New Liftoff transfer from deposited tenuis REC annotation; not original breviuscula annotation or published GENESPACE evidence',
                      donor_identity_note='Deposited hap1 filename contains Chr1_h2 and Chr2_h1; names preserved unchanged',
                      liftoff_thresholds={'coverage': .5, 'sequence_identity': .5},
                      strict_audit_thresholds={'coverage': .9, 'sequence_identity': .9})
    (a.run / 'manifest.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(f'PREPARED annotation donor: {len(donor_genes)} genes; ranked {len(ranked)} austrobrasiliensis three-sequence candidates.')
    print('Mapping remains unconfirmed; no chromosome evidence modified.')


def lift(a):
    out = a.run / 'annotation'; out.mkdir(exist_ok=True)
    cmd = ['liftoff', str(a.run / 'sources/target.fasta'), str(a.run / 'sources/donor.fasta'),
           '-g', str(a.run / 'sources/donor.gff'), '-o', str(out / 'lifted.gff'),
           '-u', str(out / 'unmapped.txt'), '-dir', str(out / 'intermediate'),
           '-p', str(a.threads), '-a', '0.5', '-s', '0.5', '-exclude_partial']
    (out / 'command.json').write_text(json.dumps(cmd, indent=2) + '\n')
    with (out / 'liftoff.out').open('w') as stdout, (out / 'liftoff.err').open('w') as stderr:
        subprocess.run(cmd, check=True, cwd=out, stdout=stdout, stderr=stderr)
    (out / 'LIFTOFF_COMPLETE').write_text('PASS\n')
    print('Liftoff finished; detailed alignment logs retained in annotation/liftoff.{out,err}')


def collect(a):
    out = a.run / 'annotation'
    if not (out / 'LIFTOFF_COMPLETE').exists():
        raise ValueError('Liftoff did not complete; see its task log and annotation/liftoff.err')
    donor_ids = {g['id'] for g in genes(a.run / 'sources/donor.gff')}
    from capture_rhynchospora_fig1a import inventory_sequences
    lengths = {r['sequence']: r['length'] for r in inventory_sequences(a.run / 'sources/target.fasta')}
    mapped = list(genes(out / 'lifted.gff')); counts = Counter(g['id'] for g in mapped)
    strict, excluded = [], Counter()
    for g in mapped:
        if g['sequence'] not in lengths or not 1 <= g['start'] <= g['end'] <= lengths[g['sequence']]:
            raise ValueError('Transferred annotation is outside exact target coordinates')
        attrs = g['attributes']
        reason = None
        if g['id'] not in donor_ids or counts[g['id']] != 1:
            reason = 'nonunique_or_unrecognized_donor_id'
        elif 'coverage' not in attrs or 'sequence_ID' not in attrs:
            reason = 'missing_transfer_quality'
        elif not all(math.isfinite(float(attrs[k])) and .9 <= float(attrs[k]) <= 1
                     for k in ('coverage', 'sequence_ID')):
            reason = 'below_strict_quality'
        elif attrs.get('partial_mapping', 'False').lower() == 'true':
            reason = 'partial_mapping'
        if reason:
            excluded[reason] += 1
        else:
            strict.append(g)
    (out / 'strict_gene_coordinates.json').write_text(json.dumps(strict, indent=2) + '\n')
    report = dict(donor_genes=len(donor_ids), transferred_gene_records=len(mapped),
                  unique_strict_genes=len(strict), exclusions=dict(excluded),
                  target_genes_by_sequence=dict(Counter(g['sequence'] for g in strict)),
                  original_publication_annotation=False, homology_comparison_complete=False,
                  chromosome_mapping_confirmed=False,
                  note='Transferred genes are candidate anchors, not independently established orthology. All GFF records retained; strict coordinates are a separately audited subset.')
    (a.run / 'annotation_audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print(f"ANNOTATION transferred={len(mapped)} strict_unique={len(strict)} donor={len(donor_ids)}")
    print('ANNOTATION AUDIT:', a.run / 'annotation_audit.json')
    print('CHROMOSOME CANDIDATES:', a.run / 'chromosome_mapping_candidates.json')
    print('No optimized PDF yet; chromosome identity and reconstructed homology must be validated.')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--stage', choices=('prepare', 'lift', 'collect'), required=True)
    p.add_argument('--run', type=Path, required=True)
    p.add_argument('--inputs', type=Path)
    p.add_argument('--published', type=Path)
    p.add_argument('--threads', type=int, default=8)
    args = p.parse_args()
    globals()[args.stage](args)
