"""Pegasus-only real block-projection comparison; no homology discovery rerun."""
import argparse
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

from syntangle.fixtures import load_fixture
from syntangle.layout import initial_layout_state, score_crossings
from syntangle.hybrid_experiment import run_hybrid
from syntangle.saved_layout import decode_saved_layout
from syntangle.visualize import _unambiguous_links

BENCHMARK = Path(__file__).resolve().parents[1] / 'benchmark'
sys.path.insert(0, str(BENCHMARK))
from compare_genespace import export_native_input, load_native_states, improve_flips, read_tsv, write_tsv
from visual_summary import riparian
import matplotlib.pyplot as plt


def atomic_json(path, value):
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n')
    tmp.replace(path)


def add_deposited_reference_anchors(fixture, prepared, native_input):
    """Supply directly observed nonadjacent blocks, never transitive homology."""
    source = prepared / 'published_block_coordinates.tsv'
    if not source.exists():
        return {'source': 'adjacent scored block links only', 'additional_anchors': 0}
    allowed = {(ref.species_id, ref.chromosome_id) for ref in fixture.chromosome_refs}
    adjacent = {frozenset((a, b)) for a, b in zip(fixture.species_ids, fixture.species_ids[1:])}
    bed = read_tsv(native_input / 'bed.tsv')
    audit = []
    with source.open() as f:
        for line, row in enumerate(csv.DictReader(f, delimiter='\t'), 2):
            a, b = row['genome1'], row['genome2']
            if a == b or frozenset((a, b)) in adjacent: continue
            sides = [(row[f'genome{i}'], row[f'chr{i}']) for i in (1, 2)]
            if not all(side in allowed for side in sides): continue
            hid = f'deposited_nonadjacent_row_{line}'
            for i, (sp, chrom) in enumerate(sides, 1):
                midpoint = (float(row[f'startBp{i}']) + float(row[f'endBp{i}'])) / 2
                bed.append(dict(genome=sp, chr=chrom, ord=midpoint, og=hid,
                                noAnchor='FALSE', isArrayRep='TRUE'))
            audit.append(dict(source_row=line, source_block_id=row['blkID'], genomes=[a, b]))
    write_tsv(native_input / 'bed.tsv', bed)
    prepared_provenance = json.loads((prepared / 'provenance.json').read_text())
    result = {'source': prepared_provenance.get('source_kind', 'Deposited block evidence'),
              'additional_anchors': len(audit), 'additional_blocks': audit,
              'use': 'GENESPACE ordering only; every method rescored on the same adjacent fixture links',
              'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest()}
    atomic_json(native_input / 'reference_anchor_audit.json', result)
    return result


def render(fixture, panels, output, provenance):
    fig, axes = plt.subplots(2, 2, figsize=(19, 11))
    expected = sum(len(_unambiguous_links(fixture, a, b))
                   for a, b in zip(fixture.species_ids, fixture.species_ids[1:]))
    if expected != provenance['retained_links']:
        raise ValueError('Renderer would omit imported links')
    for ax, (state, title, note) in zip(axes.flat, panels):
        riparian(ax, fixture, state, title, note)
        chroms = {c.ref: c for c in fixture.chromosomes}
        gap = max(c.length for c in fixture.chromosomes) * .08
        span = max(sum(chroms[r].length for r in state.chromosome_order[s]) +
                   gap * (len(state.chromosome_order[s]) - 1) for s in fixture.species_ids)
        for y, species in enumerate(fixture.species_ids):
            x = 0
            for ref in state.chromosome_order[species]:
                width = chroms[ref].length
                label = ref.chromosome_id + (' (-)' if state.chromosome_orientation[ref] < 0 else '')
                ax.text((x + width / 2) / span, y - .12, label, ha='center',
                        va='bottom', rotation=60, fontsize=5, zorder=6)
                x += width + gap
    for ax in list(axes.flat)[len(panels):]:
        ax.axis('off')
    pair_audit = provenance['pair_audit']
    if all('orientation_conflicts' in p for p in pair_audit):
        orientation_note = f"Source reciprocal orientation conflicts: {sum(len(p['orientation_conflicts']) for p in pair_audit)}."
    else:
        orientation_note = 'Reciprocal orientation conflicts not audited; source orientations retained in provenance.'
    source_note = provenance.get('source_kind', 'archived GENESPACE block evidence')
    extent_note = provenance.get('chromosome_extent', 'Maximum observed block endpoint; NOT assembly chromosome lengths')
    fig.suptitle(f'{fixture.fixture_id}: {source_note}', x=.12, ha='left', fontsize=15)
    fig.text(.12, .025,
             f'All {expected} imported adjacent-pair links retained in every panel. {orientation_note}\n'
             'C counts block-midpoint crossings. Input uses chromosome-name order; it is NOT the published drawing.\n'
             f'Chromosome lengths: {extent_note}.\n'
             'GS uses block proxy anchors, including source nonadjacent blocks when available (ordering only). All panels share the adjacent-link objective.\n'
             'this is NOT a gene-level GENESPACE rerun. Ribbon spans retain coordinates; source strands remain annotations in provenance.',
             fontsize=8, color='#465363')
    fig.subplots_adjust(left=.12, right=.98, top=.87, bottom=.14, hspace=.4, wspace=.25)
    temporary = output / 'comparison.partial.pdf'
    fig.savefig(temporary, format='pdf')
    temporary.replace(output / 'comparison.pdf')
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--prepared', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--seconds', type=float, default=150)
    a = p.parse_args()
    if a.seconds <= 0:
        p.error('--seconds must be positive')
    a.output.mkdir(parents=True, exist_ok=True)
    fixture = load_fixture(a.prepared / 'fixture.json')
    provenance = json.loads((a.prepared / 'provenance.json').read_text())
    initial = initial_layout_state(fixture)
    panels = [(initial, 'Input: chromosome-name order', 'block projection')]
    # A usable vector figure exists even if a subsequent solver is interrupted.
    render(fixture, panels, a.output, provenance)
    export_native_input(fixture, a.output / 'native_input')
    native_anchor_audit = add_deposited_reference_anchors(fixture, a.prepared, a.output / 'native_input')
    subprocess.run(['Rscript', str(BENCHMARK / 'genespace_native_order.R'),
                    str(a.output / 'native_input'), str(a.output),
                    '--skip-incomplete-variants'], check=True, timeout=180)
    variants = list(load_native_states(fixture, a.output / 'native_orders.tsv'))
    if not variants:
        raise ValueError('No GENESPACE ordering variants')
    score = lambda state: score_crossings(fixture, state).crossings
    native = min(variants, key=lambda v: (score(v[2]), v[0]))
    assisted = [(v, improve_flips(fixture, state, 1 + i)[0])
                for i, (v, metadata, state) in enumerate(variants)]
    flips = min(assisted, key=lambda v: (score(v[1]), v[0]))
    panels += [(native[2], 'GENESPACE ordering function', 'best complete reference/weight; block proxy anchors'),
               (flips[1], 'GENESPACE ordering + our flip assistance', 'fixed chromosome order; heuristic')]
    render(fixture, panels, a.output, provenance)
    start = min([initial, native[2], flips[1]], key=score)
    best = start
    def checkpoint(state):
        nonlocal best
        if score(state) <= score(best):
            best = state
            atomic_json(a.output / 'incumbent.json', {'crossings': score(state), 'state': state.to_dict()})
    checkpoint(start)
    result = run_hybrid(fixture, start, 'hybrid_adaptive', a.seconds,
                        progress=checkpoint, mirror_symmetry=True,
                        neighborhood_min_decisions=512, neighborhood_fraction=.15,
                        neighborhood_max_seconds=30, deduplicate_neighborhoods=True)
    final = decode_saved_layout(fixture, result['optimized_state'])
    if score(final) != result['upper_bound'] or score(final) > score(start):
        raise AssertionError('Final rescore or incumbent retention failed')
    atomic_json(a.output / 'result.json', result)
    panels += [(final, 'SynTangle', f"L={result['lower_bound']:,}; U={result['upper_bound']:,}; {result['optimality_status']}")]
    render(fixture, panels, a.output, provenance)
    atomic_json(a.output / 'comparison_audit.json', {
        'fixture_sha256': hashlib.sha256((a.prepared / 'fixture.json').read_bytes()).hexdigest(),
        'panels': [{'title': title, 'crossings': score(state), 'state': state.to_dict()} for state, title, note in panels],
        'provenance': provenance,
        'native_anchor_audit': native_anchor_audit,
        'native_variant_audit': (a.output / 'native_variant_audit.tsv').read_text(),
        'scope': provenance.get('scope', 'Real block projection; not reproduction of the published drawing or native gene-level workflow')})
    (a.output / 'COMPLETE').write_text('PASS\n')
    print(f"{fixture.fixture_id}: input={score(initial)} GS={score(native[2])} GS+flips={score(flips[1])} ST={score(final)} L={result['lower_bound']} {result['optimality_status']}", flush=True)
    print(a.output / 'comparison.pdf', flush=True)


if __name__ == '__main__':
    main()
