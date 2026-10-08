"""Reconstruct Fig. 4b chromosome order and compare already saved layouts."""
import argparse
import hashlib
import json
from pathlib import Path
from syntangle.fixtures import load_fixture
from syntangle.layout import LayoutState, score_crossings
from syntangle.saved_layout import decode_saved_layout
from plot_comparison import render, atomic_json

ORDERS = {
    'schMedS3h1': ['chr1', 'chr2', 'chr3', 'chr4'],
    'schMedS3h2': ['chr1', 'chr2', 'chr3', 'chr4'],
    'schPol2': ['chr3', 'chr4', 'chr1', 'chr2'],
    'schNov1': ['chr3', 'chr1', 'chr2'],
    'schLug1': ['chr1', 'chr2', 'chr3', 'chr4'],
}

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--run', type=Path, required=True)
    a = p.parse_args()
    prepared = a.run / 'prepared/planarians'
    saved = a.run / 'figures/planarians'
    output = a.run / 'figures/planarians_paper_order'
    fixture = load_fixture(prepared / 'fixture.json')
    refs = {(r.species_id, r.chromosome_id): r for r in fixture.chromosome_refs}
    expected = {(sp, c) for sp, cs in ORDERS.items() for c in cs}
    if set(refs) != expected:
        raise ValueError('Fixture differs from the 19 labelled chromosomes in Fig. 4b')
    paper = LayoutState(
        {sp: tuple(refs[sp, c] for c in ORDERS[sp]) for sp in fixture.species_ids},
        {r: 1 for r in fixture.chromosome_refs})
    audit = json.loads((saved / 'comparison_audit.json').read_text())
    if audit['fixture_sha256'] != hashlib.sha256((prepared / 'fixture.json').read_bytes()).hexdigest():
        raise ValueError('Saved comparison uses different evidence')
    result = json.loads((saved / 'result.json').read_text())
    panels = [(paper, 'Paper Fig. 4b: reconstructed order',
               'labels verified; native orientations assumed')]
    for item in audit['panels'][1:]:
        state = decode_saved_layout(fixture, item['state'])
        if score_crossings(fixture, state).crossings != item['crossings']:
            raise ValueError('Saved panel score disagrees with current scorer')
        panels.append((state, item['title'], 'saved comparison; no solver rerun'))
    if len(panels) != 4:
        raise ValueError('Need a completed four-panel comparison')
    provenance = json.loads((prepared / 'provenance.json').read_text())
    # Colour-only assignments use directly deposited reference blocks; never alter evidence.
    import csv
    from collections import defaultdict
    intervals = defaultdict(list)
    with (prepared / 'published_block_coordinates.tsv').open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            for i, j in ((1, 2), (2, 1)):
                if row[f'genome{i}'] == 'schMedS3h1':
                    lo, hi = sorted((float(row[f'startBp{j}']), float(row[f'endBp{j}'])))
                    intervals[row[f'genome{j}'], row[f'chr{j}']].append((lo, hi, row[f'chr{i}']))
    votes = defaultdict(lambda: defaultdict(float))
    for chromosome in fixture.chromosomes:
        for block in chromosome.blocks:
            for lo, hi, refchr in intervals[chromosome.ref.species_id, chromosome.ref.chromosome_id]:
                overlap = max(0, min(hi, block.end) - max(lo, block.start))
                votes[block.homology_id][refchr] += overlap / (block.end - block.start)
    assignments = {}
    for h, v in votes.items():
        ordered = sorted(v.items(), key=lambda item: (-item[1], item[0]))
        if ordered and (len(ordered) == 1 or ordered[0][1] > ordered[1][1]):
            assignments[h] = ordered[0][0]
    provenance['homology_colour_reference'] = assignments
    provenance['colour_reference'] = 'schMedS3h1'
    provenance['colour_note'] = 'Colours: dominant direct source-block overlap with S. mediterranea h1; tied/unassigned links grey. Colour only, no new homology.'
    provenance['paper_order_reconstruction'] = True
    provenance['display_top_to_bottom'] = list(reversed(fixture.species_ids))
    output.mkdir(parents=True, exist_ok=True)
    render(fixture, panels, output, provenance)
    scores = [score_crossings(fixture, state).crossings for state, _, _ in panels]
    report = {
        'paper': 'https://doi.org/10.1038/s41467-024-52380-9',
        'figure': '4b',
        'figure_image': 'https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41467-024-52380-9/MediaObjects/41467_2024_52380_Fig4_HTML.png',
        'author_script': 'https://github.com/Jeremias-Brand/PlanarianGenomeAnalysis/blob/main/scripts/genespace_plotting.R',
        'chromosome_order_evidence': 'Direct reading of labelled chromosomes in Fig. 4b',
        'orientation_evidence': 'Native coordinate directions assumed; author plotting call has no explicit flips, but Fig. 4b has no direction arrows',
        'scope': 'Paper chromosome order reconstructed on identical 585 block-midpoint links. Not an exact recreation of gene-level ribbons, chromosome lengths, gaps or alignment.',
        'fixture_sha256': audit['fixture_sha256'],
        'display_top_to_bottom': provenance['display_top_to_bottom'],
        'panels': [{'title': title, 'crossings': score, 'state': state.to_dict()}
                   for (state, title, _), score in zip(panels, scores)],
        'lower_bound': result['lower_bound'],
        'upper_bound': result['upper_bound'],
        'optimality_status': result['optimality_status'],
        'excess_interval': [max(0, scores[0] - result['upper_bound']),
                            max(0, scores[0] - result['lower_bound'])],
    }
    atomic_json(output / 'paper_order_audit.json', report)
    print('PAPER_ORDER_C=%s GS=%s GS_PLUS_FLIPS=%s ST=%s L=%s STATUS=%s' %
          (*scores, result['lower_bound'], result['optimality_status']), flush=True)
    print(output / 'comparison.pdf', flush=True)
    import zipfile
    archive = a.run.parent / 'SynTangle_planarian_paper_order.zip'
    with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        for f in output.iterdir():
            if f.is_file(): z.write(f, f.name)
        z.write(Path(__file__), 'plot_planarian_paper_order.py')
    print(archive, flush=True)

if __name__ == '__main__':
    main()
