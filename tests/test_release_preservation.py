"""Check actual rendered endpoints against unchanged real-data input tables."""
import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from syntangle.cli import main
from syntangle.genespace import import_layout
from syntangle.layout import initial_layout_state
from syntangle.public_plot import plot_comparison
from syntangle.table_input import load_tables

EXAMPLE = Path(__file__).resolve().parents[1] / 'examples/rhynchospora'


class ReleasePreservationTests(unittest.TestCase):
    def test_all_public_modes_preserve_evidence_and_rendered_endpoints(self):
        fixture, evidence = load_tables(EXAMPLE/'chromosomes.tsv', EXAMPLE/'links.tsv')
        with (EXAMPLE/'links.tsv').open() as handle:
            links = list(csv.DictReader(handle, delimiter='\t'))
        chromosomes = {c.ref: c for c in fixture.chromosomes}
        refs = {(r.species_id, r.chromosome_id): r for r in fixture.chromosome_refs}
        gap = max(sum(c.length for c in fixture.chromosomes if c.ref.species_id == sp)
                  for sp in fixture.species_ids) * .015
        start = initial_layout_state(fixture)
        with tempfile.TemporaryDirectory() as tmp:
            for mode in ('flip-assist', 'none', 'crossing-cap', 'combined'):
                with self.subTest(mode=mode):
                    out = Path(tmp)/mode
                    args = [mode if mode == 'flip-assist' else 'run',
                            '--chromosomes', str(EXAMPLE/'chromosomes.tsv'),
                            '--links', str(EXAMPLE/'links.tsv'), '--output', str(out),
                            '--seconds', '2']
                    if mode != 'flip-assist':
                        args += ['--ribbon-objective', mode, '--ribbon-seconds', '2']
                    self.assertEqual(main(args), 0)
                    self.assertEqual(json.loads((out/'fixture.json').read_text()), evidence)
                    final = import_layout(fixture, out/'layout.tsv')
                    if mode == 'flip-assist':
                        self.assertEqual(final.chromosome_order, start.chromosome_order)
                    stats = json.loads((out/'stats.json').read_text())
                    self.assertEqual(stats['output']['link_count'], 37)
                    self.assertEqual(stats['chromosome_count'], 10)
                    self.assertFalse(stats['independent_certificate'])
                    if mode == 'crossing-cap':
                        ribbon = json.loads((out/'ribbon_solver.json').read_text())
                        self.assertLessEqual(ribbon['final_metrics']['crossings'],
                                             ribbon['starting_metrics']['crossings'])
                    expected = []
                    for state in (start, final):
                        offsets = {}
                        for sp in fixture.species_ids:
                            x = 0
                            for ref in state.chromosome_order[sp]:
                                offsets[ref] = x
                                x += chromosomes[ref].length + gap
                        for row in links:
                            ends = []
                            for side in (1, 2):
                                ref = refs[row[f'species{side}'], row[f'chromosome{side}']]
                                y = -fixture.species_ids.index(ref.species_id)
                                for edge in ('start', 'end'):
                                    x = float(row[f'{edge}{side}'])
                                    if state.chromosome_orientation[ref] == -1:
                                        x = chromosomes[ref].length - x
                                    ends.append((x + offsets[ref], y))
                            expected.append(tuple(ends))
                    # Inspect paths handed to matplotlib, not merely the helper function.
                    from matplotlib.patches import PathPatch
                    captured = []
                    def capture(path, *args, **kwargs):
                        captured.append(tuple(tuple(path.vertices[i]) for i in (0, 7, 3, 4)))
                        return PathPatch(path, *args, **kwargs)
                    with patch('syntangle.public_plot.PathPatch', side_effect=capture):
                        plot_comparison(fixture, start, final, out/'checked.pdf')
                    self.assertCountEqual(captured, expected)
                    self.assertEqual(len(captured), 74)
