import itertools
import random
import unittest
from dataclasses import replace
from unittest.mock import patch

from syntangle.model import Fixture, Chromosome, ChromosomeRef, BlockOccurrence, OrientationConstraint
from syntangle.layout import LayoutState, initial_layout_state, score_crossings
from syntangle.component_space import optimization_component_map
from syntangle.component_presolve import presolve_components
from syntangle.hybrid_experiment import run_hybrid
from syntangle.saved_layout import decode_saved_layout


def strand(n, seed=1):
    rng = random.Random(seed)
    chromosomes = []
    for i in range(n):
        order = list(range(6))
        rng.shuffle(order)
        blocks = tuple(BlockOccurrence(f'{i}:{h}', str(h), j*3, j*3+1, '+') for j,h in enumerate(order))
        chromosomes.append(Chromosome(ChromosomeRef(str(i), 'Z'), 20, 1, 1, blocks))
    return Fixture(1, 'strand', 'test', 'test', tuple(chromosomes), (), {})


class PresolveTests(unittest.TestCase):
    def presolve(self, f, start=None):
        return presolve_components(f, start or initial_layout_state(f), *optimization_component_map(f))

    def test_matches_exhaustive_flips_with_intrinsic_crossings(self):
        for n in (2, 3, 4, 5):
            for seed in range(5):
                f = strand(n, seed)
                start = initial_layout_state(f)
                exact = min(score_crossings(f, LayoutState(start.chromosome_order, dict(zip(f.chromosome_refs, signs)))).crossings
                            for signs in itertools.product((-1, 1), repeat=n))
                state, records = self.presolve(f)
                self.assertEqual(score_crossings(f, state).crossings, exact)
                self.assertEqual(records[0]['component_lower'], exact)
                self.assertTrue(records[0]['frozen'])

    def test_31_rows_never_build_general_model(self):
        f = strand(31)
        with patch('syntangle.hybrid_experiment.JointModel', side_effect=AssertionError('Frozen strand entered MILP')):
            result = run_hybrid(f, initial_layout_state(f), 'milp_reclaim', 1, component_presolve=True)
        self.assertEqual(result['upper_bound'], result['lower_bound'])
        self.assertEqual(len(result['presolve_diagnostics'][0]['transition_tables']), 30)
        self.assertEqual(result['presolve_diagnostics'][0]['status'], 'exact strand dynamic program')

    def test_adjacent_hard_equation_and_nonlocal_deferral(self):
        f = strand(4)
        adjacent = replace(f, orientation_constraints=(OrientationConstraint(f.chromosome_refs[0], f.chromosome_refs[1], 1),))
        start = initial_layout_state(adjacent)
        signs = dict(start.chromosome_orientation); signs[f.chromosome_refs[1]] = -1
        state, records = self.presolve(adjacent, LayoutState(start.chromosome_order, signs))
        self.assertNotEqual(state.chromosome_orientation[f.chromosome_refs[0]], state.chromosome_orientation[f.chromosome_refs[1]])
        self.assertTrue(records[0]['frozen'])
        nonlocal_f = replace(f, orientation_constraints=(OrientationConstraint(f.chromosome_refs[0], f.chromosome_refs[-1], 0),))
        _, records = self.presolve(nonlocal_f)
        self.assertFalse(records[0]['frozen'])

    def test_missing_species_does_not_create_new_adjacency(self):
        f = strand(3)
        middle = replace(f.chromosomes[1], blocks=())
        f = replace(f, chromosomes=(f.chromosomes[0], middle, f.chromosomes[2]))
        state, records = self.presolve(f)
        self.assertEqual(score_crossings(f, state).crossings, 0)
        self.assertEqual(sum(r['component_lower'] for r in records), 0)

    def test_one_bridge_is_not_filtered_to_claim_independence(self):
        f = strand(3)
        extra = Chromosome(ChromosomeRef('0', 'autosome'), 20, 2, 1,
                           (BlockOccurrence('bridge:a', 'bridge', 2, 3, '+'),))
        b = replace(f.chromosomes[1], blocks=f.chromosomes[1].blocks+(BlockOccurrence('bridge:b', 'bridge', 2, 3, '+'),))
        f = replace(f, chromosomes=(f.chromosomes[0], extra, b, f.chromosomes[2]))
        state, records = self.presolve(f)
        self.assertEqual(len(records), 1)
        self.assertFalse(records[0]['frozen'])
        self.assertEqual(records[0]['max_chromosomes_per_row'], 2)
        self.assertEqual(score_crossings(f, state).link_count, score_crossings(f, initial_layout_state(f)).link_count)

    def test_frozen_component_survives_unresolved_component_search(self):
        from syntangle.table_input import load_tables
        from pathlib import Path
        p = Path(__file__).resolve().parents[1]/'examples/quickstart'
        f, _ = load_tables(p/'chromosomes.tsv', p/'links.tsv')
        # Exercise the existing SciPy backend without requiring highspy locally.
        result = run_hybrid(f, initial_layout_state(f), 'milp_reclaim', 2,
                            backend='scipy', component_presolve=True)
        final = decode_saved_layout(f, result['optimized_state'])
        self.assertLessEqual(score_crossings(f, final).crossings, score_crossings(f, initial_layout_state(f)).crossings)
        self.assertLessEqual(result['lower_bound'], result['upper_bound'])


if __name__ == '__main__':
    unittest.main()
