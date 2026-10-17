import itertools
import random
from pathlib import Path
import unittest
from syntangle import load_fixture, optimize_auto, score_crossings
from syntangle.bounds import AnchorCell, RelaxedCrossingBound
from syntangle.model import ChromosomeRef
from syntangle.orientation_space import OrientationBasis, orientation_basis
from syntangle.saved_layout import decode_saved_layout, validate_saved_layout
from syntangle.layout import LayoutState

ROOT = Path(__file__).resolve().parents[1]


class StrongerBoundTests(unittest.TestCase):
    def test_frustrated_cycle_strict_improvement_and_all_partial_bounds(self):
        refs = tuple(ChromosomeRef('sp', str(i)) for i in range(3))
        basis = OrientationBasis(dict.fromkeys(refs, 1), tuple((r,) for r in refs))
        cells = (AnchorCell(0, refs[0], refs[1], ((0.,0.),(1.,1.))),
                 AnchorCell(0, refs[1], refs[2], ((0.,0.),(1.,1.))),
                 AnchorCell(0, refs[2], refs[0], ((0.,1.),(1.,0.))))
        groups = {ref: i for i, ref in enumerate(refs)}
        old = RelaxedCrossingBound(basis, cells, 0, groups)
        new = RelaxedCrossingBound(basis, cells, 0, groups, cluster_size=3)
        self.assertEqual(old.lower_bound((None,)*3), 0)
        self.assertEqual(new.lower_bound((None,)*3), 1)
        for partial in itertools.product((None,0,1), repeat=3):
            optimum = min(old.lower_bound(bits) for bits in itertools.product((0,1), repeat=3)
                          if all(a is None or a == b for a,b in zip(partial,bits)))
            self.assertGreaterEqual(new.lower_bound(partial), old.lower_bound(partial))
            self.assertEqual(new.lower_bound(partial), optimum)

    def test_partitioned_bounds_with_cross_cluster_factors(self):
        rng = random.Random(719)
        refs = tuple(ChromosomeRef('sp', str(i)) for i in range(5))
        basis = OrientationBasis(dict.fromkeys(refs, 1), tuple((r,) for r in refs))
        cells = tuple(AnchorCell(0, *rng.sample(refs, 2),
                      tuple((rng.choice((0., .5, 1.)), rng.choice((0., .5, 1.)))
                            for _ in range(5))) for _ in range(12))
        groups = {ref: i for i, ref in enumerate(refs)}
        old = RelaxedCrossingBound(basis, cells, 0, groups)
        full = {bits: old.lower_bound(bits) for bits in itertools.product((0,1), repeat=5)}
        for size in (1,2,3,5):
            new = RelaxedCrossingBound(basis, cells, 0, groups, cluster_size=size)
            for partial in itertools.product((None,0,1), repeat=5):
                optimum = min(value for bits, value in full.items()
                              if all(a is None or a == b for a,b in zip(partial,bits)))
                self.assertGreaterEqual(new.lower_bound(partial), old.lower_bound(partial))
                self.assertLessEqual(new.lower_bound(partial), optimum)

    def test_saved_incumbent_is_retained(self):
        fixture = load_fixture(ROOT/'examples/fixtures/fusion_chain_closed_cycle.json')
        saved = optimize_auto(fixture, local_restarts=1).layout.optimized_state
        decoded = decode_saved_layout(fixture, saved.to_dict())
        self.assertEqual(saved, decoded)
        result = optimize_auto(fixture, local_restarts=1, starting_state=decoded,
                               time_limit_seconds=.000001, bound_cluster_size=6)
        self.assertLessEqual(result.layout.optimized_score.crossings,
                             score_crossings(fixture, saved).crossings)

    def test_invalid_saved_order_rejected(self):
        fixture = load_fixture(ROOT/'examples/fixtures/fusion_chain_closed_cycle.json')
        basis = orientation_basis(fixture, fixture.chromosome_refs)
        saved = optimize_auto(fixture, local_restarts=1).layout.optimized_state
        orders = dict(saved.chromosome_order)
        species = next(iter(orders)); orders[species] = orders[species] + orders[species][:1]
        with self.assertRaisesRegex(ValueError, 'exactly once'):
            validate_saved_layout(fixture, LayoutState(orders, saved.chromosome_orientation), basis)
