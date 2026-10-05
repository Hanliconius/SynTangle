import itertools
from pathlib import Path
import random
import unittest
from unittest.mock import patch

from syntangle import load_fixture
from syntangle.bounds import build_relaxed_crossing_bound, _cell_internal_crossings, _orientation_for_ref
from syntangle.incidence import chromosome_node_id
from syntangle.layout import initial_layout_state
from syntangle.orientation_space import orientation_basis
from syntangle.pair_cost import PreparedPairCrossings, pair_component_crossings
from syntangle.heuristic import _component_map
import syntangle.branch_bound as branch

ROOT = Path(__file__).resolve().parents[1]


class BranchCacheTests(unittest.TestCase):
    def setUp(self):
        self.fixture = load_fixture(ROOT/'examples/fixtures/fusion_chain_closed_cycle.json')
        self.nodes = frozenset(chromosome_node_id(c.ref) for c in self.fixture.chromosomes)
        self.basis = orientation_basis(self.fixture, tuple(c.ref for c in self.fixture.chromosomes))

    def test_all_partial_bounds_identical(self):
        bound = build_relaxed_crossing_bound(self.fixture, self.nodes, self.basis)
        for bits in itertools.product((None, 0, 1), repeat=len(self.basis.free_flip_groups)):
            expected = bound.constant_interleaving_lower_bound
            for cell in bound.cells:
                groups = sorted({bound.group_index_by_ref[cell.left_ref], bound.group_index_by_ref[cell.right_ref]})
                unknown = [g for g in groups if bits[g] is None]
                values = []
                for assignment in itertools.product((0, 1), repeat=len(unknown)):
                    temp = dict(zip(unknown, assignment))
                    signs = [_orientation_for_ref(self.basis, bound.group_index_by_ref, ref, bits, temp)
                             for ref in (cell.left_ref, cell.right_ref)]
                    values.append(_cell_internal_crossings(cell.anchors, *signs))
                expected += min(values)
            self.assertEqual(bound.lower_bound(bits), expected)

    def test_random_orders_and_signs_identical(self):
        scorer = PreparedPairCrossings(self.fixture, self.nodes)
        rng = random.Random(34)
        initial = initial_layout_state(self.fixture)
        for _ in range(100):
            orders = {sp: tuple(rng.sample(list(order), len(order)))
                      for sp, order in initial.chromosome_order.items()}
            signs = {ref: rng.choice((1, -1)) for ref in initial.chromosome_orientation}
            for left, right in zip(self.fixture.species_ids, self.fixture.species_ids[1:]):
                self.assertEqual(scorer.score(left, right, orders[left], orders[right], signs),
                    pair_component_crossings(self.fixture, left, right, orders[left], orders[right], signs, self.nodes))

    def test_identical_branch_search_and_context_restored(self):
        components, component_of = _component_map(self.fixture)
        initial = initial_layout_state(self.fixture)
        # Compare complete result records, including all branch/reduction counters.
        for cap in (0, 5, 100):
            cached = branch._solve_branch_component(self.fixture, initial, initial, 0, components[0], component_of, cap)
            with patch.object(branch, '_prepared_edge') as context:
                context.get.return_value = None
                uncached = branch._solve_branch_component(self.fixture, initial, initial, 0, components[0], component_of, cap)
            self.assertEqual(cached, uncached)
            self.assertIsNone(branch._prepared_edge.get())
