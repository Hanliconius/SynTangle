import importlib.util
from itertools import permutations
from pathlib import Path
import random
import unittest

from syntangle import load_fixture, optimize_auto
from syntangle.heuristic import _component_map
from syntangle.order_dp import PairwiseOrderCosts
from syntangle.model import ChromosomeRef

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('probe',ROOT/'validation/benchmark/diagnose_order_relaxation.py')
probe=importlib.util.module_from_spec(spec);spec.loader.exec_module(probe)


class OrderProbeTests(unittest.TestCase):
    def test_triangle_packing_is_admissible_exhaustively(self):
        rng=random.Random(143)
        refs=tuple(ChromosomeRef('sp',str(i)) for i in range(5))
        for _ in range(30):
            costs=PairwiseOrderCosts(refs,{(a,b):rng.randrange(9) for a in refs for b in refs if a!=b})
            lower=probe.pair_minimum(costs)+probe.triangle_packing(costs)[0]
            optimum=min(sum(costs.cost(a,b) for i,a in enumerate(order) for b in order[i+1:])
                        for order in permutations(refs))
            self.assertLessEqual(lower,optimum)

    def test_triangle_cycle_strict_lift(self):
        refs=tuple(ChromosomeRef('sp',str(i)) for i in range(3))
        costs=PairwiseOrderCosts(refs,{(refs[0],refs[1]):0,(refs[1],refs[0]):1,
                                     (refs[1],refs[2]):0,(refs[2],refs[1]):1,
                                     (refs[2],refs[0]):0,(refs[0],refs[2]):1})
        self.assertEqual(probe.pair_minimum(costs),0)
        self.assertEqual(probe.triangle_packing(costs)[0],1)

    def test_conditional_probe_bounds_and_counter_accounting(self):
        fixture=load_fixture(ROOT/'examples/fixtures/fusion_chain_closed_cycle.json')
        result=optimize_auto(fixture,local_restarts=1,transition_cap_per_component=1,
                             branch_node_cap_per_component=50)
        components,_=_component_map(fixture)
        for nodes in components:
            for species in fixture.species_ids:
                row=probe.probe_row(fixture,result.layout.optimized_state,nodes,species)
                self.assertLessEqual(row['triangle_bound'],row['exact_conditional_optimum'])
        details=result.details
        self.assertEqual(result.layout.states_evaluated,
            details['factor_table_entries_evaluated']+details['branch_search_nodes_evaluated']+
            details['factor_continuation_nodes_evaluated'])
