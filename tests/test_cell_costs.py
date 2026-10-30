import itertools
import random
import time
import unittest
from unittest.mock import patch

from syntangle.component_space import optimization_component_map
from syntangle.global_experiment import JointModel
from syntangle.hybrid_experiment import run_hybrid
from syntangle.layout import initial_layout_state, score_crossings, LayoutState, canonicalize_component_order
from syntangle.orientation_space import orientation_basis
from syntangle.bounds import _cell_internal_crossings, _pairwise_interleave_minimum
from syntangle.cell_costs import iter_cell_costs, PreparationDeadline
from test_final_exhaustive import tiny, enumerate_legal
from test_component_presolve import strand
from dataclasses import replace


class CellTests(unittest.TestCase):
    def test_factor_aggregation_matches_old_model_and_canonical_scores(self):
        for seed in range(4):
            f = tiny(seed)
            components, component_of = optimization_component_map(f)
            for i, nodes in enumerate(components):
                refs = tuple(r for r in f.chromosome_refs if component_of[r] == i)
                old = JointModel(f, nodes, refs)
                new = JointModel(f, nodes, refs, cell_aggregation=True)
                for state in enumerate_legal(f):
                    expected = score_crossings(f, state, restrict_component_nodes=nodes).crossings
                    self.assertEqual(old.objective(old.encode(state)), expected)
                    self.assertEqual(new.objective(new.encode(state)), expected)

    def test_ties_and_float_reversal_kernels_match_brute_force(self):
        rng = random.Random(7)
        values = [0.1, 0.1, 0.5, 0.9, 1.0-0.9, 0.0, 1.0]
        for _ in range(30):
            anchors = tuple((rng.choice(values), rng.choice(values)) for _ in range(8))
            for x, y in itertools.product((1,-1), repeat=2):
                transformed = [(a if x==1 else 1-a, b if y==1 else 1-b) for a,b in anchors]
                expected = sum((a<c and b>d) or (a>c and b<d) for (a,b),(c,d) in itertools.combinations(transformed,2))
                self.assertEqual(_cell_internal_crossings(anchors, x, y), expected)
            a,b = anchors[:4], anchors[4:]
            for side in (0,1):
                expected = min(sum(x[side]<y[side] for x in a for y in b), sum(x[side]>y[side] for x in a for y in b))
                self.assertEqual(_pairwise_interleave_minimum(a,b,coordinate=side), expected)

    def test_general_solve_retains_strand_and_global_bound(self):
        f = tiny(2)
        z = strand(3)
        z = replace(z, chromosomes=tuple(replace(c, ref=type(c.ref)(sp, 'Z'),
            blocks=tuple(replace(b, homology_id='z'+b.homology_id, occurrence_id='z'+b.occurrence_id) for b in c.blocks))
            for sp,c in zip(f.species_ids, z.chromosomes)))
        f = replace(f, chromosomes=tuple(c for sp in f.species_ids for c in f.chromosomes+z.chromosomes if c.ref.species_id==sp))
        raw = initial_layout_state(f)
        start = LayoutState(raw.chromosome_order, orientation_basis(f,f.chromosome_refs).base_assignment)
        with patch('syntangle.hybrid_experiment.JointModel', wraps=JointModel) as constructor:
            result = run_hybrid(f,start,'milp_reclaim',3,backend='scipy',component_presolve=True)
        self.assertEqual(constructor.call_count, 1)
        self.assertTrue(all(r.chromosome_id != 'Z' for r in constructor.call_args.args[2]))
        self.assertLessEqual(result['upper_bound'], score_crossings(f,start).crossings)
        self.assertEqual(result['upper_bound'], result['lower_bound'])

    def test_preparation_deadline_never_exports_partial_model_bound(self):
        f = tiny(2)
        nodes = optimization_component_map(f)[0][0]
        with self.assertRaises(PreparationDeadline):
            list(iter_cell_costs(f,nodes,time.perf_counter()-1))
        raw = initial_layout_state(f)
        start = LayoutState(raw.chromosome_order, orientation_basis(f,f.chromosome_refs).base_assignment)
        with patch('syntangle.hybrid_experiment.JointModel', side_effect=PreparationDeadline('test')):
            result = run_hybrid(f,start,'milp_reclaim',2,backend='scipy',component_presolve=True)
        self.assertEqual(result['lower_bound'],0)
        self.assertEqual(result['upper_bound'],score_crossings(f,canonicalize_component_order(f,start)).crossings)

    def test_legacy_and_opt_in_scipy_paths_agree(self):
        f = tiny(3)
        raw = initial_layout_state(f)
        start = LayoutState(raw.chromosome_order, orientation_basis(f,f.chromosome_refs).base_assignment)
        results = [run_hybrid(f,start,'milp_reclaim',3,backend='scipy',component_presolve=enabled)
                   for enabled in (False,True)]
        self.assertEqual(results[0]['upper_bound'],results[1]['upper_bound'])
        self.assertEqual(results[0]['lower_bound'],results[1]['lower_bound'])


if __name__ == '__main__':
    unittest.main()
