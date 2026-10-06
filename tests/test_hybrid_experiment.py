import unittest
from pathlib import Path
from syntangle import load_fixture
from unittest.mock import patch
from test_coupled_bound import fixture_for
from syntangle.global_experiment import JointModel
from syntangle.hybrid_experiment import run_hybrid, improve_neighborhoods
from syntangle.layout import LayoutState,initial_layout_state,score_crossings
from syntangle.incidence import chromosome_node_id


class HybridTests(unittest.TestCase):
    def test_hint_checkpoints_and_strict_optimum_audit(self):
        f=fixture_for(19)
        m=JointModel(f,frozenset(chromosome_node_id(r) for r in f.chromosome_refs),f.chromosome_refs)
        initial=initial_layout_state(f);initial=LayoutState(initial.chromosome_order,m.basis.base_assignment)
        checkpoints=[]
        state,bound,info=m.highs(initial,5,progress=checkpoints.append)
        self.assertTrue(info['feasible_mip_start'])
        self.assertEqual(bound,score_crossings(f,state).crossings)
        self.assertTrue(m.highs(state,5,strict_proof=True)['no_better_proven'])
        self.assertIsNotNone(m.highs(initial,5,strict_proof=True)['counterexample_crossings'])
        self.assertTrue(checkpoints)
        self.assertIs(m.linear_model(),m.linear_model())

    def test_all_variants_preserve_incumbent_and_global_bound(self):
        f=fixture_for(31)
        m=JointModel(f,frozenset(chromosome_node_id(r) for r in f.chromosome_refs),f.chromosome_refs)
        initial=initial_layout_state(f);initial=LayoutState(initial.chromosome_order,m.basis.base_assignment)
        for method in ('milp_reclaim','milp_hint','hybrid_3','hybrid_adaptive'):
            result=run_hybrid(f,initial,method,3)
            self.assertEqual(result['upper_bound'],result['lower_bound'])
            self.assertLessEqual(result['upper_bound'],score_crossings(f,initial).crossings)
            self.assertEqual(result['global_search_restarts'],0)

    def test_neighborhoods_then_global_reuse_same_matrix(self):
        f=fixture_for(19)
        m=JointModel(f,frozenset(chromosome_node_id(r) for r in f.chromosome_refs),f.chromosome_refs)
        initial=initial_layout_state(f);initial=LayoutState(initial.chromosome_order,m.basis.base_assignment)
        cached=m.linear_model()
        improved,info=improve_neighborhoods(m,initial,.5,adaptive=True)
        self.assertGreater(info['subsolves'],0)
        self.assertLessEqual(score_crossings(f,improved).crossings,score_crossings(f,initial).crossings)
        final,bound,_=m.highs(improved,5)
        self.assertEqual(bound,score_crossings(f,final).crossings)
        self.assertIs(cached,m.linear_model())

    def test_proven_zero_components_do_not_reenter_search(self):
        f=load_fixture(Path(__file__).resolve().parents[1]/'examples/fixtures/perfect_1to1_30x3.json')
        state=initial_layout_state(f)
        with patch('syntangle.hybrid_experiment.JointModel',side_effect=AssertionError('Zero component reentered search')):
            result=run_hybrid(f,state,'hybrid_adaptive',5)
        self.assertEqual(result['upper_bound'],0)
        self.assertEqual(result['lower_bound'],0)
        self.assertTrue(all(d['component_upper']==0 for d in result['component_diagnostics']))

    def test_single_remaining_component_gets_all_remaining_time(self):
        f=fixture_for(31)
        m=JointModel(f,frozenset(chromosome_node_id(r) for r in f.chromosome_refs),f.chromosome_refs)
        initial=initial_layout_state(f);initial=LayoutState(initial.chromosome_order,m.basis.base_assignment)
        result=run_hybrid(f,initial,'milp_reclaim',10)
        self.assertEqual(len(result['component_diagnostics']),1)
        self.assertGreater(result['component_diagnostics'][0]['allocated_seconds'],9)
