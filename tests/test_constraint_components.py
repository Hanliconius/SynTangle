"""Regression: hard orientation equations may join distinct homology components."""
import unittest
from syntangle.model import Fixture,Chromosome,ChromosomeRef,BlockOccurrence,OrientationConstraint
from syntangle.component_space import optimization_component_map
from syntangle.incidence import build_incidence_graph
from syntangle.layout import initial_layout_state,LayoutState,score_crossings
from syntangle.orientation_space import orientation_basis
from syntangle.hybrid_experiment import run_hybrid
from test_final_exhaustive import enumerate_legal


def coupled_fixture():
    chroms=[]
    for sp in ('A','B'):
        for i in range(2):
            hs=[2*i,2*i+1]
            if sp=='B' and i==1:hs.reverse()
            blocks=tuple(BlockOccurrence(f'{sp}_{h}',f'H{h}',k*3,k*3+1,'?') for k,h in enumerate(hs))
            chroms.append(Chromosome(ChromosomeRef(sp,str(i)),10,i+1,1,blocks))
    constraints=(OrientationConstraint(chroms[0].ref,chroms[1].ref,0),OrientationConstraint(chroms[2].ref,chroms[3].ref,0))
    return Fixture(1,'hard_bridge','hard bridge','homology disconnected but hard coupled',tuple(chroms),constraints,{})

class ConstraintComponentTests(unittest.TestCase):
    def test_all_solver_partitions_use_same_constraint_closed_space(self):
        from syntangle import layout,heuristic,branch_bound,residual_solver,layer_dp
        f=coupled_fixture();self.assertEqual(len(build_incidence_graph(f).connected_components()),2)
        expected=optimization_component_map(f);self.assertEqual(len(expected[0]),1)
        for module in (layout,heuristic,branch_bound,residual_solver):self.assertEqual(module._component_map(f),expected)
        self.assertEqual(layer_dp._component_data(f)[1:],expected)
        orientation_basis(f,tuple(sorted(f.chromosome_refs)))
        from syntangle.residual import build_residual_factorization
        graph=build_residual_factorization(f)
        self.assertTrue(all(v.component_id==0 for v in graph.variables))
        self.assertEqual(len(build_incidence_graph(f).connected_components()),2)
    def test_hard_bridge_solvers_match_independent_global_optimum(self):
        from syntangle import optimize_auto
        f=coupled_fixture();best=min(score_crossings(f,s).crossings for s in enumerate_legal(f))
        self.assertEqual(best,1)
        basis=orientation_basis(f,f.chromosome_refs);raw=initial_layout_state(f);start=LayoutState(raw.chromosome_order,basis.base_assignment)
        for mirror in (False,True):
            r=run_hybrid(f,start,'milp_hint',10,mirror_symmetry=mirror)
            self.assertEqual(r['upper_bound'],best);self.assertEqual(r['lower_bound'],best)
        r=optimize_auto(f,local_restarts=1,branch_node_cap_per_component=1000)
        self.assertEqual(r.layout.optimized_score.crossings,best)
