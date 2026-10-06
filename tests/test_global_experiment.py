from itertools import product, permutations
import unittest
import random
from syntangle.model import Fixture, Chromosome, ChromosomeRef, BlockOccurrence
from test_coupled_bound import fixture_for
from syntangle.layout import LayoutState, initial_layout_state, score_crossings
from syntangle.incidence import chromosome_node_id
from syntangle.global_experiment import JointModel, run_experiment


class GlobalExperimentTests(unittest.TestCase):
    def test_model_matches_every_legal_layout_and_milp_optimum(self):
        for seed in (19,31):
            fixture=fixture_for(seed)
            nodes=frozenset(chromosome_node_id(r) for r in fixture.chromosome_refs)
            model=JointModel(fixture,nodes,fixture.chromosome_refs)
            best=float('inf');start=None
            for bits in product((0,1),repeat=len(model.basis.free_flip_groups)):
                orientation=dict(model.basis.base_assignment)
                for group,bit in zip(model.basis.free_flip_groups,bits):
                    for r in group:orientation[r]*=(-1 if bit else 1)
                for orders in product(*(permutations(row) for row in model.rows.values())):
                    state=LayoutState(dict(zip(fixture.species_ids,orders)),orientation)
                    exact=score_crossings(fixture,state).crossings
                    encoded=model.encode(state)
                    self.assertEqual(exact,model.objective(encoded))
                    self.assertEqual(state,model.decode(encoded,state))
                    best=min(best,exact);start=state
            result=run_experiment(fixture,start,'milp',10)
            self.assertEqual(best,result['upper_bound'])
            self.assertEqual(best,result['lower_bound'])

    def test_three_chromosome_transitivity_matches_exhaustive(self):
        rng=random.Random(7);chromosomes=[]
        for sp in ('A','B'):
            homologies=list(range(6));rng.shuffle(homologies)
            for i in range(3):
                ref=ChromosomeRef(sp,str(i))
                blocks=tuple(BlockOccurrence(f'{sp}_{h}',f'H{h}',j*3,j*3+1,'+')
                    for j,h in enumerate(homologies[2*i:2*i+2]))
                chromosomes.append(Chromosome(ref,20,i,1,blocks))
        fixture=Fixture(1,'three','three','transitivity',tuple(chromosomes),(),{})
        model=JointModel(fixture,frozenset(chromosome_node_id(r) for r in fixture.chromosome_refs),fixture.chromosome_refs)
        best=float('inf');start=initial_layout_state(fixture)
        for flips in product((0,1),repeat=len(model.basis.free_flip_groups)):
            signs=dict(model.basis.base_assignment)
            for group,bit in zip(model.basis.free_flip_groups,flips):
                for r in group:signs[r]*=(-1 if bit else 1)
            for orders in product(*(permutations(row) for row in model.rows.values())):
                state=LayoutState(dict(zip(fixture.species_ids,orders)),signs)
                value=score_crossings(fixture,state).crossings
                self.assertEqual(value,model.objective(model.encode(state)))
                best=min(best,value)
        solved,bound,_=model.milp(start,5)
        self.assertEqual(bound,best)
        self.assertEqual(score_crossings(fixture,solved).crossings,best)

    def test_neighborhood_fixing_is_conditional(self):
        fixture=fixture_for(19)
        model=JointModel(fixture,frozenset(chromosome_node_id(r) for r in fixture.chromosome_refs),fixture.chromosome_refs)
        start=initial_layout_state(fixture)
        # initial orientations may not satisfy hard parity; use propagated basis.
        start=LayoutState(start.chromosome_order,model.basis.base_assignment)
        bits=model.encode(start)
        state,_,info=model.milp(start,5,{i:v for i,v in enumerate(bits)})
        self.assertEqual(start,state)
        self.assertEqual(info['bound_scope'],'restricted neighborhood')
        result=run_experiment(fixture,start,'lns',.1)
        self.assertEqual(result['lower_bound'],0)
        self.assertLessEqual(result['upper_bound'],score_crossings(fixture,start).crossings)

    def test_sdp_relaxation_and_rounding(self):
        fixture=fixture_for(19)
        model=JointModel(fixture,frozenset(chromosome_node_id(r) for r in fixture.chromosome_refs),fixture.chromosome_refs)
        initial=initial_layout_state(fixture)
        start=LayoutState(initial.chromosome_order,model.basis.base_assignment)
        exact,_,_=model.milp(start,5)
        state,info=model.sdp(start,3,full_threshold=0,block_size=3)
        self.assertLessEqual(info['numerical_relaxation_objective'],score_crossings(fixture,exact).crossings+.02)
        self.assertIsNone(info['certified_lower_bound'])
        self.assertLessEqual(score_crossings(fixture,state).crossings,score_crossings(fixture,start).crossings)
