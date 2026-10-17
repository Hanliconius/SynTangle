from itertools import permutations, product
from pathlib import Path
import random
import unittest
from syntangle import load_fixture, optimize_auto, score_crossings
from syntangle.bounds import build_relaxed_crossing_bound
from syntangle.coupled_bound import CoupledCrossingBound, bucket_tables
from syntangle.incidence import chromosome_node_id
from syntangle.layout import LayoutState
from syntangle.model import Fixture, Chromosome, ChromosomeRef, BlockOccurrence, OrientationConstraint
from syntangle.orientation_space import orientation_basis
import syntangle.branch_bound as branch

ROOT=Path(__file__).resolve().parents[1]


def fixture_for(seed):
    rng=random.Random(seed);chromosomes=[]
    for species in ('A','B','C'):
        homologies=list(range(8));rng.shuffle(homologies)
        for j in range(2):
            ref=ChromosomeRef(species,str(j))
            blocks=tuple(BlockOccurrence(f'{species}_{h}',f'H{h}',k*3,k*3+1,'+')
                         for k,h in enumerate(homologies[j*4:j*4+4]))
            chromosomes.append(Chromosome(ref,20,j,1,blocks))
    constraints=(OrientationConstraint(chromosomes[0].ref,chromosomes[2].ref,1),)
    return Fixture(1,'tiny','tiny','coupled exhaustive',tuple(chromosomes),constraints,{})


class CoupledBoundTests(unittest.TestCase):
    def test_all_partial_orientations_and_fixed_orders_are_admissible(self):
        for seed in (19,31):
            fixture=fixture_for(seed)
            nodes=frozenset(chromosome_node_id(r) for r in fixture.chromosome_refs)
            basis=orientation_basis(fixture,fixture.chromosome_refs)
            independent=build_relaxed_crossing_bound(fixture,nodes,basis)
            bound=CoupledCrossingBound(fixture,nodes,independent,6)
            states=[]
            for bits in product((0,1),repeat=len(basis.free_flip_groups)):
                orientation=dict(basis.base_assignment)
                for group,bit in zip(basis.free_flip_groups,bits):
                    for ref in group:
                        orientation[ref]*=(-1 if bit else 1)
                row_orders=[tuple(permutations(tuple(r for r in fixture.chromosome_refs if r.species_id==sp)))
                            for sp in fixture.species_ids]
                for orders in product(*row_orders):
                    state=LayoutState(dict(zip(fixture.species_ids,orders)),orientation)
                    score=score_crossings(fixture,state).crossings
                    indexed=dict(enumerate(orders))
                    self.assertEqual(bound.lower_bound(bits,indexed),score)
                    states.append((bits,indexed,score))
            for partial in product((None,0,1),repeat=len(basis.free_flip_groups)):
                matching=[row for row in states if all(a is None or a==b for a,b in zip(partial,row[0]))]
                self.assertLessEqual(bound.lower_bound(partial),min(r[2] for r in matching))
                self.assertGreaterEqual(bound.lower_bound(partial),independent.lower_bound(partial))
                fixed={1:matching[0][1][1]}
                self.assertLessEqual(bound.lower_bound(partial,fixed),
                                     min(r[2] for r in matching if r[1][1]==fixed[1]))

    def test_shared_order_variable_captures_neighbor_conflict(self):
        factors=[((0,1),{(a,b):int(a!=b) for a,b in product((0,1),repeat=2)}),
                 ((1,2),{(a,b):int(a!=b) for a,b in product((0,1),repeat=2)})]
        tables=bucket_tables(factors,3)
        assignment=(0,None,1)
        self.assertEqual(sum(t[tuple(assignment[v] for v in scope)] for scope,t in tables),1)
        self.assertEqual(sum(min(value for bits,value in t.items()
                                 if all(assignment[v] is None or assignment[v]==b for v,b in zip(scope,bits)))
                             for scope,t in factors),0)

    def test_coupled_branch_result_matches_exhaustive_fixture(self):
        fixture=load_fixture(ROOT/'examples/fixtures/fusion_chain_closed_cycle.json')
        exact=optimize_auto(fixture,local_restarts=1)
        result=optimize_auto(fixture,local_restarts=1,transition_cap_per_component=1,
                             branch_node_cap_per_component=10000,coupled_bound_size=6)
        self.assertLessEqual(result.details['lower_bound'],exact.details['upper_bound'])
        self.assertEqual(result.details['upper_bound'],exact.details['upper_bound'])
        self.assertIsNone(branch._active_coupled_bound.get())
