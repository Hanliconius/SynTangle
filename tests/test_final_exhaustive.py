"""Pegasus final checks: independent brute force versus solver policies.
No model-derived enumeration or lower bound is used for the reference optimum.
"""
from itertools import permutations,product
import random,unittest
from syntangle.model import Fixture,Chromosome,ChromosomeRef,BlockOccurrence,OrientationConstraint
from syntangle.layout import LayoutState,initial_layout_state,score_crossings
from syntangle.audit import fixture_fingerprint
from syntangle.saved_layout import decode_saved_layout,validate_saved_layout
from syntangle.orientation_space import orientation_basis
from syntangle.hybrid_experiment import run_hybrid


def tiny(seed,three=False):
    rng=random.Random(seed);chroms=[]
    for sp in ('A','B') if three else ('A','B','C'):
        hs=list(range(6));rng.shuffle(hs)
        n=3 if three else 2
        for j in range(n):
            blocks=tuple(BlockOccurrence(f'{sp}_{h}',f'H{h}',k*3,k*3+1,'?') for k,h in enumerate(hs[j*(6//n):(j+1)*(6//n)]))
            chroms.append(Chromosome(ChromosomeRef(sp,str(j)),20,j+1,1,blocks))
    constraints=(OrientationConstraint(chroms[0].ref,chroms[n].ref,seed%2),)
    return Fixture(1,f'final_{seed}_{three}','final','independent exhaustive legal-layout check',tuple(chroms),constraints,{})


def enumerate_legal(f):
    refs=f.chromosome_refs;rows=[tuple(permutations(tuple(r for r in refs if r.species_id==sp))) for sp in f.species_ids]
    for signs in product((-1,1),repeat=len(refs)):
        orient=dict(zip(refs,signs))
        if any((orient[c.a]!=orient[c.b])!=bool(c.xor) for c in f.orientation_constraints):continue
        for orders in product(*rows):yield LayoutState(dict(zip(f.species_ids,orders)),orient)


class FinalExhaustiveTests(unittest.TestCase):
    def test_random_graphs_direct_mirror_and_selective_match_bruteforce(self):
        for seed in range(10):
            for three in (False,True):
                f=tiny(seed,three);fingerprint=fixture_fingerprint(f)
                best=min(score_crossings(f,s).crossings for s in enumerate_legal(f))
                basis=orientation_basis(f,f.chromosome_refs);raw=initial_layout_state(f);start=LayoutState(raw.chromosome_order,basis.base_assignment)
                for method,options in [('milp_hint',{}),('milp_hint',dict(mirror_symmetry=True)),('hybrid_adaptive',dict(mirror_symmetry=True,neighborhood_min_decisions=0,neighborhood_fraction=.15,neighborhood_max_seconds=1,deduplicate_neighborhoods=True))]:
                    checkpoints=[];result=run_hybrid(f,start,method,10,progress=lambda s:checkpoints.append(score_crossings(f,s).crossings),**options)
                    state=decode_saved_layout(f,result['optimized_state']);validate_saved_layout(f,state,basis)
                    self.assertEqual(result['upper_bound'],best,(seed,three,method,options))
                    self.assertEqual(result['lower_bound'],best)
                    self.assertEqual(score_crossings(f,state).crossings,best)
                    self.assertEqual(fixture_fingerprint(f),fingerprint)
                    self.assertTrue(all(a>=b for a,b in zip(checkpoints,checkpoints[1:])))
                    self.assertEqual(result['global_search_restarts'],0)
    def test_zero_deadline_preserves_valid_candidate_and_honest_bound(self):
        f=tiny(3);basis=orientation_basis(f,f.chromosome_refs);raw=initial_layout_state(f);start=LayoutState(raw.chromosome_order,basis.base_assignment)
        optimum=min(score_crossings(f,s).crossings for s in enumerate_legal(f))
        result=run_hybrid(f,start,'milp_hint',0)
        state=decode_saved_layout(f,result['optimized_state']);validate_saved_layout(f,state,basis)
        self.assertLessEqual(result['lower_bound'],optimum)
        self.assertGreaterEqual(result['upper_bound'],optimum)
        self.assertLessEqual(result['upper_bound'],score_crossings(f,start).crossings)
        self.assertEqual(score_crossings(f,state).crossings,result['upper_bound'])
