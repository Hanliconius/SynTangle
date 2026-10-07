"""Scientific metric contracts; synthetic references are not native GENESPACE."""
from itertools import permutations, product
import unittest
from unittest.mock import patch

from syntangle.audit import fixture_fingerprint
from syntangle.layout import LayoutState, initial_layout_state, score_crossings
from syntangle.model import Fixture, Chromosome, ChromosomeRef, BlockOccurrence
from syntangle.orientation_space import orientation_basis
from syntangle.refinement import select_refinement, refine_reference_layout
from syntangle.tangledness import tangledness_metrics
from test_coupled_bound import fixture_for


def presentation_fixture():
    chromosomes=[]
    for species in ('A','B','C'):
        for j in range(8):
            ref=ChromosomeRef(species,str(j))
            blocks=(BlockOccurrence(f'{species}_{j}',f'H{j}',2,3,'+'),)
            chromosomes.append(Chromosome(ref,10,7-j if species=='B' else j,1,blocks))
    return Fixture(1,'presentation','presentation','synthetic display artifact',tuple(chromosomes),(),{})


def untangled(fixture):
    return LayoutState({sp:tuple(r for r in fixture.chromosome_refs if r.species_id==sp)
                        for sp in fixture.species_ids},dict.fromkeys(fixture.chromosome_refs,1))


class TanglednessTests(unittest.TestCase):
    def test_high_visual_tangle_can_have_proven_zero_intrinsic_tangle(self):
        fixture=presentation_fixture(); reference=initial_layout_state(fixture)
        fingerprint=fixture_fingerprint(fixture)
        self.assertEqual(score_crossings(fixture,reference).crossings,56)
        result=select_refinement(fixture,reference,untangled(fixture),reference_name='synthetic')
        metrics=result['reference_metrics']
        self.assertEqual(metrics['intrinsic_tangledness'],0)
        self.assertEqual(metrics['excess_layout_tangledness'],56)
        self.assertTrue(metrics['optimality_established'])
        self.assertEqual(fixture_fingerprint(fixture),fingerprint)

    def test_presentation_changes_visual_score_without_changing_zero_optimum(self):
        fixture=presentation_fixture(); zero=untangled(fixture)
        for row in permutations(zero.chromosome_order['B'][:3]):
            orders=dict(zero.chromosome_order)
            orders['B']=row+zero.chromosome_order['B'][3:]
            display=LayoutState(orders,dict(zero.chromosome_orientation))
            result=select_refinement(fixture,display,zero,reference_name='synthetic')
            self.assertEqual(result['reference_metrics']['intrinsic_tangledness'],0)
            self.assertEqual(result['reference_metrics']['excess_layout_tangledness'],
                             score_crossings(fixture,display).crossings)

    def test_unresolved_optimum_reports_intervals_not_exact_excess(self):
        metrics=tangledness_metrics(100,40,20)
        self.assertIsNone(metrics['intrinsic_tangledness'])
        self.assertIsNone(metrics['excess_layout_tangledness'])
        self.assertEqual(metrics['intrinsic_tangledness_bounds'],dict(lower=20,upper=40))
        self.assertEqual(metrics['excess_layout_tangledness_bounds'],dict(lower=60,upper=80))
        self.assertEqual(metrics['best_retained_excess_bounds'],dict(lower=0,upper=20))
        self.assertEqual(metrics['demonstrated_avoidable_crossings'],60)

    def test_intervals_cover_exhaustively_known_nonzero_optimum(self):
        fixture=fixture_for(19);basis=orientation_basis(fixture,fixture.chromosome_refs)
        rows=[tuple(permutations(tuple(r for r in fixture.chromosome_refs if r.species_id==sp)))
              for sp in fixture.species_ids]
        scores=[]
        for bits in product((0,1),repeat=len(basis.free_flip_groups)):
            signs=dict(basis.base_assignment)
            for group,bit in zip(basis.free_flip_groups,bits):
                for ref in group:signs[ref]*=(-1 if bit else 1)
            for orders in product(*rows):
                scores.append(score_crossings(fixture,LayoutState(dict(zip(fixture.species_ids,orders)),signs)).crossings)
        optimum=min(scores);self.assertGreater(optimum,0)
        display=max(scores);candidate=sorted(set(scores))[1]
        metrics=tangledness_metrics(display,candidate,0)
        bounds=metrics['intrinsic_tangledness_bounds'];excess=metrics['excess_layout_tangledness_bounds']
        self.assertLessEqual(bounds['lower'],optimum);self.assertGreaterEqual(bounds['upper'],optimum)
        self.assertLessEqual(excess['lower'],display-optimum);self.assertGreaterEqual(excess['upper'],display-optimum)
        exact=tangledness_metrics(display,optimum,optimum)
        self.assertEqual(exact['excess_layout_tangledness'],display-optimum)

    def test_known_display_tightens_inferior_candidate(self):
        metrics=tangledness_metrics(20,30,5)
        self.assertEqual(metrics['intrinsic_tangledness_bounds']['upper'],20)
        self.assertEqual(metrics['demonstrated_avoidable_crossings'],0)

    def test_bad_counts_and_inconsistent_bound_rejected(self):
        for args in ((-1,0,0),(1,True,0),(1,0,0.5),(3,2,3)):
            with self.subTest(args=args),self.assertRaises(ValueError):tangledness_metrics(*args)


class ReferenceRefinementTests(unittest.TestCase):
    def setUp(self):
        self.fixture=presentation_fixture()
        self.reference=initial_layout_state(self.fixture)
        self.zero=untangled(self.fixture)

    def test_regression_and_absent_candidate_keep_reference_without_win(self):
        for candidate in (self.reference,None):
            result=select_refinement(self.fixture,self.zero,candidate,reference_name='synthetic')
            self.assertEqual(result['returned_crossings'],0)
            self.assertFalse(result['improved']);self.assertTrue(result['reference_retained'])
            self.assertEqual(result['selected_state'],self.zero.to_dict())
        self.assertIn('regression',select_refinement(self.fixture,self.zero,self.reference)['outcome'])

    def test_different_tied_layout_retains_original(self):
        tied=LayoutState({sp:tuple(reversed(row)) for sp,row in self.zero.chromosome_order.items()},
                         dict(self.zero.chromosome_orientation))
        result=select_refinement(self.fixture,self.zero,tied)
        self.assertEqual(result['raw_candidate_crossings'],0)
        self.assertEqual(result['selected_state'],self.zero.to_dict())
        self.assertIn('tie',result['outcome'])

    def test_illegal_order_and_orientation_are_not_comparisons(self):
        orders=dict(self.zero.chromosome_order);orders['A']=(orders['A'][0],)*8
        bad=LayoutState(orders,dict(self.zero.chromosome_orientation))
        with self.assertRaises(ValueError):select_refinement(self.fixture,self.reference,bad)
        signs=dict(self.zero.chromosome_orientation);signs[next(iter(signs))]=0
        with self.assertRaises(ValueError):select_refinement(self.fixture,LayoutState(self.zero.chromosome_order,signs))
        constrained=fixture_for(19)
        with self.assertRaises(ValueError):select_refinement(constrained,initial_layout_state(constrained))

    def test_real_solver_keeps_baseline_evidence_and_consistent_metrics(self):
        before=self.reference.to_dict();fingerprint=fixture_fingerprint(self.fixture)
        result=refine_reference_layout(self.fixture,self.reference,seconds=2,method='milp_reclaim',reference_name='synthetic')
        self.assertEqual(result['returned_crossings'],0)
        self.assertEqual(self.reference.to_dict(),before)
        self.assertEqual(fixture_fingerprint(self.fixture),fingerprint)
        self.assertEqual(result['solver_result']['tangledness']['intrinsic_tangledness'],0)

    def test_optimizer_receives_a_copy_of_the_reference(self):
        saved=self.reference.to_dict()
        def mutating_optimizer(fixture,start,*args,**kwargs):
            start.chromosome_order.clear()
            return dict(optimized_state=self.zero.to_dict(),upper_bound=0,lower_bound=0)
        with patch('syntangle.hybrid_experiment.run_hybrid',side_effect=mutating_optimizer):
            result=refine_reference_layout(self.fixture,self.reference,reference_name='synthetic')
        self.assertEqual(self.reference.to_dict(),saved)
        self.assertEqual(result['reference_crossings'],56)

    def test_adapter_checks_raw_score_and_retains_worse_solver_result(self):
        raw=dict(optimized_state=self.reference.to_dict(),upper_bound=56,lower_bound=0)
        with patch('syntangle.hybrid_experiment.run_hybrid',return_value=raw):
            result=refine_reference_layout(self.fixture,self.zero,reference_name='synthetic')
        self.assertEqual(result['raw_candidate_crossings'],56)
        self.assertEqual(result['returned_crossings'],0)
        raw['upper_bound']=1
        with patch('syntangle.hybrid_experiment.run_hybrid',return_value=raw),self.assertRaises(AssertionError):
            refine_reference_layout(self.fixture,self.zero)


if __name__=='__main__':unittest.main()
