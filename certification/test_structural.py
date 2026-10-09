import itertools
import random
import unittest
from proof import Evidence
from structural import lower_bound, report, construct_thin
from test_proof import fixture


def square(two_per_bundle=False):
    rows = []
    for sp in ['A','B']:
        chroms = []
        for c in range(2):
            blocks = []
            for other in range(2):
                a,b = (c,other) if sp=='A' else (other,c)
                for k in range(2 if two_per_bundle else 1):
                    # Exactly one bundle prefers opposite orientations.
                    slot = 1-k if sp=='B' and a==b==1 and two_per_bundle else k
                    start = other*4 + slot
                    h = f'h{a}{b}{k}'
                    blocks.append(dict(occurrence_id=sp+h,homology_id=h,start=start,end=start+1,strand='+'))
            chroms.append(dict(id=str(c),length=10,blocks=blocks))
        rows.append(dict(id=sp,chromosomes=chroms))
    return dict(fixture_version=1,id='square',species=rows)


class StructuralTests(unittest.TestCase):
    def test_positive_structural_optimum_and_construction(self):
        evidence = Evidence(fixture())
        self.assertEqual(lower_bound(evidence)['lower'],1)
        result = construct_thin(evidence)
        self.assertEqual(result['status'],'EXACT_THIN_CONSTRUCTION')
        self.assertEqual(report(evidence,result['layout'])['status'],'EXACT_STRUCTURAL')

    def test_balance_is_not_sufficient_for_zero(self):
        evidence = Evidence(square())
        self.assertEqual(lower_bound(evidence)['lower'],0)
        optimum = min(evidence.score(state) for state in evidence.states())
        self.assertEqual(optimum,1)
        self.assertEqual(construct_thin(evidence)['status'],'OUTSIDE_THIN_CLASS')

    def test_negative_cycle_supplies_positive_penalty(self):
        evidence = Evidence(square(True))
        result = lower_bound(evidence)
        self.assertEqual(result['base'],0)
        self.assertEqual(result['lower'],1)
        self.assertEqual(len(result['cycles']),1)
        for state in evidence.states():
            self.assertGreaterEqual(evidence.score(state),result['lower'])

    def test_hard_forced_parity_and_inconsistent_hard_cycle(self):
        raw = fixture()
        raw['orientation_constraints'] = [dict(a='A:1',b='B:1',xor=1)]
        evidence = Evidence(raw)
        self.assertEqual(lower_bound(evidence)['lower'],2)
        for state in evidence.states():
            self.assertEqual(report(evidence,state)['status'],'EXACT_STRUCTURAL')
        raw['orientation_constraints'].append(dict(a='A:1',b='B:1',xor=0))
        with self.assertRaises(ValueError):
            lower_bound(Evidence(raw))

    def test_random_bounds_against_every_legal_small_layout(self):
        rng = random.Random(414)
        for repeat in range(20):
            rows = []
            for sp in ['A','B','C']:
                blocks = [[],[]]
                for h in range(5):
                    chrom = rng.randrange(2)
                    start = rng.randrange(15)
                    blocks[chrom].append(dict(occurrence_id=f'{sp}{h}',homology_id=str(h),start=start,end=start+1,strand='+'))
                rows.append(dict(id=sp,chromosomes=[dict(id=str(c),length=20,blocks=blocks[c]) for c in range(2)]))
            raw = dict(fixture_version=1,id=f'random{repeat}',species=rows)
            evidence = Evidence(raw)
            bound = lower_bound(evidence)['lower']
            for state in evidence.states():
                self.assertGreaterEqual(evidence.score(state),bound)


if __name__ == '__main__':
    unittest.main()
