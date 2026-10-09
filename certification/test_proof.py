import itertools
import unittest
import tempfile
import subprocess
import json
from pathlib import Path
from fractions import Fraction
from proof import Evidence, IntegerModel, exhaustive


def fixture():
    def row(sp, homologies):
        return {'id':sp, 'chromosomes':[{'id':'1', 'length':10,
            'blocks':[{'occurrence_id':sp+h, 'homology_id':h,
                       'start':i*2, 'end':i*2+1, 'strand':'+'}
                      for i,h in enumerate(homologies)]}]}
    return {'fixture_version':1, 'id':'nonzero',
            'species':[row('A',['a','b','c']), row('B',['a','c','b'])]}


class ExactProofTests(unittest.TestCase):
    def test_certificate_replay_rejects_tampering(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = fixture()
            (root/'fixture.json').write_text(json.dumps(raw))
            state = next(Evidence(raw).states())
            (root/'layout.json').write_text(json.dumps(state))
            common = [str(root/'fixture.json'), str(root/'layout.json'), str(root/'certificate.json')]
            script = str(Path(__file__).with_name('proof.py'))
            import sys
            subprocess.run([sys.executable,script,'check',*common], check=True, capture_output=True)
            subprocess.run([sys.executable,script,'verify',*common], check=True, capture_output=True)
            certificate = json.loads((root/'certificate.json').read_text())
            certificate['checked'] += 1
            (root/'certificate.json').write_text(json.dumps(certificate))
            self.assertNotEqual(subprocess.run([sys.executable,script,'verify',*common], capture_output=True).returncode,0)
            subprocess.run([sys.executable,script,'export',*common[:2],str(root/'model')], check=True, capture_output=True)
            self.assertEqual(json.loads((root/'model/manifest.json').read_text())['status'],'EXPORTED_NOT_CERTIFIED')

    def test_nonzero_optimum_and_counterexample(self):
        evidence = Evidence(fixture())
        state = next(evidence.states())
        result = exhaustive(evidence, state, 100, 10)
        self.assertEqual((result['status'],result['upper']), ('EXACT_EXHAUSTIVE',1))
        worse = {**state, 'chromosome_orientation':{'A:1':1, 'B:1':-1}}
        self.assertEqual(exhaustive(evidence,worse,100,10)['status'], 'COUNTEREXAMPLE')

    def test_model_matches_geometry_over_every_small_layout(self):
        raw = fixture()
        # Split selected chromosomes: exercises order/flip, order/order and flip/flip.
        for row in raw['species']:
            old = row['chromosomes'][0]
            row['chromosomes'].append({'id':'2','length':10,'blocks':[old['blocks'].pop()]})
            row['chromosomes'].append({'id':'3','length':10,'blocks':[]})
        evidence = Evidence(raw)
        model = IntegerModel(evidence)
        for state in evidence.states():
            values = model.assignment(state)
            self.assertEqual(model.evaluate(values), evidence.score(state))
            for terms,sense,rhs in model.rows:
                lhs = sum(c*values[v] for v,c in terms.items())
                self.assertTrue({'=':lhs==rhs, '<=':lhs<=rhs, '>=':lhs>=rhs}[sense])

    def test_triangle_constraints_allow_exactly_permutations(self):
        raw = fixture()
        raw['species'] = [{'id':'A','chromosomes':[
            {'id':str(i),'length':1,'blocks':[]} for i in range(4)]}]
        model = IntegerModel(Evidence(raw))
        variables = list(model.orders.values())
        accepted = 0
        for bits in itertools.product((0,1),repeat=len(variables)):
            values = dict(zip(variables,bits))
            if all((sum(c*values[v] for v,c in terms.items()) <= rhs if sense=='<='
                    else sum(c*values[v] for v,c in terms.items()) >= rhs)
                   for terms,sense,rhs in model.rows):
                accepted += 1
        self.assertEqual(accepted,24)

    def test_constraints_and_refusal_to_claim_after_limits(self):
        raw = fixture()
        raw['orientation_constraints'] = [{'a':'A:1','b':'B:1','xor':0}]
        evidence = Evidence(raw)
        self.assertEqual(len(list(evidence.states())),2)
        state = next(evidence.states())
        self.assertEqual(exhaustive(evidence,state,0,10)['status'],'UNRESOLVED_SPACE_LIMIT')
        self.assertEqual(exhaustive(evidence,state,100,0)['status'],'UNRESOLVED_DEADLINE')
        model = IntegerModel(evidence)
        self.assertEqual(model.evaluate(model.assignment(state)),1)

    def test_zero_ties_and_exact_coordinates(self):
        raw = fixture()
        raw['species'][1]['chromosomes'][0]['blocks'] = [
            {**b,'start':0,'end':1} for b in raw['species'][1]['chromosomes'][0]['blocks']]
        evidence = Evidence(raw)
        state = next(evidence.states())
        self.assertEqual(exhaustive(evidence,state,0,0)['status'],'EXACT_ZERO')
        self.assertEqual(IntegerModel(evidence).evaluate(IntegerModel(evidence).assignment(state)),0)
        raw = fixture()
        for row in raw['species']:
            for i,b in enumerate(row['chromosomes'][0]['blocks']):
                b['start'] = Fraction(i,10**30)
                b['end'] = b['start'] + Fraction(1,10**31)
        evidence = Evidence(raw)
        self.assertEqual(evidence.score(next(evidence.states())),1)

    def test_ambiguous_or_extra_constraints_fail_closed(self):
        raw = fixture()
        raw['species'][0]['chromosomes'][0]['blocks'].append(
            {**raw['species'][0]['chromosomes'][0]['blocks'][0],'occurrence_id':'extra'})
        with self.assertRaises(ValueError):
            Evidence(raw)
        raw = fixture()
        raw['precedence_constraints'] = []
        with self.assertRaises(ValueError):
            Evidence(raw)


if __name__ == '__main__':
    unittest.main()
