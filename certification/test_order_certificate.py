import copy
from fractions import Fraction
import unittest
import tempfile
import subprocess
import sys
import json
from pathlib import Path
from proof import Evidence, IntegerModel
from test_proof import fixture
from order_certificate import bound, canonical_rows, replay

class OrderCertificateTests(unittest.TestCase):
    def setUp(self):
        self.evidence = Evidence(fixture())
        self.model = IntegerModel(self.evidence)
        self.rows = canonical_rows(self.model)

    def node(self,fixed):
        values = ['0']*len(self.rows)
        return dict(fixed=fixed,multipliers=values,lower=bound(self.model,self.rows,fixed,values))

    def test_arbitrary_rational_multipliers_never_exclude_a_legal_layout(self):
        multipliers = [Fraction(i+1,7) if sense=='<=' else Fraction(-2,3)
                       for i,(_,sense,_) in enumerate(self.rows)]
        for state in self.evidence.states():
            assignment = self.model.assignment(state)
            score = self.evidence.score(state)
            for fixed in ({},{self.model.binary[0]:assignment[self.model.binary[0]]}):
                self.assertLessEqual(bound(self.model,self.rows,fixed,multipliers),score)

    def test_negative_inequality_multiplier_rejected(self):
        multipliers = ['0']*len(self.rows)
        multipliers[next(i for i,r in enumerate(self.rows) if r[1]=='<=')] = '-1'
        with self.assertRaises(ValueError):
            bound(self.model,self.rows,{},multipliers)

    def test_tree_partition_and_tampering(self):
        variable = self.model.binary[0]
        nodes = [self.node({}),self.node({variable:0}),self.node({variable:1})]
        nodes[0].update(branch=variable,children=[1,2])
        certificate = dict(nodes=nodes)
        result = replay(self.model,1,certificate)
        self.assertEqual(result['status'],'UNRESOLVED_ORDER_GAP')
        for mutate in ('inheritance','shared','bound','unreachable'):
            altered = copy.deepcopy(certificate)
            if mutate=='inheritance': altered['nodes'][2]['fixed'][variable]=0
            if mutate=='shared': altered['nodes'][0]['children']=[1,1]
            if mutate=='bound': altered['nodes'][1]['lower']=999
            if mutate=='unreachable': altered['nodes'].append(self.node({}))
            with self.assertRaises(ValueError): replay(self.model,1,altered)

    def test_generation_and_fresh_process_replay(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            f,l,c = [root/name for name in ('fixture.json','layout.json','certificate.json')]
            f.write_text(json.dumps(fixture()))
            l.write_text(json.dumps(next(self.evidence.states())))
            script = str(Path(__file__).with_name('order_certificate.py'))
            args = [str(f),str(l),str(c)]
            subprocess.run([sys.executable,script,'check',*args,'--seconds','5','--max-nodes','31'],
                           check=True,capture_output=True)
            saved = json.loads(c.read_text())
            self.assertEqual(saved['summary']['status'],'EXACT_ORDER_CERTIFICATE')
            subprocess.run([sys.executable,script,'verify',*args],check=True,capture_output=True)
            saved['summary']['lower'] += 1
            c.write_text(json.dumps(saved))
            self.assertNotEqual(subprocess.run([sys.executable,script,'verify',*args],
                                              capture_output=True).returncode,0)

    def test_integer_rounding_can_certify_fractional_lower_bound(self):
        # Min x over binary x with 2x>=1 is exactly 1, LP lower bound 1/2.
        class Model:
            binary=['x']; products={}; objective={'x':1}; constant=0
            rows=[({'x':2},'>=',1)]
        model=Model()
        certificate=dict(nodes=[dict(fixed={},multipliers=['1/2'],lower=1)])
        self.assertEqual(replay(model,1,certificate)['status'],'EXACT_ORDER_CERTIFICATE')

if __name__=='__main__': unittest.main()
