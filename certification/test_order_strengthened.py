import copy
import itertools
import json
from pathlib import Path
import tempfile
import unittest
from proof import Evidence,IntegerModel
from test_proof import fixture
from order_certificate import bound,canonical_rows
from order_strengthened import cycle_cuts,run

class StrengthenedTests(unittest.TestCase):
    def test_four_cycle_validity_and_fractional_obstruction(self):
        class Model:
            products={('a','b'):'ab',('b','c'):'bc',('c','d'):'cd',('a','d'):'ad'};objective={}
        rows,count=cycle_cuts(Model(),1)
        self.assertEqual((len(rows),count),(8,1))
        for a,b,c,d in itertools.product((0,1),repeat=4):
            values=dict(a=a,b=b,c=c,d=d,ab=a*b,bc=b*c,cd=c*d,ad=a*d)
            for terms,sense,rhs in rows:self.assertLessEqual(sum(values[v]*coef for v,coef in terms.items()),rhs)
        values=dict(a=.5,b=.5,c=.5,d=.5,ab=0,bc=.5,cd=.5,ad=.5)
        self.assertTrue(any(sum(values[v]*coef for v,coef in terms.items())>rhs for terms,sense,rhs in rows))
        self.assertEqual(cycle_cuts(Model(),0),([],0))

    def test_real_model_cycle_cuts_cover_every_small_legal_layout(self):
        def species(sp,groups):
            return dict(id=sp,chromosomes=[dict(id=str(i),length=10,blocks=[dict(occurrence_id=sp+h,homology_id=h,start=j*2,end=j*2+1,strand='+') for j,h in enumerate(group)]) for i,group in enumerate(groups)])
        raw=dict(fixture_version=1,id='cycle_model',species=[species('A',[['a','e','b','f'],['c','g','d','h']]),species('B',[['a','e','c','g'],['b','f','d','h']])])
        evidence=Evidence(raw);model=IntegerModel(evidence);extra,count=cycle_cuts(model,20)
        self.assertGreater(len(extra),0)
        states=list(evidence.states());best=min(states,key=evidence.score)
        for state in states:
            values=model.assignment(state)
            for terms,sense,rhs in extra:self.assertLessEqual(sum(values[v]*coef for v,coef in terms.items()),rhs)
        with tempfile.TemporaryDirectory() as d:
            f,l,c=[Path(d)/n for n in ('fixture.json','layout.json','certificate.json')]
            f.write_text(json.dumps(raw));l.write_text(json.dumps(best))
            summary=run('check',f,l,c,5,31,20)
            self.assertLessEqual(summary['lower'],evidence.score(best))
            self.assertEqual(run('verify',f,l,c),summary)

    def test_inherited_witness_is_monotone_on_both_children(self):
        model=IntegerModel(Evidence(fixture()));rows=canonical_rows(model)
        weights=[str(i+1)+'/7' if sense=='<=' else '-2/3' for i,(_,sense,_) in enumerate(rows)]
        parent=bound(model,rows,{},weights)
        for bit in (0,1):self.assertGreaterEqual(bound(model,rows,{model.binary[0]:bit},weights),parent)

    def test_generate_replay_and_tamper_rejection(self):
        with tempfile.TemporaryDirectory() as d:
            f,l,c=[Path(d)/n for n in ('fixture.json','layout.json','certificate.json')]
            raw=fixture();f.write_text(json.dumps(raw));l.write_text(json.dumps(next(Evidence(raw).states())))
            self.assertEqual(run('check',f,l,c,5,31,10)['gap'],0)
            run('verify',f,l,c)
            saved=json.loads(c.read_text());saved['nodes'][0]['lower']+=1;c.write_text(json.dumps(saved))
            with self.assertRaises(ValueError):run('verify',f,l,c)

if __name__=='__main__':unittest.main()
