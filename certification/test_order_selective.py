import json,tempfile,unittest
from pathlib import Path
from proof import Evidence
from test_proof import fixture
from order_selective import run,violations
class SelectiveTests(unittest.TestCase):
    def test_only_violated_unused_rows_selected(self):
        pool=[({'a':1},'<=',0),({'a':-1},'<=',0),({'a':2},'<=',0)]
        self.assertEqual(violations(pool,{'a':.5},[],1),[2]);self.assertEqual(violations(pool,{'a':.5},[2],5),[0])
    def test_exact_replay_and_tamper(self):
        with tempfile.TemporaryDirectory() as d:
            f,l,c=[Path(d)/n for n in ('f.json','l.json','c.json')]
            raw=fixture();f.write_text(json.dumps(raw));l.write_text(json.dumps(next(Evidence(raw).states())))
            summary=run('check',f,l,c,5,2,3,10);self.assertEqual(summary,run('verify',f,l,c))
            saved=json.loads(c.read_text());saved['selected']=[-1];c.write_text(json.dumps(saved))
            with self.assertRaises(ValueError):run('verify',f,l,c)
    def test_real_cycle_model_exhaustive_and_monotone(self):
        def species(sp,groups):
            return dict(id=sp,chromosomes=[dict(id=str(i),length=10,blocks=[dict(occurrence_id=sp+h,homology_id=h,start=j*2,end=j*2+1,strand='+') for j,h in enumerate(group)]) for i,group in enumerate(groups)])
        raw=dict(fixture_version=1,id='cycle_model',species=[species('A',[['a','e','b','f'],['c','g','d','h']]),species('B',[['a','e','c','g'],['b','f','d','h']])])
        evidence=Evidence(raw);best=min(evidence.states(),key=evidence.score)
        with tempfile.TemporaryDirectory() as d:
            f,l,c=[Path(d)/n for n in ('f.json','l.json','c.json')]
            f.write_text(json.dumps(raw));l.write_text(json.dumps(best))
            summary=run('check',f,l,c,5,1,20,20);saved=json.loads(c.read_text())
            history=[r['lower'] for r in saved['history']];self.assertEqual(history,sorted(history))
            self.assertEqual(summary['lower'],evidence.score(best));self.assertEqual(summary,run('verify',f,l,c))
    def test_zero_budget_valid_bound(self):
        with tempfile.TemporaryDirectory() as d:
            f,l,c=[Path(d)/n for n in ('f.json','l.json','c.json')]
            raw=fixture();f.write_text(json.dumps(raw));l.write_text(json.dumps(next(Evidence(raw).states())))
            summary=run('check',f,l,c,0,2,3,10);self.assertGreaterEqual(summary['gap'],0);self.assertEqual(summary,run('verify',f,l,c))
if __name__=='__main__':unittest.main()
