import json,tempfile,unittest
from pathlib import Path
from test_proof import fixture
from proof import Evidence
from order_selective import run as selective
from order_resume import run
class ResumeTests(unittest.TestCase):
    def test_zero_budget_keeps_complete_tree_and_bound(self):
        with tempfile.TemporaryDirectory() as d:
            f,l,c,o=[Path(d)/n for n in ('f','l','c','o')];raw=fixture();f.write_text(json.dumps(raw));l.write_text(json.dumps(next(Evidence(raw).states())))
            old=selective('check',f,l,c,0,16,10,20,21)
            before=json.loads(c.read_text());new=run(f,l,c,o,0,101);after=json.loads(o.read_text())
            self.assertEqual(new,old);self.assertEqual(before['nodes'],after['nodes']);self.assertEqual(before['selected'],after['selected'])
    def test_replay_before_resume_rejects_changed_tree(self):
        with tempfile.TemporaryDirectory() as d:
            f,l,c,o=[Path(d)/n for n in ('f','l','c','o')];raw=fixture();f.write_text(json.dumps(raw));l.write_text(json.dumps(next(Evidence(raw).states())))
            selective('check',f,l,c,0,16,10,20,21);s=json.loads(c.read_text());s['nodes'][0]['lower']+=1;c.write_text(json.dumps(s))
            with self.assertRaises(ValueError):run(f,l,c,o,0,101)
            self.assertFalse(o.exists())
    def test_additional_search_keeps_old_bound(self):
        with tempfile.TemporaryDirectory() as d:
            f,l,c,o=[Path(d)/n for n in ('f','l','c','o')];raw=fixture();f.write_text(json.dumps(raw));l.write_text(json.dumps(next(Evidence(raw).states())))
            old=selective('check',f,l,c,0,16,10,20,21);new=run(f,l,c,o,2,101)
            self.assertGreaterEqual(new['lower'],old['lower']);self.assertEqual(new['upper'],old['upper'])
if __name__=='__main__':unittest.main()
