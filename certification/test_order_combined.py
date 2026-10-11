import json,tempfile,unittest
from pathlib import Path
from test_proof import fixture
from proof import Evidence
from order_selective import run
class CombinedTests(unittest.TestCase):
    def test_frontier_keeps_exact_root_witness_and_rejects_tamper(self):
        with tempfile.TemporaryDirectory() as d:
            f,l,c=[Path(d)/n for n in ('f','l','c')];raw=fixture();f.write_text(json.dumps(raw));l.write_text(json.dumps(next(Evidence(raw).states())))
            summary=run('check',f,l,c,3,16,10,20,21);saved=json.loads(c.read_text())
            self.assertGreaterEqual(summary['lower'],saved['root_lower'])
            self.assertEqual(summary,run('verify',f,l,c))
            saved['root_lower']+=1;c.write_text(json.dumps(saved))
            with self.assertRaises(ValueError):run('verify',f,l,c)
    def test_zero_time_carries_valid_witness(self):
        with tempfile.TemporaryDirectory() as d:
            f,l,c=[Path(d)/n for n in ('f','l','c')];raw=fixture();f.write_text(json.dumps(raw));l.write_text(json.dumps(next(Evidence(raw).states())))
            self.assertEqual(run('check',f,l,c,0,16,10,20,21),run('verify',f,l,c))
if __name__=='__main__':unittest.main()
