import hashlib,json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import simple_papers as p
class SourceTests(unittest.TestCase):
    def test_bad_download_never_committed(self):
        with tempfile.TemporaryDirectory() as d:
            dest=Path(d)/'data'
            def fake(cmd,**kw):Path(cmd[cmd.index('--output')+1]).write_bytes(b'wrong')
            with patch.object(p.subprocess,'run',side_effect=fake):
                with self.assertRaises(ValueError):p.fetch('https://example.invalid/data',dest,hashlib.md5(b'correct').hexdigest())
            self.assertFalse(dest.exists());self.assertFalse(dest.with_suffix('.partial').exists())
    def test_verified_cache_reused(self):
        with tempfile.TemporaryDirectory() as d:
            dest=Path(d)/'data';dest.write_bytes(b'correct')
            with patch.object(p.subprocess,'run') as run:
                p.fetch('https://example.invalid/data',dest,hashlib.md5(b'correct').hexdigest());run.assert_not_called()
    def test_collector_does_not_claim_missing_comparisons(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)/'run';root.mkdir()
            p.collect(root,[{'id':'example'}])
            rows=json.loads((root/'summary.json').read_text());self.assertEqual([r['status'] for r in rows],['NOT_RUN']*3)
            import zipfile
            with zipfile.ZipFile(root/'SynTangle_simple_papers.zip') as z:self.assertNotIn('SynTangle_simple_papers.zip',z.namelist())
if __name__=='__main__':unittest.main()
