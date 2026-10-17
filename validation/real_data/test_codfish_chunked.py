import tempfile, unittest, json
from pathlib import Path
from unittest.mock import patch
import codfish_chunked as c
class ChunkTests(unittest.TestCase):
    def test_partition_and_checkpoint(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d); work=run/'work'; case={'species':[{'id':'a'},{'id':'b'}]}
            for s in ['a','b']:
                root=run/'sources'/s; root.mkdir(parents=True)
                (root/'COMPLETE').write_text('PASS'); (root/'proteins.faa').write_text(''.join(f'>{s}{i}\nAAAA\n' for i in range(9)))
            with patch.object(c.subprocess,'run'): c.setup(run,work,case,4)
            self.assertEqual(len(json.loads((work/'tasks.json').read_text())),8)
            ids=[r.splitlines()[0] for i in range(4) for r in c.records(work/f'a.{i}.faa')]
            self.assertEqual(len(ids),len(set(ids))); self.assertEqual(len(ids),9)
            calls=[]
            def fake(cmd,**kw):
                calls.append(cmd); Path(cmd[cmd.index('-out')+1]).write_text('hit\n')
            with patch.object(c.subprocess,'check_output',return_value='test'),patch.object(c.subprocess,'run',side_effect=fake):
                c.search(work,0); c.search(work,0)
                self.assertEqual(len(calls),1)
                (work/'0.tsv').write_text('tampered')
                c.search(work,0); self.assertEqual(len(calls),2)
    def test_failed_search_no_checkpoint(self):
        with tempfile.TemporaryDirectory() as d:
            w=Path(d); (w/'tasks.json').write_text(json.dumps([dict(query='a',target='b',chunk=0)]))
            with patch.object(c.subprocess,'check_output',return_value='test'),patch.object(c.subprocess,'run',side_effect=RuntimeError('killed')):
                with self.assertRaises(RuntimeError): c.search(w,0)
            self.assertFalse((w/'0.complete.json').exists())
if __name__=='__main__': unittest.main()
