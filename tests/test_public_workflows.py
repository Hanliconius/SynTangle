import json
import tempfile
import unittest
from pathlib import Path
from syntangle.cli import main
from syntangle.table_input import load_tables
from syntangle.public_plot import directed_endpoints
from syntangle.fixtures import fixture_from_dict

EXAMPLE = Path(__file__).resolve().parents[1]/'examples'/'quickstart'


class PublicWorkflowTests(unittest.TestCase):
    def test_two_workflows_and_statistics(self):
        with tempfile.TemporaryDirectory() as tmp:
            for command,expected in [('flip-assist',8),('run',0)]:
                out=Path(tmp)/command
                args=[command,'--chromosomes',str(EXAMPLE/'chromosomes.tsv'),'--links',str(EXAMPLE/'links.tsv'),'--output',str(out),'--seconds','5']
                self.assertEqual(main(args),0)
                stats=json.loads((out/'stats.json').read_text())
                self.assertEqual(stats['input']['crossings'],10)
                self.assertEqual(stats['output']['crossings'],expected)
                self.assertEqual(stats['input']['link_count'],stats['output']['link_count'])
                self.assertEqual(stats['output']['link_count'],8)
                self.assertFalse(stats['independent_certificate'])
                self.assertEqual(len((out/'layout.tsv').read_text().splitlines()),7)
                if command=='flip-assist':
                    rows=lambda p:[line.rsplit('\t',1)[0] for line in p.read_text().splitlines()]
                    self.assertEqual(rows(out/'layout.tsv'),rows(out/'input_layout.tsv'))
                with self.assertRaises(ValueError):main(args)

    def test_invalid_evidence_is_not_filtered(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'links.tsv'
            source=(EXAMPLE/'links.tsv').read_text()
            for text in [source.replace('ab1\tA\t1','ab1\tA\tmissing'),source.replace('ab2\t','ab1\t'),source.replace('ab1\tA\t1\t10','ab1\tA\t1\tnan')]:
                path.write_text(text)
                with self.assertRaises(ValueError):load_tables(EXAMPLE/'chromosomes.tsv',path)

    def test_ribbon_cap_and_vector_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'result'
            main(['run','--chromosomes',str(EXAMPLE/'chromosomes.tsv'),'--links',str(EXAMPLE/'links.tsv'),'--output',str(out),'--seconds','3','--ribbon-objective','crossing-cap','--ribbon-seconds','3','--figure'])
            stats=json.loads((out/'stats.json').read_text())
            self.assertEqual(stats['output']['crossings'],0)
            self.assertEqual(stats['output']['twists'],0)
            self.assertTrue((out/'run_comparison.pdf').read_bytes().startswith(b'%PDF'))

    def test_directed_endpoint_whole_flip(self):
        fixture,_=load_tables(EXAMPLE/'chromosomes.tsv',EXAMPLE/'links.tsv')
        for c in fixture.chromosomes:
            for b in c.blocks:
                a=directed_endpoints(c,b,1);z=directed_endpoints(c,b,-1)
                self.assertEqual(tuple(c.length-v for v in z),a)
                self.assertEqual(a[0]<a[1],b.strand=='+')

    def test_invalid_budget_does_not_write_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'result'
            with self.assertRaises(ValueError):main(['run','--seconds','nan','--output',str(out)])
            self.assertFalse(out.exists())
