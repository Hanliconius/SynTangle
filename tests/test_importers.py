import csv
import json
import tempfile
import unittest
from pathlib import Path
from syntangle.cli import main
from syntangle.table_input import load_tables


class ImporterTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.q=self.write('q.sizes','q1\t100\nq_unused\t50\n')
        self.t=self.write('t.sizes','t1\t100\n')
    def tearDown(self):self.tmp.cleanup()
    def write(self,name,text):
        p=self.root/name;p.write_text(text);return str(p)
    def args(self,command,out):
        return [command,'--query-species','query','--target-species','target','--query-sizes',self.q,'--target-sizes',self.t,'--output',str(self.root/out)]
    def test_paf_direction_and_explicit_filter_audit(self):
        paf=self.write('test.paf','q1\t100\t10\t20\t-\tt1\t100\t30\t40\t9\t10\t60\nq1\t100\t60\t70\t+\tt1\t100\t80\t90\t8\t10\t255\n')
        main(self.args('import-paf','all')+['--paf',paf])
        f,_=load_tables(self.root/'all/chromosomes.tsv',self.root/'all/links.tsv')
        self.assertEqual(len(f.chromosomes),3)
        with (self.root/'all/links.tsv').open() as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
        self.assertEqual((rows[0]['start2'],rows[0]['end2']),('40','30'))
        main(self.args('import-paf','filtered')+['--paf',paf,'--min-mapq','30'])
        audit=json.loads((self.root/'filtered/import_audit.json').read_text())
        self.assertEqual((audit['source_records'],audit['retained_links'],audit['excluded_records']),(2,1,1))
        self.assertIn('unknown',(self.root/'filtered/excluded.tsv').read_text())
    def test_paf_wrong_assembly_fails_before_output(self):
        paf=self.write('bad.paf','q1\t101\t10\t20\t+\tt1\t100\t30\t40\t9\t10\t60\n')
        with self.assertRaises(ValueError):main(self.args('import-paf','bad')+['--paf',paf])
        self.assertFalse((self.root/'bad').exists())
    def test_jcvi_envelopes_orientation_and_duplicate_rows(self):
        qbed=self.write('q.bed','q1\t10\t20\tqa\nq1\t60\t70\tqb\n')
        tbed=self.write('t.bed','t1\t20\t30\tta\nt1\t80\t90\ttb\n')
        simple=self.write('test.simple','StartGeneA EndGeneA StartGeneB EndGeneB Score Orientation\nqa qb ta tb 2 -\nred*qa qb ta tb 2 +\n')
        main(self.args('import-jcvi','jcvi')+['--simple',simple,'--query-bed',qbed,'--target-bed',tbed])
        with (self.root/'jcvi/links.tsv').open() as h:rows=list(csv.DictReader(h,delimiter='\t'))
        self.assertEqual(len(rows),2)
        self.assertEqual((rows[0]['start1'],rows[0]['end1'],rows[0]['start2'],rows[0]['end2']),('10','70','90','20'))
        self.assertEqual((rows[1]['start2'],rows[1]['end2']),('20','90'))
        self.assertNotEqual(rows[0]['link_id'],rows[1]['link_id'])
    def test_jcvi_missing_gene_is_not_dropped(self):
        qbed=self.write('q.bed','q1\t10\t20\tqa\n')
        tbed=self.write('t.bed','t1\t20\t30\tta\n')
        simple=self.write('test.simple','qa missing ta ta 1 +\n')
        with self.assertRaises(ValueError):main(self.args('import-jcvi','bad')+['--simple',simple,'--query-bed',qbed,'--target-bed',tbed])
        self.assertFalse((self.root/'bad').exists())
