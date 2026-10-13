import json,tempfile,unittest
from pathlib import Path
from simple_reconstruction import choose,gene_records,paf_rows,lengths

class ReconstructionTests(unittest.TestCase):
    def test_no_longest_n_substitution(self):
        spec={'id':'s','expected_chromosomes':2,'chromosome_regex':r'chr[12]'}
        with self.assertRaises(ValueError):choose({'chr1':10,'unplaced':1000000},spec)
    def test_parent_graph_and_fixed_coordinates(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a.gff';p.write_text('chr1\ts\tgene\t2\t90\t.\t-\t.\tID=g\nchr1\ts\tmRNA\t2\t90\t.\t-\t.\tID=t;Parent=g\nchr1\ts\tCDS\t3\t20\t.\t-\t0\tID=c;Parent=t;protein_id=p\n')
            records,mapping,_=gene_records(p,{'chr1':100})
            self.assertEqual(mapping['p'],'g');self.assertEqual(mapping['t'],'g');self.assertEqual(records['g'][0][3:7],['2','90','.','-'])
    def test_signed_alignment_no_internal_rearrangement(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a.paf';p.write_text('q\t100\t10\t40\t-\tt\t200\t80\t110\t30\t30\t60\ttp:A:P\nq\t100\t50\t60\t+\tt\t200\t140\t150\t10\t10\t60\ttp:A:S\n')
            rows=paf_rows(p,'a','b',{'q':{'length':100}},{'t':{'length':200}},0,20,20)
            self.assertEqual(len(rows),1);self.assertEqual(rows[0]['orient'],'-');self.assertEqual(rows[0]['startBp1'],11);self.assertEqual(rows[0]['endBp2'],110)
    def test_bad_report_length(self):
        report=[dict(role='assembled-molecule',assignedMoleculeLocationType='Chromosome',chrName='1',genbankAccession='a',length=20)]
        with self.assertRaises(ValueError):choose({'a':10},dict(id='s',expected_chromosomes=1),report)
    def test_scopes(self):
        m=json.loads(Path(__file__).with_name('simple_reconstruction_manifest.json').read_text());cs={c['id']:c for c in m['cases']}
        self.assertEqual(len(cs['faba']['species']),6);self.assertEqual(len(cs['rhodotorula']['species']),11)
        self.assertEqual(cs['ficus']['pairs'],[['Ficus_benjamina',s] for s in ['Ficus_microcarpa','Ficus_hispida','Ficus_hirta','Ficus_pumila']])
        self.assertTrue(cs['ficus']['species'][-1]['blocked'])
if __name__=='__main__':unittest.main()
