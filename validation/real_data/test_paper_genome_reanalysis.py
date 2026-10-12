import unittest,tempfile,json
from pathlib import Path
from paper_genome_reanalysis import select_chromosomes,collect
from syntangle.fixtures import load_fixture
class PaperInputsTests(unittest.TestCase):
    def test_no_longest_n_fallback(self):
        with self.assertRaises(ValueError):select_chromosomes({'chr1':20,'chr2':5},dict(id='x',min_chromosome_length=10,expected_chromosomes=2))
    def test_published_littledalei_set(self):
        ids=[f'CM{i:06d}.1' for i in range(22079,22108)];lengths={i:100 for i in ids};lengths['unplaced']=1000
        rows=[dict(genbankAccession=i,length=100) for i in ids]
        chosen,_=select_chromosomes(lengths,dict(id='C',resolve_from_nucleotide='CM022079',expected_chromosomes=29),rows)
        self.assertEqual(set(chosen),set(ids));self.assertNotIn('unplaced',chosen)
    def test_signed_blocks_and_unlinked_chromosomes_retained(self):
        with tempfile.TemporaryDirectory() as d:
            run=Path(d);case=dict(id='test',species=[dict(id='a'),dict(id='b')],paper='test',figure='1',size_order=False,colour_reference='a',input_layout_note='test')
            for sp in ('a','b'):
                out=run/'sources'/sp;out.mkdir(parents=True)
                (out/'chromosomes.json').write_text(json.dumps({'chr1':dict(id='chr1',length=100),'chr2':dict(id='chr2',length=120)}))
                (out/'input_audit.json').write_text('{}')
            out=run/'pairs/0';out.mkdir(parents=True);(out/'COMPLETE').write_text('PASS')
            (out/'blocks.json').write_text(json.dumps([dict(blkID='link',orient='-',genome1='a',genome2='b',chr1='chr1',chr2='chr1',startBp1=11,endBp1=40,startBp2=51,endBp2=80)]))
            collect(run,case);f=load_fixture(run/'prepared/test/fixture.json')
            self.assertEqual(len(f.chromosomes),4)
            blocks=[b for c in f.chromosomes for b in c.blocks]
            self.assertEqual([(b.start,b.end,b.strand) for b in blocks],[(10,40,'+'),(50,80,'-')])
if __name__=='__main__':unittest.main()
