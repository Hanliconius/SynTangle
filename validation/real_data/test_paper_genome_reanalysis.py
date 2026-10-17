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

class AnnotationMappingTests(unittest.TestCase):
    def test_zero_padding_keeps_signed_coordinates_and_attributes(self):
        from paper_genome_reanalysis import normalize_annotation
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);src=root/'gff';dst=root/'mapped'
            line='01\tmaker\tCDS\t10\t20\t.\t-\t2\tID=c;Parent=t\n';src.write_text(line)
            report=[dict(chrName='1',refseqAccession='NC_044048.1',genbankAccession='LR633943.1',sequenceName='Super_Scaffold_20')]
            audit=normalize_annotation(src,dst,{'NC_044048.1':100},report,True)
            self.assertEqual(dst.read_text(),'NC_044048.1'+line[2:]);self.assertTrue(audit['coordinates_strands_attributes_unchanged'])
            self.assertEqual(audit['sequence_mapping'],{'01':'NC_044048.1'})
    def test_ambiguous_numeric_report_is_rejected(self):
        from paper_genome_reanalysis import normalize_annotation
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);src=root/'gff';src.write_text('01\tx\tCDS\t1\t2\t.\t+\t0\tParent=t\n')
            report=[dict(chrName='1',refseqAccession='A'),dict(chrName='01',refseqAccession='B')]
            with self.assertRaisesRegex(ValueError,'Ambiguous sequence alias'):
                normalize_annotation(src,root/'mapped',{'A':100,'B':100},report,True)
    def test_out_of_bounds_alias_is_rejected(self):
        from paper_genome_reanalysis import normalize_annotation
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);src=root/'gff';src.write_text('01\tx\tCDS\t1\t101\t.\t+\t0\tParent=t\n')
            with self.assertRaisesRegex(ValueError,'outside mapped'):
                normalize_annotation(src,root/'mapped',{'A':100},[dict(chrName='1',refseqAccession='A')],True)

if __name__=='__main__':unittest.main()
