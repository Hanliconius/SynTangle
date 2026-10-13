import json,tempfile,unittest
from pathlib import Path
from simple_reconstruction import choose,gene_records,paf_rows,lengths,asset,fixture,sha,dump
from unittest.mock import patch

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
    def test_directory_cache_is_separate_per_url(self):
        with tempfile.TemporaryDirectory() as d:
            pages=[]
            def fake(url,path,md5=None):
                if url.endswith('/'):
                    pages.append(path)
                    path.write_text('<a href="'+url+('genes.gff.gz' if 'annotation' in url else 'dna.fna.gz')+'">file</a>')
                else:path.write_text('data')
                return {'url':url}
            with patch('simple_reconstruction.fetch',side_effect=fake):
                asset(dict(directory='https://x/genome/',filename_suffix='.fna.gz'),Path(d),None)
                p,_=asset(dict(directory='https://x/annotation/',filename_suffix='.gff.gz'),Path(d),None)
            self.assertEqual(p.name,'genes.gff.gz');self.assertNotEqual(pages[0],pages[1])
    def test_distinct_duplicate_is_not_silently_merged(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a.gff';line='chr1\ts\tgene\t2\t20\t.\t+\t.\tID=g\n'
            p.write_text(line+line)
            records,_,_=gene_records(p,{'chr1':100});self.assertEqual(len(records),1)
            p.write_text(line+line.replace('2\t20','30\t50'))
            records,_,_=gene_records(p,{'chr1':100});self.assertEqual(len(records),2)
    def test_repeated_gene_parent_resolved_by_placement(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a.gff'
            p.write_text('chr1\ts\tgene\t2\t20\t.\t+\t.\tID=g\nchr1\ts\tgene\t30\t50\t.\t-\t.\tID=g\nchr1\ts\tmRNA\t31\t49\t.\t-\t.\tID=t;Parent=g\nchr1\ts\tfive_prime_UTR\t31\t32\t.\t-\t.\tID=u;Parent=t\nchr1\ts\tfive_prime_UTR\t35\t36\t.\t-\t.\tID=u;Parent=t\n')
            records,mapping,_=gene_records(p,{'chr1':100})
            self.assertEqual(records[mapping['t']][0][3:5],['30','50'])
    def test_nonadjacent_anchors_do_not_change_scored_links(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);case=dict(id='test',species=[dict(id=x) for x in 'abc'],pairs=[['a','b'],['b','c'],['a','c']],paper='p',figure='f',method='minimap2',reconstruction_note='test',colour_reference='a')
            for sp in 'abc':dump(root/'sources'/sp/'chromosomes.json',{'chr1':dict(length=100,display_name='1')})
            for i,(a,b) in enumerate(case['pairs']):
                out=root/'pairs'/str(i);dump(out/'blocks.json',[dict(blkID=str(i),genome1=a,genome2=b,chr1='chr1',chr2='chr1',startBp1=11,endBp1=20,startBp2=21,endBp2=30,orient='+')]);(out/'COMPLETE').write_text(sha(out/'blocks.json'))
            prepared=fixture(root,case)
            provenance=json.loads((prepared/'provenance.json').read_text())
            self.assertEqual(provenance['retained_links'],2)
            self.assertEqual(len((prepared/'published_block_coordinates.tsv').read_text().splitlines()),4)
    def test_scopes(self):
        m=json.loads(Path(__file__).with_name('simple_reconstruction_manifest.json').read_text());cs={c['id']:c for c in m['cases']}
        self.assertEqual(len(cs['faba']['species']),6);self.assertEqual(len(cs['rhodotorula']['species']),11)
        self.assertEqual(cs['ficus']['pairs'],[['Ficus_benjamina',s] for s in ['Ficus_microcarpa','Ficus_hispida','Ficus_hirta','Ficus_pumila']])
        self.assertTrue(cs['ficus']['species'][-1]['blocked'])
if __name__=='__main__':unittest.main()
