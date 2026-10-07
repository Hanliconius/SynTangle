import importlib.util
from pathlib import Path
import unittest
spec=importlib.util.spec_from_file_location('block_import',Path(__file__).resolve().parents[1]/'validation/real_data/prepare_genespace_blocks.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def row(a='A',b='B',c1='1',c2='2',lo1='1',hi1='10',lo2='20',hi2='11'):
    return dict(gen1=a,gen2=b,chr1=c1,chr2=c2,startBp1=lo1,endBp1=hi1,startBp2=lo2,endBp2=hi2,orient='-',blkID='x')

class ImportTests(unittest.TestCase):
    def test_reciprocal_reverse_and_duplicate_links_preserved(self):
        direct=row();reverse=row('B','A','2','1','20','11','1','10')
        data,p=m.convert([direct,reverse,direct,reverse],['A','B'],'test')
        self.assertEqual(p['retained_links'],2)
        self.assertEqual(len(data['species'][0]['chromosomes'][0]['blocks']),2)
        self.assertEqual(len({x['homology_id'] for x in p['edges']}),2)
        self.assertEqual(data['species'][1]['chromosomes'][0]['blocks'][0]['start'],11)
        self.assertEqual(data['orientation_constraints'],[])
    def test_reciprocal_difference_fails_closed(self):
        with self.assertRaises(ValueError):m.convert([row(),row('B','A','2','1','20','11','1','9')],['A','B'],'test')
    def test_missing_pair_fails(self):
        with self.assertRaises(ValueError):m.convert([row()],['A','B'],'test')
    def test_invalid_coordinates_fail(self):
        with self.assertRaises(ValueError):m.signature(row(lo1='nan'))
