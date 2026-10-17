import tempfile
from pathlib import Path
import unittest
from rhynchospora_native_regions import annotations,import_regions,IDS

class NativeTests(unittest.TestCase):
    def test_parent_mapping_and_repeat_cds_id(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a.gff'
            p.write_text('c\t.\tgene\t1\t100\t.\t+\t.\tID=g\nc\t.\tmRNA\t1\t100\t.\t+\t.\tID=t;Parent=g\nc\t.\tCDS\t1\t10\t.\t+\t0\tID=cds;Parent=t\nc\t.\tCDS\t20\t30\t.\t+\t0\tID=cds;Parent=t\n')
            rec,tx,rows=annotations(p,{'c'})
            self.assertEqual(tx,{'t':'g'});self.assertEqual(len(rows),4)

    def test_exact_native_row_import_retains_lengths_and_reference(self):
        species=[dict(id=sp,chromosomes=[dict(id='Chr1_h1',length=1000,blocks=[])]) for sp in IDS]
        row=dict(genome1=IDS[0],genome2=IDS[1],chr1='Chr1_h1',chr2='Chr1_h1',startBp1='10',endBp1='20',startBp2='90',endBp2='70',refChr='Chr2_h1',color='#F9AC60')
        data,colours,refs,audit=import_regions([row],species)
        self.assertEqual(data['species'][0]['chromosomes'][0]['blocks'][0]['start'],9)
        self.assertEqual(data['species'][1]['chromosomes'][0]['blocks'][0]['start'],69)
        self.assertEqual(data['species'][2]['chromosomes'][0]['length'],1000)
        from syntangle.fixtures import fixture_from_dict
        fixture_from_dict(data)
        self.assertEqual(data['species'][1]['chromosomes'][0]['blocks'][0]['strand'],'-')
        self.assertEqual(len(audit),1);self.assertEqual(list(refs.values()),['Chr2_h1'])
        bad=dict(row,startBp1='0')
        with self.assertRaises(ValueError):import_regions([bad],species)
        with self.assertRaises(ValueError):import_regions([dict(row,genome2=IDS[2])],species)

if __name__=='__main__':unittest.main()
