from pathlib import Path
import tempfile
import unittest
from rhynchospora_te_audit import coding_genes,overlap,read_te


class TEAuditTests(unittest.TestCase):
    def test_duplicate_repeats_not_double_counted_and_introns_not_counted(self):
        self.assertEqual(overlap([(0,10),(90,100)],[(5,95),(5,95)]),10)
        self.assertEqual(overlap([(0,10),(90,100)],[(20,80)]),0)

    def test_cds_uses_parent_gene_and_representative_transcript(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'genes.gff';p.write_text('c\tx\tgene\t1\t100\t.\t+\t.\tID=g\n'
                'c\tx\tmRNA\t1\t100\t.\t+\t.\tID=t;Parent=g\n'
                'c\tx\tCDS\t1\t10\t.\t+\t0\tParent=t\n'
                'c\tx\tCDS\t91\t100\t.\t+\t0\tParent=t\n')
            result=coding_genes(p)
            self.assertEqual(result['g']['cds'],[(0,10),(90,100)])

    def test_satellite_and_unknown_not_called_te_and_names_fail_closed(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'repeats.bed';p.write_text('c\t0\t10\tLTR/Gypsy\nc\t20\t30\tSatellite/Tyba\nc\t40\t50\tUnknown\n')
            spans,counts=read_te(p,{'c':100})
            self.assertEqual(spans,{'c':[(0,10)]})
            self.assertEqual(counts['unclassified_rows'],1)
            with self.assertRaises(ValueError):read_te(p,{'other':100})


if __name__=='__main__':unittest.main()
