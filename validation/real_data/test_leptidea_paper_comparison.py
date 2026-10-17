import unittest
from leptidea_paper_comparison import best_proteins
class TranslationTests(unittest.TestCase):
    def test_longest_valid_isoform_and_rejected_internal_stop(self):
        best,rejected=best_proteins([('a','MA*'),('b','MAAA'),('c','MA*AA'),('d','MAAAAA')],{'a':'g','b':'g','c':'h','d':'excluded'},{'g','h'})
        self.assertEqual(best,{'g':('b','MAAA')});self.assertEqual(rejected,1)
    def test_unknown_and_duplicate_transcript_rejected(self):
        for peptides in [[('unknown','MA')],[('a','MA'),('a','MA')]]:
            with self.assertRaises(ValueError):best_proteins(peptides,{'a':'g'},{'g'})
if __name__=='__main__':unittest.main()
