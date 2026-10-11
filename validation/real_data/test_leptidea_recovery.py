import copy,unittest
from unittest.mock import patch
from recover_leptidea_chromosomes import select,LEPTIDEA
class RecoveryTests(unittest.TestCase):
    def rows(self):
        name='LrealiM'
        return [dict(assemblyAccession=LEPTIDEA[name],sequenceName=f'HiC_scaffold_{i}',genbankAccession=f'A{i}',length=n) for i,n in enumerate((5,4,1),1)]
    def test_known_paper_count_fraction_and_descending_series(self):
        with patch('recover_leptidea_chromosomes.TABLE',{'LrealiM':(2,3,.9)}):
            result=select(self.rows(),'LrealiM');self.assertEqual(len(result['chromosomes']),2);self.assertEqual(result['selected_fraction'],.9)
    def test_numbering_need_not_follow_length_order(self):
        with patch('recover_leptidea_chromosomes.TABLE',{'LrealiM':(2,3,.9)}):
            rows=self.rows();rows[0]['length']=4;rows[1]['length']=5
            selected=select(rows,'LrealiM');self.assertFalse(selected['numbered_lengths_descending'])
            self.assertEqual([r['sequence_name'] for r in selected['chromosomes']],['HiC_scaffold_1','HiC_scaffold_2'])
    def test_wrong_accession_names_lengths_and_counts_rejected(self):
        with patch('recover_leptidea_chromosomes.TABLE',{'LrealiM':(2,3,.9)}):
            for key,value in [('assemblyAccession','wrong'),('sequenceName','other'),('length',20)]:
                rows=self.rows();rows[0][key]=value
                with self.assertRaises(ValueError):select(rows,'LrealiM')
            with self.assertRaises(ValueError):select(self.rows()[:-1],'LrealiM')
if __name__=='__main__':unittest.main()
