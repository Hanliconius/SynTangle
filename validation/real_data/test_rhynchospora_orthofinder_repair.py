import tempfile,unittest
from pathlib import Path
from repair_rhynchospora_orthofinder import IDS,assignment,read_ids,validate_results
class RepairTests(unittest.TestCase):
    def data(self,root):
        peptides=root/'peptide';peptides.mkdir()
        for i,sp in enumerate(IDS):(peptides/(sp+'.fa')).write_text(f'>{sp}_g1\nACDE\n')
        og=root/'result/Orthogroups';og.mkdir(parents=True)
        head='Orthogroup\t'+'\t'.join(IDS)+'\n'
        (og/'Orthogroups.tsv').write_text(head+'OG0001\t'+'\t'.join(sp+'_g1' for sp in IDS)+'\n')
        (og/'Orthogroups_UnassignedGenes.tsv').write_text(head)
        tree=root/'result/Species_Tree';tree.mkdir();(tree/'SpeciesTree_rooted.txt').write_text('('+','.join(IDS)+');\n')
        hog=root/'result/Phylogenetic_Hierarchical_Orthogroups';hog.mkdir();(hog/'N0.tsv').write_text('HOG\tOG\tGene Tree Parent Clade\t'+'\t'.join(IDS)+'\nN0.HOG000\tOG0001\tN0\t'+'\t'.join(sp+'_g1' for sp in IDS)+'\n')
        return root/'result',peptides
    def test_complete_assignment_and_final_outputs_required(self):
        with tempfile.TemporaryDirectory() as d:
            result,peps=self.data(Path(d));self.assertEqual(validate_results(result,peps)['groups'],1)
            (result/'Phylogenetic_Hierarchical_Orthogroups/N0.tsv').unlink()
            with self.assertRaises(ValueError):validate_results(result,peps)
    def test_duplicate_or_unknown_peptides_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            result,peps=self.data(Path(d));p=result/'Orthogroups/Orthogroups.tsv'
            p.write_text(p.read_text().replace('OG0001','OG0001')+'OG0002\t'+IDS[0]+'_g1\t\t\n')
            with self.assertRaises(ValueError):validate_results(result,peps)
    def test_species_cache_exact_match(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'SpeciesIDs.txt';p.write_text(''.join(f'{i}: {sp}.fa\n' for i,sp in enumerate(IDS)))
            self.assertEqual(set(read_ids(p).values()),set(IDS));p.write_text('0: other.fa\n')
            with self.assertRaises(ValueError):read_ids(p)
if __name__=='__main__':unittest.main()
