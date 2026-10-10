import tempfile
from pathlib import Path
import unittest
from rhynchospora_comparison import strict

class TransferImportTests(unittest.TestCase):
    def test_ambiguous_low_quality_unknown_and_nonfinite_are_excluded(self):
        records=[('valid','coverage=0.9;sequence_ID=0.99'),
                 ('duplicate','coverage=1;sequence_ID=1'),
                 ('duplicate','coverage=1;sequence_ID=1'),
                 ('poor','coverage=0.89;sequence_ID=1'),
                 ('nan','coverage=nan;sequence_ID=1'),
                 ('unknown','coverage=1;sequence_ID=1')]
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'genes.gff'
            path.write_text(''.join(f'chr1\tLiftoff\tgene\t1\t10\t.\t+\t.\tID={name};{attrs}\n' for name,attrs in records))
            retained,excluded=strict(path,{'valid','duplicate','poor','nan'},{'chr1':100})
            self.assertEqual(set(retained),{'valid'})
            self.assertEqual(sum(excluded.values()),5)

    def test_outside_target_coordinates_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'genes.gff'
            path.write_text('chr1\tLiftoff\tgene\t1\t101\t.\t+\t.\tID=g;coverage=1;sequence_ID=1\n')
            with self.assertRaises(ValueError):strict(path,{'g'},{'chr1':100})

if __name__=='__main__':unittest.main()
