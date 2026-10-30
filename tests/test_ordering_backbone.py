import unittest
from syntangle.ordering_backbone import select_links

def link(i,chrom='B',anchor=None):
 return dict(link_id=str(i),buscoID=str(i) if anchor is None else anchor,species1='s1',species2='s2',chromosome1='A',chromosome2=chrom,start1=i,end1=i+.1,start2=i,end2=i+.1)
class Tests(unittest.TestCase):
 def test_singleton_breaks_runs_and_all_rows_audited(self):
  rows=[link(i,'C' if i==3 else 'B') for i in range(7)]
  kept,audit=select_links(rows)
  self.assertEqual([r['link_id'] for r in kept],['0','1','2','4','5','6'])
  self.assertEqual(len(audit),7);self.assertTrue(all(r['retained_for_final_plot'] for r in audit))
 def test_copies_do_not_inflate_support(self):
  rows=[link(i,anchor='one') for i in range(5)]
  self.assertEqual(select_links(rows)[0],[])
 def test_short_fragment_is_not_rescued_by_total_pair_support(self):
  rows=[link(i,'C' if i in (2,5) else 'B') for i in range(8)]
  self.assertEqual(select_links(rows)[0],[])
 def test_input_unchanged_and_threshold_one(self):
  rows=[link(0)];before=repr(rows)
  self.assertEqual(select_links(rows,1,1)[0],rows);self.assertEqual(repr(rows),before)
if __name__=='__main__':unittest.main()
