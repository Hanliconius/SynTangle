import unittest
from rhynchospora_blocks import build


class BlockTests(unittest.TestCase):
    def data(self,order):
        species=[]
        for sp,ids in [('a',list(range(len(order)))),('b',order)]:
            species.append(dict(id=sp,chromosomes=[dict(id='c',length=1000,blocks=[dict(homology_id=str(h),start=i*10,end=i*10+5,strand='+') for i,h in enumerate(ids)])]))
        return dict(id='test',species=species)

    def test_forward_and_reverse_runs(self):
        for order in ([0,1,2,3],[3,2,1,0]):
            out,audit=build(self.data(order))
            self.assertEqual(len(audit),1)
            self.assertEqual(set(audit[0]['members']),{'0','1','2','3'})

    def test_direction_change_and_singletons_preserved(self):
        out,audit=build(self.data([0,2,1,3]))
        members=[h for r in audit for h in r['members']]
        self.assertEqual(sorted(members),['0','1','2','3'])
        self.assertEqual(len(members),len(set(members)))
        self.assertGreater(len(audit),1)
        self.assertEqual(len(self.data([0,2,1,3])['species'][0]['chromosomes'][0]['blocks']),4)


if __name__=='__main__':unittest.main()
