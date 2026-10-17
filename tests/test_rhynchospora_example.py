"""Regression checks for the shipped real-data example and archived layouts."""
import json
from pathlib import Path
import unittest
from syntangle.table_input import load_tables
from syntangle.layout import initial_layout_state,score_crossings
from syntangle.saved_layout import decode_saved_layout,validate_saved_layout
from syntangle.orientation_space import orientation_basis
from syntangle.genespace import flip_assist

BASE=Path(__file__).resolve().parents[1]/'examples/rhynchospora'
class RhynchosporaExampleTests(unittest.TestCase):
    def test_saved_layouts_and_fixed_order_helper(self):
        fixture,_=load_tables(BASE/'chromosomes.tsv',BASE/'links.tsv')
        audit=json.loads((BASE/'reference_layouts.json').read_text())
        basis=orientation_basis(fixture,fixture.chromosome_refs)
        self.assertEqual(len(fixture.chromosomes),10)
        for item,expected in zip(audit['panels'],[55,13,9,9]):
            state=decode_saved_layout(fixture,item['state'])
            validate_saved_layout(fixture,state,basis)
            score=score_crossings(fixture,state)
            self.assertEqual(score.crossings,expected)
            self.assertEqual(score.link_count,37)
        start=initial_layout_state(fixture)
        result=flip_assist(fixture,start,seconds=5)
        self.assertEqual(result.chromosome_order,start.chromosome_order)
        self.assertEqual(score_crossings(fixture,result).crossings,9)
        changed=[r.label for r in fixture.chromosome_refs if result.chromosome_orientation[r]!=start.chromosome_orientation[r]]
        self.assertEqual(changed,['R_breviuscula:Chr2_h1'])
if __name__=='__main__':unittest.main()
