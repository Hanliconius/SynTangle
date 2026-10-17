import csv,json,tempfile,unittest
from pathlib import Path
from syntangle.fixtures import fixture_from_dict
from syntangle.layout import initial_layout_state,score_crossings
from syntangle.genespace import import_layout,flip_assist

class RefinementTests(unittest.TestCase):
    def fixture(self):
        return fixture_from_dict(dict(fixture_version=1,id='flip',title='test',purpose='test',species=[dict(id=sp,chromosomes=[dict(id='1',length=100,blocks=[dict(occurrence_id=sp+h,homology_id=h,start=x,end=x+5,strand='+') for h,x in [('h',10 if sp=='A' else 80),('j',80 if sp=='A' else 10)]])]) for sp in ('A','B')]))
    def test_flip_improves_and_preserves_order_and_evidence(self):
        f=self.fixture();before=repr(f);start=initial_layout_state(f);out=flip_assist(f,start,1)
        self.assertEqual(out.chromosome_order,start.chromosome_order)
        self.assertEqual(score_crossings(f,out).crossings,0)
        self.assertEqual(repr(f),before)
    def test_import_missing_duplicate_or_invalid_sign_fails(self):
        f=self.fixture()
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'layout.tsv';header='species\tchromosome\trank\torientation\n'
            for body in ['A\t1\t1\t1\n','A\t1\t1\t1\nA\t1\t2\t1\nB\t1\t1\t1\n','A\t1\t1\t0\nB\t1\t1\t1\n']:
                p.write_text(header+body)
                with self.assertRaises(ValueError):import_layout(f,p)
            p.write_text(header+'A\t1\t1\t1\nB\t1\t1\t-1\n')
            self.assertEqual(len(import_layout(f,p).chromosome_orientation),2)
    def test_coupled_flips_respect_equations(self):
        from dataclasses import replace
        from syntangle.model import OrientationConstraint
        from syntangle.orientation_space import orientation_basis
        from syntangle.saved_layout import validate_saved_layout
        f=self.fixture();a,b=f.chromosome_refs
        f=replace(f,orientation_constraints=(OrientationConstraint(a,b,0),))
        start=initial_layout_state(f);out=flip_assist(f,start,1)
        validate_saved_layout(f,out,orientation_basis(f,f.chromosome_refs))
        self.assertEqual(out.chromosome_orientation[a],out.chromosome_orientation[b])
if __name__=='__main__':unittest.main()
