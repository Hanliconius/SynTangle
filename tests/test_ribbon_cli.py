import csv
import json
import tempfile
import unittest
from pathlib import Path
from syntangle.cli import main
from syntangle.genespace import import_layout
from syntangle.fixtures import load_fixture

class RibbonCliTests(unittest.TestCase):
    def test_exported_layout_and_audit_round_trip(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = root / 'fixture.json'
            fixture.write_text(json.dumps(dict(fixture_version=1,id='cli',title='t',purpose='t',species=[dict(id=sp,chromosomes=[dict(id='1',length=100,blocks=[dict(occurrence_id=sp,homology_id='h',start=10,end=80,strand='+')])]) for sp in ('a','b')])))
            layout = root / 'input.tsv'
            layout.write_text('species\tchromosome\trank\torientation\na\t1\t1\t-1\nb\t1\t1\t1\n')
            before = fixture.read_bytes(), layout.read_bytes()
            for objective in ('crossing-cap', 'combined'):
                out = root / objective
                self.assertEqual(main(['refine-ribbons',str(fixture),'--layout',str(layout),'--output',str(out),'--objective',objective,'--seconds','2']),0)
                audit = json.loads((out/'audit.json').read_text())
                self.assertEqual(audit['final_metrics']['twists'],0)
                self.assertEqual(audit['objective'],objective)
                self.assertEqual(len(audit['fixture_sha256']),64)
                self.assertEqual(len(import_layout(load_fixture(fixture),out/'layout.tsv').chromosome_orientation),2)
            self.assertEqual(before,(fixture.read_bytes(),layout.read_bytes()))
