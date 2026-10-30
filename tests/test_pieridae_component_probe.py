"""Small runner/submission check; no Slurm or remote filesystem required."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

from syntangle.table_input import load_tables
from syntangle.layout import initial_layout_state
from syntangle.workflow import write_layout


class ProbeTests(unittest.TestCase):
    def test_audit_solve_and_frozen_worker(self):
        repo = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            source = base/'source'
            (source/'tables').mkdir(parents=True)
            (source/'optimization').mkdir()
            for name in ('chromosomes.tsv','links.tsv'):
                shutil.copy2(repo/'examples/quickstart'/name,source/'tables'/name)
            fixture, _ = load_tables(source/'tables/chromosomes.tsv',source/'tables/links.tsv')
            for name in ('input_layout.tsv','layout.tsv'):
                write_layout(source/'optimization'/name,initial_layout_state(fixture))
            def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
            (source/'optimization/COMPLETE.json').write_text(json.dumps(dict(
                tables_sha256={n:sha(source/'tables'/n) for n in ('chromosomes.tsv','links.tsv')},
                layout_sha256=sha(source/'optimization/layout.tsv'))))
            env = dict(os.environ,PYTHONPATH=str(repo/'src'))
            for mode in ('audit','solve','recover'):
                cmd = [sys.executable,str(repo/'validation/real_data/probe_pieridae_components.py'),
                       str(source),str(base/mode),'--seconds','2']
                if mode=='audit':cmd.append('--audit-only')
                if mode=='recover':cmd.append('--recover')
                subprocess.run(cmd,check=True,env=env,capture_output=True,text=True)
                audit = json.loads((base/mode/'audit.json').read_text())
                self.assertEqual(audit['presolved']['crossings'],0)
            binary = base/'bin'; binary.mkdir()
            (binary/'sbatch').write_text('#!/bin/bash\nprintf "12345\\n"\n')
            (binary/'git').write_text('#!/bin/bash\nif [[ "$*" == *rev-parse* ]]; then echo test-provenance; fi\n')
            for p in binary.iterdir():p.chmod(0o755)
            env.update(PATH=str(binary)+':'+env['PATH'],SCT_PYTHON=sys.executable)
            subprocess.run(['bash',str(repo/'validation/real_data/submit_pieridae_components.sh'),
                            str(source),'2','--audit-only'],check=True,env=env,capture_output=True,text=True)
            out = next(source.glob('component_test.*'))
            subprocess.run(['bash',str(out/'code/worker.sh'),str(source),str(out),sys.executable,
                            '2','--audit-only','no'],check=True,env=env,capture_output=True,text=True)
            self.assertTrue((out/'audit.json').exists())
            self.assertTrue((out/'code/src/syntangle/component_presolve.py').exists())


if __name__ == '__main__':
    unittest.main()
