"""Timed block comparison plus full gene rescoring of every saved panel."""
import argparse
import faulthandler
import json
from pathlib import Path
import sys
import time
import plot_comparison as plot
from syntangle.fixtures import load_fixture
from syntangle.layout import score_crossings
from syntangle.saved_layout import decode_saved_layout


def timed(name,function):
    def call(*args,**kwargs):
        start=time.monotonic();print(f'START {name}',flush=True)
        result=function(*args,**kwargs)
        print(f'DONE {name} seconds={time.monotonic()-start:.3f}',flush=True)
        return result
    return call


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--prepared',type=Path,required=True);p.add_argument('--genes',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--seconds',type=float,default=300)
    a=p.parse_args()
    faulthandler.dump_traceback_later(120,repeat=True)
    for name in ('render','export_native_input','load_native_states','improve_flips','run_hybrid'):
        setattr(plot,name,timed(name,getattr(plot,name)))
    sys.argv=[sys.argv[0],'--prepared',str(a.prepared),'--output',str(a.output),'--seconds',str(a.seconds)]
    plot.main()
    (a.output/'COMPLETE').unlink()
    genes=load_fixture(a.genes/'fixture.json')
    audit=json.loads((a.output/'comparison_audit.json').read_text())
    scores=[]
    for panel in audit['panels']:
        state=decode_saved_layout(genes,panel['state'])
        c=score_crossings(genes,state).crossings
        scores.append(dict(title=panel['title'],run_crossings=panel['crossings'],full_gene_crossings=c))
        print(f'RESCORE {panel["title"]}: runs={panel["crossings"]} genes={c}',flush=True)
    plot.atomic_json(a.output/'full_gene_rescore.json',dict(panels=scores,
        warning='Run optimum, if certified, is not a gene-crossing optimum. Both objectives reported separately.'))
    (a.output/'COMPLETE').write_text('PASS\n')
    faulthandler.cancel_dump_traceback_later()


if __name__=='__main__':main()
