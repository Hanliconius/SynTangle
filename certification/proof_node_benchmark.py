"""Sequential single-core proof comparison; reads saved real fixtures only."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);p.add_argument('--seconds',type=float,default=90);a=p.parse_args()
    code=Path(__file__).resolve().parent;report=[]
    variants=[('previous', 'order_certificate.py',101,[]),('frontier','order_strengthened.py',101,['--cycles','0']),('cycle_root','order_strengthened.py',1,['--cycles','500']),('cycle_frontier','order_strengthened.py',101,['--cycles','500'])]
    for task,case in ((0,'planarians'),(1,'flatworms_fig4f')):
        root=a.run/f'task_{task}'
        for label,script,nodes,extra in variants:
            cert=root/(label+'.json');common=[str(root/'fixture.json'),str(root/'layout.json'),str(cert)]
            log=root/(label+'.log');started=time.monotonic()
            print(f'START {case} {label} cap={nodes}',flush=True)
            with log.open('w') as out:
                subprocess.run([sys.executable,str(code/script),'check',*common,'--seconds',str(a.seconds),'--max-nodes',str(nodes),*extra],stdout=out,stderr=subprocess.STDOUT,check=True)
                search_elapsed=time.monotonic()-started
                subprocess.run([sys.executable,str(code/script),'verify',*common],stdout=out,stderr=subprocess.STDOUT,check=True)
            saved=json.loads(cert.read_text());summary=saved['summary']
            row=dict(case=case,variant=label,seconds=search_elapsed,replay_seconds=time.monotonic()-started-search_elapsed,**summary,lp_calls=saved['lp_calls'],cuts=saved.get('cycle_cuts',0),rows=saved['rows'],verified=True)
            report.append(row);(a.run/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
            print(f"DONE {case} {label} L={row['lower']} U={row['upper']} gap={row['gap']} nodes={row['nodes']} LP={row['lp_calls']} cuts={row['cuts']} {search_elapsed:.3f}s replay=PASS",flush=True)
    (a.run/'COMPLETE').write_text('PASS all certificates independently replayed\n')

if __name__=='__main__':main()
