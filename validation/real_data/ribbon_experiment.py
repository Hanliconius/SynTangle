"""Run opt-in twist objectives from an archived complete midpoint layout."""
import argparse,json,shutil,hashlib
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from syntangle.fixtures import load_fixture
from syntangle.saved_layout import decode_saved_layout
from syntangle.ribbon_objective import optimize_ribbons,ribbon_metrics
CONFIGS=[('crossing_cap_count','lex',False),('crossing_cap_width','lex',True),('combined_count','weighted',False),('combined_width','weighted',True)]
def main():
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['run','collect']);p.add_argument('--run',type=Path,required=True);p.add_argument('--index',type=int,default=0);a=p.parse_args()
 if a.stage=='collect':
  results=[]
  for name,_,_ in CONFIGS:
   d=a.run/name;result=json.loads((d/'result.json').read_text());results.append(dict(case=name,**{k:result[k] for k in ['starting_metrics','final_metrics','solver_status','seconds']}))
  (a.run/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
  with ZipFile(a.run/'SynTangle_ribbon_objectives.zip','w',ZIP_DEFLATED) as z:
   for d in [a.run/n for n,_,_ in CONFIGS]:
    for f in d.iterdir():
     if f.is_file():z.write(f,str(f.relative_to(a.run)))
   z.write(a.run/'summary.json','summary.json')
  print(json.dumps(results,indent=2));print('ARCHIVE',a.run/'SynTangle_ribbon_objectives.zip');return
 from plot_comparison import render
 from preservation_audit import require_audit
 source=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/SynTangle_pipeline_test/local_results/simple_reconstructed_QI1cVs/rhodotorula')
 prepared=source/'prepared/rhodotorula_600';saved=source/'figures/rhodotorula_600_600/comparison_audit.json'
 fixture=load_fixture(prepared/'fixture.json');audit=json.loads(saved.read_text());provenance=json.loads((prepared/'provenance.json').read_text())
 baseline=decode_saved_layout(fixture,audit['panels'][3]['state']);name,mode,width=CONFIGS[a.index];out=a.run/name;out.mkdir(parents=True,exist_ok=True)
 result=optimize_ribbons(fixture,baseline,mode=mode,width_weighted=width,weight=1,seconds=600)
 result['input_fixture_sha256']=hashlib.sha256((prepared/'fixture.json').read_bytes()).hexdigest()
 result['baseline_audit_sha256']=hashlib.sha256(saved.read_bytes()).hexdigest()
 final=decode_saved_layout(fixture,result['optimized_state']);(out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
 # Retain the actual four-panel baseline evidence, not reconstructed display states.
 panels=[]
 for index,title in [(0,'Input'),(2,'Supplemented GENESPACE + flip assistance'),(3,'Saved midpoint-optimal SynTangle')]:
  state=decode_saved_layout(fixture,audit['panels'][index]['state']);m=ribbon_metrics(fixture,state)
  panels.append((state,title,f"C={m['crossings']}; twists={m['twists']}; width-twists={m['width_weighted_twists']:.2f}"))
 m=result['final_metrics'];panels.append((final,name,f"C={m['crossings']}; twists={m['twists']}; width-twists={m['width_weighted_twists']:.2f}"))
 provenance['workflow_note']='Opt-in ribbon objective experiment; numerical solver status is not an independent optimality certificate.'
 render(fixture,panels,out,provenance);require_audit(out)
 shutil.copy2(prepared/'fixture.json',out/'fixture.json')
 print(name,json.dumps(result['starting_metrics']),'->',json.dumps(m),flush=True);print('PDF',out/'comparison.pdf')
if __name__=='__main__':main()
