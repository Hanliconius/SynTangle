import argparse,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--manifest',type=Path,required=True);a=p.parse_args()
lines=['# Real block-link comparisons','','GENESPACE ordering-function comparison on matched block projections; not reproduced published layouts.','','| Dataset | Budget | GS order | GS + our flips | Raw SynTangle | Lower | Returned | Status |','|---|---:|---:|---:|---:|---:|---:|---|'];records=[]
for d in json.loads(a.manifest.read_text())['datasets']:
 for seconds in (150,600):
  path=a.root/d['id']/f'benchmark_{seconds}'/'result.json'
  if not path.exists():
   lines.append(f"| {d['id']} | {seconds} | — | — | — | — | — | incomplete; inspect Slurm logs |");continue
  r=json.loads(path.read_text());s=r['raw_solver'];v=r['validated_selection'];records.append(r)
  status=s['optimality_status']+'; '+('raw improvement' if r['strict_raw_improvement_over_assisted'] else 'raw tie' if s['upper_bound']==r['assisted_crossings'] else 'raw regression; baseline retained')
  lines.append(f"| {d['id']} | {seconds} | {r['native_order_crossings']} | {r['assisted_crossings']} | {s['upper_bound']} | {s['lower_bound']} | {v['returned_crossings']} | {status} |")
lines+=['','Only saved complete results are compared. Guarded baseline retention is not counted as a raw solver win.','Six archived cases come from one study; the two grass panels overlap in biology.']
(a.root/'comparison_summary.md').write_text('\n'.join(lines)+'\n');(a.root/'comparison_records.json').write_text(json.dumps(records,indent=2)+'\n');print('\n'.join(lines))
