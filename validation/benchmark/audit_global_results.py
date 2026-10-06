"""Re-score previous reported MILP optima and test strict improvement infeasible.

This is a separate verification problem, never a restart of optimization's
excluded frontier. Numerical HiGHS infeasibility remains a solver certificate,
not an independently checked rational proof.
"""
import argparse
import csv
import json
from pathlib import Path
from syntangle import load_validation_bundle, score_crossings
from syntangle.heuristic import _component_map
from syntangle.saved_layout import decode_saved_layout,validate_saved_layout
from syntangle.orientation_space import orientation_basis
from syntangle.global_experiment import JointModel
from benchmark_pipeline import atomic


def audit(root,history,output,seconds=30,case_index=None):
    records=[]
    with (root/'benchmark_manifest.tsv').open() as handle:
        entries=list(csv.DictReader(handle,delimiter='\t'))
    if case_index is not None:entries=[entries[case_index]]
    for entry in entries:
        case_id=entry['case_id'];fixture=load_validation_bundle(root/entry['case_dir'])
        candidates=[]
        runs=sorted(list(history.glob('global_methods_*'))+list(history.glob('hybrid_methods_*')))
        for run in runs:
            for path in (run/case_id).glob('*/result.json'):
                result=json.loads(path.read_text())
                state=decode_saved_layout(fixture,result['optimized_state'])
                validate_saved_layout(fixture,state,orientation_basis(fixture,fixture.chromosome_refs))
                score=score_crossings(fixture,state).crossings
                if score!=result['upper_bound']:raise AssertionError('Recorded crossing score changed: '+str(path))
                if result['lower_bound']==score:candidates.append((score,str(path),state))
        if not candidates:
            records.append(dict(case_id=case_id,status='no prior reported optimum to audit'))
            print('AUDIT',case_id,'no prior optimum',flush=True);continue
        score,source,state=min(candidates,key=lambda x:(x[0],x[1]))
        components,component_of=_component_map(fixture);checks=[];component_sum=0
        for i,nodes in enumerate(components):
            value=score_crossings(fixture,state,restrict_component_nodes=nodes).crossings
            component_sum+=value
            if value==0:
                checks.append(dict(component=i,no_better_proven=True,reason='zero crossings; nonnegative objective'))
                continue
            model=JointModel(fixture,nodes,tuple(r for r in fixture.chromosome_refs if component_of[r]==i))
            if abs(model.objective(model.encode(state))-value)>1e-6:raise AssertionError('Model/canonical score mismatch')
            check=model.highs(state,seconds,strict_proof=True);checks.append(dict(component=i,**check))
            if check['counterexample_crossings'] is not None:
                raise AssertionError('Counterexample to prior optimum: '+source)
        if component_sum!=score:raise AssertionError('Component sum differs from canonical whole-fixture objective')
        verified=all(r['no_better_proven'] for r in checks)
        records.append(dict(case_id=case_id,canonical_score=score,source=source,
            status='strict-improvement infeasibility verified' if verified else 're-score passed; proof unresolved within audit budget',
            components=checks))
        atomic(output,records)
        print('AUDIT',case_id,'C='+str(score),'STRICT_NO_BETTER='+str(verified),flush=True)
    atomic(output,records)
    return records

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--root',required=True);p.add_argument('--history',required=True);p.add_argument('--output',required=True)
    p.add_argument('--seconds',type=float,default=30);p.add_argument('--case-index',type=int)
    args=p.parse_args();audit(Path(args.root),Path(args.history),Path(args.output),args.seconds,args.case_index)
