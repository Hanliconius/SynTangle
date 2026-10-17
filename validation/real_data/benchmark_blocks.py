"""Same-evidence block-midpoint comparison; never a reproduced paper figure."""
import argparse,json,os,subprocess,sys,time
from pathlib import Path
from syntangle import load_fixture,score_crossings
from syntangle.layout import initial_layout_state
from syntangle.saved_layout import decode_saved_layout
from syntangle.refinement import select_refinement
from syntangle.hybrid_experiment import run_hybrid
from syntangle.audit import fixture_fingerprint
from syntangle.structural import build_structural_projection
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'benchmark'))
from compare_genespace import export_native_input,load_native_states,improve_flips


def write(path,data):
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(data,indent=2)+'\n');tmp.replace(path)


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--manifest',type=Path,required=True);p.add_argument('--index',type=int,required=True);p.add_argument('--seconds',type=int,default=150);a=p.parse_args()
    d=json.loads(a.manifest.read_text())['datasets'][a.index];f=load_fixture(a.root/d['id']/'fixture.json')
    out=a.root/d['id']/f'benchmark_{a.seconds}';out.mkdir(parents=True,exist_ok=True)
    fingerprint=fixture_fingerprint(f);initial=initial_layout_state(f)
    result=dict(dataset=d['id'],fingerprint=fingerprint,budget=a.seconds,scope='GENESPACE 1.3.1 native chromosome ordering applied to the SAME block-link projection; NOT the original published riparian layout or a complete native gene-level GENESPACE pipeline',input_crossings=score_crossings(f,initial).crossings,graph=build_structural_projection(f).to_dict())
    export_native_input(f,out/'native_input')
    helper=Path(__file__).resolve().parents[1]/'benchmark/genespace_native_order.R'
    subprocess.run(['bash','-c','set -euo pipefail; eval "$("$1" shell hook --shell bash)"; micromamba activate "$2"; exec Rscript "$3" "$4" "$5"','native',os.environ['SCT_MAMBA'],os.environ.get('SCT_GENESPACE_ENV','lep_busco_painter_clean'),str(helper),str(out/'native_input'),str(out)],check=True,timeout=180)
    variants=list(load_native_states(f,out/'native_orders.tsv'))
    native=min(variants,key=lambda v:(score_crossings(f,v[2]).crossings,v[0]));assisted=[];variant_rows=[]
    started=time.perf_counter()
    for j,(name,meta,state) in enumerate(variants):
        flipped,n=improve_flips(f,state,100+j)
        assisted.append((name,flipped));variant_rows.append(dict(variant=name,native_crossings=score_crossings(f,state).crossings,assisted_crossings=score_crossings(f,flipped).crossings,candidates=n))
    best=min(assisted,key=lambda v:(score_crossings(f,v[1]).crossings,v[0]))
    result.update(native_order_crossings=score_crossings(f,native[2]).crossings,native_variant=native[0],assisted_crossings=score_crossings(f,best[1]).crossings,assisted_variant=best[0],flip_seconds=time.perf_counter()-started,variants=variant_rows)
    write(out/'GENESPACE.layout.json',native[2].to_dict());write(out/'GENESPACE_plus_flips.layout.json',best[1].to_dict());write(out/'baseline.json',result)
    checkpoints=[]
    def progress(state):
        c=score_crossings(f,state).crossings;checkpoints.append(c)
        write(out/'incumbent.json',dict(crossings=c,state=state.to_dict(),proof_pending=True))
    # Independent public display start. Baseline retention is reported separately.
    raw=run_hybrid(f,initial,'hybrid_adaptive',a.seconds,progress=progress,mirror_symmetry=True,neighborhood_min_decisions=512,neighborhood_fraction=.15,neighborhood_max_seconds=30,deduplicate_neighborhoods=True)
    if fixture_fingerprint(f)!=fingerprint:raise AssertionError('Evidence changed')
    decision=select_refinement(f,best[1],decode_saved_layout(f,raw['optimized_state']),lower_bound=raw['lower_bound'],reference_name='native GENESPACE ordering + our flips, block projection')
    result.update(raw_solver=raw,validated_selection=decision,checkpoint_crossings=checkpoints,strict_raw_improvement_over_native=raw['upper_bound']<result['native_order_crossings'],strict_raw_improvement_over_assisted=raw['upper_bound']<result['assisted_crossings'])
    write(out/'result.json',result)
    print(f"{d['id']} budget={a.seconds} GS_order={result['native_order_crossings']} GS_plus_flips={result['assisted_crossings']} raw_ST={raw['upper_bound']} lower={raw['lower_bound']} returned={decision['returned_crossings']} {raw['optimality_status']}",flush=True)

if __name__=='__main__':main()
