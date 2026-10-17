"""Size-filtered search guidance; full evidence is always restored and rescored."""
import argparse
from dataclasses import replace
import hashlib
import json
import math
from pathlib import Path
import time
from syntangle.heuristic import _component_map
from syntangle.fixtures import load_fixture
from syntangle.hybrid_experiment import run_hybrid
from syntangle.layout import score_crossings, _anchor_key
from syntangle.saved_layout import decode_saved_layout
from syntangle.visualize import _unambiguous_links

FRACTIONS = (0, .10, .25, .50, .75)

def atomic(path, data):
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(data, indent=2) + '\n')
    tmp.replace(path)

def links(fixture):
    return [link for a,b in zip(fixture.species_ids, fixture.species_ids[1:])
            for link in _unambiguous_links(fixture,a,b)]

def sizes(fixture):
    spans = {}
    for hid,(a,x),(b,y) in links(fixture):
        # Require substantial span on both assemblies; raw bp, not BLAST scores.
        spans.setdefault(hid, []).extend((x.end-x.start,y.end-y.start))
    return {h:math.exp(sum(math.log(max(v,1e-300)) for v in values)/len(values))
            for h,values in spans.items()}

def filtered(fixture, importance, fraction):
    ordered = sorted(importance, key=lambda h:(importance[h],h))
    removed = set(ordered[:int(len(ordered)*fraction)])
    # Retain ALL chromosomes and hard equations, including now-unlinked chromosomes.
    return replace(fixture, chromosomes=tuple(replace(c, blocks=tuple(
        b for b in c.blocks if b.homology_id not in removed)) for c in fixture.chromosomes)), removed

def metrics(fixture, state):
    ranks = {r:i for order in state.chromosome_order.values() for i,r in enumerate(order)}
    count=0; weighted=0.0
    for a,b in zip(fixture.species_ids,fixture.species_ids[1:]):
        entries=[]
        for hid,(c,x),(d,y) in _unambiguous_links(fixture,a,b):
            left=_anchor_key(c,x,ranks[c.ref],state.chromosome_orientation[c.ref])
            right=_anchor_key(d,y,ranks[d.ref],state.chromosome_orientation[d.ref])
            w=math.sqrt((x.end-x.start)/c.length * (y.end-y.start)/d.length)
            entries.append((left,right,w))
        for i,(l,r,w) in enumerate(entries):
            for ll,rr,ww in entries[i+1:]:
                if (l<ll and r>rr) or (l>ll and r<rr):
                    count+=1; weighted+=math.sqrt(w*ww)
    assert count == score_crossings(fixture,state).crossings
    return {'crossings':count,'span_weighted_crossings':weighted}

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--index',type=int,required=True)
    p.add_argument('--case',default='annotated_lep_reanalysis')
    p.add_argument('--start-panel',type=int,default=None)
    args=p.parse_args(); fraction=FRACTIONS[args.index]
    prepared=args.source/'prepared'/args.case
    figure=args.source/'figures'/args.case
    fixture=load_fixture(prepared/'fixture.json')
    audit=json.loads((figure/'comparison_audit.json').read_text())
    sha=hashlib.sha256((prepared/'fixture.json').read_bytes()).hexdigest()
    assert audit['fixture_sha256']==sha
    baseline=(min(audit['panels'],key=lambda x:x['crossings']) if args.start_panel is None
              else audit['panels'][args.start_panel])
    start=decode_saved_layout(fixture,baseline['state'])
    assert score_crossings(fixture,start).crossings==baseline['crossings']
    destination=args.output/f'task_{args.index}';destination.mkdir(parents=True,exist_ok=False)
    importance=sizes(fixture);best=start;visual_best=start;started=time.perf_counter();events=[]
    importance_note='geometric mean link endpoint spans; globally ranked, lexical ties'
    provenance=json.loads((prepared/'provenance.json').read_text())
    if 'original_coordinates' in provenance:
        importance={h:math.exp(sum(math.log(r['end']-r['start']) for r in rows)/len(rows))
                    for h,rows in provenance['original_coordinates'].items()}
        if set(importance)!=set(sizes(fixture)):raise ValueError('Source gene spans differ from retained evidence')
        importance_note='geometric mean deposited BUSCO gene span in bp; remove whole shared BUSCO across rows; lexical ties'
    def checkpoint(state):
        nonlocal best,visual_best
        measured=metrics(fixture,state)
        if measured['crossings']<metrics(fixture,best)['crossings']:
            best=state;events.append(dict(seconds=time.perf_counter()-started,**measured))
            atomic(destination/'incumbent.json',dict(state=state.to_dict(),**measured))
        if measured['span_weighted_crossings']<metrics(fixture,visual_best)['span_weighted_crossings']:
            visual_best=state
    checkpoint(start);state=start;stages=[]
    # Equal total requested solver budget: 180s for every task.
    schedule=[(0,180)] if fraction==0 else [(fraction,60),(fraction/2,60),(0,60)]
    full_lower=0
    for removed_fraction,seconds in schedule:
        reduced,removed=filtered(fixture,importance,removed_fraction)
        result=run_hybrid(reduced,state,'hybrid_adaptive',seconds,progress=checkpoint,
                          mirror_symmetry=True,deduplicate_neighborhoods=True)
        state=decode_saved_layout(fixture,result['optimized_state']);checkpoint(state)
        if removed_fraction==0:full_lower=max(full_lower,result['lower_bound'])
        components,_=_component_map(reduced)
        stages.append(dict(component_sizes=sorted([len(c) for c in components],reverse=True),removed_fraction=removed_fraction,removed_ids=sorted(removed),
                           retained_links=len(links(reduced)),full_data_metrics=metrics(fixture,state),
                           result=result,bound_scope='full evidence' if not removed else 'filtered evidence only'))
        print(f'filter={fraction} stage_removed={removed_fraction} full_C={metrics(fixture,state)["crossings"]}',flush=True)
    final=metrics(fixture,best)
    assert final['crossings']<=baseline['crossings']
    atomic(destination/'report.json',dict(fixture_sha256=sha,filter_fraction=fraction,
        initial=metrics(fixture,start),final=final,best_state=best.to_dict(),
        visual_best=metrics(fixture,visual_best),visual_best_state=visual_best.to_dict(),
        full_lower_bound=full_lower,full_gap=final['crossings']-full_lower,
        proof_scope='full-data bound only; filtered optima are not full-data certificates',
        importance=importance_note,
        start_panel=args.start_panel,case=args.case,
        visual_metric='crossing pair weight sqrt(w_i*w_j); w=geometric mean fractional chromosome span',
        events=events,stages=stages,wall_seconds=time.perf_counter()-started))
    (destination/'COMPLETE').write_text('PASS\n')
    print(f'FINAL filter={fraction} C={final["crossings"]} L={full_lower}',flush=True)

if __name__=='__main__':main()
