"""Refine a complete exported GENESPACE display against fixed supplied evidence."""
import csv
import hashlib
import json
import math
import time
from pathlib import Path
from .fixtures import load_fixture
from .layout import LayoutState, score_crossings
from .orientation_space import orientation_basis
from .saved_layout import validate_saved_layout


def import_layout(fixture, path):
    """TSV columns: species, chromosome, rank, orientation (+1/-1)."""
    refs={(r.species_id,r.chromosome_id):r for r in fixture.chromosome_refs}
    rows={sp:[] for sp in fixture.species_ids};signs={}
    with Path(path).open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            ref=refs[row['species'],row['chromosome']]
            if ref in signs:raise ValueError('Duplicate chromosome in imported layout')
            rank=int(row['rank']);sign=int(row['orientation'])
            rows[ref.species_id].append((rank,ref));signs[ref]=sign
    for values in rows.values():
        if len({rank for rank,_ in values})!=len(values):raise ValueError('Duplicate chromosome rank')
    state=LayoutState({sp:tuple(r for _,r in sorted(values)) for sp,values in rows.items()},signs)
    validate_saved_layout(fixture,state,orientation_basis(fixture,fixture.chromosome_refs))
    return state


def flip_assist(fixture, start, seconds=60):
    """Strictly improving group flips; fixed order and hard orientation equations."""
    basis=orientation_basis(fixture,fixture.chromosome_refs)
    validate_saved_layout(fixture,start,basis)
    current=start;value=score_crossings(fixture,current).crossings;deadline=time.monotonic()+seconds
    while time.monotonic()<deadline:
        best=current;best_value=value
        for group in basis.free_flip_groups:
            if time.monotonic()>=deadline:break
            signs=dict(current.chromosome_orientation)
            for ref in group:signs[ref]*=-1
            candidate=LayoutState(current.chromosome_order,signs)
            cost=score_crossings(fixture,candidate).crossings
            if cost<best_value:best,best_value=candidate,cost
        if best_value>=value:break
        current,value=best,best_value
    validate_saved_layout(fixture,current,basis)
    if current.chromosome_order!=start.chromosome_order:raise AssertionError('Flip assistance changed order')
    return current


def refine(args):
    if not math.isfinite(args.seconds) or args.seconds<=0:raise ValueError('seconds must be finite and positive')
    fixture=load_fixture(args.fixture);start=import_layout(fixture,args.layout)
    evidence=json.dumps(fixture.to_dict(),sort_keys=True) if hasattr(fixture,'to_dict') else Path(args.fixture).read_text()
    before=hashlib.sha256(evidence.encode()).hexdigest()
    if args.mode=='flips':
        final=flip_assist(fixture,start,args.seconds);status='heuristic; no optimality certificate'
    else:
        from .hybrid_experiment import run_hybrid
        from .saved_layout import decode_saved_layout
        result=run_hybrid(fixture,start,'hybrid_adaptive',args.seconds)
        final=decode_saved_layout(fixture,result['optimized_state']);status=result['optimality_status']
    validate_saved_layout(fixture,final,orientation_basis(fixture,fixture.chromosome_refs))
    initial_score=score_crossings(fixture,start).crossings;final_score=score_crossings(fixture,final).crossings
    if final_score>initial_score:raise AssertionError('Refinement regressed')
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    with (out/'layout.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['species','chromosome','rank','orientation'])
        for sp,order in final.chromosome_order.items():
            for rank,ref in enumerate(order,1):w.writerow([sp,ref.chromosome_id,rank,final.chromosome_orientation[ref]])
    (out/'audit.json').write_text(json.dumps(dict(mode=args.mode,input_crossings=initial_score,output_crossings=final_score,optimality_status=status,fixture_sha256=before,layout_input_sha256=hashlib.sha256(Path(args.layout).read_bytes()).hexdigest(),allowed_moves='whole-chromosome flips' if args.mode=='flips' else 'whole-chromosome order and flips',state=final.to_dict()),indent=2)+'\n')
    print(f'{args.mode}: {initial_score} -> {final_score}; {out / "layout.tsv"}')
    return 0
