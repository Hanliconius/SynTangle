"""User workflows with explicit statistical and preservation reporting."""
import csv
import hashlib
import json
import math
import time
from collections import Counter
from importlib.metadata import version, PackageNotFoundError
from pathlib import Path
from .fixtures import fixture_from_dict
from .genespace import flip_assist, import_layout
from .layout import initial_layout_state, score_crossings
from .orientation_space import orientation_basis
from .saved_layout import decode_saved_layout, validate_saved_layout
from .table_input import load_tables
from .visualize import _unambiguous_links


def metrics(fixture, state):
    result = score_crossings(fixture, state).to_dict()
    terms = []
    for a, b in zip(fixture.species_ids, fixture.species_ids[1:]):
        for _, (c1,b1), (c2,b2) in _unambiguous_links(fixture,a,b):
            twist = ((b1.strand == '-') ^ (b2.strand == '-') ^
                     (state.chromosome_orientation[c1.ref] < 0) ^
                     (state.chromosome_orientation[c2.ref] < 0))
            terms.append((twist, min(b1.end-b1.start,b2.end-b2.start)))
    mean = sum(w for _,w in terms)/len(terms) if terms else 1
    result.update(twists=sum(t for t,_ in terms),
                  width_weighted_twists=sum(w/mean for t,w in terms if t))
    return result


def write_layout(path, state):
    with Path(path).open('w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['species','chromosome','rank','orientation'])
        for sp, order in state.chromosome_order.items():
            for rank, ref in enumerate(order,1):
                writer.writerow([sp,ref.chromosome_id,rank,state.chromosome_orientation[ref]])


def run(args):
    for name in ('seconds','ribbon_seconds'):
        value = getattr(args,name,60)
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f'--{name.replace("_","-")} must be finite and positive')
    if not math.isfinite(getattr(args,'weight',1)) or getattr(args,'weight',1) < 0:
        raise ValueError('--weight must be finite and nonnegative')
    out = Path(args.output)
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError('Output must be a new or empty directory')
    if args.fixture:
        if args.chromosomes or args.links:
            raise ValueError('Use --fixture OR --chromosomes and --links')
        data = json.loads(Path(args.fixture).read_text())
        fixture = fixture_from_dict(data)
        inputs = [args.fixture]
    else:
        if not args.chromosomes or not args.links:
            raise ValueError('Provide --chromosomes and --links, or --fixture')
        fixture, data = load_tables(args.chromosomes,args.links)
        inputs = [args.chromosomes,args.links]
    for c in fixture.chromosomes:
        if not math.isfinite(c.length) or any(not math.isfinite(x) for b in c.blocks for x in (b.start,b.end)):
            raise ValueError('All coordinates and lengths must be finite')
        if any(b.strand == '?' for b in c.blocks):
            raise ValueError('Directed ribbon workflows require known + or - block strands')
    for sp in fixture.species_ids:
        counts = Counter(b.homology_id for c in fixture.chromosomes if c.ref.species_id==sp for b in c.blocks)
        if any(n>1 for n in counts.values()):
            raise ValueError('Repeated homology within a species is ambiguous; supply explicit pairwise links')
    start = import_layout(fixture,args.layout) if args.layout else initial_layout_state(fixture)
    if args.layout: inputs.append(args.layout)
    basis = orientation_basis(fixture,fixture.chromosome_refs)
    validate_saved_layout(fixture,start,basis)
    before = metrics(fixture,start)
    if not before['link_count']: raise ValueError('No links between adjacent species')
    if args.figure:
        from .public_plot import plot_comparison
    begin = time.monotonic()
    solver = None
    ribbon = None
    lower = None
    if args.command == 'flip-assist':
        final = flip_assist(fixture,start,args.seconds)
        assert final.chromosome_order == start.chromosome_order
    else:
        from .hybrid_experiment import run_hybrid
        solver = run_hybrid(fixture,start,'hybrid_adaptive',args.seconds,seed=args.seed)
        final = decode_saved_layout(fixture,solver['optimized_state'])
        lower = solver['lower_bound']
        if args.ribbon_objective != 'none':
            from .ribbon_objective import optimize_ribbons
            ribbon = optimize_ribbons(fixture,final,
                mode='lex' if args.ribbon_objective=='crossing-cap' else 'weighted',
                weight=args.weight,width_weighted=args.twist_metric=='width',seconds=args.ribbon_seconds)
            final = decode_saved_layout(fixture,ribbon['optimized_state'])
    validate_saved_layout(fixture,final,basis)
    after = metrics(fixture,final)
    assert before['link_count'] == after['link_count']
    status = 'heuristic_no_certificate' if lower is None else ('solver_bound_closed' if lower==after['crossings'] else 'bounded_best_known')
    try: software_version = version('syntangle')
    except PackageNotFoundError: software_version = 'source-checkout'
    stats = dict(schema_version=1,software_version=software_version,mode=args.command,
        crossing_objective='Strict inversions of link midpoints between adjacent species; ties excluded',
        ribbon_objective=getattr(args,'ribbon_objective','none'),input=before,output=after,
        crossing_lower_bound=lower,crossing_gap=None if lower is None else after['crossings']-lower,
        crossing_status=status,independent_certificate=False,seconds=time.monotonic()-begin,
        species_count=len(fixture.species_ids),chromosome_count=len(fixture.chromosomes),
        input_sha256={str(p):hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in inputs},
        allowed_moves='whole-chromosome flips' if args.command=='flip-assist' else 'whole-chromosome permutations and flips',
        preservation=dict(all_chromosomes=True,all_links=True,internal_coordinates=True),
        state=final.to_dict())
    out.mkdir(parents=True,exist_ok=True)
    write_layout(out/'layout.tsv',final)
    write_layout(out/'input_layout.tsv',start)
    for name, value in [('fixture.json',data),('stats.json',stats),('solver.json',solver),('ribbon_solver.json',ribbon)]:
        if value is not None: (out/name).write_text(json.dumps(value,indent=2)+'\n')
    with (out/'flips.tsv').open('w') as handle:
        writer=csv.writer(handle,delimiter='\t');writer.writerow(['genome','chr'])
        for ref in fixture.chromosome_refs:
            if final.chromosome_orientation[ref] != start.chromosome_orientation[ref]:
                writer.writerow([ref.species_id,ref.chromosome_id])
    with (out/'pair_stats.tsv').open('w') as handle:
        writer=csv.writer(handle,delimiter='\t');writer.writerow(['species1','species2','input_crossings','output_crossings'])
        for a,b in zip(before['pair_crossings'],after['pair_crossings']):
            writer.writerow([a['species1'],a['species2'],a['crossings'],b['crossings']])
    with (out/'summary.tsv').open('w') as handle:
        writer=csv.writer(handle,delimiter='\t')
        writer.writerow(['input_crossings','output_crossings','links','chromosomes','lower_bound','gap','status'])
        writer.writerow([before['crossings'],after['crossings'],after['link_count'],len(fixture.chromosomes),lower,stats['crossing_gap'],status])
    if args.figure: plot_comparison(fixture,start,final,out/(args.command.replace('-','_')+'_comparison.pdf'))
    print(f'{args.command}: {before["crossings"]} -> {after["crossings"]} crossings; {status}; {out.resolve()}')
    return 0
