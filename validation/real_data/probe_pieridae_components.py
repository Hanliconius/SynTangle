#!/usr/bin/env python3
"""Audit and test an existing Pieridae view using its full, verified tables."""
import argparse
import csv
import hashlib
import importlib.util
import json
from collections import Counter
from pathlib import Path
import shutil
import sys
import time

from syntangle.component_space import optimization_component_map
from syntangle.component_presolve import presolve_components
from syntangle.genespace import import_layout
from syntangle.hybrid_experiment import run_hybrid
from syntangle.layout import score_crossings
from syntangle.saved_layout import decode_saved_layout, validate_saved_layout
from syntangle.orientation_space import orientation_basis
from syntangle.table_input import load_tables
from syntangle.workflow import write_layout


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def dump(path, value):
    temp = path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value, indent=2)+'\n')
    temp.replace(path)


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write(path, rows, columns):
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--seconds', type=float, default=600)
    parser.add_argument('--audit-only', action='store_true')
    parser.add_argument('--render', action='store_true')
    args = parser.parse_args()
    source, out = args.source.resolve(), args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if (out/'audit.json').exists():
        raise ValueError('Audit already exists; do not overwrite a scientific run')
    receipt = json.loads((source/'optimization/COMPLETE.json').read_text())
    hashes = {name: sha(source/'tables'/name) for name in ('chromosomes.tsv', 'links.tsv')}
    if any(receipt['tables_sha256'].get(name) != digest for name,digest in hashes.items()):
        raise ValueError('Source table hash differs from its completion receipt')
    if sha(source/'optimization/layout.tsv') != receipt['layout_sha256']:
        raise ValueError('Source layout hash differs from receipt')
    f, _ = load_tables(source/'tables/chromosomes.tsv', source/'tables/links.tsv')
    basis = orientation_basis(f, f.chromosome_refs)
    candidates = []
    for name in ('input_layout.tsv', 'layout.tsv'):
        state = import_layout(f, source/'optimization'/name)
        validate_saved_layout(f, state, basis)
        candidates.append((score_crossings(f,state).crossings, name, state))
    _, selected, start = min(candidates, key=lambda x:(x[0], x[1]))
    before = score_crossings(f,start)
    components, component_of = optimization_component_map(f)
    begun = time.perf_counter()
    state, records = presolve_components(f,start,components,component_of)
    normalized = score_crossings(f,state)
    if normalized.crossings > before.crossings or normalized.link_count != before.link_count:
        raise ValueError('Presolve regression or evidence loss')
    labels = {(r['species'],r['chromosome']):r.get('label',r['chromosome']) for r in read(source/'tables/chromosomes.tsv')}
    rows = [dict(species=r.species_id, chromosome=r.chromosome_id,
                 label=labels[r.species_id,r.chromosome_id], component=component_of[r]) for r in f.chromosome_refs]
    write(out/'chromosome_components.tsv', rows, ['species','chromosome','label','component'])
    support = Counter()
    interface_links = Counter()
    cells = {}
    for r in read(source/'tables/links.tsv'):
        key = tuple(r[k] for k in ('species1','chromosome1','species2','chromosome2'))
        support[key] += 1
        interface = r['species1'],r['species2']
        interface_links[interface] += 1
        cells.setdefault(interface,set()).add((r['chromosome1'],r['chromosome2']))
    write(out/'chromosome_pair_support.tsv', [dict(species1=a,chromosome1=b,species2=c,chromosome2=d,links=n)
        for (a,b,c,d),n in sorted(support.items())], ['species1','chromosome1','species2','chromosome2','links'])
    z_components = sorted({component_of[r] for r in f.chromosome_refs
                           if labels[r.species_id,r.chromosome_id].lower() in ('z','chrz','chr_z')})
    audit = dict(source=str(source), input_sha256=hashes, selected_start=selected,
        candidate_scores=[dict(file=n,crossings=v) for v,n,_ in candidates],
        input=before.to_dict(), presolved=normalized.to_dict(), components=records,
        component_count=len(components), frozen_components=sum(r['frozen'] for r in records),
        explicit_z_label_components=z_components,
        z_note='Includes only explicit Z labels; chromosome 1 is never assumed to be Z. Inspect all members and partner counts.',
        interfaces=[dict(species1=a,species2=b,links=n,old_link_pairs=n*(n-1)//2,
                         chromosome_pair_cells=len(cells[a,b])) for (a,b),n in interface_links.items()],
        presolve_seconds=time.perf_counter()-begun,
        species_order=list(f.species_ids), species_order_changed=False,
        all_links_retained=True, independently_verified_certificate=False)
    dump(out/'audit.json', audit)
    write_layout(out/'presolved_layout.tsv', state)
    print('PRESOLVE',len(components),'components;',audit['frozen_components'],'frozen;',before.crossings,'->',normalized.crossings,flush=True)
    if args.audit_only:
        return
    dest = out/'optimization'; dest.mkdir()
    def checkpoint(candidate):
        write_layout(dest/'checkpoint.tmp.tsv',candidate)
        (dest/'checkpoint.tmp.tsv').replace(dest/'checkpoint_layout.tsv')
        value = score_crossings(f,candidate)
        dump(dest/'checkpoint.json',dict(crossings=value.crossings,layout_sha256=sha(dest/'checkpoint_layout.tsv'),status='feasible incumbent; search incomplete'))
    checkpoint(state)
    # Start from the original incumbent so all exact presolve bounds are retained
    # in the package's final result. This is a fresh experiment, not a resumed frontier.
    result = run_hybrid(f,start,'hybrid_adaptive',args.seconds,component_presolve=True,progress=checkpoint)
    final = decode_saved_layout(f,result['optimized_state'])
    validate_saved_layout(f,final,basis)
    after = score_crossings(f,final)
    if after.crossings > before.crossings or after.link_count != before.link_count:
        raise ValueError('Final full-evidence regression or evidence loss')
    write_layout(dest/'input_layout.tsv',start); write_layout(dest/'layout.tsv',final)
    dump(dest/'solver.json',result)
    dump(dest/'COMPLETE.json',dict(input=before.to_dict(),output=after.to_dict(),
        tables_sha256=hashes,layout_sha256=sha(dest/'layout.tsv'),
        status='bound closed under canonical scorer' if result['optimality_gap']==0 else 'bounded best known',
        lower_bound=result['lower_bound'],independent_certificate=False))
    if args.render:
        # Reuse the source view's frozen plotter/style and immutable homology.
        for name in ('tables','homology','envs'):
            (out/name).symlink_to(source/name,target_is_directory=True)
        for name in ('tasks.tsv','species_order.tsv','tree_source.json'):
            if (source/name).exists():shutil.copy2(source/name,out/name)
        config = json.loads((source/'config.json').read_text())
        config.update(comparison_input_label='Best supplied all-link layout',
                      comparison_output_label='Component-first SynTangle; all links')
        dump(out/'config.json',config)
        plotter = out/'code/plotter/analyse.py'
        sys.path.insert(0,str(plotter.parent))
        spec = importlib.util.spec_from_file_location('pieridae_plotter',plotter)
        module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
        module.plot(out)
    print('COMPONENT_TEST_COMPLETE',after.crossings,'L',result['lower_bound'],str(out),flush=True)


if __name__ == '__main__':
    main()
