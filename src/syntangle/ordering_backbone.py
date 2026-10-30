"""Experimental ordering-only evidence mask, never an evidence deletion policy."""
from collections import defaultdict
from itertools import groupby


def select_links(rows, min_support=3, min_run=3):
    """Require distinct-anchor chromosome-pair support and runs on both sides.

    Runs are consecutive anchors with one chromosome partner in native coordinate
    order. They need not be collinear. Ambiguous partner assignments break runs.
    All copies of an accepted anchor retain their original link identities.
    """
    if min_support < 1 or min_run < 1:
        raise ValueError('Positive support and run thresholds required')
    if len({r['link_id'] for r in rows}) != len(rows):
        raise ValueError('Duplicate link IDs')
    support = defaultdict(set)
    for r in rows:
        support[tuple(r[k] for k in ('species1','species2','chromosome1','chromosome2'))].add(r['buscoID'])
    accepted = []
    for side, other in ((1,2),(2,1)):
        tracks = defaultdict(lambda: defaultdict(list))
        for r in rows:
            key = (r['species1'],r['species2'],r['chromosome'+str(side)])
            point = (float(r['start'+str(side)])+float(r['end'+str(side)]))/2
            tracks[key][(point,r['buscoID'])].append(r)
        keep = set()
        for points in tracks.values():
            sequence=[]
            for point, links in sorted(points.items()):
                partners={r['chromosome'+str(other)] for r in links}
                partner=next(iter(partners)) if len(partners)==1 else None
                sequence.append((partner,point[1],links))
            for partner, run in groupby(sequence,key=lambda x:x[0]):
                run=list(run)
                if partner is not None and len({x[1] for x in run})>=min_run:
                    keep.update(r['link_id'] for x in run for r in x[2])
        accepted.append(keep)
    kept=[];audit=[]
    for r in rows:
        key=tuple(r[k] for k in ('species1','species2','chromosome1','chromosome2'))
        n=len(support[key]);runs=all(r['link_id'] in mask for mask in accepted)
        yes=n>=min_support and runs
        if yes:kept.append(dict(r))
        audit.append(dict(link_id=r['link_id'],buscoID=r['buscoID'],distinct_pair_support=n,
                          supported_runs_on_both_sides=runs,used_for_ordering=yes,
                          retained_for_final_plot=True))
    return kept,audit
