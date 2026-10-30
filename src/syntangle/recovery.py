"""Full-evidence recovery proposals; no contraction or optimality claim."""
from collections import Counter, defaultdict
import time

from .component_space import optimization_component_map
from .crossing_cache import CrossingCostCache
from .layout import LayoutState, score_crossings
from .orientation_space import orientation_basis
from .saved_layout import validate_saved_layout


def connection_audit(rows, labels):
    """Label cell edges and find bridges; reciprocal best chains are proposals only."""
    support = Counter()
    for row in rows:
        a = (row['species1'], row['chromosome1'])
        b = (row['species2'], row['chromosome2'])
        support[tuple(sorted((a, b)))] += 1
    adjacency = defaultdict(set)
    neighbors = defaultdict(dict)
    for (a, b), count in support.items():
        adjacency[a].add(b)
        adjacency[b].add(a)
        neighbors[a, b[0]][b] = count
        neighbors[b, a[0]][a] = count
    best = {}
    for key, values in neighbors.items():
        maximum = max(values.values())
        winners = [r for r, value in values.items() if value == maximum]
        if len(winners) == 1:
            best[key] = winners[0]
    chains = defaultdict(set)
    edges = []
    for (a, b), count in sorted(support.items()):
        reciprocal = best.get((a, b[0])) == b and best.get((b, a[0])) == a
        if reciprocal:
            chains[a].add(b)
            chains[b].add(a)
        # Connectivity audit only: deleting a cell here never edits the fixture.
        seen = {a}
        stack = [a]
        while stack and b not in seen:
            u = stack.pop()
            for v in adjacency[u]:
                if (u == a and v == b) or (u == b and v == a):
                    continue
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
        edges.append(dict(species1=a[0], chromosome1=a[1], label1=labels.get(a,a[1]),
                          species2=b[0], chromosome2=b[1], label2=labels.get(b,b[1]),
                          links=count, fraction1=count/sum(neighbors[a,b[0]].values()),
                          fraction2=count/sum(neighbors[b,a[0]].values()),
                          reciprocal_best=reciprocal, graph_bridge=b not in seen,
                          explicit_z_incident=any(labels.get(x,x[1]).lower() in ('z','chrz','chr_z') for x in (a,b))))
    groups, visited = [], set()
    for root in sorted(chains):
        if root in visited:
            continue
        group, stack = set(), [root]
        while stack:
            u = stack.pop()
            if u in group:
                continue
            group.add(u)
            stack.extend(chains[u]-group)
        visited.update(group)
        # Never use a chain containing two chromosomes from one species.
        if len({r[0] for r in group}) == len(group):
            groups.append(tuple(sorted(group)))
    return edges, groups


def recover_layout(fixture, start, groups, seconds, progress=None, active_refs=None):
    """Monotone chain moves, whole-chromosome insertions and legal group flips.

    Time includes factor preparation. Each accepted move is independently
    rescored by the canonical scorer before checkpointing. A timed-out scan is
    not reported as a local optimum. No general MILP is constructed.
    """
    if seconds <= 0:
        raise ValueError('Recovery seconds must be positive')
    begun = time.perf_counter()
    deadline = begun + seconds
    basis = orientation_basis(fixture, fixture.chromosome_refs)
    validate_saved_layout(fixture, start, basis)
    cache = CrossingCostCache(fixture, cell_aggregation=True, preparation_deadline=deadline)
    state = start
    current = cache.prepare(state)
    if current != score_crossings(fixture, state).crossings:
        raise ValueError('Factor score differs from canonical score')
    refs = {(r.species_id,r.chromosome_id):r for r in fixture.chromosome_refs}
    active = set(fixture.chromosome_refs) if active_refs is None else set(active_refs)
    _, component_of = optimization_component_map(fixture)
    chains = [tuple(refs[x] for x in g) for g in groups if all(refs[x] in active for x in g)]
    history, evaluations = [], 0

    def candidates(state):
        # Coordinated end placement can cross a barrier to single-row moves.
        for group in chains:
            for end in (0, -1):
                order = dict(state.chromosome_order)
                for ref in group:
                    row = order[ref.species_id]
                    slots = [i for i,x in enumerate(row) if component_of[x] == component_of[ref]]
                    remaining = tuple(row[i] for i in slots if row[i] != ref)
                    moved = (ref,)+remaining if end == 0 else remaining+(ref,)
                    rebuilt = list(row)
                    for i,x in zip(slots,moved):
                        rebuilt[i] = x
                    order[ref.species_id] = tuple(rebuilt)
                yield 'chain_end', LayoutState(order, dict(state.chromosome_orientation))
            # Flip the entire candidate chain, including its hard-coupled partners.
            flip_refs = {r for free in basis.free_flip_groups if set(free) & set(group) for r in free}
            if flip_refs <= active:
                signs = dict(state.chromosome_orientation)
                for ref in flip_refs:
                    signs[ref] *= -1
                yield 'chain_flip', LayoutState(dict(state.chromosome_order), signs)
        # Propagate partner order in both directions, retaining every chromosome.
        # A neighboring row supplies a proposal, never a hard equality constraint.
        for left,right in zip(fixture.species_ids,fixture.species_ids[1:]):
            for target,neighbor in ((left,right),(right,left)):
                ranks = {r:i for i,r in enumerate(state.chromosome_order[neighbor])}
                partner = {}
                for chain in chains:
                    members = {r.species_id:r for r in chain}
                    if target in members and neighbor in members:
                        partner[members[target]] = members[neighbor]
                row = state.chromosome_order[target]
                for component in sorted({component_of[r] for r in partner}):
                    slots = [i for i,r in enumerate(row) if r in partner and component_of[r] == component]
                    ordered = sorted((row[i] for i in slots),key=lambda r:ranks[partner[r]])
                    rebuilt = list(row)
                    for i,r in zip(slots,ordered):
                        rebuilt[i] = r
                    if tuple(rebuilt) != row:
                        yield 'neighbor_order', LayoutState({**state.chromosome_order,target:tuple(rebuilt)},
                                                           dict(state.chromosome_orientation))
        for group in basis.free_flip_groups:
            if not set(group) <= active:
                continue
            signs = dict(state.chromosome_orientation)
            for ref in group:
                signs[ref] *= -1
            yield 'legal_flip', LayoutState(dict(state.chromosome_order), signs)
        # Insertion crosses several adjacent-swap barriers in one legal move.
        for species, row in state.chromosome_order.items():
            for ref in row:
                if ref not in active:
                    continue
                remaining = tuple(x for x in row if x != ref)
                for slot in range(len(row)):
                    new = remaining[:slot]+(ref,)+remaining[slot:]
                    if new != row and all(component_of[a] == component_of[b] for a,b in zip(new,row)):
                        yield 'insertion', LayoutState({**state.chromosome_order, species:new}, dict(state.chromosome_orientation))

    exhausted = False
    while time.perf_counter() < deadline:
        improved = False
        for kind, candidate in candidates(state):
            if time.perf_counter() >= deadline:
                break
            evaluations += 1
            value = cache.score_candidate(candidate)
            if value < current:
                validate_saved_layout(fixture,candidate,basis)
                actual = score_crossings(fixture,candidate)
                if actual.crossings != value:
                    raise ValueError('Recovery delta disagrees with canonical score')
                history.append(dict(move=kind,before=current,after=value,seconds=time.perf_counter()-begun,
                    changed_orders={sp:[r.label for r in order] for sp,order in candidate.chromosome_order.items()
                                    if order != state.chromosome_order[sp]},
                    changed_orientations={r.label:s for r,s in candidate.chromosome_orientation.items()
                                          if s != state.chromosome_orientation[r]}))
                state, current = candidate, value
                cache.prepare(state)
                if progress:
                    progress(state)
                improved = True
                break  # Regenerate all moves against the new incumbent.
        if not improved:
            exhausted = time.perf_counter() < deadline
            break
    return state, dict(status='neighborhood exhausted' if exhausted else 'time budget reached',
                       candidate_evaluations=evaluations, accepted_moves=history,
                       factor_lower_bound=cache.constant+sum(min(c) for a,b,c in cache.terms),
                       seconds=time.perf_counter()-begun, independent_certificate=False,
                       scope='full-evidence heuristic; no search-space pruning or new bound')
