"""Exact easy-component presolve on immutable, full evidence.

The two-state chain recurrence is exact for the existing midpoint objective.
It is not an independent rational-arithmetic certificate or a homology filter.
"""
from collections import Counter

from .layout import LayoutState, canonicalize_component_order, score_crossings
from .layout import _occurrences_by_species_homology
from .pair_cost import pair_component_crossings


def presolve_components(fixture, state, components, component_of, progress=None):
    """Return a canonical incumbent and records for every decision component.

    A strand has at most one chromosome per displayed species. Its only free
    decisions are signs, with interactions on adjacent rows. Nonlocal hard
    equations are conservatively deferred; none are dropped or approximated.
    """
    state = canonicalize_component_order(fixture, state)
    index = _occurrences_by_species_homology(fixture)
    species = fixture.species_ids
    positions = {sp: i for i, sp in enumerate(species)}
    records = []
    for component_id, nodes in enumerate(components):
        refs = tuple(r for r in fixture.chromosome_refs if component_of[r] == component_id)
        counts = Counter(r.species_id for r in refs)
        record = dict(component=component_id, chromosomes=[r.label for r in refs],
                      max_chromosomes_per_row=max(counts.values(), default=0),
                      component_lower=0, status='unresolved', frozen=False)
        strand = all(n <= 1 for n in counts.values())
        constraints = [c for c in fixture.orientation_constraints if c.a in refs]
        nonlocal_constraint = any(abs(positions[c.a.species_id] - positions[c.b.species_id]) > 1
                                  for c in constraints)
        if strand and not nonlocal_constraint:
            by_species = {r.species_id: r for r in refs}
            domains = [((1, -1) if sp in by_species else (0,)) for sp in species]
            # Tuples (cost, flips from incumbent, sign path) break ties stably.
            first = by_species.get(species[0])
            dp = {s: (0, int(first is not None and s != state.chromosome_orientation[first]), (s,))
                  for s in domains[0]}
            tables = []
            for i in range(1, len(species)):
                a, b = species[i-1:i+1]
                left, right = by_species.get(a), by_species.get(b)
                table = {}
                following = {}
                for t in domains[i]:
                    candidates = []
                    for s, (cost, flips, path) in dp.items():
                        signs = {r: v for r, v in ((left, s), (right, t)) if r is not None}
                        if any(c.a in signs and c.b in signs and
                               ((signs[c.a] != signs[c.b]) != bool(c.xor)) for c in constraints):
                            continue
                        value = pair_component_crossings(
                            fixture, a, b, (left,) if left else (), (right,) if right else (),
                            signs, nodes, occurrence_index=index)
                        table[f'{s},{t}'] = value
                        candidates.append((cost+value,
                            flips+int(right is not None and t != state.chromosome_orientation[right]), path+(t,)))
                    if candidates:
                        following[t] = min(candidates)
                if not following:
                    raise ValueError('No legal orientation for strand component')
                dp = following
                tables.append(dict(species1=a, species2=b, costs=table))
            optimum, _, path = min(dp.values())
            signs = dict(state.chromosome_orientation)
            signs.update({by_species[sp]: s for sp, s in zip(species, path) if sp in by_species})
            state = LayoutState(dict(state.chromosome_order), signs)
            checked = score_crossings(fixture, state, restrict_component_nodes=nodes).crossings
            if checked != optimum:
                raise ValueError('Strand DP reconstruction disagrees with canonical score')
            record.update(status='exact strand dynamic program', frozen=True,
                          component_lower=optimum, component_upper=optimum,
                          transition_tables=tables, selected_signs=list(path),
                          proof_scope='exhaustive two-state recurrence; canonical midpoint arithmetic')
        else:
            value = score_crossings(fixture, state, restrict_component_nodes=nodes).crossings
            record['component_upper'] = value
            if value == 0:
                record.update(status='zero bound attained', frozen=True,
                              proof_scope='nonnegative crossing objective and feasible zero layout')
            elif nonlocal_constraint:
                record['deferred_reason'] = 'nonadjacent hard orientation equation'
        records.append(record)
        if progress is not None and record['frozen']:
            progress(state)
    return state, records
