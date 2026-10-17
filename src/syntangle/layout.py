from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations, product
from math import factorial

from .incidence import build_incidence_graph, chromosome_node_id
from .model import BlockOccurrence, Chromosome, ChromosomeRef, Fixture


class AmbiguousHomologyError(ValueError):
    """Raised when exact crossing scoring would require resolving ambiguous copies."""


class InfeasibleOrientationConstraints(ValueError):
    """Raised when no whole-chromosome orientation satisfies hard XOR constraints."""


class SearchSpaceTooLarge(ValueError):
    """Raised when exact enumeration exceeds the configured component state cap."""


@dataclass(frozen=True)
class CrossingScore:
    crossings: int
    pair_crossings: tuple[tuple[str, str, int], ...]
    link_count: int

    def to_dict(self) -> dict[str, object]:
        return {
            "crossings": self.crossings,
            "pair_crossings": [
                {"species1": a, "species2": b, "crossings": count}
                for a, b, count in self.pair_crossings
            ],
            "link_count": self.link_count,
        }


@dataclass(frozen=True)
class LayoutState:
    chromosome_order: dict[str, tuple[ChromosomeRef, ...]]
    chromosome_orientation: dict[ChromosomeRef, int]

    def to_dict(self) -> dict[str, object]:
        return {
            "chromosome_order": {
                species: [ref.chromosome_id for ref in refs]
                for species, refs in self.chromosome_order.items()
            },
            "chromosome_orientation": {
                ref.label: orientation
                for ref, orientation in sorted(self.chromosome_orientation.items())
            },
        }


@dataclass(frozen=True)
class ExactLayoutResult:
    initial_state: LayoutState
    normalized_state: LayoutState
    optimized_state: LayoutState
    initial_score: CrossingScore
    normalized_score: CrossingScore
    optimized_score: CrossingScore
    states_evaluated: int
    optimality_status: str

    @property
    def excess_layout_crossings_initial(self) -> int:
        return self.initial_score.crossings - self.optimized_score.crossings

    def to_dict(self) -> dict[str, object]:
        return {
            "initial_state": self.initial_state.to_dict(),
            "normalized_state": self.normalized_state.to_dict(),
            "optimized_state": self.optimized_state.to_dict(),
            "initial_score": self.initial_score.to_dict(),
            "normalized_score": self.normalized_score.to_dict(),
            "optimized_score": self.optimized_score.to_dict(),
            "states_evaluated": self.states_evaluated,
            "optimality_status": self.optimality_status,
            "excess_layout_crossings_initial": self.excess_layout_crossings_initial,
        }


def initial_layout_state(fixture: Fixture) -> LayoutState:
    order: dict[str, tuple[ChromosomeRef, ...]] = {}
    orientation = {
        chromosome.ref: chromosome.display_orientation
        for chromosome in fixture.chromosomes
    }

    for species in fixture.species_ids:
        chroms = [
            chromosome
            for chromosome in fixture.chromosomes
            if chromosome.ref.species_id == species
        ]
        indexed = list(enumerate(chroms))
        indexed.sort(
            key=lambda item: (
                item[1].display_rank is None,
                item[1].display_rank if item[1].display_rank is not None else item[0],
                item[0],
            )
        )
        order[species] = tuple(chrom.ref for _, chrom in indexed)

    return LayoutState(
        chromosome_order=order,
        chromosome_orientation=orientation,
    )


def _component_map(fixture: Fixture) -> tuple[tuple[frozenset[str], ...], dict[ChromosomeRef, int]]:
    from .component_space import optimization_component_map
    return optimization_component_map(fixture)


def canonicalize_component_order(fixture: Fixture, state: LayoutState) -> LayoutState:
    """Make disconnected components contiguous in one common deterministic order.

    Under R7 this removes representational layout noise without consuming an
    optimization degree of freedom. Order inside each component is preserved.
    """

    components, component_of = _component_map(fixture)
    order: dict[str, tuple[ChromosomeRef, ...]] = {}

    for species in fixture.species_ids:
        current = state.chromosome_order[species]
        by_component: dict[int, list[ChromosomeRef]] = {}
        for ref in current:
            by_component.setdefault(component_of[ref], []).append(ref)
        canonical: list[ChromosomeRef] = []
        for component_id in range(len(components)):
            canonical.extend(by_component.get(component_id, ()))
        order[species] = tuple(canonical)

    return LayoutState(
        chromosome_order=order,
        chromosome_orientation=dict(state.chromosome_orientation),
    )


def _occurrences_by_species_homology(
    fixture: Fixture,
) -> dict[str, dict[str, list[tuple[Chromosome, BlockOccurrence]]]]:
    index: dict[str, dict[str, list[tuple[Chromosome, BlockOccurrence]]]] = {}
    for chromosome in fixture.chromosomes:
        species = chromosome.ref.species_id
        for block in chromosome.blocks:
            index.setdefault(species, {}).setdefault(block.homology_id, []).append(
                (chromosome, block)
            )
    return index


def _anchor_key(
    chromosome: Chromosome,
    block: BlockOccurrence,
    rank: int,
    orientation: int,
) -> tuple[int, float]:
    midpoint = (block.start + block.end) / 2.0
    fraction = midpoint / chromosome.length
    if orientation == -1:
        fraction = 1.0 - fraction
    return rank, fraction


class _Fenwick:
    def __init__(self, size: int) -> None:
        self.tree = [0] * (size + 1)

    def add(self, index: int, value: int = 1) -> None:
        i = index + 1
        while i < len(self.tree):
            self.tree[i] += value
            i += i & -i

    def prefix(self, index: int) -> int:
        total = 0
        i = index + 1
        while i > 0:
            total += self.tree[i]
            i -= i & -i
        return total


def _strict_inversion_count(
    links: list[tuple[tuple[int, float], tuple[int, float]]]
) -> int:
    if len(links) < 2:
        return 0

    right_values = sorted({right for _, right in links})
    right_rank = {value: idx for idx, value in enumerate(right_values)}
    links.sort(key=lambda item: (item[0], item[1]))

    fenwick = _Fenwick(len(right_values))
    seen = 0
    crossings = 0
    i = 0

    # Links sharing an identical left endpoint do not cross one another, so
    # query an equal-left group before inserting that group's right endpoints.
    while i < len(links):
        j = i + 1
        while j < len(links) and links[j][0] == links[i][0]:
            j += 1

        for _, right in links[i:j]:
            rank = right_rank[right]
            leq = fenwick.prefix(rank)
            crossings += seen - leq

        for _, right in links[i:j]:
            fenwick.add(right_rank[right], 1)
            seen += 1

        i = j

    return crossings


def score_crossings(
    fixture: Fixture,
    state: LayoutState,
    *,
    restrict_component_nodes: frozenset[str] | None = None,
) -> CrossingScore:
    """Count unweighted anchor crossings between adjacent species layers.

    Exact scoring currently requires each homology ID to occur at most once per
    species. Ambiguous/duplicated homology is rejected rather than coerced into
    an arbitrary one-to-one mapping.
    """

    occurrence_index = _occurrences_by_species_homology(fixture)
    chromosome_lookup = {chrom.ref: chrom for chrom in fixture.chromosomes}
    rank_lookup: dict[ChromosomeRef, int] = {}
    for species, refs in state.chromosome_order.items():
        for rank, ref in enumerate(refs):
            rank_lookup[ref] = rank

    pair_crossings: list[tuple[str, str, int]] = []
    total_links = 0

    for species1, species2 in zip(fixture.species_ids, fixture.species_ids[1:]):
        homologies = sorted(
            set(occurrence_index.get(species1, {}))
            & set(occurrence_index.get(species2, {}))
        )
        links: list[tuple[tuple[int, float], tuple[int, float]]] = []

        for homology_id in homologies:
            occ1 = occurrence_index[species1][homology_id]
            occ2 = occurrence_index[species2][homology_id]
            if len(occ1) != 1 or len(occ2) != 1:
                raise AmbiguousHomologyError(
                    f"Homology {homology_id!r} is not one-to-one between "
                    f"{species1} and {species2}; exact crossing scoring will not guess"
                )

            chromosome1, block1 = occ1[0]
            chromosome2, block2 = occ2[0]

            if restrict_component_nodes is not None:
                if (
                    chromosome_node_id(chromosome1.ref) not in restrict_component_nodes
                    or chromosome_node_id(chromosome2.ref) not in restrict_component_nodes
                ):
                    continue

            left = _anchor_key(
                chromosome1,
                block1,
                rank_lookup[chromosome1.ref],
                state.chromosome_orientation[chromosome1.ref],
            )
            right = _anchor_key(
                chromosome2,
                block2,
                rank_lookup[chromosome2.ref],
                state.chromosome_orientation[chromosome2.ref],
            )
            links.append((left, right))

        pair_count = _strict_inversion_count(links)
        pair_crossings.append((species1, species2, pair_count))
        total_links += len(links)

    return CrossingScore(
        crossings=sum(count for _, _, count in pair_crossings),
        pair_crossings=tuple(pair_crossings),
        link_count=total_links,
    )


def _component_orientation_assignments(
    fixture: Fixture,
    refs: tuple[ChromosomeRef, ...],
) -> tuple[dict[ChromosomeRef, int], ...]:
    ref_set = set(refs)
    relevant = tuple(
        constraint
        for constraint in fixture.orientation_constraints
        if constraint.a in ref_set or constraint.b in ref_set
    )
    if any(constraint.a not in ref_set or constraint.b not in ref_set for constraint in relevant):
        raise ValueError(
            "Orientation constraint crosses disconnected incidence components; "
            "component-wise exact layout optimization is not yet valid for this case"
        )

    assignments: list[dict[ChromosomeRef, int]] = []
    for bits in product((1, -1), repeat=len(refs)):
        assignment = dict(zip(refs, bits))
        ok = True
        for constraint in relevant:
            abit = 0 if assignment[constraint.a] == 1 else 1
            bbit = 0 if assignment[constraint.b] == 1 else 1
            if (abit ^ bbit) != constraint.xor:
                ok = False
                break
        if ok:
            assignments.append(assignment)

    if not assignments:
        raise InfeasibleOrientationConstraints(
            "No whole-chromosome orientation satisfies all hard orientation constraints"
        )
    return tuple(assignments)


def exact_optimize_small(
    fixture: Fixture,
    *,
    state_cap_per_component: int = 250_000,
) -> ExactLayoutResult:
    """Prove a minimum crossing layout by exact component-wise enumeration.

    Disconnected components are first placed in one common canonical order.
    The solver then enumerates only whole-chromosome permutations and complete
    chromosome reversals *within* each incidence component. Chromosome interiors
    are never altered.
    """

    initial = initial_layout_state(fixture)
    normalized = canonicalize_component_order(fixture, initial)
    initial_score = score_crossings(fixture, initial)
    normalized_score = score_crossings(fixture, normalized)

    components, component_of = _component_map(fixture)

    chosen_orders: dict[int, dict[str, tuple[ChromosomeRef, ...]]] = {}
    chosen_orientation: dict[int, dict[ChromosomeRef, int]] = {}
    states_evaluated = 0

    for component_id, component_nodes in enumerate(components):
        refs = tuple(
            sorted(
                ref for ref in fixture.chromosome_refs
                if component_of[ref] == component_id
            )
        )
        refs_by_species: dict[str, tuple[ChromosomeRef, ...]] = {}
        for species in fixture.species_ids:
            species_refs = tuple(ref for ref in refs if ref.species_id == species)
            if species_refs:
                refs_by_species[species] = species_refs

        rough_states = 2 ** len(refs)
        for species_refs in refs_by_species.values():
            rough_states *= factorial(len(species_refs))
        if rough_states > state_cap_per_component:
            raise SearchSpaceTooLarge(
                f"Component {component_id} has up to {rough_states} exact states; "
                f"cap is {state_cap_per_component}"
            )

        orientation_assignments = _component_orientation_assignments(fixture, refs)
        species = tuple(refs_by_species)
        permutation_options = tuple(
            tuple(permutations(refs_by_species[sp]))
            for sp in species
        )

        best_key: tuple[object, ...] | None = None
        best_orders: dict[str, tuple[ChromosomeRef, ...]] | None = None
        best_orientation: dict[ChromosomeRef, int] | None = None

        for perm_choice in product(*permutation_options):
            local_orders = dict(zip(species, perm_choice))

            # Build one full legal state with all other components left in their
            # normalized order. Cross-component links do not exist; common
            # component ordering prevents inter-component layout crossings.
            full_orders: dict[str, tuple[ChromosomeRef, ...]] = {}
            for sp in fixture.species_ids:
                base = list(normalized.chromosome_order[sp])
                component_positions = [
                    idx for idx, ref in enumerate(base)
                    if component_of[ref] == component_id
                ]
                if sp in local_orders:
                    for idx, ref in zip(component_positions, local_orders[sp]):
                        base[idx] = ref
                full_orders[sp] = tuple(base)

            for local_orientation in orientation_assignments:
                full_orientation = dict(normalized.chromosome_orientation)
                full_orientation.update(local_orientation)
                candidate = LayoutState(
                    chromosome_order=full_orders,
                    chromosome_orientation=full_orientation,
                )
                score = score_crossings(
                    fixture,
                    candidate,
                    restrict_component_nodes=component_nodes,
                ).crossings
                states_evaluated += 1

                order_key = tuple(
                    ref.label
                    for sp in fixture.species_ids
                    for ref in local_orders.get(sp, ())
                )
                orientation_key = tuple(
                    local_orientation[ref] for ref in sorted(local_orientation)
                )
                key = (score, order_key, orientation_key)
                if best_key is None or key < best_key:
                    best_key = key
                    best_orders = {sp: tuple(order) for sp, order in local_orders.items()}
                    best_orientation = dict(local_orientation)

        assert best_orders is not None
        assert best_orientation is not None
        chosen_orders[component_id] = best_orders
        chosen_orientation[component_id] = best_orientation

    final_order: dict[str, tuple[ChromosomeRef, ...]] = {}
    for species in fixture.species_ids:
        refs: list[ChromosomeRef] = []
        for component_id in range(len(components)):
            refs.extend(chosen_orders[component_id].get(species, ()))
        final_order[species] = tuple(refs)

    final_orientation: dict[ChromosomeRef, int] = {}
    for component_id in range(len(components)):
        final_orientation.update(chosen_orientation[component_id])

    optimized = LayoutState(
        chromosome_order=final_order,
        chromosome_orientation=final_orientation,
    )
    optimized_score = score_crossings(fixture, optimized)

    return ExactLayoutResult(
        initial_state=initial,
        normalized_state=normalized,
        optimized_state=optimized,
        initial_score=initial_score,
        normalized_score=normalized_score,
        optimized_score=optimized_score,
        states_evaluated=states_evaluated,
        optimality_status="proven optimum",
    )
