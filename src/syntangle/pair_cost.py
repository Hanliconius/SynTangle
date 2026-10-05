from __future__ import annotations

from .incidence import chromosome_node_id
from .layout import (
    AmbiguousHomologyError,
    _anchor_key,
    _occurrences_by_species_homology,
    _strict_inversion_count,
)
from .model import ChromosomeRef, Fixture


def pair_component_crossings(
    fixture: Fixture,
    species1: str,
    species2: str,
    order1: tuple[ChromosomeRef, ...],
    order2: tuple[ChromosomeRef, ...],
    orientation: dict[ChromosomeRef, int],
    component_nodes: frozenset[str],
    *,
    occurrence_index=None,
) -> int:
    """Crossings contributed by one incidence component between two layers."""

    if not order1 or not order2:
        return 0

    if occurrence_index is None:
        occurrence_index = _occurrences_by_species_homology(fixture)
    rank1 = {ref: rank for rank, ref in enumerate(order1)}
    rank2 = {ref: rank for rank, ref in enumerate(order2)}
    homologies = sorted(
        set(occurrence_index.get(species1, {}))
        & set(occurrence_index.get(species2, {}))
    )
    links = []

    for homology_id in homologies:
        occ1 = occurrence_index[species1][homology_id]
        occ2 = occurrence_index[species2][homology_id]
        if len(occ1) != 1 or len(occ2) != 1:
            raise AmbiguousHomologyError(
                f"Homology {homology_id!r} is not one-to-one between "
                f"{species1} and {species2}; exact crossing scoring will not guess"
            )

        chrom1, block1 = occ1[0]
        chrom2, block2 = occ2[0]
        if (
            chromosome_node_id(chrom1.ref) not in component_nodes
            or chromosome_node_id(chrom2.ref) not in component_nodes
        ):
            continue

        links.append(
            (
                _anchor_key(
                    chrom1, block1, rank1[chrom1.ref], orientation[chrom1.ref]
                ),
                _anchor_key(
                    chrom2, block2, rank2[chrom2.ref], orientation[chrom2.ref]
                ),
            )
        )

    return _strict_inversion_count(links)


class PreparedPairCrossings:
    """Worker-local immutable anchor data; no search states or exclusions stored."""

    def __init__(self, fixture, component_nodes):
        index = _occurrences_by_species_homology(fixture)
        self.links = {}
        for species1, species2 in zip(fixture.species_ids, fixture.species_ids[1:]):
            links = []
            for homology in sorted(set(index.get(species1, {})) &
                                   set(index.get(species2, {}))):
                left, right = index[species1][homology], index[species2][homology]
                if len(left) != 1 or len(right) != 1:
                    raise AmbiguousHomologyError(f"Homology {homology!r} is not one-to-one")
                chrom1, block1 = left[0]
                chrom2, block2 = right[0]
                if any(chromosome_node_id(chrom.ref) not in component_nodes
                       for chrom in (chrom1, chrom2)):
                    continue
                # Precompute both orientations independently, including strict
                # floating-point ties; do not infer reverse crossings by subtraction.
                endpoints = tuple((chrom.ref, {
                    sign: _anchor_key(chrom, block, 0, sign)[1]
                    for sign in (1, -1)})
                    for chrom, block in ((chrom1, block1), (chrom2, block2)))
                links.append(endpoints)
            self.links[species1, species2] = tuple(links)

    def score(self, species1, species2, order1, order2, orientation):
        if not order1 or not order2:
            return 0
        rank1 = {ref: rank for rank, ref in enumerate(order1)}
        rank2 = {ref: rank for rank, ref in enumerate(order2)}
        return _strict_inversion_count([
            ((rank1[left], positions1[orientation[left]]),
             (rank2[right], positions2[orientation[right]]))
            for (left, positions1), (right, positions2)
            in self.links[species1, species2]])
