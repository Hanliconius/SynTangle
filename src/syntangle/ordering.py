from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

from .incidence import build_incidence_graph
from .model import ChromosomeRef, Fixture


@dataclass(frozen=True)
class ReversibleChain:
    chromosome: ChromosomeRef
    occurrence_ids: tuple[str, ...]
    homology_ids: tuple[str, ...]

    @property
    def reversed_occurrence_ids(self) -> tuple[str, ...]:
        return tuple(reversed(self.occurrence_ids))

    @property
    def reversed_homology_ids(self) -> tuple[str, ...]:
        return tuple(reversed(self.homology_ids))

    def to_dict(self) -> dict[str, object]:
        return {
            "chromosome": self.chromosome.label,
            "occurrence_ids": list(self.occurrence_ids),
            "homology_ids": list(self.homology_ids),
            "reverse_homology_ids": list(self.reversed_homology_ids),
        }


@dataclass(frozen=True)
class UnresolvedChromosomeOrder:
    species_id: str
    a: ChromosomeRef
    b: ChromosomeRef
    component_id: int

    def to_dict(self) -> dict[str, object]:
        return {
            "species": self.species_id,
            "a": self.a.label,
            "b": self.b.label,
            "component_id": self.component_id,
        }


@dataclass(frozen=True)
class ComponentOrderEquivalence:
    component_id: int
    chromosomes_by_species: tuple[tuple[str, tuple[ChromosomeRef, ...]], ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "component_id": self.component_id,
            "chromosomes_by_species": {
                species: [ref.label for ref in refs]
                for species, refs in self.chromosomes_by_species
            },
        }


@dataclass(frozen=True)
class OrderingConstraintState:
    reversible_chains: tuple[ReversibleChain, ...]
    equivalent_components: tuple[ComponentOrderEquivalence, ...]
    unresolved_chromosome_orders: tuple[UnresolvedChromosomeOrder, ...]
    forced_internal_chain_count: int

    def to_dict(self) -> dict[str, object]:
        return {
            "forced_internal_chain_count": self.forced_internal_chain_count,
            "reversible_chains": [chain.to_dict() for chain in self.reversible_chains],
            "equivalent_components": [
                component.to_dict() for component in self.equivalent_components
            ],
            "unresolved_chromosome_orders": [
                relation.to_dict() for relation in self.unresolved_chromosome_orders
            ],
        }


def derive_ordering_constraints(fixture: Fixture) -> OrderingConstraintState:
    """Derive only order facts that are already implied by the input.

    Two distinct ideas are kept separate:

    * Within a chromosome, block occurrence order is a HARD reversible chain:
      native order or complete reverse, never an internal permutation.
    * Between whole chromosomes in the same connected homology component,
      relative display order is UNRESOLVED unless another hard constraint says
      otherwise.

    Disconnected incidence components are EQUIVALENT at the global display
    level under R7 when a common component order is used consistently across
    species. This function deliberately does *not* turn same-chromosome block
    adjacency into a cross-species chromosome-adjacency rule; that would add
    biological assumptions not present in the fixture.
    """

    graph = build_incidence_graph(fixture)
    components = graph.connected_components()

    node_to_component: dict[str, int] = {}
    for component_id, nodes in enumerate(components):
        for node in nodes:
            node_to_component[node] = component_id

    reversible_chains: list[ReversibleChain] = []
    for chromosome in sorted(fixture.chromosomes, key=lambda c: c.ref):
        reversible_chains.append(
            ReversibleChain(
                chromosome=chromosome.ref,
                occurrence_ids=tuple(block.occurrence_id for block in chromosome.blocks),
                homology_ids=tuple(block.homology_id for block in chromosome.blocks),
            )
        )

    equivalent_components: list[ComponentOrderEquivalence] = []
    for component_id, nodes in enumerate(components):
        by_species: dict[str, list[ChromosomeRef]] = {}
        for chromosome in fixture.chromosomes:
            cnode = f"chrom::{chromosome.ref.species_id}::{chromosome.ref.chromosome_id}"
            if cnode not in nodes:
                continue
            by_species.setdefault(chromosome.ref.species_id, []).append(chromosome.ref)
        equivalent_components.append(
            ComponentOrderEquivalence(
                component_id=component_id,
                chromosomes_by_species=tuple(
                    (species, tuple(sorted(refs)))
                    for species, refs in sorted(by_species.items())
                ),
            )
        )

    unresolved: list[UnresolvedChromosomeOrder] = []
    species_to_refs: dict[str, list[ChromosomeRef]] = {}
    for chromosome in fixture.chromosomes:
        species_to_refs.setdefault(chromosome.ref.species_id, []).append(chromosome.ref)

    for species_id, refs in sorted(species_to_refs.items()):
        for a, b in combinations(sorted(refs), 2):
            a_component = node_to_component[f"chrom::{a.species_id}::{a.chromosome_id}"]
            b_component = node_to_component[f"chrom::{b.species_id}::{b.chromosome_id}"]
            if a_component == b_component:
                unresolved.append(
                    UnresolvedChromosomeOrder(
                        species_id=species_id,
                        a=a,
                        b=b,
                        component_id=a_component,
                    )
                )

    return OrderingConstraintState(
        reversible_chains=tuple(reversible_chains),
        equivalent_components=tuple(equivalent_components),
        unresolved_chromosome_orders=tuple(unresolved),
        forced_internal_chain_count=sum(
            len(chain.occurrence_ids) >= 2 for chain in reversible_chains
        ),
    )
