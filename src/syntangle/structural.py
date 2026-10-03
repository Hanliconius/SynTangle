from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from .decomposition import GraphDecomposition, decompose_incidence_graph
from .incidence import IncidenceEdge, IncidenceGraph
from .model import ChromosomeRef, Fixture


@dataclass(frozen=True)
class StructuralBundle:
    bundle_id: str
    signature: tuple[tuple[ChromosomeRef, int], ...]
    homology_ids: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "bundle_id": self.bundle_id,
            "signature": [
                {"chromosome": ref.label, "occurrence_count": count}
                for ref, count in self.signature
            ],
            "homology_ids": list(self.homology_ids),
            "homology_count": len(self.homology_ids),
        }


@dataclass(frozen=True)
class StructuralProjection:
    graph: IncidenceGraph
    bundles: tuple[StructuralBundle, ...]
    decomposition: GraphDecomposition

    @property
    def cycle_rank_total(self) -> int:
        return sum(
            self.graph.summarize_component(component).cycle_rank
            for component in self.graph.connected_components()
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "bundle_count": len(self.bundles),
            "cycle_rank_total": self.cycle_rank_total,
            "hard_kernel_count": len(self.decomposition.hard_kernels),
            "core_vertex_count": len(self.decomposition.core_node_ids),
            "core_edge_count": len(self.decomposition.core_edge_ids),
            "bundles": [bundle.to_dict() for bundle in self.bundles],
        }


def _homology_signatures(
    fixture: Fixture,
) -> dict[str, tuple[tuple[ChromosomeRef, int], ...]]:
    occurrences: dict[str, Counter[ChromosomeRef]] = {}
    for chromosome in fixture.chromosomes:
        for block in chromosome.blocks:
            occurrences.setdefault(block.homology_id, Counter())[
                chromosome.ref
            ] += 1

    return {
        homology_id: tuple(sorted(counter.items()))
        for homology_id, counter in occurrences.items()
    }


def build_structural_projection(fixture: Fixture) -> StructuralProjection:
    """Collapse redundant homology groups with identical chromosome incidence.

    This projection is for chromosome-structural decomposition only. It does
    not replace the full block/anchor data used for crossing scoring, internal
    order, or inversion evidence.

    Multiple homology groups that touch the exact same chromosome multiset are
    bundled into one structural relationship. Copy multiplicity *within* one
    homology group is retained in the signature and in graph edge multiplicity.
    """

    signatures = _homology_signatures(fixture)
    grouped: dict[
        tuple[tuple[ChromosomeRef, int], ...], list[str]
    ] = {}
    for homology_id, signature in signatures.items():
        grouped.setdefault(signature, []).append(homology_id)

    ordered_signatures = sorted(
        grouped,
        key=lambda signature: tuple(
            (ref.species_id, ref.chromosome_id, count)
            for ref, count in signature
        ),
    )

    bundles: list[StructuralBundle] = []
    edges: list[IncidenceEdge] = []

    for index, signature in enumerate(ordered_signatures, start=1):
        bundle_id = f"SB{index:04d}"
        homology_ids = tuple(sorted(grouped[signature]))
        bundle = StructuralBundle(
            bundle_id=bundle_id,
            signature=signature,
            homology_ids=homology_ids,
        )
        bundles.append(bundle)

        for ref, count in signature:
            for copy_index in range(1, count + 1):
                edges.append(
                    IncidenceEdge(
                        occurrence_id=(
                            f"{bundle_id}:{ref.label}:copy{copy_index}"
                        ),
                        chromosome=ref,
                        homology_id=bundle_id,
                    )
                )

    graph = IncidenceGraph(
        chromosome_nodes=fixture.chromosome_refs,
        homology_ids=(bundle.bundle_id for bundle in bundles),
        edges=edges,
    )
    decomposition = decompose_incidence_graph(graph)
    return StructuralProjection(
        graph=graph,
        bundles=tuple(bundles),
        decomposition=decomposition,
    )
