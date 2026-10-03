from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, product

from .incidence import chromosome_node_id
from .layout import AmbiguousHomologyError
from .model import ChromosomeRef, Fixture
from .orientation_space import OrientationBasis


@dataclass(frozen=True)
class AnchorCell:
    layer_index: int
    left_ref: ChromosomeRef
    right_ref: ChromosomeRef
    anchors: tuple[tuple[float, float], ...]


@dataclass(frozen=True)
class RelaxedCrossingBound:
    """Safe crossing lower bound after relaxing global chromosome-order consistency.

    Every pair of links between adjacent species belongs to one of four classes:

    1. same chromosome on both sides: exact internal crossing contribution;
    2. same left chromosome, different right chromosomes: choose the better
       pairwise right-chromosome precedence independently;
    3. different left chromosomes, same right chromosome: choose the better
       pairwise left-chromosome precedence independently;
    4. different chromosomes on both sides: relaxed contribution zero.

    Classes 2 and 3 deliberately ignore transitivity among chromosome-order
    decisions, so the result can underestimate the achievable score but can
    never overestimate it. That makes it a valid branch-and-bound lower bound.
    """

    basis: OrientationBasis
    cells: tuple[AnchorCell, ...]
    constant_interleaving_lower_bound: int
    group_index_by_ref: dict[ChromosomeRef, int]

    def lower_bound(self, bits: tuple[int | None, ...]) -> int:
        if len(bits) != len(self.basis.free_flip_groups):
            raise ValueError("Orientation bit vector does not match orientation basis")

        total = self.constant_interleaving_lower_bound

        for cell in self.cells:
            group_indices = tuple(
                sorted(
                    {
                        self.group_index_by_ref[cell.left_ref],
                        self.group_index_by_ref[cell.right_ref],
                    }
                )
            )
            unresolved = tuple(
                index for index in group_indices if bits[index] is None
            )

            best: int | None = None
            for values in product((0, 1), repeat=len(unresolved)):
                assignment = dict(zip(unresolved, values))
                left_orientation = _orientation_for_ref(
                    self.basis,
                    self.group_index_by_ref,
                    cell.left_ref,
                    bits,
                    assignment,
                )
                right_orientation = _orientation_for_ref(
                    self.basis,
                    self.group_index_by_ref,
                    cell.right_ref,
                    bits,
                    assignment,
                )
                value = _cell_internal_crossings(
                    cell.anchors,
                    left_orientation,
                    right_orientation,
                )
                if best is None or value < best:
                    best = value

            total += 0 if best is None else best

        return total


def _orientation_for_ref(
    basis: OrientationBasis,
    group_index_by_ref: dict[ChromosomeRef, int],
    ref: ChromosomeRef,
    bits: tuple[int | None, ...],
    temporary: dict[int, int],
) -> int:
    group_index = group_index_by_ref[ref]
    bit = bits[group_index]
    if bit is None:
        bit = temporary[group_index]
    orientation = basis.base_assignment[ref]
    return -orientation if bit else orientation


def _oriented_fraction(value: float, orientation: int) -> float:
    return value if orientation == 1 else 1.0 - value


def _cell_internal_crossings(
    anchors: tuple[tuple[float, float], ...],
    left_orientation: int,
    right_orientation: int,
) -> int:
    count = 0
    transformed = tuple(
        (
            _oriented_fraction(left, left_orientation),
            _oriented_fraction(right, right_orientation),
        )
        for left, right in anchors
    )
    for (left1, right1), (left2, right2) in combinations(transformed, 2):
        if (left1 < left2 and right1 > right2) or (
            left1 > left2 and right1 < right2
        ):
            count += 1
    return count


def _pairwise_interleave_minimum(
    first: tuple[tuple[float, float], ...],
    second: tuple[tuple[float, float], ...],
    *,
    coordinate: int,
) -> int:
    first_before_second = 0
    second_before_first = 0

    for anchor_a in first:
        for anchor_b in second:
            value_a = anchor_a[coordinate]
            value_b = anchor_b[coordinate]
            if value_a > value_b:
                first_before_second += 1
            elif value_b > value_a:
                second_before_first += 1

    return min(first_before_second, second_before_first)


def build_relaxed_crossing_bound(
    fixture: Fixture,
    component_nodes: frozenset[str],
    basis: OrientationBasis,
) -> RelaxedCrossingBound:
    chromosome_lookup = {
        chromosome.ref: chromosome for chromosome in fixture.chromosomes
    }
    occurrence_index: dict[
        str,
        dict[str, list[tuple[ChromosomeRef, object]]],
    ] = {}

    for chromosome in fixture.chromosomes:
        species = chromosome.ref.species_id
        for block in chromosome.blocks:
            occurrence_index.setdefault(species, {}).setdefault(
                block.homology_id, []
            ).append((chromosome.ref, block))

    cell_map: dict[
        tuple[int, ChromosomeRef, ChromosomeRef],
        list[tuple[float, float]],
    ] = {}

    for layer_index, (left_species, right_species) in enumerate(
        zip(fixture.species_ids, fixture.species_ids[1:])
    ):
        shared = sorted(
            set(occurrence_index.get(left_species, {}))
            & set(occurrence_index.get(right_species, {}))
        )

        for homology_id in shared:
            left_occ = occurrence_index[left_species][homology_id]
            right_occ = occurrence_index[right_species][homology_id]
            if len(left_occ) != 1 or len(right_occ) != 1:
                raise AmbiguousHomologyError(
                    f"Homology {homology_id!r} is not one-to-one between "
                    f"{left_species} and {right_species}"
                )

            left_ref, left_block = left_occ[0]
            right_ref, right_block = right_occ[0]
            if (
                chromosome_node_id(left_ref) not in component_nodes
                or chromosome_node_id(right_ref) not in component_nodes
            ):
                continue

            left_chromosome = chromosome_lookup[left_ref]
            right_chromosome = chromosome_lookup[right_ref]
            left_fraction = (
                (left_block.start + left_block.end)
                / 2.0
                / left_chromosome.length
            )
            right_fraction = (
                (right_block.start + right_block.end)
                / 2.0
                / right_chromosome.length
            )
            cell_map.setdefault(
                (layer_index, left_ref, right_ref), []
            ).append((left_fraction, right_fraction))

    cells = tuple(
        AnchorCell(
            layer_index=layer_index,
            left_ref=left_ref,
            right_ref=right_ref,
            anchors=tuple(anchors),
        )
        for (layer_index, left_ref, right_ref), anchors
        in sorted(
            cell_map.items(),
            key=lambda item: (
                item[0][0],
                item[0][1].label,
                item[0][2].label,
            ),
        )
    )

    by_layer_left: dict[
        tuple[int, ChromosomeRef],
        list[AnchorCell],
    ] = {}
    by_layer_right: dict[
        tuple[int, ChromosomeRef],
        list[AnchorCell],
    ] = {}

    for cell in cells:
        by_layer_left.setdefault(
            (cell.layer_index, cell.left_ref), []
        ).append(cell)
        by_layer_right.setdefault(
            (cell.layer_index, cell.right_ref), []
        ).append(cell)

    constant_interleaving = 0

    # Same left chromosome, different right chromosomes. Only the left
    # within-chromosome coordinate matters; the right relative chromosome
    # precedence is relaxed independently for each chromosome pair.
    for grouped_cells in by_layer_left.values():
        for first, second in combinations(grouped_cells, 2):
            constant_interleaving += _pairwise_interleave_minimum(
                first.anchors,
                second.anchors,
                coordinate=0,
            )

    # Different left chromosomes, same right chromosome. Symmetric case.
    for grouped_cells in by_layer_right.values():
        for first, second in combinations(grouped_cells, 2):
            constant_interleaving += _pairwise_interleave_minimum(
                first.anchors,
                second.anchors,
                coordinate=1,
            )

    group_index_by_ref: dict[ChromosomeRef, int] = {}
    for group_index, group in enumerate(basis.free_flip_groups):
        for ref in group:
            group_index_by_ref[ref] = group_index

    return RelaxedCrossingBound(
        basis=basis,
        cells=cells,
        constant_interleaving_lower_bound=constant_interleaving,
        group_index_by_ref=group_index_by_ref,
    )
