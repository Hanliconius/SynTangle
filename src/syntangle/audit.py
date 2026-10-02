from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json

from .layout import ExactLayoutResult
from .model import Fixture


@dataclass(frozen=True)
class ChromosomeMove:
    species_id: str
    chromosome_id: str
    from_rank: int
    to_rank: int

    def to_dict(self) -> dict[str, object]:
        return {
            "species": self.species_id,
            "chromosome": self.chromosome_id,
            "from_rank": self.from_rank,
            "to_rank": self.to_rank,
        }


@dataclass(frozen=True)
class LayoutAudit:
    fixture_id: str
    input_fingerprint: str
    whole_chromosome_moves: tuple[ChromosomeMove, ...]
    whole_chromosome_flips: tuple[str, ...]
    orientation_constraints: tuple[dict[str, object], ...]
    objective_before: int
    objective_after: int
    excess_layout_crossings_initial: int
    optimality_status: str
    states_evaluated: int

    def to_dict(self) -> dict[str, object]:
        return {
            "fixture_id": self.fixture_id,
            "input_fingerprint": self.input_fingerprint,
            "whole_chromosome_moves": [move.to_dict() for move in self.whole_chromosome_moves],
            "whole_chromosome_flips": list(self.whole_chromosome_flips),
            "orientation_constraints": list(self.orientation_constraints),
            "objective": {
                "name": "unweighted_adjacent_layer_crossings",
                "before": self.objective_before,
                "after": self.objective_after,
                "excess_initial": self.excess_layout_crossings_initial,
            },
            "optimality_status": self.optimality_status,
            "states_evaluated": self.states_evaluated,
        }


def fixture_fingerprint(fixture: Fixture) -> str:
    payload = {
        "fixture_version": fixture.fixture_version,
        "fixture_id": fixture.fixture_id,
        "chromosomes": [
            {
                "species": chromosome.ref.species_id,
                "chromosome": chromosome.ref.chromosome_id,
                "length": chromosome.length,
                "display_rank": chromosome.display_rank,
                "blocks": [
                    {
                        "occurrence_id": block.occurrence_id,
                        "homology_id": block.homology_id,
                        "start": block.start,
                        "end": block.end,
                        "strand": block.strand,
                    }
                    for block in chromosome.blocks
                ],
            }
            for chromosome in fixture.chromosomes
        ],
        "orientation_constraints": [
            {
                "a": constraint.a.label,
                "b": constraint.b.label,
                "xor": constraint.xor,
            }
            for constraint in fixture.orientation_constraints
        ],
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return sha256(canonical.encode("utf-8")).hexdigest()


def build_layout_audit(fixture: Fixture, result: ExactLayoutResult) -> LayoutAudit:
    moves: list[ChromosomeMove] = []

    for species in fixture.species_ids:
        initial_order = result.initial_state.chromosome_order[species]
        final_order = result.optimized_state.chromosome_order[species]
        initial_rank = {ref: index + 1 for index, ref in enumerate(initial_order)}
        final_rank = {ref: index + 1 for index, ref in enumerate(final_order)}
        for ref in sorted(initial_rank):
            if initial_rank[ref] != final_rank[ref]:
                moves.append(
                    ChromosomeMove(
                        species_id=species,
                        chromosome_id=ref.chromosome_id,
                        from_rank=initial_rank[ref],
                        to_rank=final_rank[ref],
                    )
                )

    flips = tuple(
        ref.label
        for ref in sorted(result.optimized_state.chromosome_orientation)
        if result.initial_state.chromosome_orientation[ref]
        != result.optimized_state.chromosome_orientation[ref]
    )

    constraints = tuple(
        {
            "a": constraint.a.label,
            "b": constraint.b.label,
            "xor": constraint.xor,
        }
        for constraint in fixture.orientation_constraints
    )

    return LayoutAudit(
        fixture_id=fixture.fixture_id,
        input_fingerprint=fixture_fingerprint(fixture),
        whole_chromosome_moves=tuple(moves),
        whole_chromosome_flips=flips,
        orientation_constraints=constraints,
        objective_before=result.initial_score.crossings,
        objective_after=result.optimized_score.crossings,
        excess_layout_crossings_initial=result.excess_layout_crossings_initial,
        optimality_status=result.optimality_status,
        states_evaluated=result.states_evaluated,
    )
