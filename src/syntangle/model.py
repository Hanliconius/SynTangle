from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping


class FixtureValidationError(ValueError):
    """Raised when a synthetic fixture violates the canonical data model."""


@dataclass(frozen=True, order=True)
class ChromosomeRef:
    species_id: str
    chromosome_id: str

    @property
    def label(self) -> str:
        return f"{self.species_id}:{self.chromosome_id}"


@dataclass(frozen=True)
class BlockOccurrence:
    occurrence_id: str
    homology_id: str
    start: float
    end: float
    strand: str


@dataclass(frozen=True)
class Chromosome:
    ref: ChromosomeRef
    length: float
    display_rank: int | None
    blocks: tuple[BlockOccurrence, ...]


@dataclass(frozen=True)
class OrientationConstraint:
    a: ChromosomeRef
    b: ChromosomeRef
    xor: int


@dataclass(frozen=True)
class Fixture:
    fixture_version: int
    fixture_id: str
    title: str
    purpose: str
    chromosomes: tuple[Chromosome, ...]
    orientation_constraints: tuple[OrientationConstraint, ...]
    expected: Mapping[str, Any]

    @property
    def species_ids(self) -> tuple[str, ...]:
        seen: set[str] = set()
        ordered: list[str] = []
        for chromosome in self.chromosomes:
            if chromosome.ref.species_id not in seen:
                seen.add(chromosome.ref.species_id)
                ordered.append(chromosome.ref.species_id)
        return tuple(ordered)

    @property
    def chromosome_refs(self) -> tuple[ChromosomeRef, ...]:
        return tuple(chrom.ref for chrom in self.chromosomes)

    @property
    def homology_ids(self) -> tuple[str, ...]:
        return tuple(sorted({block.homology_id for chrom in self.chromosomes for block in chrom.blocks}))


def parse_chromosome_ref(value: str) -> ChromosomeRef:
    if not isinstance(value, str) or ":" not in value:
        raise FixtureValidationError(f"Invalid chromosome reference: {value!r}")
    species_id, chromosome_id = value.split(":", 1)
    if not species_id or not chromosome_id:
        raise FixtureValidationError(f"Invalid chromosome reference: {value!r}")
    return ChromosomeRef(species_id=species_id, chromosome_id=chromosome_id)
