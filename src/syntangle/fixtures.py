from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

from .model import (
    BlockOccurrence,
    Chromosome,
    ChromosomeRef,
    Fixture,
    FixtureValidationError,
    OrientationConstraint,
    parse_chromosome_ref,
)


def _require(mapping: Mapping[str, Any], key: str, context: str) -> Any:
    if key not in mapping:
        raise FixtureValidationError(f"Missing {key!r} in {context}")
    return mapping[key]


def fixture_from_dict(data: Mapping[str, Any]) -> Fixture:
    version = _require(data, "fixture_version", "fixture")
    if version != 1:
        raise FixtureValidationError(f"Unsupported fixture_version: {version!r}")

    fixture_id = str(_require(data, "id", "fixture"))
    title = str(_require(data, "title", "fixture"))
    purpose = str(_require(data, "purpose", "fixture"))
    species_rows = _require(data, "species", "fixture")
    if not isinstance(species_rows, list) or not species_rows:
        raise FixtureValidationError("fixture.species must be a non-empty list")

    chromosomes: list[Chromosome] = []
    seen_species: set[str] = set()
    seen_chromosomes: set[ChromosomeRef] = set()
    seen_occurrences: set[str] = set()

    for species in species_rows:
        species_id = str(_require(species, "id", "species"))
        if species_id in seen_species:
            raise FixtureValidationError(f"Duplicate species id: {species_id}")
        seen_species.add(species_id)

        chrom_rows = _require(species, "chromosomes", f"species {species_id}")
        if not isinstance(chrom_rows, list):
            raise FixtureValidationError(f"species {species_id}.chromosomes must be a list")

        seen_display_ranks: set[int] = set()
        for chrom_row in chrom_rows:
            chromosome_id = str(_require(chrom_row, "id", f"species {species_id} chromosome"))
            ref = ChromosomeRef(species_id=species_id, chromosome_id=chromosome_id)
            if ref in seen_chromosomes:
                raise FixtureValidationError(f"Duplicate chromosome: {ref.label}")
            seen_chromosomes.add(ref)

            length = float(_require(chrom_row, "length", ref.label))
            if length <= 0:
                raise FixtureValidationError(f"Chromosome length must be > 0: {ref.label}")

            display_rank_raw = chrom_row.get("display_rank")
            display_rank = None if display_rank_raw is None else int(display_rank_raw)
            if display_rank is not None:
                if display_rank < 1:
                    raise FixtureValidationError(f"display_rank must be >= 1: {ref.label}")
                if display_rank in seen_display_ranks:
                    raise FixtureValidationError(
                        f"Duplicate display_rank {display_rank} within species {species_id}"
                    )
                seen_display_ranks.add(display_rank)

            block_rows = _require(chrom_row, "blocks", ref.label)
            if not isinstance(block_rows, list):
                raise FixtureValidationError(f"blocks must be a list: {ref.label}")

            blocks: list[BlockOccurrence] = []
            for block_row in block_rows:
                occurrence_id = str(_require(block_row, "occurrence_id", ref.label))
                if occurrence_id in seen_occurrences:
                    raise FixtureValidationError(f"Duplicate occurrence_id: {occurrence_id}")
                seen_occurrences.add(occurrence_id)

                homology_id = str(_require(block_row, "homology_id", occurrence_id))
                start = float(_require(block_row, "start", occurrence_id))
                end = float(_require(block_row, "end", occurrence_id))
                strand = str(_require(block_row, "strand", occurrence_id))

                if not homology_id:
                    raise FixtureValidationError(f"Empty homology_id: {occurrence_id}")
                if strand not in {"+", "-", "?"}:
                    raise FixtureValidationError(f"Invalid strand {strand!r}: {occurrence_id}")
                if start < 0 or end <= start or end > length:
                    raise FixtureValidationError(
                        f"Invalid coordinates {start}-{end} on {ref.label} length {length}"
                    )

                blocks.append(
                    BlockOccurrence(
                        occurrence_id=occurrence_id,
                        homology_id=homology_id,
                        start=start,
                        end=end,
                        strand=strand,
                    )
                )

            # Genomic coordinate order is canonical and immutable downstream.
            blocks.sort(key=lambda b: (b.start, b.end, b.occurrence_id))
            chromosomes.append(
                Chromosome(
                    ref=ref,
                    length=length,
                    display_rank=display_rank,
                    blocks=tuple(blocks),
                )
            )

    constraints: list[OrientationConstraint] = []
    for raw in data.get("orientation_constraints", []):
        a = parse_chromosome_ref(_require(raw, "a", "orientation constraint"))
        b = parse_chromosome_ref(_require(raw, "b", "orientation constraint"))
        xor = int(_require(raw, "xor", "orientation constraint"))
        if xor not in {0, 1}:
            raise FixtureValidationError(f"Orientation xor must be 0/1: {xor!r}")
        if a not in seen_chromosomes or b not in seen_chromosomes:
            raise FixtureValidationError(
                f"Orientation constraint references unknown chromosome: {a.label}, {b.label}"
            )
        constraints.append(OrientationConstraint(a=a, b=b, xor=xor))

    expected = data.get("expected", {})
    if not isinstance(expected, Mapping):
        raise FixtureValidationError("expected must be an object")

    return Fixture(
        fixture_version=version,
        fixture_id=fixture_id,
        title=title,
        purpose=purpose,
        chromosomes=tuple(chromosomes),
        orientation_constraints=tuple(constraints),
        expected=expected,
    )


def load_fixture(path: str | Path) -> Fixture:
    fixture_path = Path(path)
    with fixture_path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, Mapping):
        raise FixtureValidationError("Fixture JSON root must be an object")
    return fixture_from_dict(data)
