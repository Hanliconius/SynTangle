from __future__ import annotations

import csv
from pathlib import Path

from .fixtures import fixture_from_dict
from .model import Fixture, FixtureValidationError


REQUIRED_FILES = ("species.tsv", "chromosomes.tsv", "occurrences.tsv")


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_validation_bundle(directory: str | Path) -> Fixture:
    """Load the public/input side of a simulator validation bundle.

    Hidden evolutionary and display-tangle logs may coexist in the same
    directory, but this importer deliberately ignores them.
    """

    root = Path(directory)
    for filename in REQUIRED_FILES:
        if not (root / filename).is_file():
            raise FixtureValidationError(
                f"Validation bundle is missing required file: {filename}"
            )

    species_rows = _read_tsv(root / "species.tsv")
    chromosome_rows = _read_tsv(root / "chromosomes.tsv")
    occurrence_rows = _read_tsv(root / "occurrences.tsv")

    if not species_rows:
        raise FixtureValidationError("Validation bundle has no species rows")

    try:
        species_rows.sort(key=lambda row: int(row["species_rank"]))
    except (KeyError, TypeError, ValueError) as exc:
        raise FixtureValidationError(
            "species.tsv requires integer species_rank"
        ) from exc

    species_ids = [row.get("species_id", "") for row in species_rows]
    if any(not species_id for species_id in species_ids):
        raise FixtureValidationError("species.tsv contains an empty species_id")
    if len(species_ids) != len(set(species_ids)):
        raise FixtureValidationError("species.tsv contains duplicate species_id")

    blocks_by_chromosome: dict[tuple[str, str], list[dict[str, object]]] = {}
    for row in occurrence_rows:
        try:
            species_id = row["species_id"]
            chromosome_id = row["chromosome_id"]
            block = {
                "occurrence_id": row["occurrence_id"],
                "homology_id": row["homology_id"],
                "start": float(row["start"]),
                "end": float(row["end"]),
                "strand": row["strand"],
            }
        except (KeyError, TypeError, ValueError) as exc:
            raise FixtureValidationError(
                "occurrences.tsv has an invalid or missing required field"
            ) from exc
        blocks_by_chromosome.setdefault(
            (species_id, chromosome_id), []
        ).append(block)

    chromosomes_by_species: dict[str, list[dict[str, object]]] = {
        species_id: [] for species_id in species_ids
    }
    for row in chromosome_rows:
        try:
            species_id = row["species_id"]
            chromosome_id = row["chromosome_id"]
            length = float(row["length"])
            display_rank = int(row["display_rank"])
        except (KeyError, TypeError, ValueError) as exc:
            raise FixtureValidationError(
                "chromosomes.tsv has an invalid or missing required field"
            ) from exc

        if species_id not in chromosomes_by_species:
            raise FixtureValidationError(
                f"chromosomes.tsv references unknown species: {species_id}"
            )

        chromosomes_by_species[species_id].append(
            {
                "id": chromosome_id,
                "length": length,
                "display_rank": display_rank,
                "blocks": blocks_by_chromosome.get(
                    (species_id, chromosome_id), []
                ),
            }
        )

    for species_id in species_ids:
        chromosomes_by_species[species_id].sort(
            key=lambda chromosome: int(chromosome["display_rank"])
        )

    metadata_path = root / "metadata.tsv"
    fixture_id = root.name
    title = f"Simulator validation bundle: {fixture_id}"
    purpose = "Forward-simulated extant chromosomes with hidden truth withheld"

    if metadata_path.is_file():
        metadata_rows = _read_tsv(metadata_path)
        metadata = {
            row.get("key", ""): row.get("value", "")
            for row in metadata_rows
            if row.get("key")
        }
        fixture_id = metadata.get("fixture_id", fixture_id)
        title = metadata.get("title", title)
        purpose = metadata.get("purpose", purpose)

    data = {
        "fixture_version": 1,
        "id": fixture_id,
        "title": title,
        "purpose": purpose,
        "species": [
            {
                "id": species_id,
                "chromosomes": chromosomes_by_species[species_id],
            }
            for species_id in species_ids
        ],
        "orientation_constraints": [],
        "expected": {},
    }
    return fixture_from_dict(data)
