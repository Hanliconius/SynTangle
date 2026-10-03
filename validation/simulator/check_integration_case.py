from __future__ import annotations

import json
import sys
from pathlib import Path

from syntangle import (
    build_incidence_graph,
    build_layout_audit,
    build_structural_projection,
    exact_optimize_layer_dp,
    fixture_fingerprint,
    load_validation_bundle,
)


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit(
            "Usage: python validation/simulator/check_integration_case.py BUNDLE_DIR"
        )

    bundle = Path(sys.argv[1])
    fixture = load_validation_bundle(bundle)

    if fixture.species_ids != ("speciesA", "speciesB", "speciesC"):
        raise AssertionError(f"Unexpected species order: {fixture.species_ids}")

    public_fingerprint = fixture_fingerprint(fixture)

    state_rows = (bundle / "input_display_state.tsv").read_text(
        encoding="utf-8"
    ).splitlines()
    state_header = state_rows[0].split("\t")
    state_index = {name: idx for idx, name in enumerate(state_header)}
    expected_orientation = {}
    for line in state_rows[1:]:
        fields = line.split("\t")
        expected_orientation[
            (fields[state_index["species"]], fields[state_index["chrom"]])
        ] = int(fields[state_index["orientation"]])

    observed_orientation = {
        (chromosome.ref.species_id, chromosome.ref.chromosome_id):
            chromosome.display_orientation
        for chromosome in fixture.chromosomes
    }
    if observed_orientation != expected_orientation:
        raise AssertionError("Public input display orientation was not preserved")

    hidden_evolution = bundle / "hidden_evolution_log.tsv"
    hidden_tangle = bundle / "hidden_tangle_log.tsv"
    if not hidden_evolution.is_file() or not hidden_tangle.is_file():
        raise AssertionError("Expected hidden-truth files were not written")

    raw_graph = build_incidence_graph(fixture)
    raw_cycle_rank = sum(
        raw_graph.summarize_component(component).cycle_rank
        for component in raw_graph.connected_components()
    )
    structural = build_structural_projection(fixture)
    if len(structural.bundles) >= len(fixture.homology_ids):
        raise AssertionError(
            "Simulator anchors were not compressed into chromosome-signature bundles"
        )
    if structural.cycle_rank_total > raw_cycle_rank:
        raise AssertionError("Structural projection increased cycle rank")

    result = exact_optimize_layer_dp(
        fixture,
        orientation_cap_per_component=4096,
        permutation_cap_per_species=720,
        transition_cap_per_component=2000000,
    )
    audit = build_layout_audit(fixture, result.layout)

    if result.layout.optimized_score.crossings > result.layout.initial_score.crossings:
        raise AssertionError("Optimizer made the deliberately tangled layout worse")
    if result.layout.optimality_status != "proven optimum":
        raise AssertionError("Integration case did not retain exact optimality")
    if audit.input_fingerprint != public_fingerprint:
        raise AssertionError("Audit fingerprint differs from public solver input")

    # Hidden truth is inspected only after the solve, for benchmark validation.
    hidden_text = hidden_evolution.read_text(encoding="utf-8")
    if "fission" not in hidden_text or "fusion" not in hidden_text or "inversion" not in hidden_text:
        raise AssertionError("Hidden simulator history is missing a planned event type")

    summary = {
        "fixture_id": fixture.fixture_id,
        "species": list(fixture.species_ids),
        "chromosomes": len(fixture.chromosomes),
        "homology_groups": len(fixture.homology_ids),
        "structural_bundles": len(structural.bundles),
        "raw_cycle_rank": raw_cycle_rank,
        "structural_cycle_rank": structural.cycle_rank_total,
        "crossings_initial": result.layout.initial_score.crossings,
        "crossings_optimized": result.layout.optimized_score.crossings,
        "optimality_status": result.layout.optimality_status,
        "dp_transitions": result.layout.states_evaluated,
        "input_fingerprint": public_fingerprint,
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
