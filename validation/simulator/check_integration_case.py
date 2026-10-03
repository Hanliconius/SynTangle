from __future__ import annotations

import json
import sys
from pathlib import Path

from syntangle import (
    build_layout_audit,
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

    hidden_evolution = bundle / "hidden_evolution_log.tsv"
    hidden_tangle = bundle / "hidden_tangle_log.tsv"
    if not hidden_evolution.is_file() or not hidden_tangle.is_file():
        raise AssertionError("Expected hidden-truth files were not written")

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
