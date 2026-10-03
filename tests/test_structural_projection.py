from pathlib import Path
import unittest

from syntangle import build_incidence_graph, build_structural_projection, fixture_from_dict, load_fixture

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "examples" / "fixtures"


class StructuralProjectionTests(unittest.TestCase):
    def test_stage0_cycle_rank_is_unchanged(self):
        for path in sorted(FIXTURES.glob("*.json")):
            if path.name == "fixture.schema.json":
                continue
            fixture = load_fixture(path)
            if "incidence_cycle_rank_total" in fixture.expected:
                self.assertEqual(
                    build_structural_projection(fixture).cycle_rank_total,
                    fixture.expected["incidence_cycle_rank_total"],
                )

    def test_repeated_anchor_evidence_is_compressed(self):
        species = []
        for sp in ("sp1", "sp2", "sp3"):
            blocks = []
            for i, h in enumerate(("H1", "H2", "H3")):
                blocks.append({
                    "occurrence_id": f"{sp}:{h}",
                    "homology_id": h,
                    "start": i * 10,
                    "end": i * 10 + 5,
                    "strand": "+",
                })
            species.append({
                "id": sp,
                "chromosomes": [{"id": "A", "length": 100, "blocks": blocks}],
            })

        fixture = fixture_from_dict({
            "fixture_version": 1,
            "id": "repeated_evidence",
            "title": "Repeated evidence",
            "purpose": "Test structural compression",
            "species": species,
            "orientation_constraints": [],
            "expected": {},
        })
        raw = build_incidence_graph(fixture)
        raw_rank = sum(
            raw.summarize_component(c).cycle_rank
            for c in raw.connected_components()
        )
        structural = build_structural_projection(fixture)
        self.assertEqual(raw_rank, 4)
        self.assertEqual(structural.cycle_rank_total, 0)
        self.assertEqual(len(structural.bundles), 1)


if __name__ == "__main__":
    unittest.main()
