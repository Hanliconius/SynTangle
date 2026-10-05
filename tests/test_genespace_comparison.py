import importlib.util
import tempfile
import unittest
from pathlib import Path

from syntangle.fixtures import fixture_from_dict
from syntangle.layout import initial_layout_state, score_crossings

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "genespace_comparison", ROOT / "validation/benchmark/compare_genespace.py"
)
comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparison)


def fixture():
    species = []
    for name, sign in (("A", 1), ("B", -1)):
        species.append(dict(id=name, chromosomes=[dict(
            id="chr1", length=100, display_rank=1, display_orientation=sign,
            blocks=[dict(occurrence_id=f"{name}:{i}", homology_id=f"h{i}",
                         start=start, end=start + 10, strand="+")
                    for i, start in enumerate((10, 60))]
        )]))
    return fixture_from_dict(dict(fixture_version=1, id="flipped",
                                  title="flipped", purpose="comparison test",
                                  species=species))


class GenespaceComparisonTests(unittest.TestCase):
    def test_flip_assistance_preserves_order_and_fixes_reversed_pair(self):
        data = fixture()
        initial = initial_layout_state(data)
        self.assertEqual(score_crossings(data, initial).crossings, 1)
        result, evaluations = comparison.improve_flips(data, initial, 1)
        self.assertEqual(score_crossings(data, result).crossings, 0)
        self.assertEqual(result.chromosome_order, initial.chromosome_order)
        self.assertGreater(evaluations, 0)
        comparison.validate_state(data, result)

    def test_native_adapter_reflects_public_presentation_only(self):
        data = fixture()
        with tempfile.TemporaryDirectory() as directory:
            comparison.export_native_input(data, Path(directory))
            rows = comparison.read_tsv(Path(directory) / "bed.tsv")
            coords = {(r["genome"], r["og"]): float(r["ord"]) for r in rows}
            self.assertEqual(coords[("A", "h0")], 15)
            self.assertEqual(coords[("B", "h0")], 85)
            self.assertEqual(coords[("B", "h1")], 35)
            # Original biological coordinates must remain unchanged.
            self.assertEqual(data.chromosomes[1].blocks[0].start, 10)

    def test_native_output_must_retain_every_chromosome(self):
        data = fixture()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "orders.tsv"
            comparison.write_tsv(path, [dict(variant="A_w1", genome="A",
                chr="chr1", plotOrd=1, ordering_seconds=0.1)])
            with self.assertRaises(ValueError):
                list(comparison.load_native_states(data, path))


if __name__ == "__main__":
    unittest.main()
