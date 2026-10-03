from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from syntangle import fixture_fingerprint, load_validation_bundle


class ValidationBundleTests(unittest.TestCase):
    def test_bundle_import_ignores_hidden_truth(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "species.tsv").write_text(
                "species_id\tspecies_rank\nsp1\t1\nsp2\t2\n",
                encoding="utf-8",
            )
            (root / "chromosomes.tsv").write_text(
                "species_id\tchromosome_id\tlength\tdisplay_rank\n"
                "sp1\tA\t100\t1\n"
                "sp2\tB\t100\t1\n",
                encoding="utf-8",
            )
            (root / "occurrences.tsv").write_text(
                "occurrence_id\thomology_id\tspecies_id\tchromosome_id\tstart\tend\tstrand\n"
                "sp1:g1\tg1\tsp1\tA\t0\t50\t+\n"
                "sp2:g1\tg1\tsp2\tB\t0\t50\t+\n",
                encoding="utf-8",
            )
            (root / "metadata.tsv").write_text(
                "key\tvalue\nfixture_id\tbundle_test\n",
                encoding="utf-8",
            )
            (root / "input_display_state.tsv").write_text(
                "species\tchrom\tsource_rank\tdisplay_rank\torientation\n"
                "sp1\tA\t1\t1\t-1\n"
                "sp2\tB\t1\t1\t1\n",
                encoding="utf-8",
            )
            hidden = root / "hidden_evolution_log.tsv"
            hidden.write_text(
                "event_id\ttype\n1\tfusion\n",
                encoding="utf-8",
            )

            fixture1 = load_validation_bundle(root)
            fingerprint1 = fixture_fingerprint(fixture1)

            hidden.write_text(
                "event_id\ttype\n1\tfission\n2\tinversion\n",
                encoding="utf-8",
            )
            fixture2 = load_validation_bundle(root)
            fingerprint2 = fixture_fingerprint(fixture2)

            self.assertEqual(fixture1.fixture_id, "bundle_test")
            self.assertEqual(fixture1.species_ids, ("sp1", "sp2"))
            orientations = {
                chromosome.ref.label: chromosome.display_orientation
                for chromosome in fixture1.chromosomes
            }
            self.assertEqual(orientations["sp1:A"], -1)
            self.assertEqual(orientations["sp2:B"], 1)
            self.assertEqual(fingerprint1, fingerprint2)

    def test_species_rank_controls_layer_order(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "species.tsv").write_text(
                "species_id\tspecies_rank\nlate\t2\nearly\t1\n",
                encoding="utf-8",
            )
            (root / "chromosomes.tsv").write_text(
                "species_id\tchromosome_id\tlength\tdisplay_rank\n"
                "early\tA\t100\t1\n"
                "late\tB\t100\t1\n",
                encoding="utf-8",
            )
            (root / "occurrences.tsv").write_text(
                "occurrence_id\thomology_id\tspecies_id\tchromosome_id\tstart\tend\tstrand\n"
                "early:g\tg\tearly\tA\t0\t50\t+\n"
                "late:g\tg\tlate\tB\t0\t50\t+\n",
                encoding="utf-8",
            )

            fixture = load_validation_bundle(root)
            self.assertEqual(fixture.species_ids, ("early", "late"))


if __name__ == "__main__":
    unittest.main()
