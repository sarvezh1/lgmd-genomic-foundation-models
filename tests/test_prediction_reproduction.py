import csv
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "reproduction"))

from reproduce_metrics import load_and_validate, reproduce, verify_bundle_hashes, verify_error_intersections


class PredictionScientificRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bundle = verify_bundle_hashes()
        cls.canonical, cls.registry, cls.by_config = load_and_validate()

    def test_canonical_contract(self):
        self.assertEqual(len(self.canonical), 6960)
        self.assertEqual(sum(row["split"] == "train" for row in self.canonical), 4866)
        test = [row for row in self.canonical if row["split"] == "test"]
        self.assertEqual(len(test), 2094)
        self.assertEqual(sum(int(row["binary_label"]) for row in test), 509)
        self.assertEqual(len({row["locus_id"] for row in self.canonical}), 6094)

    def test_configuration_coverage(self):
        self.assertEqual(len(self.registry), 28)
        self.assertEqual(sum(row["primary"] == "true" for row in self.registry), 26)
        self.assertEqual(sum(row["sensitivity"] == "true" for row in self.registry), 2)
        self.assertTrue(all(len(rows) == 2094 for rows in self.by_config.values()))

    def test_metrics_reproduce(self):
        rows, failures = reproduce(self.registry, self.by_config, 1e-12)
        self.assertEqual(failures, [])
        self.assertEqual(len(rows), 28)
        self.assertLessEqual(max(float(row["max_abs_difference"]) for row in rows), 1e-12)

    def test_error_landscape_reproduces(self):
        result, failures = verify_error_intersections(self.registry, self.by_config)
        self.assertEqual(failures, [])
        self.assertEqual((result["all_12_101bp_wrong"], result["all_24_sequence_wrong"], result["all_26_primary_wrong"]), (5, 2, 1))
        self.assertEqual(result["universal_residual_variation_id"], "4538462")

    def test_evidence_manifest_and_source_map(self):
        manifest = json.loads(
            (ROOT / "metadata/evidence_manifest.json").read_text(encoding="utf-8")
        )
        self.assertEqual(manifest["files"]["primary_predictions"]["records"], 54444)
        with (ROOT / "metadata/source_to_result.tsv").open(
            encoding="utf-8", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 452)
        self.assertFalse(any(row["verification_status"] == "FAIL" for row in rows))


if __name__ == "__main__":
    unittest.main()
