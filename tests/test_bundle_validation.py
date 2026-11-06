import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "validation"))

from validate_canonical_bundle import proximity


class BundleValidationTests(unittest.TestCase):
    def test_nearest_locus_proximity(self):
        train = {("1", 100), ("1", 200), ("2", 50)}
        test = {("1", 102), ("1", 160), ("2", 80)}
        result = proximity(train, test)
        self.assertEqual(result["unique_test_loci"], 3)
        self.assertEqual(result["median_bp"], 30)
        self.assertEqual(result["within_50_bp"], 3)

    def test_missing_training_chromosome_fails(self):
        with self.assertRaises(ValueError):
            proximity({("1", 100)}, {("2", 100)})


if __name__ == "__main__":
    unittest.main()
