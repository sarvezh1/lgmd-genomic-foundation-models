import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "validation"))

from validate_repository import validate


class RepositoryTests(unittest.TestCase):
    def test_repository_integrity_and_hygiene(self):
        errors, checks = validate(ROOT)
        self.assertEqual(errors, [])
        self.assertEqual(checks["scientific_counts"]["variants"], 6960)
        self.assertEqual(checks["scientific_counts"]["fixed_predictions"], 58632)
        self.assertEqual(checks["source_to_result"], 452)
        self.assertEqual(checks["approved_assets"], 14)
        self.assertEqual(checks["markdown_files"], 4)


if __name__ == "__main__":
    unittest.main()
