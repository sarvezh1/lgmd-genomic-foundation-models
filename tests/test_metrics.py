import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "validation"))

from metrics import auroc, average_precision, classification_metrics


class MetricTests(unittest.TestCase):
    def test_perfect_metrics(self):
        metrics = classification_metrics([0, 0, 1, 1], [0, 0, 1, 1])
        self.assertTrue(math.isclose(metrics["mcc"], 1.0))
        self.assertTrue(math.isclose(auroc([0, 0, 1, 1], [0.1, 0.2, 0.8, 0.9]), 1.0))
        self.assertTrue(math.isclose(average_precision([0, 0, 1, 1], [0.1, 0.2, 0.8, 0.9]), 1.0))

    def test_tied_auc(self):
        self.assertTrue(math.isclose(auroc([0, 1], [0.5, 0.5]), 0.5))

    def test_known_confusion_counts(self):
        metrics = classification_metrics([0, 0, 1, 1], [0, 1, 0, 1])
        self.assertEqual((metrics["tn"], metrics["fp"], metrics["fn"], metrics["tp"]), (1, 1, 1, 1))
        self.assertTrue(math.isclose(metrics["mcc"], 0.0))


if __name__ == "__main__":
    unittest.main()
