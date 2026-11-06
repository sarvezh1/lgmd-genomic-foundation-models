import copy
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "validation"))

from prediction_integrity import validate_prediction_records


def fixtures():
    test = [
        {"VariationID": "1", "GeneSymbol": "GENE1", "binary_label": "0"},
        {"VariationID": "2", "GeneSymbol": "GENE2", "binary_label": "1"},
    ]
    registry = [
        {"configuration_id": "primary", "primary": "true", "sensitivity": "false", "probability_score": "true", "score_direction": "higher_is_more_pathogenic_like"},
        {"configuration_id": "sensitivity", "primary": "false", "sensitivity": "true", "probability_score": "false", "score_direction": "higher_is_more_pathogenic_like"},
    ]
    predictions = []
    for config, scores in (("primary", ("0.2", "0.8")), ("sensitivity", ("-1.0", "2.0"))):
        for row, score in zip(test, scores):
            predictions.append({**row, "configuration_id": config, "score": score, "threshold": "0.5", "prediction": str(int(float(score) >= 0.5)), "score_direction": "higher_is_more_pathogenic_like"})
    return test, registry, predictions


class PredictionIntegrityTests(unittest.TestCase):
    def test_valid_fixture(self):
        test, registry, predictions = fixtures()
        self.assertEqual(set(validate_prediction_records(predictions, registry, test, 2)), {"primary", "sensitivity"})

    def test_duplicate_identifier_fails(self):
        test, registry, predictions = fixtures()
        with self.assertRaisesRegex(ValueError, "duplicate prediction"):
            validate_prediction_records(predictions + [copy.deepcopy(predictions[0])], registry, test, 2)

    def test_missing_prediction_fails(self):
        test, registry, predictions = fixtures()
        with self.assertRaisesRegex(ValueError, "coverage failure"):
            validate_prediction_records(predictions[:-1], registry, test, 2)

    def test_threshold_inconsistency_fails(self):
        test, registry, predictions = fixtures()
        predictions[0]["prediction"] = "1"
        with self.assertRaisesRegex(ValueError, "threshold inconsistency"):
            validate_prediction_records(predictions, registry, test, 2)

    def test_score_direction_fails(self):
        test, registry, predictions = fixtures()
        predictions[0]["score_direction"] = "lower_is_more_pathogenic_like"
        with self.assertRaisesRegex(ValueError, "score direction"):
            validate_prediction_records(predictions, registry, test, 2)

    def test_primary_sensitivity_separation_fails(self):
        test, registry, predictions = fixtures()
        registry[0]["sensitivity"] = "true"
        with self.assertRaisesRegex(ValueError, "exactly primary or sensitivity"):
            validate_prediction_records(predictions, registry, test, 2)


if __name__ == "__main__":
    unittest.main()
