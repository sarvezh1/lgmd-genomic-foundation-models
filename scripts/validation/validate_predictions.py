"""Validate alignment and calculate metrics for a frozen prediction file."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from metrics import auroc, average_precision, classification_metrics


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(8192)
        handle.seek(0)
        delimiter = "\t" if sample.count("\t") >= sample.count(",") else ","
        return list(csv.DictReader(handle, delimiter=delimiter))


def find_column(rows: list[dict[str, str]], aliases: tuple[str, ...]) -> str:
    if not rows:
        raise ValueError("table is empty")
    match = next((name for name in aliases if name in rows[0]), None)
    if match is None:
        raise ValueError(f"missing column; accepted aliases: {aliases}")
    return match


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--test", type=Path, required=True)
    parser.add_argument("--expect-canonical", action="store_true")
    args = parser.parse_args()
    try:
        predictions, test = read_table(args.predictions), read_table(args.test)
        pred_id = find_column(predictions, ("variant_id", "VariationID", "variation_id", "id"))
        test_id = find_column(test, ("variant_id", "VariationID", "variation_id", "id"))
        score_col = find_column(predictions, ("score", "probability", "y_score", "pathogenic_probability"))
        pred_col = find_column(predictions, ("prediction", "pred", "y_pred", "binary_prediction"))
        label_col = find_column(predictions, ("label", "binary_label", "y_true"))
        test_label = find_column(test, ("label", "binary_label", "Label", "y"))
        pred_ids = [row[pred_id].strip() for row in predictions]
        test_ids = [row[test_id].strip() for row in test]
        if len(set(pred_ids)) != len(pred_ids):
            raise ValueError("prediction identifiers are duplicated")
        if pred_ids != test_ids:
            raise ValueError("prediction rows are not in exact test-table order")
        labels = [int(float(row[label_col])) for row in predictions]
        expected_labels = [int(float(row[test_label])) for row in test]
        if labels != expected_labels:
            raise ValueError("prediction labels do not align with test labels")
        predictions_binary = [int(float(row[pred_col])) for row in predictions]
        scores = [float(row[score_col]) for row in predictions]
        if args.expect_canonical and len(predictions) != 2094:
            raise ValueError(f"prediction row count is {len(predictions)}, expected 2,094")
        metrics = classification_metrics(labels, predictions_binary)
        metrics["auroc"] = auroc(labels, scores)
        metrics["average_precision"] = average_precision(labels, scores)
    except (OSError, ValueError) as exc:
        print(f"Prediction validation FAILED: {exc}", file=sys.stderr)
        return 1
    ordered = ("mcc", "auroc", "average_precision", "accuracy", "balanced_accuracy", "precision", "recall", "specificity", "f1")
    print(f"Prediction validation passed for {len(predictions)} aligned rows.")
    print(" ".join(f"{name}={metrics[name]:.12g}" for name in ordered))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
