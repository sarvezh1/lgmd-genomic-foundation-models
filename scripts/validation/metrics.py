"""Dependency-free binary classification metrics used by repository checks."""

from __future__ import annotations

import math
from collections.abc import Iterable


def _binary(values: Iterable[int], name: str) -> list[int]:
    result = [int(value) for value in values]
    if not result or any(value not in (0, 1) for value in result):
        raise ValueError(f"{name} must be a non-empty binary vector")
    return result


def confusion_counts(y_true: Iterable[int], y_pred: Iterable[int]) -> dict[str, int]:
    truth = _binary(y_true, "y_true")
    pred = _binary(y_pred, "y_pred")
    if len(truth) != len(pred):
        raise ValueError("y_true and y_pred lengths differ")
    return {
        "tn": sum(t == 0 and p == 0 for t, p in zip(truth, pred)),
        "fp": sum(t == 0 and p == 1 for t, p in zip(truth, pred)),
        "fn": sum(t == 1 and p == 0 for t, p in zip(truth, pred)),
        "tp": sum(t == 1 and p == 1 for t, p in zip(truth, pred)),
    }


def classification_metrics(y_true: Iterable[int], y_pred: Iterable[int]) -> dict[str, float | int]:
    c = confusion_counts(y_true, y_pred)
    tn, fp, fn, tp = (c[key] for key in ("tn", "fp", "fn", "tp"))
    total = tn + fp + fn + tp
    recall = tp / (tp + fn) if tp + fn else math.nan
    specificity = tn / (tn + fp) if tn + fp else math.nan
    precision = tp / (tp + fp) if tp + fp else 0.0
    denominator = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = ((tp * tn - fp * fn) / denominator) if denominator else 0.0
    return {
        **c,
        "accuracy": (tp + tn) / total,
        "balanced_accuracy": (recall + specificity) / 2,
        "precision": precision,
        "recall": recall,
        "specificity": specificity,
        "f1": 2 * precision * recall / (precision + recall) if precision + recall else 0.0,
        "mcc": mcc,
    }


def auroc(y_true: Iterable[int], scores: Iterable[float]) -> float:
    truth = _binary(y_true, "y_true")
    values = [float(value) for value in scores]
    if len(truth) != len(values):
        raise ValueError("y_true and score lengths differ")
    positives = sum(truth)
    negatives = len(truth) - positives
    if not positives or not negatives:
        raise ValueError("AUROC requires both classes")
    ordered = sorted(zip(values, truth), key=lambda item: item[0])
    rank_sum = 0.0
    i = 0
    while i < len(ordered):
        j = i + 1
        while j < len(ordered) and ordered[j][0] == ordered[i][0]:
            j += 1
        average_rank = ((i + 1) + j) / 2
        rank_sum += average_rank * sum(label for _, label in ordered[i:j])
        i = j
    return (rank_sum - positives * (positives + 1) / 2) / (positives * negatives)


def average_precision(y_true: Iterable[int], scores: Iterable[float]) -> float:
    truth = _binary(y_true, "y_true")
    values = [float(value) for value in scores]
    if len(truth) != len(values):
        raise ValueError("y_true and score lengths differ")
    positives = sum(truth)
    if not positives:
        raise ValueError("average precision requires a positive example")
    ordered = sorted(zip(values, truth), key=lambda item: item[0], reverse=True)
    tp = fp = 0
    previous_recall = result = 0.0
    i = 0
    while i < len(ordered):
        score = ordered[i][0]
        j = i
        while j < len(ordered) and ordered[j][0] == score:
            if ordered[j][1]:
                tp += 1
            else:
                fp += 1
            j += 1
        recall = tp / positives
        precision = tp / (tp + fp)
        result += (recall - previous_recall) * precision
        previous_recall = recall
        i = j
    return result
