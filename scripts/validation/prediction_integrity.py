"""Reusable integrity checks for normalized variant-level prediction records."""

from __future__ import annotations

import math
from collections import defaultdict


def validate_prediction_records(
    predictions: list[dict[str, str]],
    registry: list[dict[str, str]],
    test_rows: list[dict[str, str]],
    expected_per_configuration: int,
) -> dict[str, list[dict[str, str]]]:
    registry_by_id = {row["configuration_id"]: row for row in registry}
    if len(registry_by_id) != len(registry):
        raise ValueError("duplicate configuration identifiers")
    for config in registry:
        primary = config["primary"] == "true"
        sensitivity = config["sensitivity"] == "true"
        if primary == sensitivity:
            raise ValueError(f"configuration must be exactly primary or sensitivity: {config['configuration_id']}")
    test_by_id = {row["VariationID"]: row for row in test_rows}
    if len(test_by_id) != len(test_rows):
        raise ValueError("duplicate canonical test VariationIDs")
    by_config: dict[str, list[dict[str, str]]] = defaultdict(list)
    seen = set()
    for row in predictions:
        key = (row["configuration_id"], row["VariationID"])
        if key in seen:
            raise ValueError(f"duplicate prediction key: {key}")
        seen.add(key)
        config = registry_by_id.get(row["configuration_id"])
        if config is None:
            raise ValueError(f"unregistered configuration: {row['configuration_id']}")
        canonical = test_by_id.get(row["VariationID"])
        if canonical is None or row["GeneSymbol"] != canonical["GeneSymbol"] or row["binary_label"] != canonical["binary_label"]:
            raise ValueError(f"prediction/canonical alignment failure: {key}")
        score, threshold, prediction = float(row["score"]), float(row["threshold"]), int(row["prediction"])
        if not math.isfinite(score):
            raise ValueError(f"non-finite score: {key}")
        if prediction != int(score >= threshold):
            raise ValueError(f"threshold inconsistency: {key}")
        if config["probability_score"] == "true" and not 0.0 <= score <= 1.0:
            raise ValueError(f"probability range failure: {key}")
        if row["score_direction"] != config["score_direction"] or row["score_direction"] != "higher_is_more_pathogenic_like":
            raise ValueError(f"score direction failure: {key}")
        by_config[row["configuration_id"]].append(row)
    if set(by_config) != set(registry_by_id):
        raise ValueError("prediction configuration coverage failure")
    expected_ids = set(test_by_id)
    for config_id, rows in by_config.items():
        if len(rows) != expected_per_configuration or {row["VariationID"] for row in rows} != expected_ids:
            raise ValueError(f"prediction VariationID coverage failure: {config_id}")
    return dict(by_config)
