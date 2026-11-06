"""Verify the canonical evidence bundle and reproduce fixed-test metrics."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / "scripts" / "validation"))
from metrics import auroc, average_precision, classification_metrics  # noqa: E402
from prediction_integrity import validate_prediction_records  # noqa: E402

SEQUENCE_KEYS = {"dnabert2", "nt_v2", "hyenadna", "caduceus_ph", "moderngena", "evo2"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ensure_repo_path(path: Path) -> Path:
    resolved, root = path.resolve(), ROOT.resolve()
    if resolved != root and root not in resolved.parents:
        raise ValueError(f"output path must remain inside repository: {resolved}")
    return resolved


def write_json(path: Path, value: object) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(value, indent=2, sort_keys=True) + "\n")


def expected_metric_rows() -> dict[tuple[str, str, int | None], dict[str, str]]:
    rows = read_tsv(ROOT / "tables" / "supplementary_28_configurations.tsv")
    result = {}
    for row in rows:
        context = int(float(row["context_bp"])) if row["context_bp"] else None
        result[(row["model_key_integrated"], row["regime_integrated"], context)] = row
    if len(result) != 28:
        raise ValueError("archived metric table does not resolve to 28 unique keys")
    return result


def verify_bundle_hashes() -> dict[str, object]:
    bundle = json.loads((ROOT / "metadata" / "evidence_manifest.json").read_text(encoding="utf-8"))
    for name, record in bundle["files"].items():
        path = ROOT / record["path"]
        if not path.is_file() or sha256(path) != record["sha256"]:
            raise ValueError(f"bundle hash failure: {name}")
    return bundle


def load_and_validate() -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, list[dict[str, str]]]]:
    canonical = read_tsv(ROOT / "data" / "canonical_variants.tsv")
    registry = read_tsv(ROOT / "data" / "configurations.tsv")
    predictions = read_tsv(ROOT / "data" / "predictions_primary.tsv") + read_tsv(ROOT / "data" / "predictions_sensitivity.tsv")
    if len(canonical) != 6960 or len({row["VariationID"] for row in canonical}) != 6960:
        raise ValueError("canonical VariationID count failure")
    train = [row for row in canonical if row["split"] == "train"]
    test = [row for row in canonical if row["split"] == "test"]
    if len(train) != 4866 or len(test) != 2094:
        raise ValueError("canonical split count failure")
    if sum(int(row["binary_label"]) for row in test) != 509:
        raise ValueError("canonical test label distribution failure")
    if len({row["locus_id"] for row in canonical}) != 6094 or len({row["locus_id"] for row in test}) != 1834:
        raise ValueError("canonical locus count failure")
    if {row["locus_id"] for row in train} & {row["locus_id"] for row in test}:
        raise ValueError("exact train/test locus overlap")
    if len(registry) != 28 or sum(row["primary"] == "true" for row in registry) != 26 or sum(row["sensitivity"] == "true" for row in registry) != 2:
        raise ValueError("configuration registry primary/sensitivity failure")

    by_config = validate_prediction_records(predictions, registry, test, 2094)
    return canonical, registry, by_config


def reproduce(registry: list[dict[str, str]], by_config: dict[str, list[dict[str, str]]], tolerance: float) -> tuple[list[dict[str, object]], list[str]]:
    expected = expected_metric_rows()
    failures: list[str] = []
    output: list[dict[str, object]] = []
    field_map = {"mcc": "MCC", "auroc": "AUROC", "average_precision": "AUPRC", "accuracy": "Accuracy", "balanced_accuracy": "Balanced_Accuracy", "precision": "Precision", "recall": "Recall", "specificity": "Specificity", "f1": "F1", "tn": "TN", "fp": "FP", "fn": "FN", "tp": "TP"}
    for config in registry:
        config_id = config["configuration_id"]
        rows = sorted(by_config[config_id], key=lambda row: int(row["VariationID"]))
        labels = [int(row["binary_label"]) for row in rows]
        predictions = [int(row["prediction"]) for row in rows]
        scores = [float(row["score"]) for row in rows]
        calculated = classification_metrics(labels, predictions)
        calculated["auroc"] = auroc(labels, scores)
        calculated["average_precision"] = average_precision(labels, scores)
        regime = config["regime"] if config["model_key"] in SEQUENCE_KEYS else "native"
        context = int(config["context_bp"]) if config["context_bp"] else None
        archived = expected.get((config["model_key"], regime, context))
        if archived is None:
            failures.append(f"no archived metric row for {config_id}")
            continue
        differences = []
        compared = 0
        for calc_name, archived_name in field_map.items():
            if archived[archived_name] != "":
                compared += 1
                differences.append(abs(float(calculated[calc_name]) - float(archived[archived_name])))
        maximum = max(differences, default=0.0)
        status = "PASS" if maximum <= tolerance else "FAIL"
        if status == "FAIL":
            failures.append(f"metric mismatch {config_id}: max abs difference {maximum:.17g}")
        output.append({
            "configuration_id": config_id,
            "primary": config["primary"],
            "sensitivity": config["sensitivity"],
            "n": len(rows),
            "MCC": calculated["mcc"],
            "AUROC": calculated["auroc"],
            "AUPRC": calculated["average_precision"],
            "Accuracy": calculated["accuracy"],
            "Balanced_Accuracy": calculated["balanced_accuracy"],
            "Precision": calculated["precision"],
            "Recall": calculated["recall"],
            "Specificity": calculated["specificity"],
            "F1": calculated["f1"],
            "TN": calculated["tn"], "FP": calculated["fp"], "FN": calculated["fn"], "TP": calculated["tp"],
            "archived_MCC": archived["MCC"], "archived_AUROC": archived["AUROC"], "archived_AUPRC": archived["AUPRC"],
            "compared_fields": compared,
            "max_abs_difference": maximum,
            "status": status,
        })
    return output, failures


def verify_error_intersections(registry: list[dict[str, str]], by_config: dict[str, list[dict[str, str]]]) -> tuple[dict[str, object], list[str]]:
    failures: list[str] = []
    registry_by_id = {row["configuration_id"]: row for row in registry}
    predictions = {config_id: {row["VariationID"]: int(row["prediction"]) for row in rows} for config_id, rows in by_config.items()}
    labels = {row["VariationID"]: int(row["binary_label"]) for row in next(iter(by_config.values()))}
    genes = {row["VariationID"]: row["GeneSymbol"] for row in next(iter(by_config.values()))}
    original12 = [row["configuration_id"] for row in registry if row["model_key"] in SEQUENCE_KEYS and row["context_bp"] == "101"]
    sequence24 = [row["configuration_id"] for row in registry if row["model_key"] in SEQUENCE_KEYS]
    context201 = [row["configuration_id"] for row in registry if row["model_key"] in SEQUENCE_KEYS and row["context_bp"] == "201"]
    context501 = [row["configuration_id"] for row in registry if row["model_key"] in SEQUENCE_KEYS and row["context_bp"] == "501"]
    native = [row["configuration_id"] for row in registry if row["primary"] == "true" and row["regime"] == "native"]
    primary26 = [row["configuration_id"] for row in registry if row["primary"] == "true"]
    if not (len(original12) == 12 and len(sequence24) == 24 and len(native) == 2 and len(primary26) == 26):
        raise ValueError("configuration grouping failure during error verification")

    def all_wrong(identifier: str, configs: list[str]) -> bool:
        return all(predictions[config][identifier] != labels[identifier] for config in configs)

    all12_ids = {identifier for identifier in labels if all_wrong(identifier, original12)}
    all24_ids = {identifier for identifier in labels if all_wrong(identifier, sequence24)}
    all26_ids = {identifier for identifier in labels if all_wrong(identifier, primary26)}
    observed = (len(all12_ids), len(all24_ids), len(all26_ids))
    if observed != (5, 2, 1):
        failures.append(f"shared-error intersections differ: {observed}")
    if all26_ids != {"4538462"} or genes.get("4538462") != "ANO5":
        failures.append(f"universal residual differs: {sorted(all26_ids)}")

    recurrent_rows = read_tsv(ROOT / "data" / "recurrent_hard_variants.tsv")
    recurrent_ids = {row["VariationID"] for row in recurrent_rows}
    if len(recurrent_ids) != 787:
        failures.append(f"recurrent hard set has {len(recurrent_ids)} IDs")
    any201 = sum(not all_wrong(identifier, context201) for identifier in recurrent_ids)
    any501 = sum(not all_wrong(identifier, context501) for identifier in recurrent_ids)
    either_native = sum(not all_wrong(identifier, native) for identifier in recurrent_ids)
    later = list(dict.fromkeys(context201 + context501 + native))
    any_later = sum(not all_wrong(identifier, later) for identifier in recurrent_ids)
    recurrent_all24 = sum(all_wrong(identifier, sequence24) for identifier in recurrent_ids)
    recurrent_all26 = sum(all_wrong(identifier, primary26) for identifier in recurrent_ids)
    recurrent_observed = (any201, any501, either_native, any_later, recurrent_all24, recurrent_all26)
    if recurrent_observed != (709, 731, 772, 782, 2, 1):
        failures.append(f"recurrent-hard retrospective coverage differs: {recurrent_observed}")
    result = {
        "all_12_101bp_wrong": len(all12_ids),
        "all_24_sequence_wrong": len(all24_ids),
        "all_26_primary_wrong": len(all26_ids),
        "universal_residual_variation_id": next(iter(all26_ids)) if len(all26_ids) == 1 else None,
        "universal_residual_gene": genes.get(next(iter(all26_ids))) if len(all26_ids) == 1 else None,
        "recurrent_hard_n": len(recurrent_ids),
        "recurrent_rescued_by_any_201": any201,
        "recurrent_rescued_by_any_501": any501,
        "recurrent_rescued_by_either_native": either_native,
        "recurrent_rescued_by_any_later_strategy": any_later,
        "recurrent_wrong_all_24_sequence": recurrent_all24,
        "recurrent_wrong_all_26_primary": recurrent_all26,
        "interpretation": "retrospective oracle coverage only; no ensemble was trained or evaluated",
    }
    return result, failures


def write_reports(report_dir: Path, metric_rows: list[dict[str, object]], summary: dict[str, object]) -> None:
    report_dir = ensure_repo_path(report_dir)
    report_dir.mkdir(parents=True, exist_ok=True)
    metric_path = report_dir / "metric_verification.tsv"
    with metric_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metric_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(metric_rows)
    summary["metric_report"] = metric_path.relative_to(ROOT).as_posix()
    summary["metric_report_sha256"] = sha256(metric_path)
    write_json(report_dir / "verification_summary.json", summary)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tolerance", type=float, default=1e-12)
    parser.add_argument("--report-dir", type=Path, default=ROOT / "verification_outputs" / "prediction_metrics")
    args = parser.parse_args()
    try:
        bundle = verify_bundle_hashes()
        canonical, registry, by_config = load_and_validate()
        metric_rows, metric_failures = reproduce(registry, by_config, args.tolerance)
        intersections, intersection_failures = verify_error_intersections(registry, by_config)
        failures = metric_failures + intersection_failures
        summary = {
            "status": "PASS" if not failures else "FAIL",
            "tolerance": args.tolerance,
            "canonical_variants": len(canonical),
            "primary_configurations": sum(row["primary"] == "true" for row in registry),
            "sensitivity_configurations": sum(row["sensitivity"] == "true" for row in registry),
            "configurations_reproduced": len(metric_rows),
            "configurations_passing": sum(row["status"] == "PASS" for row in metric_rows),
            "maximum_metric_absolute_difference": max(float(row["max_abs_difference"]) for row in metric_rows),
            "error_landscape": intersections,
            "bundle_manifest_version": bundle["version"],
            "failures": failures,
        }
        write_reports(args.report_dir, metric_rows, summary)
    except (OSError, ValueError, KeyError) as exc:
        print(f"Prediction reproduction FAILED: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
