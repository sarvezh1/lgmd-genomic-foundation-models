"""Build the LGMD canonical train/test split.

The study procedure operates on the corrected 6,960-row variant table in
lexicographic VariationID order. It follows the scikit-learn 1.6.1
StratifiedGroupKFold behavior with 10 folds,
``shuffle=True``, and ``random_state=42``.  Every three-fold test combination
was evaluated.  Candidates in the inclusive 0.30 +/- 0.001 test-fraction
window were ranked by class-rate deviation, mean absolute per-gene fraction
deviation, maximum per-gene deviation, test-size deviation, and fold tuple.

Later scikit-learn releases changed shuffled-group bookkeeping. This module
implements the study's 1.6.1 fold-assignment semantics directly and records
the installed runtime versions. Existing assignments, when supplied or
embedded in the canonical table, are read only after candidate selection.

Outputs are generated train/test tables, a per-variant assignment, all 120
candidate scores, and a JSON verification report in a caller-selected output
directory.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd

EXPECTED_COLUMNS = (
    "VariationID", "GeneSymbol", "binary_label", "Chromosome", "Position",
    "REF", "ALT",
)
OPTIONAL_COMPARISON_COLUMNS = ("split", "fold_id")
TARGET_TEST_FRACTION = 0.30
TEST_FRACTION_TOLERANCE = 0.001
N_SPLITS = 10
N_TEST_FOLDS = 3
RANDOM_STATE = 42


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_benchmark(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(
        path,
        sep="\t",
        dtype={
            "VariationID": str,
            "GeneSymbol": str,
            "binary_label": int,
            "Chromosome": str,
            "Position": str,
            "REF": str,
            "ALT": str,
            "split": str,
            "fold_id": "Int64",
        },
    )
    missing = sorted(set(EXPECTED_COLUMNS) - set(frame.columns))
    if missing:
        raise ValueError(f"benchmark is missing columns: {missing}")

    retained = list(EXPECTED_COLUMNS) + [
        column for column in OPTIONAL_COMPARISON_COLUMNS if column in frame.columns
    ]
    frame = frame.loc[:, retained].copy()
    frame["VariationID"] = frame["VariationID"].str.strip()
    frame["GeneSymbol"] = frame["GeneSymbol"].str.strip()
    frame["Chromosome"] = (
        frame["Chromosome"].str.strip().str.upper().str.replace("^CHR", "", regex=True)
    )
    frame["Position"] = pd.to_numeric(frame["Position"], errors="raise").astype(int)
    frame["binary_label"] = pd.to_numeric(
        frame["binary_label"], errors="raise"
    ).astype(int)
    frame["REF"] = frame["REF"].str.strip().str.upper()
    frame["ALT"] = frame["ALT"].str.strip().str.upper()

    if len(frame) != 6960:
        raise ValueError(f"expected 6,960 records, found {len(frame):,}")
    if frame["VariationID"].duplicated().any():
        raise ValueError("VariationID values are not unique")
    if frame["VariationID"].tolist() != sorted(frame["VariationID"].tolist()):
        raise ValueError("input is not lexicographically sorted by string VariationID")
    if not frame["binary_label"].isin((0, 1)).all():
        raise ValueError("binary_label must contain only 0 and 1")
    if (frame["REF"] == frame["ALT"]).any():
        raise ValueError("corrected benchmark still contains REF=ALT records")

    frame["locus_key"] = frame["Chromosome"] + ":" + frame["Position"].astype(str)
    return frame


def _find_best_fold(
    counts_per_fold: np.ndarray,
    class_counts: np.ndarray,
    group_counts: np.ndarray,
) -> int:
    best_fold = -1
    minimum_evaluation = np.inf
    minimum_samples = np.inf
    for fold in range(N_SPLITS):
        counts_per_fold[fold] += group_counts
        evaluation = float(
            np.mean(np.std(counts_per_fold / class_counts.reshape(1, -1), axis=0))
        )
        counts_per_fold[fold] -= group_counts
        samples = float(np.sum(counts_per_fold[fold]))
        if evaluation < minimum_evaluation or (
            np.isclose(evaluation, minimum_evaluation)
            and samples < minimum_samples
        ):
            minimum_evaluation = evaluation
            minimum_samples = samples
            best_fold = fold
    return best_fold


def sklearn_1_6_1_compatible_stratified_group_folds(
    labels: np.ndarray, groups: np.ndarray
) -> np.ndarray:
    """Return folds with the shuffled-group semantics used by sklearn 1.6.1."""
    _, label_inverse, class_counts = np.unique(
        labels, return_inverse=True, return_counts=True
    )
    _, group_inverse, group_counts = np.unique(
        groups, return_inverse=True, return_counts=True
    )
    counts_per_group = np.zeros((len(group_counts), len(class_counts)))
    for class_index, group_index in zip(label_inverse, group_inverse):
        counts_per_group[group_index, class_index] += 1

    # This preserves the material scikit-learn 1.6.1 behavior. The count matrix was
    # shuffled without remapping group_inverse.  Later releases changed this
    # bookkeeping, producing a different deterministic split for the same seed.
    np.random.RandomState(RANDOM_STATE).shuffle(counts_per_group)
    order = np.argsort(-np.std(counts_per_group, axis=1), kind="mergesort")

    counts_per_fold = np.zeros((N_SPLITS, len(class_counts)))
    groups_per_fold = [set() for _ in range(N_SPLITS)]
    for group_index in order:
        fold = _find_best_fold(
            counts_per_fold, class_counts, counts_per_group[group_index]
        )
        counts_per_fold[fold] += counts_per_group[group_index]
        groups_per_fold[fold].add(int(group_index))

    folds = np.full(len(labels), -1, dtype=int)
    for fold, fold_groups in enumerate(groups_per_fold):
        folds[np.isin(group_inverse, list(fold_groups))] = fold
    if (folds < 0).any():
        raise RuntimeError("at least one record was not assigned to a fold")
    return folds


def evaluate_candidates(frame: pd.DataFrame, folds: np.ndarray) -> pd.DataFrame:
    labels = frame["binary_label"].to_numpy()
    global_rate = float(labels.mean())
    genes = sorted(frame["GeneSymbol"].unique())
    gene_counts = frame["GeneSymbol"].value_counts()
    records = []
    for test_folds in itertools.combinations(range(N_SPLITS), N_TEST_FOLDS):
        mask = np.isin(folds, test_folds)
        test = frame.loc[mask]
        test_fraction = len(test) / len(frame)
        gene_fractions = np.asarray(
            [
                int((test["GeneSymbol"] == gene).sum()) / int(gene_counts[gene])
                for gene in genes
            ]
        )
        records.append(
            {
                "test_folds": ",".join(str(value) for value in test_folds),
                "train_n": int(len(frame) - len(test)),
                "test_n": int(len(test)),
                "train_fraction": float(1 - test_fraction),
                "test_fraction": float(test_fraction),
                "test_positive_rate": float(test["binary_label"].mean()),
                "size_deviation": float(abs(test_fraction - TARGET_TEST_FRACTION)),
                "class_deviation": float(
                    abs(float(test["binary_label"].mean()) - global_rate)
                ),
                "gene_mae": float(
                    np.mean(np.abs(gene_fractions - TARGET_TEST_FRACTION))
                ),
                "gene_max_deviation": float(
                    np.max(np.abs(gene_fractions - TARGET_TEST_FRACTION))
                ),
            }
        )
    return pd.DataFrame(records)


def read_partition(path: Path) -> dict[str, str]:
    frame = pd.read_csv(path, sep="\t", dtype={"VariationID": str})
    return {value.strip(): path.name for value in frame["VariationID"]}


def build(args: argparse.Namespace) -> dict[str, object]:
    frame = load_benchmark(args.variants)
    folds = sklearn_1_6_1_compatible_stratified_group_folds(
        frame["binary_label"].to_numpy(), frame["locus_key"].to_numpy()
    )
    scores = evaluate_candidates(frame, folds)
    eligible = scores.loc[
        (scores["test_fraction"] >= TARGET_TEST_FRACTION - TEST_FRACTION_TOLERANCE)
        & (scores["test_fraction"] <= TARGET_TEST_FRACTION + TEST_FRACTION_TOLERANCE)
    ].copy()
    eligible = eligible.sort_values(
        [
            "class_deviation",
            "gene_mae",
            "gene_max_deviation",
            "size_deviation",
            "test_folds",
        ]
    ).reset_index(drop=True)
    if eligible.empty:
        raise RuntimeError("no candidate is inside the inclusive test-size window")
    selected = eligible.iloc[0]
    selected_folds = tuple(int(value) for value in selected["test_folds"].split(","))
    test_mask = np.isin(folds, selected_folds)

    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    train_path = output / "train_variants_generated.tsv"
    test_path = output / "test_variants_generated.tsv"
    assignment_path = output / "locus_disjoint_split_assignment_generated.tsv"
    scores_path = output / "canonical_split_candidate_scores.tsv"
    report_path = output / "canonical_split_verification.json"

    original = list(EXPECTED_COLUMNS)
    frame.loc[~test_mask, original].sort_values("VariationID").to_csv(
        train_path, sep="\t", index=False, lineterminator="\n"
    )
    frame.loc[test_mask, original].sort_values("VariationID").to_csv(
        test_path, sep="\t", index=False, lineterminator="\n"
    )
    assignment = frame[
        ["VariationID", "GeneSymbol", "binary_label", "Chromosome", "Position", "locus_key"]
    ].copy()
    assignment["fold_id"] = folds
    assignment["split"] = np.where(test_mask, "test", "train")
    assignment.sort_values(
        ["split", "GeneSymbol", "Chromosome", "Position", "VariationID"]
    ).to_csv(assignment_path, sep="\t", index=False, lineterminator="\n")
    scores.to_csv(scores_path, sep="\t", index=False, lineterminator="\n")

    generated = dict(zip(frame["VariationID"], np.where(test_mask, "test", "train")))
    archived = None
    archived_fold_matches = None
    if args.archived_assignment:
        archived_frame = pd.read_csv(
            args.archived_assignment,
            sep="\t",
            dtype={"VariationID": str, "fold_id": int, "split": str},
        )
        archived = dict(zip(archived_frame["VariationID"], archived_frame["split"]))
        archived_folds = dict(zip(archived_frame["VariationID"], archived_frame["fold_id"]))
        archived_fold_matches = sum(
            archived_folds.get(variation_id) == int(fold)
            for variation_id, fold in zip(frame["VariationID"], folds)
        )
    elif {"split", "fold_id"}.issubset(frame.columns):
        archived = dict(zip(frame["VariationID"], frame["split"]))
        archived_fold_matches = sum(
            int(expected) == int(observed)
            for expected, observed in zip(frame["fold_id"], folds)
        )

    generated_train = {key for key, value in generated.items() if value == "train"}
    generated_test = {key for key, value in generated.items() if value == "test"}
    archived_train = set()
    archived_test = set()
    if args.archived_train:
        archived_train = set(read_partition(args.archived_train))
    if args.archived_test:
        archived_test = set(read_partition(args.archived_test))
    if "split" in frame.columns:
        if not archived_train:
            archived_train = set(frame.loc[frame["split"] == "train", "VariationID"])
        if not archived_test:
            archived_test = set(frame.loc[frame["split"] == "test", "VariationID"])

    assignment_mismatches = None
    if archived is not None:
        assignment_mismatches = sum(
            archived.get(key) != value for key, value in generated.items()
        ) + len(set(archived) - set(generated))

    try:
        import sklearn
        runtime_sklearn = sklearn.__version__
    except ImportError:
        runtime_sklearn = None

    report = {
        "implementation_status": "portable study split-generation implementation",
        "compatibility_target": "scikit-learn 1.6.1 StratifiedGroupKFold behavior",
        "authoritative_assignment_unchanged": True,
        "input": {
            "path": str(args.variants),
            "sha256": sha256(args.variants),
            "records": int(len(frame)),
            "unique_variation_ids": int(frame["VariationID"].nunique()),
            "unique_loci": int(frame["locus_key"].nunique()),
            "lexicographic_variation_id_order": True,
            "binary_label_values": sorted(int(value) for value in frame["binary_label"].unique()),
        },
        "parameters": {
            "method": "StratifiedGroupKFold",
            "group": "Chromosome + Position",
            "n_splits": N_SPLITS,
            "shuffle": True,
            "random_state": RANDOM_STATE,
            "combined_test_folds": N_TEST_FOLDS,
            "target_test_fraction": TARGET_TEST_FRACTION,
            "inclusive_tolerance": TEST_FRACTION_TOLERANCE,
            "selection_priority": [
                "class_deviation",
                "gene_mae",
                "gene_max_deviation",
                "size_deviation",
                "test_folds",
            ],
        },
        "fold_sizes": [int((folds == value).sum()) for value in range(N_SPLITS)],
        "candidate_combinations": int(len(scores)),
        "eligible_combinations": int(len(eligible)),
        "selected": {
            **{
                key: (int(value) if key in {"train_n", "test_n"} else float(value))
                for key, value in selected.items()
                if key != "test_folds"
            },
            "test_folds": list(selected_folds),
        },
        "generated": {
            "train_variants": len(generated_train),
            "test_variants": len(generated_test),
            "train_loci": int(frame.loc[~test_mask, "locus_key"].nunique()),
            "test_loci": int(frame.loc[test_mask, "locus_key"].nunique()),
            "variation_id_overlap": len(generated_train & generated_test),
            "locus_overlap": len(
                set(frame.loc[~test_mask, "locus_key"])
                & set(frame.loc[test_mask, "locus_key"])
            ),
        },
        "comparison": {
            "train_identifier_mismatches": (
                len(generated_train ^ archived_train) if archived_train else None
            ),
            "test_identifier_mismatches": (
                len(generated_test ^ archived_test) if archived_test else None
            ),
            "assignment_mismatches": assignment_mismatches,
            "fold_id_matches": archived_fold_matches,
            "train_byte_identical": (
                sha256(train_path) == sha256(args.archived_train)
                if args.archived_train else None
            ),
            "test_byte_identical": (
                sha256(test_path) == sha256(args.archived_test)
                if args.archived_test else None
            ),
            "assignment_byte_identical": (
                sha256(assignment_path) == sha256(args.archived_assignment)
                if args.archived_assignment else None
            ),
        },
        "output_sha256": {
            train_path.name: sha256(train_path),
            test_path.name: sha256(test_path),
            assignment_path.name: sha256(assignment_path),
            scores_path.name: sha256(scores_path),
        },
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "installed_scikit_learn": runtime_sklearn,
            "complete_original_environment_equivalence_claimed": False,
        },
    }
    report["exact_assignment_reproduction"] = all(
        value == 0
        for value in (
            report["comparison"]["train_identifier_mismatches"],
            report["comparison"]["test_identifier_mismatches"],
            report["comparison"]["assignment_mismatches"],
        )
        if value is not None
    ) and any(
        value is not None
        for value in (
            report["comparison"]["train_identifier_mismatches"],
            report["comparison"]["test_identifier_mismatches"],
            report["comparison"]["assignment_mismatches"],
        )
    )
    with report_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--archived-train", type=Path)
    parser.add_argument("--archived-test", type=Path)
    parser.add_argument("--archived-assignment", type=Path)
    args = parser.parse_args()
    report = build(args)
    return 0 if report["exact_assignment_reproduction"] or not any(
        (args.archived_train, args.archived_test, args.archived_assignment)
    ) else 1


if __name__ == "__main__":
    raise SystemExit(main())
