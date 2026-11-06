"""Read-only integrity and proximity checks for a supplied canonical bundle.

This validator checks the canonical data contract and split invariants.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import statistics
import sys
from pathlib import Path

ALIASES = {
    "id": ("variant_id", "VariationID", "variation_id", "id"),
    "chrom": ("chrom", "CHROM", "Chromosome", "chr"),
    "pos": ("pos", "POS", "Position", "position"),
    "ref": ("ref", "REF", "ReferenceAllele"),
    "alt": ("alt", "ALT", "AlternateAllele"),
    "gene": ("gene", "Gene", "GeneSymbol", "SYMBOL"),
    "label": ("label", "binary_label", "Label", "y"),
}


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(8192)
        handle.seek(0)
        delimiter = "\t" if sample.count("\t") >= sample.count(",") else ","
        return list(csv.DictReader(handle, delimiter=delimiter))


def columns(rows: list[dict[str, str]]) -> dict[str, str]:
    if not rows:
        raise ValueError("table is empty")
    names = set(rows[0])
    resolved = {}
    for field, candidates in ALIASES.items():
        match = next((candidate for candidate in candidates if candidate in names), None)
        if match is None:
            raise ValueError(f"missing {field} column; accepted aliases: {candidates}")
        resolved[field] = match
    return resolved


def keys(rows: list[dict[str, str]], cols: dict[str, str]) -> tuple[list[str], set[tuple[str, int]]]:
    ids, loci = [], set()
    for index, row in enumerate(rows, 2):
        identifier = row[cols["id"]].strip()
        chrom = row[cols["chrom"]].strip().removeprefix("chr")
        try:
            pos = int(row[cols["pos"]])
            label = int(float(row[cols["label"]]))
        except ValueError as exc:
            raise ValueError(f"invalid numeric value at row {index}") from exc
        ref, alt = row[cols["ref"]].upper(), row[cols["alt"]].upper()
        if not identifier or pos < 1 or label not in (0, 1):
            raise ValueError(f"invalid identifier, position, or label at row {index}")
        if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT" or ref == alt:
            raise ValueError(f"row {index} is not a canonical SNV")
        ids.append(identifier)
        loci.add((chrom, pos))
    return ids, loci


def proximity(train_loci: set[tuple[str, int]], test_loci: set[tuple[str, int]]) -> dict[str, int | float]:
    by_chrom: dict[str, list[int]] = {}
    for chrom, pos in train_loci:
        by_chrom.setdefault(chrom, []).append(pos)
    for positions in by_chrom.values():
        positions.sort()
    distances = []
    for chrom, pos in sorted(test_loci):
        positions = by_chrom.get(chrom, [])
        if not positions:
            continue
        index = bisect.bisect_left(positions, pos)
        candidates = positions[max(0, index - 1): index + 1]
        distances.append(min(abs(pos - candidate) for candidate in candidates))
    if len(distances) != len(test_loci):
        raise ValueError("at least one test chromosome has no training locus")
    return {
        "unique_test_loci": len(test_loci),
        "median_bp": statistics.median(distances),
        **{f"within_{limit}_bp": sum(distance <= limit for distance in distances) for limit in (50, 100, 250, 500)},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", type=Path, required=True)
    parser.add_argument("--train", type=Path)
    parser.add_argument("--test", type=Path)
    parser.add_argument("--expect-canonical", action="store_true")
    args = parser.parse_args()
    try:
        variants = read_table(args.variants)
        if bool(args.train) != bool(args.test):
            raise ValueError("--train and --test must be supplied together")
        if args.train:
            train, test = read_table(args.train), read_table(args.test)
        else:
            if not variants or "split" not in variants[0]:
                raise ValueError("variant table has no split column; supply --train and --test")
            train = [row for row in variants if row["split"] == "train"]
            test = [row for row in variants if row["split"] == "test"]
        variant_cols, train_cols, test_cols = (columns(rows) for rows in (variants, train, test))
        variant_ids, _ = keys(variants, variant_cols)
        train_ids, train_loci = keys(train, train_cols)
        test_ids, test_loci = keys(test, test_cols)
        if len(set(variant_ids)) != len(variant_ids) or len(set(train_ids)) != len(train_ids) or len(set(test_ids)) != len(test_ids):
            raise ValueError("duplicate variant identifiers detected")
        if set(train_ids) & set(test_ids):
            raise ValueError("train/test variant identifiers overlap")
        if train_loci & test_loci:
            raise ValueError("train/test exact loci overlap")
        if set(train_ids) | set(test_ids) != set(variant_ids):
            raise ValueError("train/test identifiers do not exactly cover the variant table")
        stats = proximity(train_loci, test_loci)
        if args.expect_canonical:
            expected_counts = (6960, 4866, 2094)
            if (len(variants), len(train), len(test)) != expected_counts:
                raise ValueError(f"counts differ from canonical {expected_counts}")
            expected = {"unique_test_loci": 1834, "median_bp": 3, "within_50_bp": 1715, "within_100_bp": 1770, "within_250_bp": 1820, "within_500_bp": 1827}
            if stats != expected:
                raise ValueError(f"proximity statistics differ from audited canonical values: {stats}")
    except (OSError, ValueError) as exc:
        print(f"Canonical bundle validation FAILED: {exc}", file=sys.stderr)
        return 1
    print(f"Canonical bundle validation passed. Counts={len(variants)}/{len(train)}/{len(test)}; proximity={stats}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
