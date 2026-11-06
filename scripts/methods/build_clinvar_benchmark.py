"""Build the LGMD SNV benchmark from a local ClinVar variant summary.

This portable implementation applies the study filters in order: GRCh38,
the ten-gene LGMD panel, six accepted aggregate clinical-significance labels,
and ClinVar ``single nucleotide variant`` type.  Pathogenic classes map to 1
and benign classes to 0. The two identified REF=ALT records are
then excluded, and records are sorted lexicographically by string VariationID.

The input is supplied explicitly and is never modified.  Processing is
chunked so the complete raw table is not duplicated in memory or output.
Generated files are written only below ``--output-dir``.  Complete numerical
reproduction is claimed in the report only when all study stage counts
and the optional canonical comparison agree exactly.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import pandas as pd

GENES = {"CAPN3", "DYSF", "FKRP", "SGCA", "SGCB", "SGCG", "ANO5", "DNAJB6", "TCAP", "LMNA"}
PATHOGENIC = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
BENIGN = {"Benign", "Likely benign", "Benign/Likely benign"}
ACCEPTED = PATHOGENIC | BENIGN
EXCLUDED_REF_ALT_IDS = {"177801", "255604"}
EXPECTED_COUNTS = {
    "raw_clinvar": 8_981_173,
    "grch38": 4_456_834,
    "ten_gene_panel": 14_931,
    "accepted_clinical_significance": 8_482,
    "snvs": 6_962,
    "corrected_canonical": 6_960,
}
USE_COLUMNS = [
    "VariationID", "GeneSymbol", "Assembly", "ClinicalSignificance", "Type",
    "Chromosome", "Start", "ReferenceAlleleVCF", "AlternateAlleleVCF",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def filter_chunk(chunk: pd.DataFrame) -> tuple[dict[str, int], pd.DataFrame]:
    counts = {"raw_clinvar": len(chunk)}
    grch38 = chunk.loc[chunk["Assembly"] == "GRCh38"]
    counts["grch38"] = len(grch38)
    panel = grch38.loc[grch38["GeneSymbol"].isin(GENES)]
    counts["ten_gene_panel"] = len(panel)
    labelled = panel.loc[panel["ClinicalSignificance"].isin(ACCEPTED)].copy()
    counts["accepted_clinical_significance"] = len(labelled)
    labelled["binary_label"] = labelled["ClinicalSignificance"].isin(PATHOGENIC).astype(int)
    snvs = labelled.loc[labelled["Type"] == "single nucleotide variant"].copy()
    counts["snvs"] = len(snvs)
    return counts, snvs


def canonical_comparison(frame: pd.DataFrame, path: Path) -> dict[str, object]:
    archived = pd.read_csv(path, sep="\t", dtype=str)
    generated = frame.rename(
        columns={
            "Start": "Position",
            "ReferenceAlleleVCF": "REF",
            "AlternateAlleleVCF": "ALT",
        }
    )
    columns = ["VariationID", "GeneSymbol", "binary_label", "Chromosome", "Position", "REF", "ALT"]
    for table in (generated, archived):
        for column in columns:
            table[column] = table[column].astype(str).str.strip()
    left = generated.set_index("VariationID")[columns[1:]].sort_index()
    right = archived.set_index("VariationID")[columns[1:]].sort_index()
    generated_ids = set(left.index)
    archived_ids = set(right.index)
    shared = sorted(generated_ids & archived_ids)
    field_mismatches = int((left.loc[shared] != right.loc[shared]).any(axis=1).sum())
    return {
        "archived_path": str(path),
        "archived_sha256": sha256(path),
        "identifier_symmetric_difference": len(generated_ids ^ archived_ids),
        "shared_identifier_field_mismatches": field_mismatches,
        "exact_canonical_record_equivalence": generated_ids == archived_ids and field_mismatches == 0,
    }


def build(args: argparse.Namespace) -> dict[str, object]:
    totals = {key: 0 for key in EXPECTED_COUNTS if key != "corrected_canonical"}
    pieces = []
    reader = pd.read_csv(
        args.clinvar,
        sep="\t",
        compression="infer",
        dtype=str,
        usecols=USE_COLUMNS,
        chunksize=args.chunk_size,
        low_memory=False,
    )
    for chunk in reader:
        counts, snvs = filter_chunk(chunk)
        for key, value in counts.items():
            totals[key] += int(value)
        if not snvs.empty:
            pieces.append(snvs)
    benchmark = pd.concat(pieces, ignore_index=True)
    benchmark["VariationID"] = benchmark["VariationID"].astype(str).str.strip()
    benchmark["ReferenceAlleleVCF"] = benchmark["ReferenceAlleleVCF"].astype(str).str.upper()
    benchmark["AlternateAlleleVCF"] = benchmark["AlternateAlleleVCF"].astype(str).str.upper()
    observed_ref_alt = set(
        benchmark.loc[
            benchmark["ReferenceAlleleVCF"] == benchmark["AlternateAlleleVCF"],
            "VariationID",
        ]
    )
    corrected = benchmark.loc[~benchmark["VariationID"].isin(EXCLUDED_REF_ALT_IDS)].copy()
    corrected = corrected.sort_values("VariationID").reset_index(drop=True)
    totals["corrected_canonical"] = len(corrected)

    if corrected["VariationID"].duplicated().any():
        raise RuntimeError("corrected benchmark contains duplicate VariationIDs")
    if not corrected["binary_label"].isin((0, 1)).all():
        raise RuntimeError("corrected benchmark contains a non-binary label")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_path = args.output_dir / "clinvar_lgmd_benchmark_generated.tsv"
    report_path = args.output_dir / "clinvar_preprocessing_verification.json"
    corrected.to_csv(output_path, sep="\t", index=False)
    comparison = canonical_comparison(corrected, args.canonical) if args.canonical else None
    stages_match = totals == EXPECTED_COUNTS
    report = {
        "implementation_status": "portable study preprocessing implementation",
        "input": {"path": str(args.clinvar), "sha256": sha256(args.clinvar)},
        "filters": {
            "assembly": "GRCh38",
            "genes": sorted(GENES),
            "pathogenic_labels": sorted(PATHOGENIC),
            "benign_labels": sorted(BENIGN),
            "variant_type": "single nucleotide variant",
            "excluded_ref_equals_alt_variation_ids": sorted(EXCLUDED_REF_ALT_IDS),
            "ordering": "lexicographic string VariationID",
        },
        "stage_counts": totals,
        "expected_stage_counts": EXPECTED_COUNTS,
        "stage_counts_match": stages_match,
        "observed_ref_equals_alt_variation_ids": sorted(observed_ref_alt),
        "output": {"path": str(output_path), "sha256": sha256(output_path)},
        "canonical_comparison": comparison,
        "runtime": {
            "python": platform.python_version(),
            "pandas": pd.__version__,
            "historical_environment_equivalence_claimed": False,
        },
    }
    report["complete_numerical_reproduction"] = bool(
        stages_match
        and observed_ref_alt == EXCLUDED_REF_ALT_IDS
        and (comparison is None or comparison["exact_canonical_record_equivalence"])
    )
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinvar", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--canonical", type=Path)
    parser.add_argument("--chunk-size", type=int, default=250_000)
    args = parser.parse_args()
    return 0 if build(args)["complete_numerical_reproduction"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
