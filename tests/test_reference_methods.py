import argparse
import csv
import hashlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "methods"))

from build_canonical_split import build
from build_clinvar_benchmark import filter_chunk
from sequence_extraction import construct


class ReferenceMethodTests(unittest.TestCase):
    def test_evidence_manifest_hashes(self):
        manifest = json.loads(
            (ROOT / "metadata/evidence_manifest.json").read_text(encoding="utf-8")
        )
        records = list(manifest["files"].values())
        records += list(manifest["gene_heldout"]["files"].values())
        for item in records:
            path = ROOT / item["path"]
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), item["sha256"])

    def test_sequence_substitution_and_coordinate_conversion(self):
        handle = io.BytesIO(b"ACGTACGTACG")
        index = {"NC_000001.11": (11, 0, 11, 11)}
        reference, alternate, accession = construct(
            handle, index, "1", 6, "C", "T", 5
        )
        self.assertEqual(accession, "NC_000001.11")
        self.assertEqual(reference, "TACGT")
        self.assertEqual(alternate, "TATGT")

    def test_method_validation_summary(self):
        report = json.loads(
            (ROOT / "metadata/method_validation.json").read_text(encoding="utf-8")
        )
        self.assertEqual(
            report["clinvar_preprocessing"]["stage_counts"]["corrected_canonical"],
            6960,
        )
        self.assertEqual(report["canonical_split"]["assignment_mismatches"], 0)
        self.assertEqual(report["canonical_split"]["selected_test_folds"], [1, 7, 8])
        self.assertEqual(report["sequence_extraction"]["mismatches"], 0)
        self.assertTrue(
            report["gene_heldout_bootstrap"]["all_intervals_reproduced"]
        )

    def test_split_generation_matches_distributed_assignments(self):
        output_parent = ROOT / "verification_outputs"
        output_parent.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=output_parent) as temporary:
            args = argparse.Namespace(
                variants=ROOT / "data/canonical_variants.tsv",
                output_dir=Path(temporary),
                archived_train=None,
                archived_test=None,
                archived_assignment=None,
            )
            report = build(args)
        self.assertTrue(report["exact_assignment_reproduction"])
        self.assertEqual(
            report["fold_sizes"],
            [699, 698, 691, 685, 713, 694, 692, 703, 693, 692],
        )
        self.assertEqual(report["candidate_combinations"], 120)
        self.assertEqual(report["eligible_combinations"], 49)
        self.assertEqual(report["selected"]["test_folds"], [1, 7, 8])
        self.assertEqual(report["comparison"]["fold_id_matches"], 6960)
        self.assertEqual(report["comparison"]["assignment_mismatches"], 0)

    def test_clinvar_filter_labels_and_snv_selection(self):
        rows = pd.DataFrame([
            {"VariationID": "1", "GeneSymbol": "CAPN3", "Assembly": "GRCh38", "ClinicalSignificance": "Pathogenic", "Type": "single nucleotide variant"},
            {"VariationID": "2", "GeneSymbol": "CAPN3", "Assembly": "GRCh38", "ClinicalSignificance": "Benign", "Type": "single nucleotide variant"},
            {"VariationID": "3", "GeneSymbol": "CAPN3", "Assembly": "GRCh38", "ClinicalSignificance": "Uncertain significance", "Type": "single nucleotide variant"},
            {"VariationID": "4", "GeneSymbol": "OTHER", "Assembly": "GRCh38", "ClinicalSignificance": "Pathogenic", "Type": "single nucleotide variant"},
        ])
        counts, filtered = filter_chunk(rows)
        self.assertEqual(
            counts,
            {
                "raw_clinvar": 4,
                "grch38": 4,
                "ten_gene_panel": 3,
                "accepted_clinical_significance": 2,
                "snvs": 2,
            },
        )
        self.assertEqual(filtered["binary_label"].tolist(), [1, 0])

    def test_bootstrap_validation_rows(self):
        with (ROOT / "metadata/gene_heldout_bootstrap_validation.tsv").open(
            encoding="utf-8", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 24)
        self.assertTrue(all(row["status"] == "PASS" for row in rows))
        maximum = max(
            max(float(row["lower_abs_difference"]), float(row["upper_abs_difference"]))
            for row in rows
        )
        self.assertLess(maximum, 5e-15)


if __name__ == "__main__":
    unittest.main()
