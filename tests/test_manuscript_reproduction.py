import csv, hashlib, json, struct, sys, unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts/reproduction"))
from reproduce_manuscript import (
    figure_data, fixed_test_maps, reproduce_table1, reproduce_table2,
    reproduce_table3, reproduce_table4, verify_gene_heldout_bundle,
)

def read_tsv(path):
    with path.open(encoding="utf-8",newline="") as handle:
        return list(csv.DictReader(handle,delimiter="\t"))

class ManuscriptScientificRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bundle=verify_gene_heldout_bundle()
        cls.canonical,cls.registry,cls.by_config,cls.metrics=fixed_test_maps()
        cls.benchmark_bundle=json.loads((ROOT/"metadata/evidence_manifest.json").read_text(encoding="utf-8"))
        cls.source_hash=cls.benchmark_bundle["files"]["primary_predictions"]["sha256"]

    def test_gene_heldout_bundle_hashes_and_protocol_separation(self):
        self.assertEqual(self.bundle["configurations"],8)
        self.assertEqual(self.bundle["prediction_records"],16752)
        self.assertIn("distinct from fixed-test",self.bundle["protocol"])
        for item in self.bundle["files"].values():
            path=ROOT/item["path"]
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),item["sha256"])

    def test_table1_metrics_and_scope(self):
        provenance=[]
        rows=reproduce_table1(self.registry,self.metrics,provenance,self.source_hash)
        self.assertEqual(len(rows),14)
        self.assertEqual({int(r["MCC_rank"]) for r in rows},set(range(1,15)))
        native={r["model_display"]:r for r in rows}
        self.assertAlmostEqual(native["GPN-Star-V"]["Specificity"],1564/1585,15)
        self.assertAlmostEqual(native["AlphaGenome AVI"]["Specificity"],1575/1585,15)
        self.assertLessEqual(max(float(r["absolute_difference"]) for r in provenance),1e-12)

    def test_table2_signed_differences(self):
        provenance=[]
        rows=reproduce_table2(self.registry,self.metrics,provenance,self.source_hash)
        self.assertEqual(len(rows),6)
        for row in rows:
            for metric in ("MCC","AUROC","AUPRC"):
                self.assertAlmostEqual(row[f"adaptation_delta_{metric}"],row[f"adapted101_{metric}"]-row[f"frozen101_{metric}"],15)
                self.assertAlmostEqual(row[f"context501_delta_{metric}"],row[f"frozen501_{metric}"]-row[f"frozen101_{metric}"],15)
        self.assertLessEqual(max(float(r["absolute_difference"]) for r in provenance),1e-12)

    def test_table3_intersections(self):
        provenance=[]
        rows,summary=reproduce_table3(self.registry,self.by_config,provenance,self.source_hash)
        self.assertEqual(len(rows),2)
        self.assertEqual((summary["all_12_101bp_wrong"],summary["all_24_sequence_wrong"],summary["all_26_primary_wrong"]),(5,2,1))
        self.assertEqual(summary["universal_residual_variation_id"],"4538462")

    def test_table4_predictions_and_reproduced_intervals(self):
        provenance=[]
        rows=reproduce_table4(provenance,self.bundle)
        self.assertEqual(len(rows),8)
        self.assertEqual(sum(r["verification_status"]=="independently_recalculated" for r in provenance),72)
        self.assertLessEqual(max(float(r["absolute_difference"]) for r in provenance),1e-12)
        evidence=read_tsv(ROOT/"data/gene_heldout_predictions.tsv")
        self.assertEqual(len(evidence),8*2094)
        self.assertTrue(all(r["heldout_gene"]==r["GeneSymbol"] for r in evidence))

    def test_figure_source_values(self):
        provenance=[]
        f1,f2,f3,f4=figure_data(self.registry,self.metrics,self.by_config,provenance,self.source_hash)
        self.assertEqual((len(f1),len(f2),len(f3),len(f4)),(26,36,3,5))
        self.assertEqual([r["universally_misclassified"] for r in f3],[5,2,1])
        self.assertEqual([(r["variants"],r["correct_by_either_native"]) for r in f4],[(1352,1348),(355,351),(254,252),(131,124),(2,1)])
        self.assertEqual(len({r["configuration_id"] for r in f1}),26)

    def test_provenance_coverage_and_deterministic_outputs(self):
        rows=read_tsv(ROOT/"metadata/source_to_result.tsv")
        self.assertEqual(len(rows),452)
        self.assertEqual(sum(r["verification_status"]=="independently_recalculated" for r in rows),398)
        self.assertEqual(sum(r["verification_status"]=="mathematically_derived_from_verified_values" for r in rows),42)
        self.assertEqual(sum(r["verification_status"]=="verified_against_archived_result" for r in rows),0)
        self.assertEqual(sum(r["verification_status"]=="verified_source_metadata" for r in rows),12)
        self.assertFalse(any(r["verification_status"]=="FAIL" for r in rows))
        first_provenance=[]
        second_provenance=[]
        first=figure_data(self.registry,self.metrics,self.by_config,first_provenance,self.source_hash)
        second=figure_data(self.registry,self.metrics,self.by_config,second_provenance,self.source_hash)
        self.assertEqual(first,second)
        self.assertEqual(first_provenance,second_provenance)

    def test_graphic_dimensions_and_formats(self):
        workflow=(ROOT/"figures/figure1.png").read_bytes()
        self.assertEqual(workflow[:8],b"\x89PNG\r\n\x1a\n")
        width,height=struct.unpack(">II",workflow[16:24])
        self.assertEqual((width,height),(4000,2250))
        self.assertEqual((ROOT/"figures/figure1.pdf").read_bytes()[:4],b"%PDF")
        for number in range(2,6):
            png=ROOT/f"figures/figure{number}.png"
            pdf=ROOT/f"figures/figure{number}.pdf"
            self.assertTrue(png.is_file() and pdf.is_file())
            data=png.read_bytes()
            self.assertEqual(data[:8],b"\x89PNG\r\n\x1a\n")
            width,height=struct.unpack(">II",data[16:24])
            self.assertGreaterEqual(width,1800)
            self.assertGreaterEqual(height,1200)
            self.assertEqual(pdf.read_bytes()[:4],b"%PDF")

if __name__=="__main__":
    unittest.main()
