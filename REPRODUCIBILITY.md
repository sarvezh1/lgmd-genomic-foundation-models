# Reproducibility

## Supported scope

The package reproduces numerical analyses from canonical variants and frozen
variant-level predictions. It supports all 28 fixed-test configurations, four
main tables, four figure datasets, shared-error and recurrent-difficulty
analyses, pooled gene-held-out results, and the gene-stratified locus
bootstrap.

It does not provide complete end-to-end regeneration of every foundation-model
embedding, adapted checkpoint, or native provider score. Model weights,
tokenizers, embeddings, raw ClinVar data, and the GRCh38 FASTA are not
distributed.

## Installation

Python 3.10 or newer is recommended.

```bash
python -m venv .venv
python -m pip install -r requirements.txt
```

The fixed-test metric workflow uses only the Python standard library. NumPy
and pandas are required for split generation and bootstrap/reference tests;
Matplotlib is required to render reproduced figures.

## Verify the repository

```bash
python scripts/validation/validate_repository.py
python scripts/validation/validate_canonical_bundle.py \
  --variants data/canonical_variants.tsv \
  --expect-canonical
python -m unittest discover -s tests -p "test_*.py"
```

The repository validator checks evidence hashes, scientific counts, approved
asset hashes, source-to-result coverage, prohibited files, private paths, and
public-content hygiene.

## Reproduce fixed-test metrics

```bash
python scripts/reproduction/reproduce_metrics.py
```

This validates VariationID alignment, labels, configuration coverage, score
semantics, thresholds, and prediction consistency. It recalculates metrics for
26 primary and two sensitivity configurations at absolute tolerance `1e-12`,
then verifies the 5/2/1 shared-error intersections and the 787-variant
recurrent-difficulty analysis.

Outputs are written to `verification_outputs/prediction_metrics/`.

## Reproduce tables and figure data

```bash
python scripts/reproduction/reproduce_manuscript.py
```

The command recalculates Tables 1-4, derives all four figure datasets, produces
452 source-to-result records, and optionally renders reproducibility figures.
The approved publication figures in `figures/` are immutable and are not
overwritten. Outputs are written to `verification_outputs/manuscript/`.

## Reproduce gene-held-out bootstrap intervals

```bash
python scripts/methods/gene_heldout_bootstrap.py
```

This performs 1,000 paired within-gene whole-locus bootstrap replicates with
seed 1605 and compares all 24 intervals with the distributed references.
Outputs are written to `verification_outputs/gene_heldout_bootstrap/`.

## Reproduce the canonical split

```bash
python scripts/methods/build_canonical_split.py \
  --variants data/canonical_variants.tsv \
  --output-dir verification_outputs/canonical_split
```

The program performs fold generation and candidate selection before comparing
with the assignments embedded in the distributed canonical table. Expected
checkpoints are ten fold sizes
`[699, 698, 691, 685, 713, 694, 692, 703, 693, 692]`, 120 candidate
combinations, 49 eligible combinations, selected test folds `(1, 7, 8)`,
4,866 training variants, 2,094 test variants, and zero assignment or locus
overlap mismatches.

## Rebuild the ClinVar benchmark

```bash
python scripts/methods/build_clinvar_benchmark.py \
  --input /path/to/variant_summary.txt.gz \
  --output-dir verification_outputs/clinvar
```

Use `--help` for optional comparison arguments. Complete numerical comparison
requires the authenticated source matching SHA-256
`acc8e6aef2f66c6fe86130f06bfb0539b3c64a0f83f9629fb99024a6620122c6`.
The exact ClinVar release date and retrieval command are unavailable.

## Validate sequence construction

```bash
python scripts/methods/sequence_extraction.py --help
```

Sequence validation requires a user-supplied indexed NCBI RefSeq GRCh38.p14
FASTA (`GCF_000001405.40`) and study sequence tables, which are not
distributed. The implementation uses strict chromosome, boundary, and central
REF validation.

## Output and line-ending policy

All generated outputs are isolated under `verification_outputs/`, which is
ignored by Git. TSV and JSON writers use explicit LF line endings.
`.gitattributes` enforces LF checkout for text and preserves PNG/PDF bytes.
Distributed evidence, figures, and tables are never rewritten by documented
commands.

## Remaining limitations

- Complete original runtime equivalence is not established.
- Several model-training, embedding-extraction, adaptation, and native-score
  retrieval workflows require external artifacts or incomplete environments.
- Exact ClinVar acquisition metadata are unavailable.
- No external cohort validation was performed.
- Provider terms apply separately from the repository's MIT software license.
