# Methods

## Benchmark construction

The benchmark was derived from the ClinVar variant summary and restricted to
GRCh38 records in ten LGMD-associated genes: *ANO5, CAPN3, DNAJB6, DYSF, FKRP,
LMNA, SGCA, SGCB, SGCG,* and *TCAP*. Accepted clinical-significance classes
were Pathogenic, Likely pathogenic, Pathogenic/Likely pathogenic, Benign,
Likely benign, and Benign/Likely benign. The first three classes were mapped to
1 and the latter three to 0. Only `single nucleotide variant` records were
retained. VariationIDs 177801 and 255604 were excluded because REF and ALT
were identical.

The filtering stages contained 8,981,173 raw records, 4,456,834 GRCh38
records, 14,931 records in the ten-gene panel, 8,482 records with accepted
clinical significance, 6,962 SNVs, and 6,960 corrected canonical variants.
VariationID was handled as a string and records were sorted lexicographically.
The portable implementation is `scripts/methods/build_clinvar_benchmark.py`.
The raw source is not distributed, and its exact ClinVar release and retrieval
date are unavailable.

The bundled labels are aggregate ClinVar assertions. Disease-description
specificity and review status are heterogeneous; the labels must not be
interpreted as uniform proof of LGMD causality or as clinical diagnoses.

## Canonical train/test split

Variants were grouped by chromosome-position locus and assigned with the
scikit-learn 1.6.1 `StratifiedGroupKFold` behavior using ten folds,
`shuffle=True`, and random state 42. All 120 combinations of three test folds
were evaluated. Candidates were eligible when their test fraction was within
the inclusive interval 0.30 +/- 0.001. The 49 eligible combinations were
ranked by:

1. deviation of the test pathogenic rate from the complete benchmark;
2. mean absolute per-gene test-fraction deviation;
3. maximum per-gene deviation;
4. test-size deviation; and
5. test-fold tuple.

The selected folds were (1, 7, 8), with fold sizes
`[699, 698, 691, 685, 713, 694, 692, 703, 693, 692]`. The split contains
4,866 training variants at 4,260 loci and 2,094 test variants at 1,834 loci,
with no exact locus overlap. `scripts/methods/build_canonical_split.py`
implements the study fold-assignment behavior and reproduces all 6,960
distributed assignments.

The split is not region-disjoint. Of 1,834 test loci, 1,715 are within 50 bp
of a training locus, 1,770 within 100 bp, 1,820 within 250 bp, and 1,827
within 500 bp; the median nearest distance is 3 bp. Performance therefore
represents held-out variants at new exact loci within the same genes, not
independent genomic-region generalization.

## Sequence construction

Reference contexts use NCBI RefSeq GRCh38.p14 (`GCF_000001405.40`) primary
chromosome accessions. One-based ClinVar positions are converted to zero-based
half-open FASTA coordinates. Odd 101-, 201-, and 501-bp windows are extracted
in forward genomic orientation, the central reference base is validated, and
the ALT sequence is formed by replacing that central base.

`scripts/methods/sequence_extraction.py` rejects unrecognized chromosomes,
out-of-bound windows, and REF mismatches. Validation showed zero mismatches
across 6,960 variants at three context sizes: 41,760 reference and alternate
sequences in total. The indexed FASTA is not distributed.

## Sequence foundation models

Six architectures were studied:

| Model | Representation/readout | Contexts |
| --- | --- | --- |
| DNABERT-2 | Attention-mask-aware mean pooling; ALT - REF | 101, 201, 501 bp |
| Nucleotide Transformer V2 | Attention-mask-aware mean pooling; ALT - REF | 101, 201, 501 bp |
| HyenaDNA | Mean sequence-position pooling; ALT - REF | 101, 201, 501 bp |
| Caduceus-Ph | Reverse-complement-averaged mean pooling; ALT - REF | 101, 201, 501 bp |
| modernGENA | Attention-mask-aware mean pooling; ALT - REF | 101, 201, 501 bp |
| Evo2 | Mean nucleotide-position pooling from `blocks.28.mlp.l3`; ALT - REF | 101, 201, 501 bp |

Frozen representations were classified using L2-regularized logistic
regression. The 101-bp supervised experiments used full-backbone fine-tuning
for five architectures and LoRA for Evo2. Evo2 adapters targeted
`blocks.28.mlp.l1`, `l2`, and `l3` with rank 8, alpha 16, dropout 0.05,
and 362,496 trainable adapter parameters. Model selection used grouped
chromosome-position validation within the 4,866-variant training partition;
the fixed test set was evaluated after selection. Sequence-model hard
predictions use threshold 0.5.

Checkpoint identifiers, revisions, pooling details, and adaptation parameters
are provided in `metadata/model_provenance.tsv` and
`metadata/adaptation_provenance.tsv`. Frozen and adapted DNABERT-2 used
different checkpoint implementations, so their contrast is not a controlled
fine-tuning-only comparison.

## Native variant scores

GPN-Star values use the sign-inverted mutation-rate-calibrated log-likelihood
ratio, `-llr_calibrated`, so larger scores indicate stronger predicted
constraint. V100 is the primary score; M447 and P243 are sensitivity
configurations. Training-selected thresholds were 5.767, 4.889, and 4.897.

AlphaGenome Variant Impact scores were represented on the PHRED scale:

`AVI_PHRED = -10 * log10(1 - AVI_QUANTILE)`

The training-selected AVI threshold was 23.87662024921556. GPN-Star measures
evolutionary constraint and AVI represents predicted functional impact.
Neither score is a calibrated pathogenicity probability, and neither was
trained using the sequence-model readout protocol.

## Fixed-test evaluation

The primary analysis contains 26 configurations: 24 sequence-model
configurations plus GPN-Star V100 and AlphaGenome AVI. GPN-Star M447 and P243
form two sensitivity configurations. Every configuration has 2,094
VariationID-aligned test predictions.

Reported measures are Matthews correlation coefficient, AUROC, average
precision (reported as AUPRC), accuracy, balanced accuracy, precision, recall,
specificity, F1, and the confusion matrix. Native thresholds were selected
using training labels; test labels were not used for threshold selection.
Scores were not normalized across model families.

Shared-error analyses intersect hard-prediction errors by VariationID. The
recurrent-difficulty set contains 787 variants misclassified by at least four
of the twelve original 101-bp sequence configurations. Counts classified
correctly by any additional strategy are retrospective oracle-style coverage,
not evaluated ensemble performance.

## Gene-held-out evaluation and uncertainty

Eight gene-held-out configurations were evaluated. For each frozen sequence
model, a downstream classifier was fitted without labelled variants from the
evaluated gene. Native models were not refitted; only their thresholds were
selected using labelled variants from the other nine genes. Pooled values
combine the ten non-overlapping held-out-gene test subsets, comprising 16,752
prediction records.

Confidence intervals use 1,000 paired bootstrap replicates with seed 1605.
Genomic loci are sampled with replacement within each held-out gene, all
variants from a sampled locus are retained, and the same resample is used for
every configuration. The 2.5th and 97.5th percentiles form the interval.
The portable implementation reproduces all 24 intervals (48 bounds) with a
maximum absolute difference of `4.996003610813204e-16`.

These intervals are conditional on the evaluated genes and fixed fitted
models. They do not represent uncertainty from model retraining, gene
sampling, external cohorts, or multiplicity-adjusted inference.

## Figures and tables

- **Figure 1:** frozen reference–alternate sequence representations,
  alternate-minus-reference difference embedding, logistic-regression
  classification, and held-out evaluation.
- **Figure 2:** AUROC and AUPRC across the 26 primary configurations.
- **Figure 3:** within-model changes after 101-bp adaptation and after moving
  from frozen 101-bp to frozen 501-bp context.
- **Figure 4:** shared errors across the 12 original 101-bp, 24 sequence, and
  26 primary configuration sets.
- **Figure 5:** native-score correct classification within strata of
  sequence-model error burden.
- **Table 1:** primary 101-bp frozen/adapted readouts and native scores.
- **Table 2:** architecture trajectories across adaptation and context.
- **Table 3:** difficult-variant overlap and additional-strategy coverage.
- **Table 4:** pooled gene-held-out performance and bootstrap intervals.

Approved assets are stored in `figures/` and
`tables/publication_table*.tsv`. Full-precision comparisons use
`tables/reference_table*.tsv`.

## Interpretation boundaries

The study is retrospective and limited to ten genes and a single canonical
partition. Nearby train/test loci can share sequence and annotation context.
ClinVar phenotype specificity, review status, and ascertainment are
heterogeneous. Adapted models, frozen readouts, and native scores represent
different learning and information regimes. No independent external cohort,
prospective evaluation, ancestry-stratified analysis, or clinical workflow
validation was performed.
