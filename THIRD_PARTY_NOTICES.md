# Third-party data and resource notices

These notices document provenance and applicable provider conditions. They do
not relicense third-party material and are not legal advice.

## ClinVar

The canonical benchmark contains reduced fields derived from ClinVar. NCBI
states that it places no restrictions on use or distribution of molecular
data, while noting that submitters may retain rights in submitted material.
ClinVar requests source attribution. ClinVar also cautions that its records are
not intended for direct diagnostic use without professional review and that
submitted classifications are not independently verified by NIH.

- ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/
- NCBI policies: https://www.ncbi.nlm.nih.gov/home/about/policies/

The exact ClinVar release and retrieval date used upstream remain
unresolved. The public bundle therefore identifies the source and this
provenance limitation without asserting a release date.

## GPN-Star genome-wide scores

The V100, M447, and P243 values correspond to the Song Lab GPN-Star score
release:

https://huggingface.co/datasets/songlab/gpn-star-scores

The official dataset card identifies Apache License 2.0 and DOI
`10.57967/hf/9849`. A copy of Apache License 2.0 is provided at
`licenses/Apache-2.0.txt`. This repository contains only benchmark-aligned
selected score rows and derived threshold classifications, not the full score
dataset.

## AlphaGenome AVI scores

The included AlphaGenome fields are 2,094 benchmark-aligned test-variant AVI
scores retrieved through the AlphaGenome API with package version 0.10.0. They
are not a copy of the complete AlphaGenome Atlas and do not include general
prediction tracks, splicing-score downloads, or AVI feature-importance tracks.

The current AlphaGenome Services terms permit publication, sharing, and
adaptation of outputs subject to the agreement and expressly distinguish the
AVI Score from the pass-through notice applying AlphaGenome Output Terms to
other output types. The official downloads page separately classifies AVI SNV
scores as a permissive-use downloadable artifact for commercial and
non-commercial use; other artifacts, including merged splicing scores and AVI
feature-importance scores, have different conditions. API provenance does not
turn these rows into a copy of the downloadable Atlas artifact, so both the
service-output provisions and the output-type distinction are recorded.

The archived quantile was transformed to the preserved PHRED-scale score as
`-10 * log10(1 - quantile)`. Binary predictions use the documented threshold
selected from training data. These transformations are modifications by this
study, not provider-supplied calibrated pathogenicity probabilities.

Conditions retained for the included AVI rows:

- research use only in this repository; no clinical or medical use;
- no use to train another machine-learning model;
- cite AlphaGenome and the AlphaGenome Atlas;
- retain provider terms links, API provenance, and modification notice; and
- do not apply the repository MIT license to AVI values.

Governing sources reviewed 2026-11-01:

- https://deepmind.google.com/science/alphagenome/terms
- https://deepmind.google.com/science/alphagenome/output-terms
- https://google-deepmind-alphagenome.readthedocs-hosted.com/faqs.html
- https://github.com/google-deepmind/alphagenome

The exact AlphaGenome terms version effective at the recorded API retrieval
date (2026-10-13) was not archived locally. The current terms support limited
sharing of the selected AVI score records with these notices, but no claim of
unrestricted rights or license assignment is made.

Suggested citations:

- Avsec et al. AlphaGenome: advancing regulatory variant effect prediction
  with a unified DNA sequence model. *Nature* 649, 1206–1218 (2026).
  https://doi.org/10.1038/s41586-025-10014-0
- AlphaGenome Atlas preprint. https://doi.org/10.64898/2026.09.16.26363192

## Genomic foundation models

No pretrained weights, tokenizer files, adapted checkpoints, or embedding
arrays are redistributed. The repository includes numerical predictions and
checkpoint identifiers as scientific evidence.

- `multimolecule/dnabert2`: model card identifies AGPL-3.0.
- `zhihan1996/DNABERT-2-117M`: no explicit model-card license was verified.
- Nucleotide Transformer V2: model card identifies CC-BY-NC-SA-4.0.
- HyenaDNA: model card identifies BSD-3-Clause.
- Caduceus project source: Apache-2.0.
- modernGENA project source: MIT; the checkpoint card did not expose a
  separately verified weight license.
- Evo2 checkpoint: model card identifies Apache-2.0.

These upstream licenses govern the providers' artifacts. The repository does
not purport to relicense them, and model weights must be obtained separately
from their providers.
