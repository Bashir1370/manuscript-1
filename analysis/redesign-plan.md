# Computational redesign for the OIPN manuscript

Status: study design plus the author's executed GSE160543 Oxaliplatin versus Vehicle analysis. See [the primary result audit](gse160543-primary-results.md). The originally submitted DEG and cross-model overlap results have not been reproduced.

## Primary question

Which **bulk DRG transcriptional changes** follow oxaliplatin exposure, and which of those changes are also observed in an anatomically and temporally defined nerve-compression comparison?

The shared changes are neuropathy-associated candidates, not oxaliplatin-specific mechanisms. These data cannot establish neuronal localization, functional cell-cycle re-entry, senescence, temporal order or causality.

## Role of each dataset

- **GSE160543, primary OIPN discovery:** Use the four oxaliplatin and four vehicle rat DRG samples. The four paclitaxel samples are present in GEO but are **excluded from the manuscript's primary and planned secondary analyses**. Their existence is documented in the complete sample manifest; they do not have to be analyzed merely because they are in the same accession.
- **GSE246156, external injury context:** Compare compression with matched sham using a clearly defined DRG level and day. The full dataset includes L4, L5 and L6 DRG at days 3 and 7, three samples per design cell, plus sciatic nerve. Do not pool levels/days or assume replicate labels identify independent animals across tissues.
- **Author's independent OIPN rat experiment:** Behavioral phenotype and endpoint bulk DRG qPCR (five animals per group) validate phenotype and direction of tissue-level expression only.

## Stage 0 — Metadata and provenance (completed in part)

See [stage0-metadata-audit.md](stage0-metadata-audit.md) and the sample tables under `data/metadata/`. Remaining checks:
1. Recover the exact 18 GSM IDs and code/count source behind the original “9 sham versus 9 NC” analysis; if unavailable, describe the original comparison as unreproducible and define a new one transparently.
2. Verify count files, gene identifiers, biological independence, library QC and outlier handling.
3. Recompute the reported 384 OIPN DEGs and reconcile 295 up + 88 down = 383.
4. Recompute the reported 63-gene overlap and direction/effect sizes in each model.
5. Reassess WGCNA against the number of independent animals and its sensitivity to tissue level/day. Retain it only if robust.

## Stage 1 — OIPN discovery

Fit one prespecified differential-expression contrast in GSE160543: **Oxaliplatin versus Vehicle** (4 versus 4). Report effect size, standard error, adjusted p-value, filtering and gene mapping for all tested genes. Evaluate prespecified cell-cycle, DNA-damage/stress and inflammatory gene sets from full ranked statistics rather than only a DEG cutoff. Inspect Cdk1, Cdc20 and Cdkn1a within the broader response; do not select the entire story from PPI rank.

The four paclitaxel samples do not enter this model or its conclusions. No claim of oxaliplatin specificity follows from this two-group contrast.

## Stage 2 — Matched nerve-injury context

Define a valid compression-versus-sham contrast in GSE246156 only after auditing DRG level, day and animal independence. Analyze this dataset separately. Compare effect directions and sizes for measurable genes and gene sets between the two studies; account for differing treatments, tissue definitions and sampling times in interpretation.

A gene that passes a threshold in only one study is not established as specific to that condition. The cross-model intersection is descriptive evidence of a shared response, not a filter for false positives or a test of OIPN-specific biology.

## Stage 3 — Alternative cellular interpretation

If a suitable published DRG single-cell or spatial reference is verified, examine neuronal, satellite-glial, Schwann-cell and immune-cell signatures as an exploratory explanation for bulk expression. A reference-based estimate does not prove that a transcript came from a specific cell type.

## Stage 4 — Independent rat results and writing

Use each animal as a biological replicate in qPCR and use an appropriate repeated-measures analysis for longitudinal behavior if raw data permit. Revise title, abstract, discussion and conclusion around bulk DRG transcription. Remove claims of neuronal cell-cycle activation, demonstrated senescence, causal pain drivers and validated therapeutic targets.

## Outcome rules

- Robust concordant response: report an injury-associated bulk DRG transcriptional pattern observed after oxaliplatin and compression.
- Weak or discordant cell-cycle response: report the actual stable findings and do not force a Cdk1/Cdc20/Cdkn1a mechanism.
- Apparent non-neuronal signature: discuss it as a plausible interpretation requiring cell-resolved validation.
- No independent or well-matched NC contrast: focus the paper on OIPN transcription plus the independent in vivo qPCR/behavioral results, and remove cross-injury conservation claims.

## Provenance and novelty

Store code, metadata, source checksums and derived outputs in this public repository; keep unpublished manuscript drafts and confidential animal-level data out of it. A prior OIPN-versus-nerve-injury bioinformatics study (PMID 38716040) requires a sharper question and transparent limitations. An earlier oxaliplatin/paclitaxel analysis of GSE160543 (PMID 36822350) is background literature, not a reason to add paclitaxel to this manuscript.
