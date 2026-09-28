# Computational redesign for the OIPN manuscript

Status: analysis plan; no new result has been computed or verified.

## Revised research question

How do cell-cycle-related and inflammatory **bulk DRG transcriptional responses** in oxaliplatin-treated rats compare with responses to paclitaxel and nerve compression, and how much of the observed signal is compatible with non-neuronal cell-state or composition changes?

This study cannot establish neuronal cell-cycle re-entry, senescence, temporal ordering, or a causal mechanism of pain without additional experiments.

## Existing evidence and intended role

- GSE160543: rat DRG RNA-seq; vehicle, oxaliplatin, paclitaxel groups (verify all sample metadata and biological replicate IDs before analysis).
- GSE246156: rat nerve-compression/sham RNA-seq. Audit spinal level, side, harvest time, animal ID, and independence before any model fit; do not assume that every DRG sample is an independent animal.
- Author's rat OIPN experiment: n=5 animals per group; longitudinal behavior, endpoint bulk DRG qPCR of Cdk1, Cdc20, Cdkn1a and inflammatory mediators. This validates tissue-level expression and behavior, not cell identity or causation.

## Stage 0 — Metadata and reproducibility gate

1. Export a one-row-per-sample sheet with accession, animal ID if available, group, drug, dose, time, DRG level, side, sequencing batch and source.
2. Inspect original raw/processed count format, gene identifiers, library sizes, sample relationships, PCA, and outliers using prespecified rules.
3. Recompute the manuscript's reported counts and reconcile the mismatch: 295 upregulated + 88 downregulated = 383, while 384 OIPN DEGs are reported.
4. Verify the 63-gene overlap, membership of Cdk1/Cdc20/Cdkn1a, and direction, log2FC, standard error, and adjusted p-value of every shared gene in both contrasts.
5. Confirm WGCNA's inputs, biological sample count, soft threshold and module stability; if sample independence or module robustness fails, remove WGCNA from the primary evidence rather than preserving it for continuity.

**Decision gate:** If the compression contrast is not well defined or independent, do not use it for the primary claim.

## Stage 1 — Primary differential expression and ranked pathway analysis

Within GSE160543, fit a model with treatment group and estimate:
- oxaliplatin versus vehicle;
- paclitaxel versus vehicle;
- oxaliplatin versus paclitaxel.

Use the same normalization, gene filtering, annotation, and multiple-testing convention for all three contrasts. Report effect sizes and uncertainty. A gene significant in one contrast and not another is not automatically drug-specific.

For GSE246156, fit a design that respects animal, side, level, and time where identifiable; otherwise restrict to a coherent matched comparison. Analyze each study separately. Do not combine raw counts across studies or interpret cross-study differences as an unconfounded drug-versus-injury effect.

Analyze prespecified cell-cycle, DNA-damage/stress, inflammatory, and cell-type-associated gene sets using full ranked statistics, followed by sensitivity analyses. Do not infer pathway protein activity from enrichment of mRNA.

## Stage 2 — Shared and differential response

Show:
- direction concordance and effect-size scatterplots for shared measurable genes;
- ranked gene-set results for each within-study contrast;
- a table of shared, divergent and insufficiently resolved signals;
- direct oxaliplatin-versus-paclitaxel results where supported within GSE160543.

Cross-study comparison with nerve compression is descriptive or a carefully qualified meta-analysis of comparable effect estimates; time, tissue-level and protocol differences remain limitations. An absent DEG call is not evidence of absence of effect.

## Stage 3 — Alternative cellular interpretation

Use a suitable published rat/mouse DRG single-cell or spatial reference only after verifying tissue, species mapping, cell labels and coverage. Examine neuronal, satellite-glial, Schwann-cell and immune-cell marker/signature behavior. If a reference-based deconvolution is feasible, test its robustness to reference choice and composition assumptions. Treat such results as estimates/hypothesis-generating, not proof that a given transcript originated from a specific cell type. Bulk qPCR cannot resolve this.

## Stage 4 — Validation and manuscript

Place the animal phenotype and endpoint qPCR after the computational analyses as an independent check of behavioral change and bulk mRNA direction. Use animals, not PCR technical wells, as biological replicates. Reassess longitudinal behavior with a repeated-measures approach appropriate to the raw measurements.

Proposed headline: **Shared and context-dependent transcriptional responses in rat dorsal root ganglia after oxaliplatin exposure and nerve injury**. Revise the headline once results are available.

Remove from central conclusions: “in post-mitotic sensory neurons,” “cell-cycle re-entry,” “senescence-like inflammatory state,” “drivers,” and “therapeutic targets.” These may appear only as explicitly untested hypotheses where relevant.

## Outcomes and stopping rules

- **Robust shared signal:** report a conserved bulk DRG injury-associated transcriptional program.
- **Differential oxaliplatin signal:** report context-dependent expression only if the direct contrast and effect-size analyses support it; do not claim molecular specificity beyond the tested conditions.
- **Predominantly non-neuronal signature:** report this as an alternative explanation for bulk observations, without assigning cellular origin to each gene.
- **Weak/unstable cell-cycle result:** do not force the original three-gene narrative; refocus on the stable biological result or acknowledge that the current data do not support a strong cell-cycle-centered paper.

## Data and provenance

Store metadata, scripts, environment versions, figure-generating code, outputs and accession/source citations in this repository as they are actually obtained. Do not commit unpublished animal-level identifiers, credentials, or invented/raw data. Record exact URLs, download dates and checksums for public input files.

Existing related literature to assess for novelty:
- Comparative transcriptome of oxaliplatin and paclitaxel DRG using GSE160543: PMID 36822350.
- Prior OIPN-versus-nerve-injury bioinformatics study: PMID 38716040.
