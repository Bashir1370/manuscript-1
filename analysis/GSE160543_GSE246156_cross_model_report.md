# GSE160543 OIPN versus GSE246156 nerve compression: cross-model DRG report

## Scope and source

This is a comparison of **within-study** rat bulk dorsal-root-ganglion (DRG) expression results. The OIPN discovery contrast is GSE160543 oxaliplatin versus vehicle (4 versus 4). The independent injury contrast is GSE246156 L5 day-7 compression versus sham (3 versus 3). Both cohorts assay DRG; the GSE160543 DRG level and sampling time are not documented to match the locked NC L5/day-7 cell. Sciatic-nerve specimens and other NC levels/days are excluded from this primary comparison.

Source analysis scripts: `scripts/02_GSE160543_oxaliplatin_primary.R`, `scripts/03_GSE160543_pathway_analysis.R`, `scripts/04_GSE160543_GO_GSVA_analysis.R` on `main`; `scripts/05_GSE246156_prepare_NC_L5_day7.R` and `scripts/GSE246156/02_DESeq2_GSE246156.R`, `03_GSEA_GO_GSE246156.R`, `04_GSVA_GSE246156.R` on `NC-GSE246156-analysis`. OIPN output is committed in `results/GSE160543_Oxaliplatin_vs_Vehicle/Pathway_analysis/`. The NC author-run output was reviewed from `pathway_analysis(10).zip` (SHA-256: `496bb03ee041baf348d4190fbba41f0e90d444b61a929250dcc6d056a906a818`) with input/QC files supplied separately as `DESeq2.zip` and `GSE246156_input_QC.zip`.

The compact NC GSEA/GSVA tables and the 50-pathway comparison are stored under `results/` on this branch. NC sample-level count and VST matrices remain in the author's local results directory.

## Input and QC decision

NC sample labels match the locked sham GSM7863794–3796 and compression GSM7863770–3772 assignments. Each count file has 32,883 features with no missing entries or unmatched features. Raw library count sums range from 15,801,858 to 19,047,123 and match the input audit. The recorded filter (counts >= 10 in at least three samples) reproduces all 14,095 tested genes; the FDR < 0.05 DEG table reproduces 6,826 rows (3,245 up, 3,581 down). All three NC samples separate from all three sham samples on PC1 of the VST PCA. Within-group VST correlations are approximately 0.968–0.985; GSM7863772 is distinct on PC2 but correlates 0.973 with both other NC samples and is retained.

The top ten genes account for approximately 20–23% of sham raw counts but 9–13% of compression raw counts. This compositional shift and the unusually large DEG fraction merit caution about normalization, tissue composition and biological interpretation. Input and output consistency checks do not resolve these alternatives.

## Ranked Hallmark GSEA

Each study was analyzed independently with signed DESeq2 Wald statistics, rat Hallmark gene sets, fgsea and BH correction within its 50 pathways. The recorded sessions used R 4.5.2, msigdbr 26.1.1 and fgsea 1.36.2. Both scripts specify 10,000 permutations and set-size limits 15–500. Gene identifiers and duplicate symbols are handled differently by the two original pipelines, and the measured gene universes/pathway sizes differ; NES values are compared by direction and within-study FDR, rather than interpreted as matched effect sizes.

Across 50 matched pathways, 30 are FDR < 0.05 in OIPN and 33 in NC. Twenty-two are significant in both: **18 concordant** and **four discordant** in signed NES. Eight are significant only in OIPN, 11 only in NC, and nine in neither. “Only” describes this significance pattern and does **not** establish specificity or an interaction between models.

| Hallmark pathway | OIPN NES (FDR) | NC NES (FDR) |
| --- | ---: | ---: |
| E2F Targets | +2.77 (0.0007) | +3.26 (0.0014) |
| G2M Checkpoint | +2.90 (0.0007) | +2.65 (0.0014) |
| Interferon Gamma Response | +2.19 (0.0007) | +2.19 (0.0014) |
| Tnfa Signaling Via Nfkb | +1.69 (0.0007) | +2.57 (0.0014) |
| Il6 Jak Stat3 Signaling | +1.61 (0.0123) | +2.32 (0.0014) |
| Interferon Alpha Response | +2.22 (0.0007) | +1.70 (0.0033) |
| Mitotic Spindle | +2.26 (0.0007) | +1.35 (0.0224) |
| P53 Pathway | +2.04 (0.0007) | +1.68 (0.0014) |
| Oxidative Phosphorylation | -2.77 (0.0007) | -2.12 (0.0011) |
| Reactive Oxygen Species Pathway | -1.76 (0.0083) | +0.75 (0.8928) |
| Dna Repair | -0.98 (0.5243) | +1.67 (0.0014) |

The four significant pathways with opposite NES direction are: apical_junction, heme_metabolism, myogenesis, uv_response_dn. Their discordance argues against describing all injury-associated programs as conserved. In particular, the ROS pathway is significant and negative in OIPN but not significant in NC; DNA repair is significant and positive in NC but not significant in OIPN. Neither difference alone demonstrates model specificity.

The reproducible, pathway-level shared pattern is enrichment of E2F/G2M and several immune-response programs, together with decreased oxidative-phosphorylation enrichment. These are changes in bulk DRG transcription, not proof of cell-cycle re-entry in sensory neurons.

## Sample-level Hallmark GSVA

The NC script maps 14,095 VST genes to 12,424 unambiguous unique symbols (101 Ensembl IDs with ambiguous symbols and one duplicate-symbol row excluded). Fifty Hallmark sets have at least 15 mapped genes. The same eight pathways preselected by the original OIPN GSVA script were tested with two-sided Wilcoxon and BH correction within eight tests.

| Preselected pathway | OIPN p / FDR (4+4) | NC p / FDR (3+3) | NC sample separation |
| --- | ---: | ---: | --- |
| G2M Checkpoint | 0.029 / 0.057 | 0.100 / 0.160 | Compression_above_Sham |
| E2F Targets | 0.029 / 0.057 | 0.100 / 0.160 | Compression_above_Sham |
| Mitotic Spindle | 0.029 / 0.057 | 0.700 / 0.700 | overlap |
| P53 Pathway | 0.057 / 0.076 | 0.200 / 0.267 | overlap |
| Interferon Alpha Response | 0.029 / 0.057 | 0.400 / 0.457 | overlap |
| Interferon Gamma Response | 0.057 / 0.076 | 0.100 / 0.160 | Compression_above_Sham |
| Tnfa Signaling Via Nfkb | 0.114 / 0.131 | 0.100 / 0.160 | Compression_above_Sham |
| Il6 Jak Stat3 Signaling | 0.200 / 0.200 | 0.100 / 0.160 | Compression_above_Sham |

In NC, E2F, G2M, interferon-gamma, TNFA/NFKB and IL6/JAK/STAT3 have all three compression scores higher than all three sham scores. p53, mitotic spindle and interferon-alpha show group overlap. The direction of all eight NC mean-score differences agrees with the corresponding NC GSEA NES. None of the NC GSVA tests reaches FDR < 0.05 (minimum FDR 0.16). With three untied observations per group, even complete separation yields a minimum two-sided exact Wilcoxon p = 0.1. These scores are supportive descriptions of per-sample consistency, not a second significance claim. GSVA score magnitudes are not compared across cohorts.

## Interpretation for the manuscript

**Supported:** Within each independent DRG study, the ranked analysis detects several similarly directed cell-cycle-associated and immune-associated transcriptional programs. At the NC sample level, five of eight preselected pathway scores completely separate the groups. The external injury model shows that those shared programs are not unique to the oxaliplatin contrast.

**Not established:** neuronal proliferation or cell-cycle re-entry; senescence in a defined cell type; cell-specific origin of these signals; OIPN-specific mechanisms; causal mediation of neuropathic pain; equivalence of anatomical DRG level/time between cohorts. A difference in FDR or NES between independent analyses is not a formal test of differential response. The opposing pathways and broad NC transcriptional/compositional shift should remain visible in the Results and Discussion.

## Recommended reporting sentence

“Independent, within-study analyses of rat bulk DRG identified concordant enrichment of E2F/G2M and selected inflammatory transcriptional programs, together with reduced oxidative-phosphorylation enrichment, after oxaliplatin exposure and L5 nerve compression. Sample-level GSVA showed complete group separation for five preselected programs in the compression cohort, while exact tests with three samples per group did not meet an adjusted significance threshold. These observations identify shared tissue-level injury-associated signatures without establishing neuronal cell-cycle re-entry or a causal pain mechanism.”

## Next validation before stronger claims

1. Sensitivity analyses in separate NC level/day cells, without treating levels from an unidentified animal as independent replicates.
2. Examine shared gene-level effects and leading-edge overlap with explicit identifier matching; evaluate whether selected programs depend on the high-abundance compositional changes.
3. Interpret existing rat qPCR/behavior as tissue-level expression and phenotype data. Cell-resolved or functional experiments are required for neuronal or causal claims.
