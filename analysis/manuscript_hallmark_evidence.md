# Evidence archive for manuscript Hallmark figures

## Scope and provenance

Question: Does bulk DRG express a shared transcriptional injury response across neuropathy models, and where does OIPN fit?

Snapshot: `CCI-GSE212311-analysis`, input commit `0c53e61ac5ae45efd16039c703cd0fa5d98747b7` (2026-09-30 audit).
Four independent studies, three injury models, all rat bulk DRG; no pooled expression matrices.

| Study | Contrast | Samples | Ranking / enrichment |
| --- | --- | --- | --- |
| GSE160543 | Oxaliplatin minus vehicle | 4 / 4 | DESeq2 Wald / fgsea, 10,000 permutations; 15–500 genes |
| GSE126773 | Control+OX minus control | 3 / 3 | limma moderated t / clusterProfiler GSEA; 10–500 genes |
| GSE246156 | L5 compression minus sham, day 7 | 3 / 3 | DESeq2 Wald / fgsea, 10,000 permutations; 15–500 genes |
| GSE212311 | Ipsilateral L4–L6 CCI minus sham, day 11 | 3 / 3 | DESeq2 Wald / fgsea, 10,000 permutations; 15–500 genes; source-aware mapping |

Use the CCI `pathway_analysis_source_aware` table, never the older Ensembl-only table.
GSEA uses the eligible full signed statistic ranking, not a significant-DEG subset.
Hallmark sets come from human MSigDB H projected to rat using msigdbr; these are not rat-native curated Hallmarks.
Versions, gene universes, feature-to-symbol rules and overlapping set sizes should be retained from each study's provenance.
The GSE126773 script does not explicitly lock GSEA backend, seed, MSigDB release or BH method; do not describe every backend setting as verified from that script alone.
The extraction preserves exported `padj` / `p.adjust` as within-study adjusted P values and does not recalculate them over selected rows.

## Selection rules

These rules formalize a post hoc manuscript selection after examining the full heatmap; they are not prospectively prespecified hypothesis tests.

1. Shared positive: NES > 0 in all FOUR studies AND within-study FDR < 0.05 in at least THREE. Thus at least one OIPN study and both physical models support each selected pathway in the present data.
2. Shared negative: the exact symmetric rule, NES < 0 in all four AND FDR < 0.05 in at least three.
3. Directional divergence (updated after author clarification): NES > 0 in BOTH OIPN studies and NES < 0 in BOTH NC/CCI studies, or the exact reverse. There is NO significance requirement for inclusion in figure 2. Zero is excluded. Within-study FDR < 0.05 is annotated with an asterisk for each tile.
4. Separate support-pattern sensitivity: positive significant enrichment in BOTH favored studies, and no positive significant enrichment in EITHER comparator study. Positive nonsignificant NES is permitted here. This must not be called an opposite-direction pattern.

The 3-of-4 recurrence criterion has no combined FDR or formal between-model interaction P value.

## Observed shared-positive results

Cells below are NES (within-study FDR).

| Hallmark | OIPN GSE160543 | OIPN GSE126773 | NC GSE246156 | CCI GSE212311 | Significant studies |
| --- | ---: | ---: | ---: | ---: | ---: |
| E2F_TARGETS | +2.77 (0.000683) | +1.74 (0.000403) | +3.26 (0.00138) | +1.09 (0.395) | 3 |
| G2M_CHECKPOINT | +2.90 (0.000683) | +1.04 (0.602) | +2.65 (0.00138) | +1.48 (0.0259) | 3 |
| IL6_JAK_STAT3_SIGNALING | +1.61 (0.0123) | +1.29 (0.259) | +2.32 (0.00138) | +2.45 (0.00149) | 3 |
| INTERFERON_ALPHA_RESPONSE | +2.22 (0.000683) | +1.66 (0.00956) | +1.70 (0.00333) | +1.81 (0.00809) | 4 |
| INTERFERON_GAMMA_RESPONSE | +2.19 (0.000683) | +1.29 (0.139) | +2.19 (0.00138) | +2.35 (0.00149) | 3 |
| TNFA_SIGNALING_VIA_NFKB | +1.69 (0.000683) | +0.78 (0.947) | +2.57 (0.00138) | +2.47 (0.00149) | 3 |

Six sets meet the shared-positive rule; only INTERFERON_ALPHA_RESPONSE reaches FDR < 0.05 in all four. E2F is nonsignificant in CCI; G2M, IL6/JAK/STAT3, IFN-gamma and TNFA/NFKB are nonsignificant in the microarray OIPN study, while retaining positive NES.
No set meets the symmetric shared-negative rule. This means no qualifying consistently negative set across these four contrasts; it does not mean negative enrichment is absent in individual studies.
The immune-related sets and E2F/G2M suggest recurrent immune-response and cell-cycle-associated transcriptional enrichment at tissue level; they do not demonstrate neuronal division, cell-cycle re-entry, senescence or causation of pain.

## Observed divergent results

The direction-only rule selects TWO sets: HEME_METABOLISM (positive in both OIPN, negative in NC/CCI) and COMPLEMENT (the reverse). The earlier implementation additionally required both positive-side studies to be significant; this omitted HEME_METABOLISM and has been corrected to match the author's direction-based definition.

| Hallmark | OIPN GSE160543 NES (FDR) | OIPN GSE126773 NES (FDR) | NC NES (FDR) | CCI NES (FDR) |
| --- | ---: | ---: | ---: | ---: |
| HEME_METABOLISM | +1.434 (0.0154) | +0.797 (0.9465) | -1.395 (0.0280) | -1.239 (0.1408) |
| COMPLEMENT | -1.276 (0.0724) | -0.891 (0.8897) | +1.949 (0.00138) | +1.656 (0.00676) |

HEME_METABOLISM is significant only in GSE160543 and NC. COMPLEMENT is significant only in NC and CCI. Neither set is significant in all four studies. The negative OIPN COMPLEMENT results do not establish repression in OIPN. These are descriptive direction differences, not a statistically tested difference between models.
The separate support rule still selects no OIPN-favoring set; it intentionally retains its original significance requirement.
The separate support rule selects four physical-injury-favoring sets: COMPLEMENT, IL2_STAT5_SIGNALING, INFLAMMATORY_RESPONSE and KRAS_SIGNALING_UP. The last three have positive NES in GSE160543; INFLAMMATORY_RESPONSE is positive in both OIPN studies. They cannot be described as absent positive direction in OIPN.
Do not label either category biologically model-specific: significance in one contrast and nonsignificance in another is not a test of their difference.

OXIDATIVE_PHOSPHORYLATION is negative/significant in the three RNA-seq studies but positive/significant in GSE126773. FATTY_ACID_METABOLISM has the same direction disagreement, and EPITHELIAL_MESENCHYMAL_TRANSITION / MITOTIC_SPINDLE reverse in the microarray. Retain these findings in the full supplementary comparison and discuss cross-study heterogeneity.

## Figure specifications and legend notes

Supplementary: existing `Hallmark_four_dataset_heatmap` displays all 50 sets, including discordance.
Figure 1: the six shared-positive sets. Tiles show original signed NES, without row standardization; asterisk means exported within-study FDR < 0.05. Fixed columns: OIPN GSE160543, OIPN GSE126773, NC GSE246156, CCI GSE212311. Selection requires all four NES positive and at least three significant studies.
Figure 2 default: direction-only selection includes HEME_METABOLISM and COMPLEMENT in separate sign-pattern panels. Legend: selected for opposite NES signs in both studies per side, independent of significance; asterisk denotes within-study FDR < 0.05. Sign discordance is descriptive and does not establish model specificity. No row scaling is applied. The full figure-2 table is `results/manuscript_hallmark/Fig2_direction_selection.csv`.
Figure 2 optional support mode: four physical-injury-favoring sets, defined by positive significant support rather than opposite direction. Label it accordingly; no OIPN-favoring set qualifies.
Both outputs use the same symmetric color scale derived from all 50 sets, numerical NES labels, and PNG (600 dpi) plus vector PDF. Clustering is disabled so the mechanistic grouping/order remains explicit.
Empty selections produce header-only CSVs and an explanatory text file, not fictitious heatmap rows.

## Interpretation limits and next scientific checks

NES is relative enrichment of a gene set in a ranking, not direct biochemical pathway activity. Cross-platform NES differences are descriptive rather than calibrated differences of effect size.
Bulk DRG includes neurons, glia and immune/stromal cells; expression shifts and cell-composition changes cannot be separated here. Hallmark names do not imply literal allograft rejection, epithelial transition or a cell-cycle event in neurons.
Study time, DRG level, platform, sample size and mapping differ. Harmonized signed-ranking GSEA with locked ortholog sets/backend is a useful sensitivity analysis before making stronger cross-study claims.
Leading-edge genes may differ even for the same positively enriched set; a subsequent gene-level intersection/audit would strengthen the claim of a shared program. Overlapping Hallmark membership also means the six sets are not six independent mechanisms.
Existing OIPN missing-gene exclusions, CCI stringent MSTRG mapping and NC count-composition shifts remain relevant; this selection does not repair them.

## Execution

From repository root after pulling the branch:

```r
source("scripts/cross_model/03_manuscript_Hallmark_heatmaps.R")
```

This creates exactly two heatmaps by default, plus all classification and audit CSVs. It does not rerun GSEA or overwrite earlier comparisons.
For the support definition instead of direction:

```r
Sys.setenv(HALLMARK_DIVERGENCE_MODE = "support")
source("scripts/cross_model/03_manuscript_Hallmark_heatmaps.R")
Sys.unsetenv("HALLMARK_DIVERGENCE_MODE")
```

For table extraction only, use `Sys.setenv(HALLMARK_TABLES_ONLY = "true")` before sourcing, then unset it.
Dependency for drawing: ggplot2. CSV classification itself uses base R.

## Validation status

The archived classifications were independently recomputed from all four source CSVs with Python; every input NES/FDR matches the existing four-dataset comparison to 1e-12. No figure was generated in this task, as requested.
The author reported successful execution of the previous version. R and ggplot2 are unavailable in this execution environment, so the updated R script has not been executed here. The direction selection was independently checked against all 50 original study-level CSV entries; significance markers are verified separately. The six shared-positive selections and the optional support rules are unchanged. No figures were generated in this update.

## Methodological references

- GSEA documentation: https://docs.gsea-msigdb.org/GSEA/GSEA_User_Guide/
- GSEA FAQ: https://docs.gsea-msigdb.org/GSEA/GSEA_FAQ/
- Gelman & Stern (2006), significance patterns are not tests of differences: https://doi.org/10.1198/000313006X152649
