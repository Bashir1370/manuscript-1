# Three-study neuropathy analysis

## Scope and preserved archive

Branch: `three-dataset-neuropathy-analysis`, copied from `CCI-GSE212311-analysis` at commit `4c1823502a2b702f7f510cac5ad071320328ef85` before removing the fourth study and its mixed-study derived outputs from the new branch's working tree. Original study inputs and earlier analyses for the retained studies are preserved. The original branch remains the complete four-study archive.

The author requested excluding GSE126773 after observing marked divergence in pathway/sample profiles. This is a **post hoc restriction**, not an exclusion justified by established sample mislabeling or technical failure. No labels were changed. This branch asks about convergence among the three included RNA-seq cohorts; it cannot establish consistency across all OIPN studies. Retain the [four-study archive](https://github.com/Bashir1370/manuscript-1/tree/CCI-GSE212311-analysis) as the record of the excluded study and discordance when reporting the restricted analysis.

| Study | Independent contrast | Samples |
|---|---|---:|
| OIPN GSE160543 | Oxaliplatin minus Vehicle | 4 + 4 |
| NC GSE246156 | Compression minus Sham, L5 day 7 | 3 + 3 |
| CCI GSE212311 | CCI minus Sham, ipsilateral L4–L6 day 11 | 3 + 3 |

The central question remains whether DRG transcriptional programs converge across these neuropathy models and where OIPN lies within that spectrum. No cohort pooling, cross-study batch correction, sample removal or group relabeling is performed. All original per-study GSEA/GSVA scoring parameters and gene/probe representatives for the retained studies are preserved. Their original per-study outputs are reused: deleting an independently analyzed study does not change another study's fits or ranks. Cross-study classifications, display selections and leading-edge intersections are recalculated.

## GSEA rules and findings

Canonical script: `scripts/cross_model/03_manuscript_Hallmark_heatmaps.R`.

Analyze all 50 Hallmarks using original NES and within-study FDR. A shared positive/negative pathway requires the same strict NES sign in all three studies and FDR <0.05 in **at least three**, retaining the original threshold. With three included studies this means significance in every study. The threshold has not been relaxed to two.

- Shared positive: **10** pathways.
- Shared negative: **2**, FATTY ACID METABOLISM and OXIDATIVE PHOSPHORYLATION. The earlier four-study statement of no shared negative enrichment does not apply here.
- Opposite-direction selection: **7** pathways, using strict opposite NES signs between the one OIPN study and both physical-injury studies; no FDR selection filter. These are APICAL JUNCTION, COMPLEMENT, HEME METABOLISM, MYC TARGETS V1, MYOGENESIS, NOTCH SIGNALING and REACTIVE OXYGEN SPECIES PATHWAY.

FDR is shown separately for every cell. Direction discordance is descriptive, not a formal treatment-by-model interaction or evidence of OIPN specificity. OIPN now has one study versus two physical-injury studies, so the replication on each side is unequal. The separate positive-support pattern remains available with `HALLMARK_DIVERGENCE_MODE=support`; nonsignificance is not evidence of no activity.

E2F TARGETS remains positive in the three studies but is excluded from the shared-significant selection because CCI FDR is approximately 0.3954. Keep it in the full 50-pathway table. Do not retain the previous six-pathway selection as if it were the result of the new selection rule.

Outputs: `results/manuscript_hallmark_three_dataset/`. Local execution produces the full-50 supplementary NES heatmap, shared-positive and shared-negative heatmaps, and an opposite-direction heatmap. NES is not row standardized.

## GSVA rules and findings

Canonical script: `scripts/cross_model/04_three_dataset_GSVA_comparison.R`; helper `GSVA_helpers.R`.

Keep all three existing complete 50-Hallmark score matrices. Analyze each study separately using limma with Control as reference, `eBayes(trend=FALSE, robust=FALSE)`, and BH over 50 pathways within that study. Report original-score Neuropathy minus Control differences, moderated intervals and p/FDR. Do not select genes or pathways to improve significance.

Display the union of the newly selected 10 shared-positive, two shared-negative and seven opposite-direction pathways: **19** pathways. Selection comes from GSEA on the same samples, so this is complementary evidence, not independent validation. Heatmap rows are standardized within each study only; dotplots display individual samples with separate panel y scales. Historical Wilcoxon outputs remain unchanged.

The archive contains 150 all-pathway study fits, 57 selected pathway–study comparisons and 380 selected sample–pathway records. The three original uploaded limma/BH fits are reused unchanged apart from display flags because they were fitted independently over the same 50 pathways and identical samples. Original score means/contrasts and BH50 values were independently checked. Running the new R script refits these same models and exports runtime provenance.

Outputs: `results/GSVA_three_dataset/`. Original scoring, feature-mapping and historical gene-set-release differences still limit calibrated cross-study magnitude comparisons.

## Leading-edge rules and findings

Canonical script: `scripts/cross_model/05_three_dataset_leading_edge.R`.

Use the new **10 shared-positive** pathways, not a hard-coded list of six. Extract their original `leadingEdge` lists from each retained GSEA. Reconstruct the exact original representative used for each symbol: highest signed statistic in GSE160543; unambiguous Ensembl mapping/largest absolute statistic with Ensembl-ID tie break in NC; original source-aware feature audit in CCI. Every reconstructed representative statistic matches the archived original rank.

Membership is counted for the **same pathway** across studies. The unchanged >=3 membership rule now means all three studies. Keep the full union of contributing genes, original feature IDs, signed statistics, log2FC, gene FDR, pathway NES/FDR and rank-availability flags. A gene absent from a GSEA rank is missing, not a membership zero and not proof of absent expression. The original gene-level DE FDR is retained, with no recalculation over selected genes.

| Shared-positive Hallmark | Leading-edge genes in all three |
|---|---:|
| ALLOGRAFT REJECTION | 20 |
| APOPTOSIS | 6 |
| COAGULATION | 4 |
| EPITHELIAL MESENCHYMAL TRANSITION | 8 |
| G2M CHECKPOINT | 19 |
| IL6 JAK STAT3 SIGNALING | 8 |
| INTERFERON ALPHA RESPONSE | 5 |
| INTERFERON GAMMA RESPONSE | 18 |
| P53 PATHWAY | 6 |
| TNFA SIGNALING VIA NFKB | 15 |

These are **109 pathway–gene memberships representing 78 distinct genes**, not 109 independent genes. All shared genes have positive log2FC in the three retained studies; individual gene FDR is recorded separately. Hallmark membership and positive enrichment do not demonstrate cell-cycle re-entry, epithelial transition, apoptosis or other causal mechanisms in a specific DRG cell type.

Outputs: `results/shared_Hallmark_leading_edge_three_dataset/`: 1,025 union pathway–gene rows and 3,075 study evidence rows; all-three membership tables, overlap summaries and distinct gene/pathway counts. Local R execution creates per-pathway log2FC heatmaps for the shared genes, with gene FDR asterisks and display clipping at ±2 only; raw numerical estimates remain unchanged.

## Local checkout and execution

A separate worktree preserves any generated files or uncommitted work in the existing local checkout:

```powershell
cd D:\manuscript_1\manuscript-1-four-dataset
git fetch origin
git worktree add ..\manuscript-1-three-dataset -b three-dataset-neuropathy-analysis origin/three-dataset-neuropathy-analysis
```

Then in R:

```r
setwd("D:/manuscript_1/manuscript-1-three-dataset")
Sys.unsetenv(c("HALLMARK_TABLES_ONLY", "HALLMARK_DIVERGENCE_MODE", "GSVA_TABLES_ONLY", "LE_TABLES_ONLY"))
source("scripts/cross_model/03_manuscript_Hallmark_heatmaps.R")
source("scripts/cross_model/04_three_dataset_GSVA_comparison.R")
source("scripts/cross_model/05_three_dataset_leading_edge.R")
```

Tables-only mode uses the respective `*_TABLES_ONLY=true` variables. GSEA and leading-edge table extraction use base R; plots need ggplot2; GSVA comparison needs limma. Existing retained-study scoring pipelines remain available if a full within-study rerun is needed. The compatibility cross-dataset entry point delegates to the new canonical GSEA script. Older three-model reports may use different criteria or Ensembl-only CCI; use the canonical scripts and new output directories above for this branch.

## Provenance and validation status

`analysis/three_dataset_input_provenance.csv` records original source blobs and source commit, including the prior per-study GSVA statistics copied from the source branch before its four-study output directory is removed here. All retained GSEA/rank/DE blobs match the source snapshot.

Three-study classifications and leading-edge tables were calculated from archived CSVs in JavaScript and independently checked in Python for membership, selection, sign, effect/FDR joins, missing-data handling, unique rows and overlap totals. GSVA contrasts were checked against original score matrices; BH50 was independently recalculated from archived p values; selected score z standardization was verified. All validations passed. R scripts were reviewed for references and lexical balance. R is unavailable in this environment: new R execution and figure rendering remain local. No new figures are generated or represented as reviewed.

## OIPN-required follow-up prioritization (2026-09-30)
See [three_dataset_gene_prioritization.md](three_dataset_gene_prioritization.md). Keeping the 78 all-three leading-edge genes and positive effects in all three, requiring gene FDR <0.05 in OIPN and at least one of NC/CCI selects 20 unique genes (8 strict all-three and 12 OIPN+NC). Five NC+CCI-only significant genes are excluded from this priority list. This post hoc descriptive selection retains the full 78-gene evidence and original within-study gene FDR.
