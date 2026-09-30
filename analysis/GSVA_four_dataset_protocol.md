# Four-study Hallmark GSVA protocol

## Purpose and scope

Complement the manuscript GSEA evidence with sample-level Hallmark GSVA in four independent rat bulk DRG contrasts. This analysis uses the same samples as GSEA and is not independent validation. Group-level changes describe bulk-tissue transcriptional programs; they do not establish neuronal specificity, cell composition, or causal mechanisms.

| Study | Locked contrast | Samples | GSVA input |
|---|---|---:|---|
| GSE160543 | Oxaliplatin minus Vehicle | 4 + 4 | Existing complete RNA-seq score matrix |
| GSE126773 | Control+OX minus Control | 3 + 3 | Existing log2 RMA matrix and probe-to-symbol annotation |
| GSE246156 | Compression minus Sham, L5 day 7 | 3 + 3 | Existing complete NC score matrix |
| GSE212311 | CCI minus Sham, ipsilateral L4–L6 day 11 | 3 + 3 | Existing source-aware CCI score matrix |

Paclitaxel samples are excluded from GSE160543. Studies are not pooled or batch-corrected together.

## Implementation

1. `scripts/GSE126773_OIPN/06_GSE126773_GSVA_Hallmark.R` replaces the earlier placeholder. It reads the existing RMA matrix and only the `probe_id`/`SYMBOL` columns of the annotation table. It excludes ambiguous/multiple-symbol probes, averages log2 RMA values across probes for a symbol, and removes constant genes. No differential-expression statistic or significance threshold selects probes. It uses human MSigDB Hallmarks projected to rat by `msigdbr`, requires all 50 sets to have at least 15 measured genes, and records membership, overlaps, mapping, input checksums and package versions. GSVA parameters are explicit: Gaussian kernel, minSize=1, maxSize=Inf, tau=1, maxDiff=TRUE, absRanking=FALSE; serial execution.
2. `scripts/cross_model/04_four_dataset_GSVA_comparison.R` reuses the three historical full score matrices and the new GSE126773 matrix. It validates locked GSM identifiers, group labels, all 50 pathway identities and matched GSEA tables before writing four-study results.
3. Both entry scripts use `scripts/cross_model/GSVA_helpers.R`. Keep this helper in its repository location.

For every study separately, fit `limma::lmFit` followed by `eBayes(trend=FALSE, robust=FALSE)` to the original, unscaled scores, with Control as reference. Test all 50 Hallmarks and apply Benjamini–Hochberg correction over all 50 in that study. The reported coefficient is **delta_GSVA = mean(Neuropathy) − mean(Control)**, not a log2 expression fold change. Export group means, the coefficient, moderated 95% confidence interval, moderated t, raw p, FDR, group sizes and complete-separation flag. Use FDR < 0.05 for exploratory statistical support; small sample sizes and model assumptions limit interpretation. The 95% intervals are individual model-based intervals, not simultaneous intervals across 50 pathways.

Historical Wilcoxon outputs remain unchanged. This is an explicitly documented alternative analysis applied to every pathway, not a switch selected because a particular pathway becomes significant. Report disagreements between methods if they affect manuscript claims. With 3 versus 3 observations, the smallest possible exact two-sided Wilcoxon p value without ties is 0.1.

## Display pathways and interpretation

Display the six shared positive GSEA sets (E2F TARGETS, G2M CHECKPOINT, IL6 JAK STAT3 SIGNALING, INTERFERON ALPHA RESPONSE, INTERFERON GAMMA RESPONSE, TNFA SIGNALING VIA NFKB) plus HEME METABOLISM and COMPLEMENT from the GSEA opposite-direction selection. These eight are selected using the same data and are exploratory displays; their GSVA significance is not a selection criterion. Retain the complete 50-pathway statistical tables as supplementary evidence.

- Heatmap: eight rows with four study panels; each pathway is standardized across samples **within that study only**. Sample ordering is locked Control then Neuropathy. Colors show relative sample scores, not FDR or cross-study calibrated effect size.
- Dotplots: all individual samples and group means, with a separate y scale for each pathway-study panel.
- GSVA–GSEA table: original GSVA difference and FDR alongside GSEA NES/FDR and direction agreement. A positive raw GSVA score alone does not imply increased activity relative to control.

The three RNA-seq score matrices are reused rather than recomputed. Assay gene universes, symbol mappings and historical gene-set releases can differ. The fourth score matrix uses a newly recorded gene-set release. This is therefore **not a harmonized rescoring of four datasets**. Interpret directions and within-study statistical support; do not rank study severity by raw score differences. Record this limitation in the manuscript.

## Run locally

From the repository root after pulling `CCI-GSE212311-analysis`:

```r
setwd("D:/manuscript_1/manuscript-1-four-dataset")
Sys.unsetenv("GSVA_TABLES_ONLY")
source("scripts/GSE126773_OIPN/06_GSE126773_GSVA_Hallmark.R")
source("scripts/cross_model/04_four_dataset_GSVA_comparison.R")
```

Required packages: GSVA, msigdbr, limma, BiocParallel and ggplot2. The scoring script uses the current parameter-object GSVA API and `msigdbr(collection="H", db_species="HS")`; older package versions using older APIs require updating. If needed:

```r
install.packages(c("BiocManager", "msigdbr", "ggplot2"))
BiocManager::install(c("GSVA", "limma", "BiocParallel"), ask = FALSE)
```

For tables only, set `Sys.setenv(GSVA_TABLES_ONLY="true")` before the comparison script. The environment setting does not suppress scoring in step 06.

Outputs are written to `results/GSE126773_OIPN/GSVA_Hallmark/` and `results/GSVA_four_dataset/`. The comparison exports 200 pathway-study rows, 32 selected pathway-study rows and 208 selected pathway-sample rows, plus audits, summaries and session information. Normal execution additionally writes PDF and PNG heatmap/dotplots. Preserve generated CSVs, session files and gene-set membership with the final results; rerunning can overwrite these new output locations.

## Input verification and execution status

Verified against branch commit `8eb6d79262a545bf3c084e6de4628bf925155b35`:

- All three historical matrices have 50 distinct Hallmarks, finite numeric values and exactly the locked 8/6/6 sample identifiers. Pathway names agree across them and their GSEA tables.
- All eight display pathways are present. Independently calculated original-score group differences agree with GSEA NES signs for all eight in each of these three studies. This is direction agreement only, not a significance finding.
- GSE126773 has 31,099 unique RMA probes and six finite expression columns matching the manifest. All annotation probe IDs match the matrix. The mapping rule retains 21,901 probes representing 14,658 unique symbols before constant-gene filtering, with no ambiguous probe IDs. Actual gene-set overlaps and scores will be established during R execution.
- Input and numerical alignment checks were performed independently in Python/JavaScript; scripts were reviewed for sample order, coefficient orientation, FDR family and output separation. R is unavailable in the development environment: GSVA, limma and figure rendering have **not been executed here**. No newly calculated GSVA statistics or figures are represented as validated results. The scripts stop on missing/misaligned inputs and export runtime provenance after successful execution.

Primary method references: [GSVA package](https://bioconductor.org/packages/GSVA/), [limma user guide](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf), [msigdbr documentation](https://igordot.github.io/msigdbr/).
