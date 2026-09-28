# GSE160543: primary Oxaliplatin versus Vehicle result audit

Date: 2026-09-28. These are results supplied by the author after running `scripts/02_GSE160543_oxaliplatin_primary.R` on Windows 11, R 4.5.2, DESeq2 1.50.2. The exact supplied output tables and figures are stored under [`results/GSE160543_Oxaliplatin_vs_Vehicle/`](../results/GSE160543_Oxaliplatin_vs_Vehicle/). This is an audit of the output, not an independent rerun of DESeq2.

## Input and quality checks

- The output identifies four Vehicle GSMs (GSM4875003–GSM4875006) and four Oxaliplatin GSMs (GSM4875011–GSM4875014). Paclitaxel is excluded.
- The downloaded GEO tar archive MD5 recorded by the run is `6bc2e7967dafe30561e574b5b5291103`.
- The union of input gene IDs contains 19,848 genes, but only 16,635 are found in every selected file. Of the 3,213 excluded IDs, 66 have a count of at least 100 in one or more files where present. This is a material completeness limitation; absence was not imputed as zero. After filtering by total rounded count ≥10, 16,632 genes enter DESeq2.
- PCA separates the groups along PC1 (32.4% variance), with within-group variation, especially GSM4875011 on PC1 and GSM4875014 on PC2. This supports a detectable tissue-level contrast but does not establish a specific pathway or cell of origin.
- Of 16,632 output rows, 972 have no adjusted p-value and 4 have no raw p-value. The reported FDR threshold is applied only to nonmissing adjusted p-values.

## Descriptive findings

| Result | Value |
| --- | ---: |
| Genes with FDR < 0.05 | 830 (418 positive, 412 negative log2 fold changes) |
| Genes with FDR < 0.05 and absolute log2FC > 1 | 125 (114 positive, 11 negative) |
| Cdkn1a | log2FC 2.366; FDR 1.35 × 10⁻³⁴ |
| Cdc20 | log2FC 1.210; FDR 7.34 × 10⁻⁵ |
| Cdk1 | log2FC 1.121; FDR 7.34 × 10⁻⁵ |

Mki67 and Top2a also have positive log2 fold changes with FDR < 0.05. This is a bulk DRG transcriptional pattern involving proliferation-associated genes. It does not localize the signal to post-mitotic neurons or establish cell-cycle activity. Alox15 is strongly induced (log2FC 4.366, FDR 2.29 × 10⁻⁷¹); its interpretation requires pathway and cellular-context analysis rather than treating one gene as a mechanism.

The prior independent rat qPCR observations should be reported as endpoint, bulk-tissue findings in that separate experiment. In this GEO contrast, Il1b is not significant (log2FC −0.543, FDR 0.594). Tnf has no adjusted p-value; Il6 is absent from two source files at very low counts and was excluded by the complete-case rule. Therefore, this GEO analysis does not replicate a coordinated increase in these three cytokine transcripts. The two experiments may differ in sampling or design; their raw protocols and time points should be compared before offering an explanation.

## Interpretation and next analysis

The current defensible result is **an oxaliplatin-associated bulk DRG expression response** in this dataset, including altered cell-cycle-related transcripts. It cannot establish cell type, neuronal cell-cycle re-entry, senescence, protein activity, causal mediation of pain, or oxaliplatin specificity among neuropathic injuries.

Next, analyze prespecified pathways with a full ranked gene list and explicit background/ID mapping; inspect neuronal, glial and immune reference signatures without treating deconvolution as localization. Separately define a matched GSE246156 compression-versus-sham contrast after resolving DRG level, time and animal independence. Compare effect sizes and directions across models, not only overlap of thresholded DEG lists. Consider re-quantification from raw sequencing reads if the missing high-count genes could affect key conclusions.
