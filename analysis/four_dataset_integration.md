# Four-dataset integration status

The `CCI-GSE212311-analysis` branch combines the previous three-model rat DRG
analysis with the scripts from `GSE126773-OIPN-validation`. Its four independent
contrasts are GSE160543 (oxaliplatin RNA-seq), GSE126773 (oxaliplatin
microarray), GSE246156 (nerve compression RNA-seq), and GSE212311 (CCI RNA-seq).
The central question remains whether bulk DRG shows a shared injury-associated
transcriptional program and where OIPN lies within that spectrum.

The GSE126773 Hallmark CSV was generated in the earlier validation workflow,
but is **not committed to this public branch**. The corresponding CEL files,
normalized expression, ranked statistics, and other local result files are also
not part of this merge. Reproduction from raw CEL files requires running the
`scripts/GSE126773_OIPN/` workflow and checking its sample manifest and audit.

Run `source("scripts/cross_model/02_four_dataset_Hallmark_comparison.R")` from
the repository root after generating or placing the local GSE126773 Hallmark
CSV at `results/GSE126773_OIPN/pathway_analysis/GSE126773_Hallmark_ranked_GSEA_results.csv`.
It requires `ggplot2`, reads four Hallmark tables, checks for identical
50-pathway sets, and writes a comparison table,
PNG heatmap, and summary into `results/cross_model_four_dataset/`. The older
`scripts/cross_dataset/04_cross_dataset_pathway_comparison.R` delegates to this
four-dataset script. The established three-model scripts and results remain
available for provenance.

GSE126773 uses a limma moderated-statistic ranking with clusterProfiler GSEA;
the other three use DESeq2 Wald ranks and fgsea. Gene universes and assay
platforms also differ. Compare enrichment direction and each study's FDR; a
numerical NES difference is not a formal test of model specificity. Differences
in direction should be retained and interpreted rather than hidden by a
claim of universal agreement.
