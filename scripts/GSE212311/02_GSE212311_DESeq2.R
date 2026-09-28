#!/usr/bin/env Rscript
# Run from the manuscript-1 root after reviewing step 01 QC:
# Rscript scripts/GSE212311/02_GSE212311_DESeq2.R

suppressPackageStartupMessages(library(DESeq2))
root <- file.path("results", "GSE212311_CCI_L4L6_day11")
if (!file.exists("README.md") || !file.exists(file.path(root, "raw_count_matrix.csv")) ||
    !file.exists(file.path(root, "PCA_6samples", "PCA_coordinates.csv"))) {
  stop("Run step 01 from the repository root, then review its PCA before step 02.")
}
out <- file.path(root, "DESeq2")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
meta <- read.csv(file.path(root, "locked_primary_samples.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
ids <- c("GSM6523751", "GSM6523752", "GSM6523753",
         "GSM6523748", "GSM6523749", "GSM6523750")
if (!identical(meta$gsm, ids) ||
    !identical(meta$condition, rep(c("Sham", "CCI"), each = 3L)) ||
    any(meta$day != 11L) || any(meta$tissue != "ipsilateral_L4-L6_DRG")) {
  stop("Sample metadata disagree with the locked six-sample contrast.")
}
counts <- read.csv(file.path(root, "raw_count_matrix.csv"), row.names = 1L,
                   check.names = FALSE)
if (!identical(colnames(counts), ids) || anyDuplicated(rownames(counts)) ||
    nrow(counts) != 43350L) stop("Count matrix differs from step 01 output.")
counts <- as.matrix(counts)
if (!is.numeric(counts) || anyNA(counts) || any(!is.finite(counts)) ||
    any(counts < 0) || any(counts != round(counts)) ||
    any(counts > .Machine$integer.max)) stop("Invalid integer count matrix.")
storage.mode(counts) <- "integer"
rownames(meta) <- meta$gsm
meta$condition <- factor(meta$condition, levels = c("Sham", "CCI"))

# The six columns are independent biological libraries, not repeated levels.
# MSTRG features remain in the all-gene table but need annotation before
# gene-symbol pathway analysis. Filtering precedes multiple testing.
keep <- rowSums(counts >= 10L) >= 3L
if (sum(keep) < 1000L) stop("Too few features after count filtering.")
dds <- DESeqDataSetFromMatrix(countData = counts[keep, , drop = FALSE],
                              colData = meta, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "CCI", "Sham"), alpha = 0.05)
tbl <- as.data.frame(res)
tbl$gene_id <- rownames(tbl)
tbl$gene_class <- ifelse(startsWith(tbl$gene_id, "ENSRNOG"),
                         "annotated_Ensembl", "MSTRG_novel_locus")
tbl <- tbl[, c("gene_id", "gene_class", "baseMean", "log2FoldChange",
               "lfcSE", "stat", "pvalue", "padj")]
tbl <- tbl[order(tbl$padj, na.last = TRUE), , drop = FALSE]
write.csv(tbl, file.path(out, "GSE212311_CCI_vs_Sham_DE_all_genes.csv"),
          row.names = FALSE)
sig <- tbl[!is.na(tbl$padj) & tbl$padj < 0.05, , drop = FALSE]
write.csv(sig, file.path(out, "GSE212311_significant_DEGs_FDR_lt_0.05.csv"),
          row.names = FALSE)

# Keep a gene-ID based pathway-ready ranking. The stable ENSRNOG IDs may be
# mapped to rat gene symbols in the subsequent enrichment stage.
ranked <- tbl[!is.na(tbl$stat) & tbl$gene_class == "annotated_Ensembl",
              c("gene_id", "stat", "log2FoldChange", "pvalue", "padj")]
write.csv(ranked, file.path(out, "ranked_annotated_genes.csv"),
          row.names = FALSE)
write.csv(as.matrix(counts(dds, normalized = TRUE)),
          file.path(out, "normalized_counts_filtered.csv"))

writeLines(c(
  "GSE212311: independent within-study contrast, CCI versus Sham.",
  "Three male rat ipsilateral L4-L6 DRG biological samples per condition; day 11.",
  "Input: GEO integer count columns; FPKM was not entered into DESeq2.",
  paste("Retained features (>=10 counts in >=3 samples):", sum(keep)),
  paste("Significant features at BH adjusted p < 0.05:", nrow(sig)),
  paste("Annotated ENSRNOG significant features:",
        sum(sig$gene_class == "annotated_Ensembl")),
  "A positive log2FoldChange means higher expression in CCI than Sham.",
  "No absolute fold-change threshold or sample exclusion was applied.",
  "A 3-vs-3 comparison has limited power; inspect the PCA and count audit."
), file.path(out, "DE_analysis_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "DESeq2_sessionInfo.txt"))
message("DESeq2 complete: ", sum(keep), " tested; ", nrow(sig),
        " with FDR < 0.05. Output: ", out)
