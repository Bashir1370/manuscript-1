# GSE246156 NC L5 day7
# DESeq2 analysis: Compression vs Sham

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
})

project_dir <- "D:/manuscript_1/manuscript-1"
setwd(project_dir)

input_dir <- "results/GSE246156_NC_L5_day7"
outdir <- file.path(input_dir, "DESeq2")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

counts <- read.csv(
  file.path(input_dir, "raw_count_matrix.csv"),
  row.names = 1,
  check.names = FALSE
)

meta <- read.csv(
  file.path(input_dir, "locked_primary_samples.csv"),
  stringsAsFactors = FALSE
)

rownames(meta) <- meta$gsm

counts <- counts[, meta$gsm]

stopifnot(all(colnames(counts) == rownames(meta)))

meta$condition <- factor(
  meta$condition,
  levels = c("Sham", "Compression")
)

message("Samples:")
print(meta[, c("gsm", "condition")])

# DESeq2 object

dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData = meta,
  design = ~ condition
)

# filtering
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

message("Genes retained: ", nrow(dds))

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("condition", "Compression", "Sham"),
  alpha = 0.05
)

res <- as.data.frame(res)
res$gene <- rownames(res)
res <- res[, c("gene", setdiff(colnames(res), "gene"))]

write.csv(
  res,
  file.path(outdir, "GSE246156_Compression_vs_Sham_DE_all_genes.csv"),
  row.names = FALSE
)

sig <- res %>%
  filter(!is.na(padj), padj < 0.05)

write.csv(
  sig,
  file.path(outdir, "GSE246156_significant_DEGs.csv"),
  row.names = FALSE
)

# normalized counts
vst_mat <- assay(vst(dds, blind = FALSE))
write.csv(
  vst_mat,
  file.path(outdir, "GSE246156_vst_expression_matrix.csv")
)

# QC PCA
pca <- plotPCA(
  vst(dds, blind = FALSE),
  intgroup = "condition",
  returnData = TRUE
)

write.csv(
  pca,
  file.path(outdir, "PCA_coordinates.csv"),
  row.names = TRUE
)

p <- ggplot(pca, aes(PC1, PC2, label = name, color = condition)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8) +
  theme_classic()

ggsave(
  file.path(outdir, "PCA.png"),
  p,
  width = 7,
  height = 5,
  dpi = 300
)

writeLines(
  capture.output(sessionInfo()),
  file.path(outdir, "DESeq2_sessionInfo.txt")
)

message("GSE246156 DESeq2 analysis completed")
