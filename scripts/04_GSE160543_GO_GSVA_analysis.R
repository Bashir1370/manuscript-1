# GSE160543 downstream analysis
# GO Biological Process enrichment + GSVA pathway scoring

library(data.table)
library(dplyr)
library(clusterProfiler)
library(org.Rn.eg.db)
library(enrichplot)
library(GSVA)
library(GSEABase)
library(msigdbr)
library(ggplot2)

outdir <- "results/GSE160543_Oxaliplatin_vs_Vehicle/pathway_analysis"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Differential expression table
res <- fread("results/GSE160543_Oxaliplatin_vs_Vehicle/DE_all_genes.csv")

# Convert symbols to Entrez IDs for GO
ids <- bitr(
  res$symbol,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Rn.eg.db
)

# Significant upregulated genes
up <- res %>%
  filter(log2FoldChange > 0, padj < 0.05) %>%
  inner_join(ids, by = c("symbol" = "SYMBOL"))

go_up <- enrichGO(
  gene = up$ENTREZID,
  OrgDb = org.Rn.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

write.csv(
  as.data.frame(go_up),
  file.path(outdir, "GO_BP_upregulated.csv"),
  row.names = FALSE
)

# Save GO plot
p <- dotplot(go_up, showCategory = 20) +
  ggtitle("GSE160543 Oxaliplatin: GO Biological Processes")

ggsave(
  file.path(outdir, "GO_BP_dotplot.png"),
  p,
  width = 8,
  height = 6,
  dpi = 300
)

# GSVA requires normalized expression matrix generated previously
# Expected input:
# results/GSE160543_Oxaliplatin_vs_Vehicle/vst_expression_matrix.csv

expr_file <- "results/GSE160543_Oxaliplatin_vs_Vehicle/vst_expression_matrix.csv"

if (file.exists(expr_file)) {
  expr <- fread(expr_file) %>% as.data.frame()
  rownames(expr) <- expr$symbol
  expr$symbol <- NULL

  pathways <- msigdbr(
    species = "Rattus norvegicus",
    collection = "H"
  ) %>%
    split(x = .$gene_symbol, f = .$gs_name)

  scores <- gsva(
    as.matrix(expr),
    pathways,
    method = "gsva"
  )

  write.csv(
    scores,
    file.path(outdir, "GSVA_Hallmark_scores.csv")
  )
}

writeLines(
  capture.output(sessionInfo()),
  file.path(outdir, "GO_GSVA_sessionInfo.txt")
)
