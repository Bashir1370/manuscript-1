# GSE160543 downstream analysis
# GO Biological Process enrichment + GSVA pathway scoring + statistics + heatmap

library(data.table)
library(dplyr)
library(clusterProfiler)
library(org.Rn.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(GSVA)
library(msigdbr)
library(ggplot2)
library(pheatmap)

outdir <- "results/GSE160543_Oxaliplatin_vs_Vehicle/pathway_analysis"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# GO Biological Process analysis
# -----------------------------

res <- fread("results/GSE160543_Oxaliplatin_vs_Vehicle/DE_all_genes.csv")

ids <- bitr(
  res$symbol,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Rn.eg.db
)

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

ggsave(
  file.path(outdir, "GO_BP_dotplot.png"),
  dotplot(go_up, showCategory = 20) +
    ggtitle("GSE160543 Oxaliplatin: GO Biological Processes"),
  width = 8,
  height = 6,
  dpi = 300
)

# -----------------------------
# GSVA Hallmark analysis
# -----------------------------

expr_file <- "results/GSE160543_Oxaliplatin_vs_Vehicle/vst_expression_matrix.csv"

expr <- fread(expr_file) %>% as.data.frame()
rownames(expr) <- colnames(expr)[1] %>% {expr[[.]]}
expr[[1]] <- NULL
expr <- as.matrix(expr)
storage.mode(expr) <- "numeric"

# Convert rat Entrez IDs to gene symbols
symbols <- mapIds(
  org.Rn.eg.db,
  keys = rownames(expr),
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

expr_symbol <- expr[!is.na(symbols), ]
rownames(expr_symbol) <- symbols[!is.na(symbols)]
expr_symbol <- expr_symbol[!duplicated(rownames(expr_symbol)), ]

hallmark <- msigdbr(
  species = "Rattus norvegicus",
  collection = "H"
) %>%
  split(x = .$gene_symbol, f = .$gs_name)

param <- gsvaParam(
  exprData = expr_symbol,
  geneSets = hallmark
)

gsva_result <- gsva(param)

scores <- as.data.frame(gsva_result)
scores$pathway <- rownames(scores)

write.csv(
  scores,
  file.path(outdir, "GSVA_Hallmark_scores.csv"),
  row.names = FALSE
)

# -----------------------------
# Selected biological programs
# -----------------------------

selected_pathways <- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MITOTIC_SPINDLE",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"
)

selected <- scores %>%
  filter(pathway %in% selected_pathways)

write.csv(
  selected,
  file.path(outdir, "selected_GSVA_program_scores.csv"),
  row.names = FALSE
)

# -----------------------------
# Statistical comparison
# -----------------------------

sample_group <- ifelse(
  grepl("487500[3-6]", colnames(gsva_result)),
  "Vehicle",
  "Oxaliplatin"
)

stats <- data.frame()

for (p in selected_pathways) {
  x <- as.numeric(gsva_result[p, ])
  test <- wilcox.test(x ~ sample_group)

  stats <- rbind(
    stats,
    data.frame(
      pathway = p,
      Vehicle_mean = mean(x[sample_group == "Vehicle"]),
      Oxaliplatin_mean = mean(x[sample_group == "Oxaliplatin"]),
      p_value = test$p.value
    )
  )
}

stats$FDR <- p.adjust(stats$p_value, method = "BH")

write.csv(
  stats,
  file.path(outdir, "GSVA_program_statistics.csv"),
  row.names = FALSE
)

# -----------------------------
# Heatmap
# -----------------------------

heat_data <- gsva_result[selected_pathways, ]

pheatmap(
  heat_data,
  scale = "row",
  filename = file.path(outdir, "GSVA_selected_programs_heatmap.png"),
  width = 8,
  height = 6
)

writeLines(
  capture.output(sessionInfo()),
  file.path(outdir, "GO_GSVA_sessionInfo.txt")
)
