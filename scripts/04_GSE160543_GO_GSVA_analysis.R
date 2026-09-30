#!/usr/bin/env Rscript
# All 50 scores and BH50 fits; historical score/DE tables remain untouched.
source("scripts/cross_model/OIPN_reproducibility_helpers.R")
source("scripts/cross_model/GSVA_helpers.R")
for (pkg in c("GSVA", "limma", "AnnotationDbi", "org.Rn.eg.db", "clusterProfiler"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install package: ", pkg)
paths <- oipn_rerun_paths()
de_file <- file.path(paths$primary, "DE_all_genes.csv")
expr_file <- file.path(paths$primary, "vst_expression_matrix.csv")
if (!file.exists(expr_file)) stop("Missing VST matrix; run scripts/cross_model/09_OIPN_reproducibility_rerun.R to regenerate it in a separate directory.")
res <- read.csv(de_file, check.names = FALSE, stringsAsFactors = FALSE)
expr <- as.matrix(read.csv(expr_file, row.names = 1, check.names = FALSE))
storage.mode(expr) <- "numeric"
expected <- c(paste0("GSM487500", 3:6), paste0("GSM487501", 1:4))
if (!setequal(colnames(expr), expected) || anyDuplicated(rownames(expr)) || any(!is.finite(expr)))
  stop("Invalid VST sample IDs, features or values.")
expr <- expr[, expected, drop = FALSE]
# Keep historical Entrez-to-symbol/first-feature GSVA policy; archive the mapping.
symbols <- AnnotationDbi::mapIds(org.Rn.eg.db::org.Rn.eg.db, keys = rownames(expr),
  column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
keep <- !is.na(symbols) & nzchar(symbols)
keep[keep] <- !duplicated(unname(symbols[keep]))
write.csv(data.frame(gene_id = rownames(expr), symbol = unname(symbols), selected = keep),
  file.path(paths$out, "GSVA_feature_to_symbol_audit.csv"), row.names = FALSE, na = "")
expr_symbol <- expr[keep, , drop = FALSE]
rownames(expr_symbol) <- unname(symbols[keep])
hallmark <- oipn_locked_hallmark()
set.seed(160543L)
if ("gsvaParam" %in% getNamespaceExports("GSVA")) {
  scores <- GSVA::gsva(GSVA::gsvaParam(exprData = expr_symbol, geneSets = hallmark))
  engine <- "gsvaParam_default_parameters"
} else {
  scores <- GSVA::gsva(expr_symbol, hallmark, method = "gsva", kcdf = "Gaussian", parallel.sz = 1)
  engine <- "legacy_gsva_Gaussian"
}
scores <- gsva_check_scores(scores, expected)
if (nrow(scores) != 50L) stop("Expected 50 GSVA pathways.")
write.csv(data.frame(pathway = rownames(scores), scores, check.names = FALSE),
  file.path(paths$out, "GSVA_Hallmark_scores.csv"), row.names = FALSE)
metadata <- data.frame(sample = expected, condition = rep(c("Control", "Neuropathy"), each = 4L))
write.csv(gsva_fit_scores(scores, metadata), file.path(paths$out, "GSVA_limma_all_50.csv"), row.names = FALSE)
comparison <- oipn_compare_gsva(scores, file.path(paths$root, "Pathway_analysis", "GSVA_Hallmark_scores.csv"), paths$out)
# Preserve historical GO settings; a tested-gene background needs a separate sensitivity run.
ids <- clusterProfiler::bitr(unique(res$symbol[!is.na(res$symbol) & nzchar(res$symbol)]),
  fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Rn.eg.db::org.Rn.eg.db)
up <- unique(ids$ENTREZID[ids$SYMBOL %in% res$symbol[!is.na(res$padj) & res$padj < .05 & res$log2FoldChange > 0]])
go <- if (length(up)) as.data.frame(clusterProfiler::enrichGO(gene = up,
  OrgDb = org.Rn.eg.db::org.Rn.eg.db, ont = "BP", pAdjustMethod = "BH", readable = TRUE)) else data.frame()
write.csv(go, file.path(paths$out, "GO_BP_upregulated.csv"), row.names = FALSE)
oipn_write_provenance(paths$out, c(de_file, expr_file, "scripts/04_GSE160543_GO_GSVA_analysis.R",
  "scripts/cross_model/GSVA_helpers.R", "scripts/cross_model/OIPN_reproducibility_helpers.R"),
  data.frame(parameter = c("seed", "engine", "duplicate_symbol_policy", "statistics", "GO_background"),
    value = c(160543, engine, "historical_first_mapped_Entrez_feature", "limma_BH_over_50", "historical_package_default")), "GSVA_GO")
message("OIPN GSVA/GO rerun and archived comparison: ", paths$out)
