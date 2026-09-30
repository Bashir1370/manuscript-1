#!/usr/bin/env Rscript
# Rerun into a separate directory using a frozen current Hallmark membership.
source("scripts/cross_model/OIPN_reproducibility_helpers.R")
for (pkg in c("fgsea", "ggplot2"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install package: ", pkg)
paths <- oipn_rerun_paths()
input <- file.path(paths$primary, "DE_all_genes.csv")
if (!file.exists(input)) stop("Missing primary DE table: ", input)
res <- read.csv(input, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("stat", "symbol", "gene_id") %in% names(res))) stop("Invalid DE schema.")
# Preserve the historical highest-signed-statistic representative and source tie order.
rank_table <- res[is.finite(res$stat) & !is.na(res$symbol) & nzchar(res$symbol), ]
rank_table <- rank_table[order(-rank_table$stat, method = "radix"), ]
rank_table <- rank_table[!duplicated(rank_table$symbol), ]
ranks <- stats::setNames(rank_table$stat, rank_table$symbol)
hallmark <- oipn_locked_hallmark()
set.seed(160543L)
fg <- as.data.frame(fgsea::fgseaSimple(pathways = hallmark, stats = ranks,
  minSize = 15, maxSize = 500, nperm = 10000, nproc = 1))
fg <- fg[order(fg$padj, fg$pathway), ]
if (nrow(fg) != 50L || anyDuplicated(fg$pathway)) stop("Not all 50 Hallmarks are testable; inspect the lock.")
fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";", character(1))
write.csv(data.frame(symbol = names(ranks), statistic = ranks),
  file.path(paths$out, "ranked_statistics.csv"), row.names = FALSE)
write.csv(rank_table, file.path(paths$out, "rank_feature_audit.csv"), row.names = FALSE)
write.csv(fg, file.path(paths$out, "GSEA_Hallmark_results.csv"), row.names = FALSE)
selected <- fg[fg$pathway %in% c("HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_P53_PATHWAY", "HALLMARK_DNA_REPAIR", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_G2M_CHECKPOINT", "HALLMARK_E2F_TARGETS", "HALLMARK_OXIDATIVE_PHOSPHORYLATION"), ]
write.csv(selected, file.path(paths$out, "selected_biological_programs.csv"), row.names = FALSE)
if (tolower(Sys.getenv("OIPN_TABLES_ONLY", "false")) != "true") {
  p <- ggplot2::ggplot(selected, ggplot2::aes(reorder(pathway, NES), NES)) +
    ggplot2::geom_col() + ggplot2::coord_flip() + ggplot2::theme_classic() +
    ggplot2::labs(title = "OIPN reproducibility rerun: selected programs", x = NULL, y = "NES")
  ggplot2::ggsave(file.path(paths$out, "selected_programs_NES.png"), p, width = 8, height = 5, dpi = 300)
}
comparison <- oipn_compare_gsea(fg, file.path(paths$root, "Pathway_analysis", "GSEA_Hallmark_results.csv"), paths$out)
oipn_write_provenance(paths$out, c(input, "scripts/03_GSE160543_pathway_analysis.R",
  "scripts/cross_model/OIPN_reproducibility_helpers.R"),
  data.frame(parameter = c("seed", "engine", "nperm", "nproc", "minSize", "maxSize", "representative"),
    value = c(160543, "fgseaSimple", 10000, 1, 15, 500, "highest_signed_statistic_source_tie_order")), "GSEA")
message("OIPN GSEA rerun and archived comparison: ", paths$out)
