#!/usr/bin/env Rscript

# GSE246156: NC L5 DRG, day 7, Compression versus Sham.
# Run from the repository root after scripts/GSE246156/02_DESeq2_GSE246156.R:
# source("scripts/GSE246156/03_GSEA_GO_GSE246156.R")
# Positive DESeq2 statistic / NES means higher expression after compression.

options(stringsAsFactors = FALSE)

required_packages <- c("AnnotationDbi", "org.Rn.eg.db", "fgsea", "msigdbr",
                       "clusterProfiler", "enrichplot", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                             logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install these R packages before running the analysis: ",
       paste(missing_packages, collapse = ", "))
}

input_file <- file.path("results", "GSE246156_NC_L5_day7", "DESeq2",
                        "GSE246156_Compression_vs_Sham_DE_all_genes.csv")
outdir <- file.path("results", "GSE246156_NC_L5_day7", "pathway_analysis")
if (!file.exists(input_file)) {
  stop("Missing DESeq2 output. Run step 02 from the repository root: ", input_file)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

de <- read.csv(input_file, check.names = FALSE)
required_columns <- c("gene", "stat", "log2FoldChange", "padj")
if (!all(required_columns %in% names(de))) {
  stop("DESeq2 output lacks: ",
       paste(setdiff(required_columns, names(de)), collapse = ", "))
}
if (nrow(de) < 1000L || anyNA(de$gene) || any(!nzchar(de$gene)) ||
    anyDuplicated(de$gene)) {
  stop("Unexpected or duplicated gene IDs in DESeq2 output.")
}
de$ensembl <- sub("\\.[0-9]+$", "", de$gene)
if (anyDuplicated(de$ensembl) || !all(grepl("^ENSRNOG[0-9]+$", de$ensembl))) {
  stop("Gene IDs are not unique rat Ensembl gene IDs after version removal.")
}

# AnnotationDbi::select exposes one-to-many mappings. Exclude ambiguous
# Ensembl-to-symbol assignments, then choose one Ensembl row per symbol by
# largest absolute signed DESeq2 statistic (ties: Ensembl ID alphabetically).
mapped <- suppressMessages(AnnotationDbi::select(
  org.Rn.eg.db::org.Rn.eg.db, keys = unique(de$ensembl),
  keytype = "ENSEMBL", columns = c("SYMBOL", "ENTREZID")
))
mapped <- unique(mapped[, c("ENSEMBL", "SYMBOL", "ENTREZID")])
write.csv(mapped, file.path(outdir, "Ensembl_annotation_all_mappings.csv"),
          row.names = FALSE, na = "")

symbol_pairs <- unique(mapped[!is.na(mapped$SYMBOL) & nzchar(mapped$SYMBOL),
                              c("ENSEMBL", "SYMBOL")])
symbol_count <- table(symbol_pairs$ENSEMBL)
ambiguous_symbol_ids <- names(symbol_count)[symbol_count > 1L]
symbol_pairs <- symbol_pairs[!symbol_pairs$ENSEMBL %in% ambiguous_symbol_ids,
                             , drop = FALSE]
rank_rows <- merge(de[, c("gene", "ensembl", "stat")], symbol_pairs,
                   by.x = "ensembl", by.y = "ENSEMBL", all = FALSE,
                   sort = FALSE)
rank_rows <- rank_rows[is.finite(rank_rows$stat), , drop = FALSE]
rank_rows <- rank_rows[order(-abs(rank_rows$stat), rank_rows$ensembl),
                       , drop = FALSE]
duplicate_symbol_rows <- sum(duplicated(rank_rows$SYMBOL))
rank_rows <- rank_rows[!duplicated(rank_rows$SYMBOL), , drop = FALSE]
if (nrow(rank_rows) < 5000L || anyDuplicated(rank_rows$SYMBOL)) {
  stop("Too few unique mapped rat symbols for Hallmark GSEA: ", nrow(rank_rows))
}
ranks <- setNames(rank_rows$stat, rank_rows$SYMBOL)
ranks <- sort(ranks, decreasing = TRUE)
write.csv(data.frame(symbol = names(ranks), statistic = unname(ranks)),
          file.path(outdir, "ranked_statistics.csv"), row.names = FALSE)

audit <- data.frame(
  measured_ensembl_ids = nrow(de),
  finite_DESeq2_statistics = sum(is.finite(de$stat)),
  ensembl_ids_with_ambiguous_symbols = length(ambiguous_symbol_ids),
  duplicate_symbol_rows_discarded = duplicate_symbol_rows,
  unique_symbols_ranked = length(ranks)
)
write.csv(audit, file.path(outdir, "mapping_audit.csv"), row.names = FALSE)

hallmark_table <- msigdbr::msigdbr(species = "Rattus norvegicus",
                                   collection = "H")
if (!all(c("gs_name", "gene_symbol") %in% names(hallmark_table))) {
  stop("Unexpected msigdbr Hallmark table format.")
}
hallmark <- split(hallmark_table$gene_symbol, hallmark_table$gs_name)
hallmark <- lapply(hallmark, unique)
hallmark <- lapply(hallmark, function(x) intersect(x, names(ranks)))
hallmark <- hallmark[lengths(hallmark) >= 15L & lengths(hallmark) <= 500L]
if (length(hallmark) < 20L) {
  stop("Too few Hallmark sets overlap the ranked rat symbols: ", length(hallmark))
}

# Match the GSE160543 analysis: signed DESeq2 stat, 15-500 genes, 10,000
# permutations. Each dataset is analyzed independently; never pool counts.
set.seed(246156)
fg <- fgsea::fgsea(pathways = hallmark, stats = ranks,
                   minSize = 15L, maxSize = 500L, nperm = 10000L)
fg <- as.data.frame(fg)
if (!nrow(fg) || !all(c("pathway", "NES", "padj", "leadingEdge") %in% names(fg))) {
  stop("fgsea did not return the expected Hallmark result columns.")
}
fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";", FUN.VALUE = "")
fg <- fg[order(is.na(fg$padj), fg$padj, fg$pathway), , drop = FALSE]
write.csv(fg, file.path(outdir, "GSEA_Hallmark_results.csv"), row.names = FALSE)

selected_pathways <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_P53_PATHWAY",
  "HALLMARK_DNA_REPAIR", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_G2M_CHECKPOINT", "HALLMARK_E2F_TARGETS",
  "HALLMARK_MITOTIC_SPINDLE", "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
)
selected <- fg[fg$pathway %in% selected_pathways, , drop = FALSE]
write.csv(selected, file.path(outdir, "selected_biological_programs.csv"),
          row.names = FALSE)
if (nrow(selected)) {
  plot_data <- selected[is.finite(selected$NES), , drop = FALSE]
  if (nrow(plot_data)) {
    plot_data$pathway <- factor(plot_data$pathway,
                                levels = plot_data$pathway[order(plot_data$NES)])
    plot <- ggplot2::ggplot(plot_data,
                            ggplot2::aes(x = pathway, y = NES, fill = NES > 0)) +
      ggplot2::geom_col() + ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = c("#2878A8", "#D26038"),
                                guide = "none") +
      ggplot2::theme_classic() +
      ggplot2::labs(title = "GSE246156 L5 day 7: selected Hallmark programs",
                    x = NULL, y = "Normalized enrichment score")
    ggplot2::ggsave(file.path(outdir, "selected_programs_NES.png"),
                    plot, width = 9, height = 6, dpi = 300)
  }
}

# GO ORA uses the measurable and Entrez-mappable genes as its universe.
# Restrict to unambiguous Ensembl-to-Entrez assignments. A GO hit is an
# enrichment among DE genes, not evidence of mechanism or cell of origin.
entrez_pairs <- unique(mapped[!is.na(mapped$ENTREZID) &
                                nzchar(mapped$ENTREZID),
                              c("ENSEMBL", "ENTREZID")])
entrez_count <- table(entrez_pairs$ENSEMBL)
ambiguous_entrez_ids <- names(entrez_count)[entrez_count > 1L]
entrez_pairs <- entrez_pairs[!entrez_pairs$ENSEMBL %in% ambiguous_entrez_ids,
                             , drop = FALSE]
entrez_ids <- function(ensembl_ids) {
  unique(entrez_pairs$ENTREZID[entrez_pairs$ENSEMBL %in% ensembl_ids])
}
universe <- entrez_ids(de$ensembl[is.finite(de$stat)])
up <- entrez_ids(de$ensembl[!is.na(de$padj) & de$padj < 0.05 &
                              is.finite(de$log2FoldChange) & de$log2FoldChange > 0])
down <- entrez_ids(de$ensembl[!is.na(de$padj) & de$padj < 0.05 &
                                is.finite(de$log2FoldChange) & de$log2FoldChange < 0])
if (length(universe) < 5000L) {
  stop("Too few Entrez-mapped tested genes for GO: ", length(universe))
}
up <- intersect(up, universe)
down <- intersect(down, universe)

run_go <- function(genes, stem) {
  if (length(genes) < 10L) {
    write.csv(data.frame(), file.path(outdir, paste0(stem, ".csv")),
              row.names = FALSE)
    return(invisible(NULL))
  }
  result <- clusterProfiler::enrichGO(
    gene = genes, universe = universe, OrgDb = org.Rn.eg.db::org.Rn.eg.db,
    keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
    pvalueCutoff = 0.05, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500,
    readable = TRUE
  )
  tab <- as.data.frame(result)
  write.csv(tab, file.path(outdir, paste0(stem, ".csv")), row.names = FALSE)
  if (nrow(tab)) {
    plot <- enrichplot::dotplot(result, showCategory = min(20L, nrow(tab))) +
      ggplot2::ggtitle(paste("GSE246156", stem))
    ggplot2::ggsave(file.path(outdir, paste0(stem, "_dotplot.png")),
                    plot, width = 9, height = 7, dpi = 300)
  }
  invisible(tab)
}
run_go(up, "GO_BP_upregulated")
run_go(down, "GO_BP_downregulated")

writeLines(c(
  "GSE246156 NC L5 DRG day 7: Compression versus Sham (3 versus 3).",
  "Positive stat/NES: greater expression after compression.",
  "GSEA uses signed DESeq2 Wald statistic and rat Hallmark gene sets.",
  "GO uses FDR < 0.05 DE genes and tested Entrez-mappable genes as universe.",
  paste("Ranked symbols:", length(ranks)),
  paste("Entrez universe:", length(universe)),
  paste("Upregulated Entrez IDs:", length(up)),
  paste("Downregulated Entrez IDs:", length(down)),
  paste("Hallmark sets tested:", nrow(fg)),
  "Bulk DRG enrichment does not localize signals to neurons or establish causality."
), file.path(outdir, "analysis_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
message("GSE246156 Hallmark GSEA and GO BP completed. Results: ", outdir)
