#!/usr/bin/env Rscript
# GSE212311: rat CCI versus Sham pathway analysis.
# Run from the repository root after steps 01-03:
# source("scripts/GSE212311/04_GSE212311_pathway_analysis.R")
# Positive Wald statistic / NES means higher expression in CCI.

required_packages <- c("AnnotationDbi", "org.Rn.eg.db", "fgsea", "msigdbr",
                       "clusterProfiler", "enrichplot", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                             logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install these R packages before pathway analysis: ",
       paste(missing_packages, collapse = ", "))
}

root <- file.path("results", "GSE212311_CCI_L4L6_day11")
input_file <- file.path(root, "DESeq2",
                        "GSE212311_CCI_vs_Sham_DE_all_genes.csv")
outdir <- file.path(root, "pathway_analysis")
if (!file.exists(input_file)) {
  stop("Missing all-gene DESeq2 table. Run step 02 from the repository root: ",
       input_file)
}
if (!file.exists(file.path(root, "PCA_6samples", "PCA_coordinates.csv"))) {
  stop("Missing step 01 PCA; inspect sample QC before pathway interpretation.")
}

de <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c("gene_id", "gene_class", "stat", "log2FoldChange", "padj")
if (!all(required_columns %in% names(de))) {
  stop("DESeq2 table lacks: ",
       paste(setdiff(required_columns, names(de)), collapse = ", "))
}
if (nrow(de) < 1000L || anyNA(de$gene_id) ||
    any(!nzchar(de$gene_id)) || anyDuplicated(de$gene_id)) {
  stop("Unexpected or duplicated feature IDs in the all-gene table.")
}
is_ensembl <- startsWith(de$gene_id, "ENSRNOG")
is_mstrg <- grepl("^MSTRG\\.[0-9]+$", de$gene_id)
if (anyNA(de$gene_class) || any(!is_ensembl & !is_mstrg) ||
    any(is_ensembl != (de$gene_class == "annotated_Ensembl")) ||
    any(!grepl("^ENSRNOG[0-9]+(\\.[0-9]+)?$",
               de$gene_id[is_ensembl]))) {
  stop("Feature IDs or gene_class differ from the locked GSE212311 input.")
}
# MSTRG loci lack stable rat Ensembl IDs; keep them in the DE table,
# but exclude them from gene-symbol gene set tests.
ens <- de[is_ensembl, , drop = FALSE]
ens$ensembl <- sub("\\.[0-9]+$", "", ens$gene_id)
if (anyDuplicated(ens$ensembl)) {
  stop("Duplicate Ensembl IDs after removing version suffixes.")
}
if (nrow(ens) < 1000L) {
  stop("Too few tested Ensembl genes for pathway analysis: ", nrow(ens))
}

mapped <- suppressMessages(AnnotationDbi::select(
  org.Rn.eg.db::org.Rn.eg.db, keys = ens$ensembl,
  keytype = "ENSEMBL", columns = c("SYMBOL", "ENTREZID")
))
mapped <- unique(mapped[, c("ENSEMBL", "SYMBOL", "ENTREZID"),
                        drop = FALSE])
if (!all(ens$ensembl %in% mapped$ENSEMBL)) {
  stop("Annotation database omitted requested Ensembl keys.")
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write.csv(mapped, file.path(outdir, "Ensembl_annotation_all_mappings.csv"),
          row.names = FALSE, na = "")

# Resolve each annotation field separately. An ID with several distinct
# symbols/Entrez IDs is excluded from that field's analysis.
unique_pairs <- function(field) {
  pairs <- unique(mapped[, c("ENSEMBL", field), drop = FALSE])
  pairs <- pairs[!is.na(pairs[[field]]) & nzchar(pairs[[field]]),
                 , drop = FALSE]
  ambiguous <- names(which(table(pairs$ENSEMBL) > 1L))
  list(pairs = pairs[!pairs$ENSEMBL %in% ambiguous, , drop = FALSE],
       ambiguous = ambiguous)
}
sym <- unique_pairs("SYMBOL")
ent <- unique_pairs("ENTREZID")

rank_rows <- merge(ens[, c("gene_id", "ensembl", "stat")], sym$pairs,
                   by.x = "ensembl", by.y = "ENSEMBL",
                   all = FALSE, sort = FALSE)
rank_rows <- rank_rows[is.finite(rank_rows$stat), , drop = FALSE]
rank_rows <- rank_rows[order(-abs(rank_rows$stat), rank_rows$ensembl),
                       , drop = FALSE]
duplicate_symbols <- sum(duplicated(rank_rows$SYMBOL))
rank_rows <- rank_rows[!duplicated(rank_rows$SYMBOL), , drop = FALSE]
if (nrow(rank_rows) < 1000L) {
  stop("Too few distinct rat symbols with finite Wald statistics: ",
       nrow(rank_rows))
}
ranks <- sort(setNames(rank_rows$stat, rank_rows$SYMBOL), decreasing = TRUE)
write.csv(data.frame(symbol = names(ranks), statistic = unname(ranks)),
          file.path(outdir, "ranked_statistics.csv"), row.names = FALSE)

hallmark_table <- msigdbr::msigdbr(species = "Rattus norvegicus",
                                   collection = "H")
if (!all(c("gs_name", "gene_symbol") %in% names(hallmark_table))) {
  stop("Unexpected msigdbr Hallmark table format.")
}
hallmark <- lapply(split(hallmark_table$gene_symbol,
                         hallmark_table$gs_name),
                   function(x) intersect(unique(x), names(ranks)))
hallmark <- hallmark[lengths(hallmark) >= 15L & lengths(hallmark) <= 500L]
if (!length(hallmark)) stop("No Hallmark sets overlap the ranked symbols.")
set.seed(212311)
fg <- as.data.frame(fgsea::fgsea(
  pathways = hallmark, stats = ranks, minSize = 15L, maxSize = 500L,
  nperm = 10000L
))
if (!nrow(fg) ||
    !all(c("pathway", "NES", "padj", "leadingEdge") %in% names(fg))) {
  stop("fgsea did not return the expected Hallmark results.")
}
fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";",
                         FUN.VALUE = "")
fg <- fg[order(is.na(fg$padj), fg$padj, fg$pathway), , drop = FALSE]
write.csv(fg, file.path(outdir, "GSEA_Hallmark_results.csv"),
          row.names = FALSE)

selected_names <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_P53_PATHWAY",
  "HALLMARK_DNA_REPAIR", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_G2M_CHECKPOINT", "HALLMARK_E2F_TARGETS",
  "HALLMARK_MITOTIC_SPINDLE", "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
)
selected <- fg[fg$pathway %in% selected_names, , drop = FALSE]
write.csv(selected, file.path(outdir, "selected_biological_programs.csv"),
          row.names = FALSE)
plot_data <- selected[is.finite(selected$NES), , drop = FALSE]
if (nrow(plot_data)) {
  plot_data$pathway <- factor(plot_data$pathway,
                              levels = plot_data$pathway[order(plot_data$NES)])
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = pathway, y = NES, fill = NES > 0
  )) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c("#2878A8", "#D26038"),
                               guide = "none") +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "GSE212311 CCI L4-L6 day 11: Hallmark programs",
                  x = NULL, y = "Normalized enrichment score")
  ggplot2::ggsave(file.path(outdir, "selected_programs_NES.png"),
                  plot, width = 9, height = 6, dpi = 300)
}

# GO ORA uses all tested, finite-statistic, Entrez-mappable Ensembl
# genes as the measured universe, with FDR < 0.05 for the target sets.
entrez_ids <- function(ensembl_ids) {
  unique(ent$pairs$ENTREZID[ent$pairs$ENSEMBL %in% ensembl_ids])
}
universe <- entrez_ids(ens$ensembl[is.finite(ens$stat)])
up <- entrez_ids(ens$ensembl[!is.na(ens$padj) & ens$padj < 0.05 &
                              is.finite(ens$log2FoldChange) &
                              ens$log2FoldChange > 0])
down <- entrez_ids(ens$ensembl[!is.na(ens$padj) & ens$padj < 0.05 &
                                is.finite(ens$log2FoldChange) &
                                ens$log2FoldChange < 0])
up <- intersect(up, universe)
down <- intersect(down, universe)
if (length(universe) < 1000L) {
  stop("Too few Entrez IDs in the measured GO universe: ", length(universe))
}
run_go <- function(genes, stem) {
  path <- file.path(outdir, paste0(stem, ".csv"))
  if (length(genes) < 10L) {
    write.csv(data.frame(), path, row.names = FALSE)
    return(invisible(NULL))
  }
  result <- clusterProfiler::enrichGO(
    gene = genes, universe = universe, OrgDb = org.Rn.eg.db::org.Rn.eg.db,
    keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
    pvalueCutoff = 0.05, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500,
    readable = TRUE
  )
  tab <- as.data.frame(result)
  write.csv(tab, path, row.names = FALSE)
  if (nrow(tab)) {
    plot <- enrichplot::dotplot(result, showCategory = min(20L, nrow(tab))) +
      ggplot2::ggtitle(paste("GSE212311", stem))
    ggplot2::ggsave(file.path(outdir, paste0(stem, "_dotplot.png")),
                    plot, width = 9, height = 7, dpi = 300)
  }
  invisible(tab)
}
run_go(up, "GO_BP_upregulated")
run_go(down, "GO_BP_downregulated")

write.csv(data.frame(
  tested_features = nrow(de),
  tested_ensembl_features = nrow(ens),
  excluded_MSTRG_features = sum(is_mstrg),
  ambiguous_symbol_ids = length(sym$ambiguous),
  ambiguous_entrez_ids = length(ent$ambiguous),
  duplicate_symbol_rows_discarded = duplicate_symbols,
  ranked_unique_symbols = length(ranks),
  measured_entrez_universe = length(universe),
  significant_up_entrez_ids = length(up),
  significant_down_entrez_ids = length(down)
), file.path(outdir, "mapping_audit.csv"), row.names = FALSE)
writeLines(c(
  "GSE212311 rat ipsilateral L4-L6 DRG day 11: CCI versus Sham (3 versus 3).",
  "Positive Wald statistic and NES mean higher expression in CCI.",
  "Hallmark GSEA uses the full signed DESeq2 rank, not only significant DEGs.",
  "GO BP uses FDR < 0.05 Ensembl DEGs and the tested Entrez universe.",
  "MSTRG features remain in the DE table but cannot enter symbol gene sets.",
  paste("Ranked symbols:", length(ranks)),
  paste("GO universe:", length(universe)),
  paste("Up/down Entrez IDs:", length(up), "/", length(down)),
  paste("Hallmark sets tested:", nrow(fg)),
  "Inspect PCA and input QC before interpreting pathway results.",
  "Bulk DRG enrichment cannot establish cell type or causal mechanism."
), file.path(outdir, "analysis_summary.txt"))
writeLines(capture.output(sessionInfo()),
           file.path(outdir, "sessionInfo.txt"))
message("GSE212311 pathway analysis completed: ", outdir)
