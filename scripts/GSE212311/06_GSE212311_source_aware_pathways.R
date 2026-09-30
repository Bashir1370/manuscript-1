#!/usr/bin/env Rscript
# GSE212311 source-aware rat gene pathway analysis.
# Run from the repository root after steps 01 and 02:
# source("scripts/GSE212311/06_GSE212311_source_aware_pathways.R")
# Positive DESeq2 Wald statistic / NES denotes higher expression in CCI.

required <- c("AnnotationDbi", "org.Rn.eg.db", "fgsea", "msigdbr",
              "clusterProfiler", "enrichplot", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1),
                             quietly = TRUE)]
if (length(missing)) stop("Install required R packages: ",
                          paste(missing, collapse = ", "))

root <- file.path("results", "GSE212311_CCI_L4L6_day11")
source_file <- file.path("data", "source-cache", "GSE212311",
                         "GSE212311_genes_fpkm_expression.txt.gz")
de_file <- file.path(root, "DESeq2",
                     "GSE212311_CCI_vs_Sham_DE_all_genes.csv")
out <- file.path(root, "pathway_analysis_source_aware")
if (!file.exists(source_file) || !file.exists(de_file)) {
  stop("Missing verified GEO source or step 02 DESeq2 output.")
}
md5 <- "71033155b073b1af600c51a20cc9b33a"
if (!identical(tolower(unname(tools::md5sum(source_file))), md5)) {
  stop("GEO source differs from the locked MD5.")
}
tab <- read.delim(gzfile(source_file), check.names = FALSE, quote = "",
                  comment.char = "", stringsAsFactors = FALSE,
                  na.strings = character())
de <- read.csv(de_file, check.names = FALSE, stringsAsFactors = FALSE)
needed_source <- c("gene_id", "gene_name", "trans_type", "Description")
needed_de <- c("gene_id", "gene_class", "stat", "padj", "log2FoldChange")
if (nrow(tab) != 43350L || anyDuplicated(tab$gene_id) ||
    !all(needed_source %in% names(tab)) ||
    !all(needed_de %in% names(de)) || anyDuplicated(de$gene_id) ||
    !all(de$gene_id %in% tab$gene_id)) {
  stop("GEO and DESeq2 feature inventories are inconsistent.")
}
if (anyNA(de$gene_class) ||
    any(startsWith(de$gene_id, "ENSRNOG") !=
        (de$gene_class == "annotated_Ensembl"))) {
  stop("Unexpected gene_class in DESeq2 table.")
}
source_row <- match(de$gene_id, tab$gene_id)
name <- trimws(as.character(tab$gene_name[source_row]))
type <- as.character(tab$trans_type[source_row])
description <- as.character(tab$Description[source_row])
is_mstrg <- grepl("^MSTRG\\.[0-9]+$", de$gene_id)
is_ens <- grepl("^ENSRNOG[0-9]+(\\.[0-9]+)?$", de$gene_id)
if (any(!is_mstrg & !is_ens)) stop("Unexpected feature ID format.")

# Admit only single-name, protein-coding MSTRG rows explicitly labeled
# with an RGD Symbol accession in the verified GEO source. All duplicate
# source gene names are excluded, even if only one copy passed DE filtering.
source_name <- trimws(as.character(tab$gene_name))
source_eligible <- grepl("^MSTRG\\.[0-9]+$", tab$gene_id) &
  tab$trans_type == "protein_coding" &
  grepl("[Source:RGD Symbol;Acc:", tab$Description, fixed = TRUE) &
  !is.na(source_name) &
  grepl("^[A-Za-z][A-Za-z0-9._-]*$", source_name)
eligible_names <- source_name[source_eligible]
duplicate_source_names <- unique(eligible_names[
  duplicated(eligible_names) | duplicated(eligible_names, fromLast = TRUE)
])
mstrg_ok <- is_mstrg & source_eligible[source_row] &
  !name %in% duplicate_source_names
mstrg <- data.frame(gene_id = de$gene_id[mstrg_ok],
                    symbol = name[mstrg_ok],
                    stat = de$stat[mstrg_ok],
                    padj = de$padj[mstrg_ok],
                    log2FoldChange = de$log2FoldChange[mstrg_ok],
                    source = "GEO_MSTRG_RGD_protein_coding")

# Map the remaining rat Ensembl IDs through the rat OrgDb. Require one
# unambiguous symbol per ID; do not turn multiple mapping rows into
# multiple copies of one tested feature.
ens <- de[is_ens, needed_de, drop = FALSE]
ens$ensembl <- sub("\\.[0-9]+$", "", ens$gene_id)
if (anyDuplicated(ens$ensembl)) stop("Duplicate Ensembl ID after version removal.")
map <- suppressMessages(AnnotationDbi::select(
  org.Rn.eg.db::org.Rn.eg.db, keys = ens$ensembl,
  keytype = "ENSEMBL", columns = "SYMBOL"
))
map <- unique(map[, c("ENSEMBL", "SYMBOL"), drop = FALSE])
map <- map[!is.na(map$SYMBOL) & nzchar(map$SYMBOL), , drop = FALSE]
ambiguous_ens <- names(which(table(map$ENSEMBL) > 1L))
map <- map[!map$ENSEMBL %in% ambiguous_ens, , drop = FALSE]
ens_symbol <- map$SYMBOL[match(ens$ensembl, map$ENSEMBL)]
ens_ok <- !is.na(ens_symbol)
ens_table <- data.frame(gene_id = ens$gene_id[ens_ok],
                        symbol = ens_symbol[ens_ok],
                        stat = ens$stat[ens_ok],
                        padj = ens$padj[ens_ok],
                        log2FoldChange = ens$log2FoldChange[ens_ok],
                        source = "OrgDb_ENS")

combined <- rbind(mstrg, ens_table)
duplicate_combined <- unique(combined$symbol[
  duplicated(combined$symbol) | duplicated(combined$symbol, fromLast = TRUE)
])
combined <- combined[!combined$symbol %in% duplicate_combined &
                       is.finite(combined$stat), , drop = FALSE]
if (nrow(combined) < 5000L || anyDuplicated(combined$symbol)) {
  stop("Insufficient unique, finite-statistic rat symbols: ", nrow(combined))
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(combined, file.path(out, "feature_to_symbol_audit.csv"),
          row.names = FALSE)
ranked <- combined[order(-combined$stat, combined$symbol), , drop = FALSE]
ranks <- setNames(ranked$stat, ranked$symbol)
write.csv(data.frame(symbol = names(ranks), statistic = unname(ranks)),
          file.path(out, "ranked_statistics.csv"), row.names = FALSE)

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
  stop("fgsea did not return expected columns.")
}
fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";",
                         FUN.VALUE = "")
fg <- fg[order(is.na(fg$padj), fg$padj, fg$pathway), , drop = FALSE]
write.csv(fg, file.path(out, "GSEA_Hallmark_results.csv"),
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
write.csv(selected, file.path(out, "selected_biological_programs.csv"),
          row.names = FALSE)
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
    ggplot2::labs(title = "GSE212311 CCI: source-aware Hallmark GSEA",
                  x = NULL, y = "Normalized enrichment score")
  ggplot2::ggsave(file.path(out, "selected_programs_NES.png"),
                  plot, width = 9, height = 6, dpi = 300)
}

# GO uses tested symbols mappable to one Entrez ID as its universe.
symbol_map <- suppressMessages(AnnotationDbi::select(
  org.Rn.eg.db::org.Rn.eg.db, keys = combined$symbol,
  keytype = "SYMBOL", columns = "ENTREZID"
))
symbol_map <- unique(symbol_map[, c("SYMBOL", "ENTREZID"), drop = FALSE])
symbol_map <- symbol_map[!is.na(symbol_map$ENTREZID) &
                           nzchar(symbol_map$ENTREZID), , drop = FALSE]
ambiguous_entrez <- names(which(table(symbol_map$SYMBOL) > 1L))
symbol_map <- symbol_map[!symbol_map$SYMBOL %in% ambiguous_entrez,
                         , drop = FALSE]
entrez_for <- function(symbols) {
  unique(symbol_map$ENTREZID[symbol_map$SYMBOL %in% symbols])
}
universe <- entrez_for(combined$symbol)
up <- entrez_for(combined$symbol[!is.na(combined$padj) &
                                  combined$padj < 0.05 &
                                  combined$log2FoldChange > 0])
down <- entrez_for(combined$symbol[!is.na(combined$padj) &
                                    combined$padj < 0.05 &
                                    combined$log2FoldChange < 0])
if (length(universe) < 5000L) {
  stop("Too few mapped genes in the GO measured universe: ",
       length(universe))
}
run_go <- function(genes, stem) {
  path <- file.path(out, paste0(stem, ".csv"))
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
  tab_go <- as.data.frame(result)
  write.csv(tab_go, path, row.names = FALSE)
  if (nrow(tab_go)) {
    plot <- enrichplot::dotplot(result, showCategory = min(20L, nrow(tab_go))) +
      ggplot2::ggtitle(paste("GSE212311", stem))
    ggplot2::ggsave(file.path(out, paste0(stem, "_dotplot.png")),
                    plot, width = 9, height = 7, dpi = 300)
  }
}
run_go(up, "GO_BP_upregulated")
run_go(down, "GO_BP_downregulated")

audit <- data.frame(
  tested_features = nrow(de),
  tested_MSTRG = sum(is_mstrg),
  eligible_unique_source_MSTRG = sum(mstrg_ok),
  duplicate_source_names_excluded = length(duplicate_source_names),
  eligible_Ensembl = nrow(ens_table),
  duplicate_cross_source_symbols_excluded = length(duplicate_combined),
  ranked_unique_symbols = length(ranks),
  hallmark_sets_tested = nrow(fg),
  mapped_GO_universe = length(universe),
  significant_up_GO_ids = length(up),
  significant_down_GO_ids = length(down)
)
write.csv(audit, file.path(out, "mapping_audit.csv"), row.names = FALSE)
writeLines(c(
  "GSE212311: independent CCI versus Sham, rat L4-L6 DRG day 11.",
  "MSTRG source names admitted only for single protein-coding RGD Symbol labels.",
  "Duplicate names and ambiguous mappings were excluded, not double counted.",
  "GSEA uses signed DESeq2 Wald statistics from all eligible tested features.",
  "GO uses FDR < 0.05 genes and the measured Entrez-mappable universe.",
  paste("Ranked unique symbols:", length(ranks)),
  paste("Hallmark sets tested:", nrow(fg)),
  paste("GO up/down Entrez:", length(up), "/", length(down)),
  "Bulk DRG cannot localize cell type or establish cell-cycle re-entry."
), file.path(out, "analysis_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
message("Source-aware GSE212311 pathway analysis complete: ", out)
