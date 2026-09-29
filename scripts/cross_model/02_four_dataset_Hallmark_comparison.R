#!/usr/bin/env Rscript
# Four independent rat DRG Hallmark GSEA contrasts. Run from repository root.
# This script compares pathway directions and within-study FDR, not pooled samples.

inputs <- c(
  OIPN_GSE160543 = "results/GSE160543_Oxaliplatin_vs_Vehicle/Pathway_analysis/GSEA_Hallmark_results.csv",
  OIPN_GSE126773 = "results/GSE126773_OIPN/pathway_analysis/GSE126773_Hallmark_ranked_GSEA_results.csv",
  NC_GSE246156 = "results/GSE246156_NC_L5_day7/pathway_analysis/GSEA_Hallmark_results.csv",
  CCI_GSE212311 = "results/GSE212311_CCI_L4L6_day11/pathway_analysis_source_aware/GSEA_Hallmark_results.csv"
)

missing <- inputs[!file.exists(inputs)]
if (length(missing)) {
  stop("Missing Hallmark result(s): ",
       paste(names(missing), missing, sep = " = ", collapse = "; "))
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package ggplot2 is required to save the four-dataset heatmap.")
}

read_hallmark <- function(path, microarray = FALSE) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  id_col <- if (microarray) "ID" else "pathway"
  fdr_col <- if (microarray) "p.adjust" else "padj"
  size_col <- if (microarray) "setSize" else "size"
  required <- c(id_col, "NES", fdr_col, size_col)
  if (!all(required %in% names(x))) {
    stop("Missing columns in ", path, ": ",
         paste(setdiff(required, names(x)), collapse = ", "))
  }
  result <- data.frame(
    pathway = x[[id_col]], NES = as.numeric(x$NES),
    FDR = as.numeric(x[[fdr_col]]), set_size = as.numeric(x[[size_col]]),
    stringsAsFactors = FALSE
  )
  if (nrow(result) != 50L || anyNA(result) ||
      anyDuplicated(result$pathway) ||
      any(!startsWith(result$pathway, "HALLMARK_")) ||
      any(!is.finite(result$NES)) ||
      any(result$FDR < 0 | result$FDR > 1) ||
      any(result$set_size <= 0)) {
    stop("Invalid 50-pathway Hallmark result: ", path)
  }
  result
}

tables <- lapply(seq_along(inputs), function(i) {
  read_hallmark(inputs[[i]], microarray = names(inputs)[i] == "OIPN_GSE126773")
})
names(tables) <- names(inputs)
pathways <- sort(tables[[1L]]$pathway)
if (!all(vapply(tables, function(x) identical(sort(x$pathway), pathways),
                logical(1)))) {
  stop("The four studies do not contain identical Hallmark pathway sets.")
}

aligned <- lapply(tables, function(x) x[match(pathways, x$pathway), ])
wide <- data.frame(pathway = pathways, stringsAsFactors = FALSE)
for (model in names(aligned)) {
  wide[[paste0(model, "_NES")]] <- aligned[[model]]$NES
  wide[[paste0(model, "_FDR")]] <- aligned[[model]]$FDR
}
significant <- sapply(aligned, function(x) x$FDR < 0.05)
signs <- sapply(aligned, function(x) sign(x$NES))
wide$significant_studies <- rowSums(significant)
wide$concordant_direction <- apply(signs, 1L, function(x) all(x > 0) || all(x < 0))
wide$significant_all_four <- wide$significant_studies == 4L

out_dir <- "results/cross_model_four_dataset"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(wide, file.path(out_dir, "Hallmark_four_dataset_comparison.csv"),
          row.names = FALSE)

plot_data <- do.call(rbind, lapply(names(aligned), function(model) {
  x <- aligned[[model]]
  data.frame(pathway = x$pathway, model = model, NES = x$NES,
             significant = x$FDR < 0.05)
}))
plot_data$model <- factor(plot_data$model, levels = names(aligned))
plot_data$pathway <- factor(plot_data$pathway, levels = rev(pathways))
limit <- max(abs(plot_data$NES))
plot <- ggplot2::ggplot(plot_data,
                        ggplot2::aes(x = model, y = pathway, fill = NES)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.2) +
  ggplot2::geom_text(ggplot2::aes(label = ifelse(significant, "*", "")),
                     size = 3) +
  ggplot2::scale_fill_gradient2(low = "#286EAA", mid = "white",
                                high = "#CE553B", midpoint = 0,
                                limits = c(-limit, limit)) +
  ggplot2::labs(
    title = "Rat DRG Hallmark GSEA across four datasets",
    subtitle = "* within-study FDR < 0.05; signed NES; methods differ by platform",
    x = NULL, y = NULL) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_text(size = 8),
                 axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
ggplot2::ggsave(file.path(out_dir, "Hallmark_four_dataset_heatmap.png"),
                plot, width = 11, height = 15, dpi = 300)

writeLines(c(
  "Independent rat bulk DRG contrasts; no count matrices or samples were pooled.",
  "GSE160543: OIPN RNA-seq, oxaliplatin versus vehicle, 4 versus 4.",
  "GSE126773: OIPN microarray, oxaliplatin versus vehicle; limma ranking and clusterProfiler GSEA.",
  "GSE246156: NC RNA-seq, L5 day 7 compression versus sham, 3 versus 3.",
  "GSE212311: CCI RNA-seq, L4-L6 day 11 CCI versus sham, 3 versus 3; source-aware mapping.",
  "The three RNA-seq studies use DESeq2 Wald ranks and fgsea; GSE126773 uses limma moderated t ranks and clusterProfiler GSEA.",
  "All four tables have the same 50 Hallmark pathway labels, but measured gene universes and enrichment algorithms differ.",
  paste("Within-study FDR < 0.05:",
        paste(names(aligned), colSums(significant), sep = "=", collapse = "; ")),
  paste("FDR < 0.05 in all four:", sum(wide$significant_all_four)),
  paste("Same NES direction in all four:", sum(wide$concordant_direction)),
  "NES differences are descriptive and are not formal between-study interaction tests.",
  "Absence of significance does not demonstrate model specificity; bulk tissue does not identify cell type."
), file.path(out_dir, "comparison_summary.txt"))
message("Four-dataset Hallmark comparison complete: ", out_dir)
