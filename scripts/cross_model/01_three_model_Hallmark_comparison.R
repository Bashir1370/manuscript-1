#!/usr/bin/env Rscript
# Compare within-study rat DRG Hallmark GSEA across OIPN, NC and CCI.
# Run from the repository root after GSE212311 step 06:
# source("scripts/cross_model/01_three_model_Hallmark_comparison.R")
# No counts or samples are pooled across studies.

input <- c(
  OIPN = file.path("results", "GSE160543_Oxaliplatin_vs_Vehicle",
                   "Pathway_analysis", "GSEA_Hallmark_results.csv"),
  NC = file.path("results", "GSE246156_NC_L5_day7",
                 "pathway_analysis", "GSEA_Hallmark_results.csv"),
  CCI = file.path("results", "GSE212311_CCI_L4L6_day11",
                  "pathway_analysis_source_aware",
                  "GSEA_Hallmark_results.csv")
)
missing <- input[!file.exists(input)]
if (length(missing)) {
  stop("Missing within-study Hallmark result(s): ",
       paste(names(missing), missing, sep = " = ", collapse = "; "))
}
read_result <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  needed <- c("pathway", "NES", "padj", "size", "leadingEdge")
  if (!all(needed %in% names(x)) || nrow(x) != 50L ||
      anyNA(x$pathway) || anyDuplicated(x$pathway) ||
      anyNA(x$NES) || any(!is.finite(x$NES)) ||
      anyNA(x$padj) || any(x$padj < 0 | x$padj > 1) ||
      anyNA(x$size) || any(x$size < 15 | x$size > 500)) {
    stop("Invalid 50-pathway Hallmark table: ", path)
  }
  x
}
tables <- lapply(input, read_result)
pathways <- sort(tables[[1L]]$pathway)
if (any(!vapply(tables, function(x) identical(sort(x$pathway), pathways),
                logical(1)))) {
  stop("The three studies do not contain the same 50 Hallmark pathways.")
}

out <- file.path("results", "cross_model_three_way")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
wide <- data.frame(pathway = pathways, stringsAsFactors = FALSE)
aligned <- lapply(tables, function(x) x[match(pathways, x$pathway), , drop = FALSE])
for (model in names(aligned)) {
  x <- aligned[[model]]
  wide[[paste0(model, "_NES")]] <- x$NES
  wide[[paste0(model, "_FDR")]] <- x$padj
  wide[[paste0(model, "_set_size")]] <- x$size
}
sig <- sapply(aligned, function(x) x$padj < 0.05)
colnames(sig) <- names(aligned)
wide$significant_models <- apply(sig, 1L, function(x) {
  if (!any(x)) "none" else paste(names(aligned)[x], collapse = "+")
})
wide$direction <- apply(sapply(aligned, function(x) sign(x$NES)),
                        1L, function(x) {
  if (all(x > 0)) "positive_all_three" else
    if (all(x < 0)) "negative_all_three" else "mixed"
})
wide$significant_all_three <- rowSums(sig) == 3L
wide$concordant_significant_all_three <-
  wide$significant_all_three & wide$direction != "mixed"

leading <- function(value) {
  if (is.na(value) || !nzchar(value)) character() else
    unique(strsplit(value, ";", fixed = TRUE)[[1L]])
}
shared <- lapply(seq_along(pathways), function(i) {
  if (!wide$significant_all_three[i]) return(character())
  Reduce(intersect, lapply(aligned, function(x) leading(x$leadingEdge[i])))
})
wide$shared_leading_edge_n <- lengths(shared)
wide$shared_leading_edge_genes <- vapply(shared, paste,
                                         collapse = ";", FUN.VALUE = "")
write.csv(wide, file.path(out, "Hallmark_three_model_comparison.csv"),
          row.names = FALSE)

selected <- c(
  "HALLMARK_G2M_CHECKPOINT", "HALLMARK_E2F_TARGETS",
  "HALLMARK_MITOTIC_SPINDLE", "HALLMARK_P53_PATHWAY",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_DNA_REPAIR"
)
if (!all(selected %in% pathways)) stop("A prespecified pathway is absent.")
write.csv(wide[match(selected, pathways), , drop = FALSE],
          file.path(out, "selected_programs_three_model.csv"),
          row.names = FALSE)
write.csv(wide[wide$significant_all_three, , drop = FALSE],
          file.path(out, "shared_significant_pathways.csv"),
          row.names = FALSE)

plot_data <- do.call(rbind, lapply(names(aligned), function(model) {
  x <- aligned[[model]]
  data.frame(pathway = x$pathway, model = model, NES = x$NES,
             significant = x$padj < 0.05)
}))
plot_data$pathway <- factor(plot_data$pathway, levels = rev(pathways))
plot_data$model <- factor(plot_data$model, levels = names(aligned))
limit <- max(abs(plot_data$NES))
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot <- ggplot2::ggplot(plot_data,
                          ggplot2::aes(x = model, y = pathway, fill = NES)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(
      label = ifelse(significant, "*", "")), size = 3) +
    ggplot2::scale_fill_gradient2(
      low = "#286EAA", mid = "white", high = "#CE553B",
      midpoint = 0, limits = c(-limit, limit)) +
    ggplot2::labs(
      title = "Rat DRG Hallmark GSEA: OIPN, NC and CCI",
      subtitle = "* FDR < 0.05 within study; color shows signed NES",
      x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 8))
  ggplot2::ggsave(file.path(out, "Hallmark_three_model_heatmap.png"),
                  plot, width = 9, height = 15, dpi = 250)
}

writeLines(c(
  "Three independent rat bulk DRG contrasts: OIPN 4 vs 4, NC L5 day 7 3 vs 3, CCI L4-L6 day 11 3 vs 3.",
  "All three source tables contain the same 50 rat Hallmark pathways.",
  paste("FDR < 0.05 within OIPN/NC/CCI:",
        paste(colSums(sig), collapse = " / ")),
  paste("FDR < 0.05 in all three:", sum(wide$significant_all_three)),
  paste("Sign concordant among these:", sum(wide$concordant_significant_all_three)),
  paste("Mixed-sign among these:",
        sum(wide$significant_all_three & wide$direction == "mixed")),
  "A non-significant pathway in one model does not establish model specificity.",
  "NES values reflect different measured gene universes; compare direction and within-study support, not NES differences as interaction tests.",
  "Shared leading-edge names are descriptive overlaps, not independent replication or causal evidence.",
  "Bulk tissue cannot establish cell type or neuronal cell-cycle re-entry."
), file.path(out, "comparison_summary.txt"))
message("Three-model Hallmark comparison complete: ", out)
