#!/usr/bin/env Rscript
# Run after three-study GSEA extraction. All three complete GSVA matrices already exist.
# No pooling of cohorts and no overwriting of historical Wilcoxon results.
run_three_dataset_GSVA <- function(draw_plots = TRUE) {
  if (!requireNamespace("limma", quietly = TRUE)) stop("Install limma.")
  if (draw_plots && !requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2.")
  helper <- "scripts/cross_model/GSVA_helpers.R"
  if (!file.exists(helper)) stop("Run from repository root after pulling the branch.")
  source(helper, local = TRUE)
  selection_file <- "results/manuscript_hallmark_three_dataset/all_50_pathways_classified.csv"
  selection <- read.csv(selection_file, check.names=FALSE, stringsAsFactors=FALSE)
  if (nrow(selection)!=50L || anyDuplicated(selection$pathway)) stop("Run three-study GSEA extraction first.")
  gsva_selected_pathways <- selection$pathway[selection$shared_positive | selection$shared_negative | selection$oipn_direction | selection$physical_direction]
  if (!length(gsva_selected_pathways)) stop("No selected display pathways.")
  specs <- list(
    OIPN_GSE160543 = list(
      score = "results/GSE160543_Oxaliplatin_vs_Vehicle/Pathway_analysis/GSVA_Hallmark_scores.csv",
      meta = "data/metadata/gse160543_analysis_selection.csv", id = "gsm", group = "group",
      expected = c(paste0("GSM487500", 3:6), paste0("GSM487501", 1:4)),
      groups = rep(c("Vehicle", "Oxaliplatin"), each = 4L),
      gsea = "results/GSE160543_Oxaliplatin_vs_Vehicle/Pathway_analysis/GSEA_Hallmark_results.csv"),
    NC_GSE246156 = list(
      score = "results/GSE246156_NC_L5_day7/pathway_analysis/GSVA_Hallmark_scores.csv",
      meta = "results/GSE246156_NC_L5_day7/locked_primary_samples.csv", id = "gsm", group = "condition",
      expected = c(paste0("GSM786379", 4:6), paste0("GSM786377", 0:2)),
      groups = rep(c("Sham", "Compression"), each = 3L),
      gsea = "results/GSE246156_NC_L5_day7/pathway_analysis/GSEA_Hallmark_results.csv"),
    CCI_GSE212311 = list(
      score = "results/GSE212311_CCI_L4L6_day11/GSVA_source_aware/GSVA_Hallmark_scores.csv",
      meta = "results/GSE212311_CCI_L4L6_day11/locked_primary_samples.csv", id = "gsm", group = "condition",
      expected = c(paste0("GSM652375", 1:3), paste0("GSM652374", 8:9), "GSM6523750"),
      groups = rep(c("Sham", "CCI"), each = 3L),
      gsea = "results/GSE212311_CCI_L4L6_day11/pathway_analysis_source_aware/GSEA_Hallmark_results.csv")
  )
  files <- unlist(lapply(specs, function(x) c(x$score, x$meta, x$gsea)), use.names = FALSE)
  if (any(!file.exists(files))) stop("Missing input(s): ", paste(files[!file.exists(files)], collapse = "; "),
                                    ". Run the original retained-study scoring workflows first.")
  out <- "results/GSVA_three_dataset"
  matrices <- results <- samples <- list()
  for (study in names(specs)) {
    spec <- specs[[study]]
    meta <- read.csv(spec$meta, check.names = FALSE, stringsAsFactors = FALSE)
    if (!all(c(spec$id, spec$group) %in% names(meta)) || anyDuplicated(meta[[spec$id]]) ||
        !all(spec$expected %in% meta[[spec$id]])) stop("Invalid metadata: ", study)
    # Select the locked contrast by GSM IDs; paclitaxel is not included.
    meta <- meta[match(spec$expected, meta[[spec$id]]), , drop = FALSE]
    if (!identical(as.character(meta[[spec$group]]), spec$groups)) stop("Wrong group assignments: ", study)
    metadata <- data.frame(sample = spec$expected, source_group = spec$groups,
                            condition = rep(c("Control", "Neuropathy"), each = length(spec$expected) / 2))
    tab <- read.csv(spec$score, check.names = FALSE, stringsAsFactors = FALSE)
    if (!"pathway" %in% names(tab) || !setequal(setdiff(names(tab), "pathway"), spec$expected)) {
      stop("Unexpected GSVA score columns: ", study)
    }
    scores <- as.matrix(tab[, spec$expected, drop = FALSE])
    rownames(scores) <- tab$pathway
    scores <- gsva_check_scores(scores, spec$expected)
    if (length(matrices) && !identical(rownames(scores), rownames(matrices[[1L]]))) {
      stop("Hallmark pathway identities differ across studies.")
    }
    stats <- gsva_fit_scores(scores, metadata)
    gsea <- read.csv(spec$gsea, check.names = FALSE, stringsAsFactors = FALSE)
    id_col <- "pathway"
    fdr_col <- "padj"
    if (!all(c(id_col, "NES", fdr_col) %in% names(gsea)) || nrow(gsea) != 50L ||
        anyDuplicated(gsea[[id_col]]) || !setequal(gsea[[id_col]], stats$pathway)) stop("Invalid GSEA input: ", study)
    gsea <- gsea[match(stats$pathway, gsea[[id_col]]), , drop = FALSE]
    if (!is.numeric(gsea$NES) || !is.numeric(gsea[[fdr_col]]) ||
        any(!is.finite(gsea$NES)) || any(!is.finite(gsea[[fdr_col]])) ||
        any(gsea[[fdr_col]] < 0 | gsea[[fdr_col]] > 1)) stop("Invalid GSEA values: ", study)
    stats$study <- study
    stats$GSEA_NES <- gsea$NES
    stats$GSEA_FDR <- gsea[[fdr_col]]
    stats$direction_agrees_with_GSEA <- sign(stats$delta_GSVA) == sign(stats$GSEA_NES)
    stats$selected_for_display <- stats$pathway %in% gsva_selected_pathways
    matrices[[study]] <- scores
    results[[study]] <- stats
    samples[[study]] <- metadata
  }
  # Do not write partial three-study results when an input fails validation.
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  combined <- do.call(rbind, results); rownames(combined) <- NULL
  selected <- combined[combined$selected_for_display, , drop = FALSE]
  write.csv(combined, file.path(out, "GSVA_limma_all_50_three_studies.csv"), row.names = FALSE)
  write.csv(selected, file.path(out, "GSVA_GSEA_selected_comparison.csv"), row.names = FALSE)
  for (study in names(results)) write.csv(results[[study]], file.path(out, paste0(study, "_limma_all_50.csv")), row.names = FALSE)
  long <- do.call(rbind, lapply(names(matrices), function(study) {
    m <- matrices[[study]][gsva_selected_pathways, , drop = FALSE]
    sds <- apply(m, 1L, stats::sd)
    # Standardize each pathway within each study only, for display.
    z <- sweep(m, 1L, rowMeans(m), "-")
    z <- sweep(z, 1L, ifelse(sds > 0, sds, 1), "/")
    do.call(rbind, lapply(seq_len(ncol(m)), function(j) {
      data.frame(study = study, pathway = rownames(m), sample = colnames(m)[j],
                  condition = samples[[study]]$condition[j],
                  source_group = samples[[study]]$source_group[j],
                  GSVA_score = m[, j], within_study_z = z[, j],
                  sample_label = paste0(ifelse(samples[[study]]$condition[j] == "Control", "Ctrl", "Inj"),
                                        ifelse(samples[[study]]$condition[j] == "Control", j, j - ncol(m) / 2)),
                  stringsAsFactors = FALSE, row.names = NULL)
    }))
  }))
  write.csv(long, file.path(out, "GSVA_selected_sample_scores.csv"), row.names = FALSE)
  files <- c(files, selection_file)
  audit <- do.call(rbind, lapply(names(specs), function(study) {
    x <- specs[[study]]
    data.frame(study = study, input_type = c("GSVA_scores", "sample_metadata", "GSEA"),
                input = c(x$score, x$meta, x$gsea),
                md5 = unname(tools::md5sum(c(x$score, x$meta, x$gsea))))
  }))
  audit <- rbind(audit, data.frame(study="selection", input_type="GSEA_selection", input=selection_file, md5=unname(tools::md5sum(selection_file))))
  write.csv(audit, file.path(out, "input_audit.csv"), row.names = FALSE)
  summary <- do.call(rbind, lapply(names(results), function(study) {
    x <- results[[study]]; y <- x[x$selected_for_display, ]
    data.frame(study = study, pathways_tested = nrow(x), significant_FDR_BH_50 = sum(x$FDR_BH_50 < 0.05),
                selected_direction_agreement = sum(y$direction_agrees_with_GSEA), selected_pathways = nrow(y))
  }))
  write.csv(summary, file.path(out, "comparison_summary.csv"), row.names = FALSE)
  if (draw_plots) {
    study_labels <- c(OIPN_GSE160543 = "OIPN GSE160543",
                       NC_GSE246156 = "NC GSE246156", CCI_GSE212311 = "CCI GSE212311")
    pretty <- function(x) gsub("_", " ", sub("^HALLMARK_", "", x))
    long$study <- factor(long$study, levels = names(specs), labels = study_labels)
    long$pathway <- factor(long$pathway, levels = rev(gsva_selected_pathways), labels = rev(pretty(gsva_selected_pathways)))
    long$sample_label <- factor(long$sample_label, levels = c(paste0("Ctrl", 1:4), paste0("Inj", 1:4)))
    long$condition <- factor(long$condition, levels = c("Control", "Neuropathy"))
    limit <- max(1, ceiling(max(abs(long$within_study_z)) * 2) / 2)
    heatmap <- ggplot2::ggplot(long, ggplot2::aes(sample_label, pathway, fill = within_study_z)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.3) +
      ggplot2::facet_grid(. ~ study, scales = "free_x", space = "free_x") +
      ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
                                    limits = c(-limit, limit), name = "Within-study\nz score") +
      ggplot2::labs(title = "Selected Hallmark GSVA scores in individual DRG samples",
                     subtitle = "Each row standardized separately within each study; Ctrl = control; Inj = neuropathy",
                     x = NULL, y = NULL) + ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
    ggplot2::ggsave(file.path(out, "GSVA_selected_sample_heatmap.pdf"), heatmap, width = 12, height = max(6, 2+.3*length(gsva_selected_pathways)))
    ggplot2::ggsave(file.path(out, "GSVA_selected_sample_heatmap.png"), heatmap, width = 12, height = max(6, 2+.3*length(gsva_selected_pathways)), dpi = 600, bg = "white")
    # Free y scales for every pathway-study panel: no cross-study score calibration.
    points <- long
    points$panel <- factor(paste(points$pathway, points$study, sep = "\n"),
                            levels = as.vector(t(outer(pretty(gsva_selected_pathways), study_labels, paste, sep = "\n"))))
    dotplot <- ggplot2::ggplot(points, ggplot2::aes(condition, GSVA_score, color = condition)) +
      ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.08, height = 0, seed = 160543), size = 1.8) +
      ggplot2::stat_summary(fun = mean, geom = "point", shape = 95, size = 7, color = "black") +
      ggplot2::facet_wrap(~ panel, ncol = 3, scales = "free_y") +
      ggplot2::scale_color_manual(values = c(Control = "#607D8B", Neuropathy = "#B2182B")) +
      ggplot2::labs(title = "Sample-level GSVA: Control versus Neuropathy",
                     subtitle = "Every point is a sample; black bar = mean; y scales vary by panel", x = NULL, y = "GSVA score") +
      ggplot2::theme_bw(base_size = 9) + ggplot2::theme(legend.position = "none", panel.grid.minor = ggplot2::element_blank())
    ggplot2::ggsave(file.path(out, "GSVA_selected_sample_dotplots.pdf"), dotplot, width = 12, height = max(10, 2*length(gsva_selected_pathways)))
    ggplot2::ggsave(file.path(out, "GSVA_selected_sample_dotplots.png"), dotplot, width = 12, height = max(10, 2*length(gsva_selected_pathways)), dpi = 300, bg = "white")
  }
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  writeLines(c(
    "Three independent rat bulk DRG contrasts; 50 Hallmarks tested separately in each study.",
    "limma eBayes trend=FALSE, robust=FALSE; BH correction over 50 pathways per study.",
    "Positive delta_GSVA = Neuropathy minus Control > 0, not a positive raw score alone.",
    "Moderated confidence intervals and FDR are exploratory model-based estimates; small groups limit inference.",
    "Display pathways selected from three-study GSEA on the same samples, not independent validation.",
    "RNA-seq scores reused from the three historical workflows.",
    "Assay gene universes, symbol mappings and original gene-set releases may differ; not a harmonized rescore.",
    "Raw delta magnitudes are descriptive within each study; no calibrated cross-study effect-size comparison.",
    "Heatmap uses within-study row z scores for display only; limma uses original scores.",
    "Historical Wilcoxon tables remain unchanged; compare the two methods transparently.",
    "No sample pooling, batch correction across studies, or mechanistic/cell-type inference."
  ), file.path(out, "analysis_notes.txt"))
  message("Three-study GSVA comparison complete: ", out)
  invisible(selected)
}
GSVA_three_dataset_result <- run_three_dataset_GSVA(
  draw_plots = tolower(Sys.getenv("GSVA_TABLES_ONLY", "false")) != "true"
)
