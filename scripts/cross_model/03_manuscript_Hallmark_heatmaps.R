#!/usr/bin/env Rscript
# Run from repository root. No GSEA is rerun and no samples are pooled.
# Default figure 2 uses direction; choose 'support' for the separate FDR pattern.
run_manuscript_hallmarks <- function(draw_plots = TRUE, divergence_mode = "direction",
                                    out_dir = "results/manuscript_hallmark") {
  if (!divergence_mode %in% c("direction", "support")) stop("Unknown divergence_mode")
  cutoff <- 0.05
  inputs <- c(
    OIPN_GSE160543 = "results/GSE160543_Oxaliplatin_vs_Vehicle/Pathway_analysis/GSEA_Hallmark_results.csv",
    OIPN_GSE126773 = "results/GSE126773_OIPN/pathway_analysis/GSE126773_Hallmark_ranked_GSEA_results.csv",
    NC_GSE246156 = "results/GSE246156_NC_L5_day7/pathway_analysis/GSEA_Hallmark_results.csv",
    CCI_GSE212311 = "results/GSE212311_CCI_L4L6_day11/pathway_analysis_source_aware/GSEA_Hallmark_results.csv"
  )
  missing <- inputs[!file.exists(inputs)]
  if (length(missing)) stop("Run from repository root. Missing: ", paste(missing, collapse = "; "))
  if (draw_plots && !requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 first.")
  read_one <- function(path, microarray) {
    x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    columns <- if (microarray) c("ID", "NES", "p.adjust", "setSize") else c("pathway", "NES", "padj", "size")
    if (!all(columns %in% names(x))) stop("Required columns missing: ", path)
    y <- data.frame(pathway = as.character(x[[columns[1]]]),
                    NES = suppressWarnings(as.numeric(x[[columns[2]]])),
                    FDR = suppressWarnings(as.numeric(x[[columns[3]]])),
                    set_size = suppressWarnings(as.numeric(x[[columns[4]]])))
    if (nrow(y) != 50L || anyNA(y) || anyDuplicated(y$pathway) ||
        any(!startsWith(y$pathway, "HALLMARK_")) ||
        any(!is.finite(y$NES)) || any(!is.finite(y$FDR)) ||
        any(y$FDR < 0 | y$FDR > 1) || any(!is.finite(y$set_size)) ||
        any(y$set_size <= 0)) stop("Invalid complete Hallmark table: ", path)
    y
  }
  tables <- lapply(seq_along(inputs), function(i) read_one(inputs[i], i == 2L))
  names(tables) <- names(inputs)
  pathways <- sort(tables[[1]]$pathway)
  if (!all(vapply(tables, function(x) identical(sort(x$pathway), pathways), logical(1)))) {
    stop("Hallmark pathway identities differ across studies.")
  }
  tables <- lapply(tables, function(x) x[match(pathways, x$pathway), , drop = FALSE])
  nes <- do.call(cbind, lapply(tables, function(x) x$NES))
  fdr <- do.call(cbind, lapply(tables, function(x) x$FDR))
  positive <- nes > 0
  positive_sig <- positive & fdr < cutoff
  negative_sig <- nes < 0 & fdr < cutoff
  wide <- data.frame(pathway = pathways)
  for (i in seq_along(inputs)) {
    wide[[paste0(names(inputs)[i], "_NES")]] <- nes[, i]
    wide[[paste0(names(inputs)[i], "_FDR")]] <- fdr[, i]
    wide[[paste0(names(inputs)[i], "_set_size")]] <- tables[[i]]$set_size
  }
  wide$n_positive_significant <- rowSums(positive_sig)
  wide$n_negative_significant <- rowSums(negative_sig)
  wide$shared_positive <- rowSums(positive) == 4L & rowSums(positive_sig) >= 3L
  wide$shared_negative <- rowSums(nes < 0) == 4L & rowSums(negative_sig) >= 3L
  # Direction rule uses NES sign alone; FDR is annotated, not used for selection.
  # Require strict opposite signs in BOTH studies per side; zero is excluded.
  wide$oipn_direction <- rowSums(nes[, 1:2, drop = FALSE] > 0) == 2L & rowSums(nes[, 3:4, drop = FALSE] < 0) == 2L
  wide$physical_direction <- rowSums(nes[, 3:4, drop = FALSE] > 0) == 2L & rowSums(nes[, 1:2, drop = FALSE] < 0) == 2L
  # Support rule: comparator has NO positive significant enrichment in either study.
  # This includes positive nonsignificant NES and is NOT evidence of absent activity.
  wide$oipn_support <- rowSums(positive_sig[, 1:2, drop = FALSE]) == 2L & rowSums(positive_sig[, 3:4, drop = FALSE]) == 0L
  wide$physical_support <- rowSums(positive_sig[, 3:4, drop = FALSE]) == 2L & rowSums(positive_sig[, 1:2, drop = FALSE]) == 0L
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  save_table <- function(x, name) write.csv(x, file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)
  save_table(wide, "all_50_pathways_classified")
  save_table(wide[wide$oipn_direction | wide$physical_direction, , drop = FALSE], "Fig2_direction_selection")
  for (flag in c("shared_positive", "shared_negative", "oipn_direction", "physical_direction", "oipn_support", "physical_support")) {
    save_table(wide[wide[[flag]], , drop = FALSE], flag)
  }
  save_table(data.frame(study = names(inputs), input = unname(inputs),
                        md5 = unname(tools::md5sum(inputs)),
                        positive_significant = colSums(positive_sig),
                        negative_significant = colSums(negative_sig)), "input_audit")
  counts <- vapply(wide[c("shared_positive", "shared_negative", "oipn_direction", "physical_direction", "oipn_support", "physical_support")], sum, numeric(1))
  writeLines(c(paste(names(counts), counts, sep = " = "),
               "Shared: same NES sign in all four and within-study FDR < 0.05 in at least three.",
               "Direction: both OIPN NES > 0 and both physical-injury NES < 0, or vice versa; no FDR selection filter.",
               "Support: both favored studies NES > 0 and FDR < 0.05; neither comparator positive and significant.",
               "Asterisk annotation: within-study FDR < 0.05 in every figure.",
               "Selection is descriptive and post hoc; recurrence has no combined FDR or interaction p-value.",
               "No row scaling; fixed columns; common symmetric NES scale across figures.",
               paste("Figure 2 mode:", divergence_mode)), file.path(out_dir, "selection_summary.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  if (draw_plots) {
    labels <- c("OIPN\nGSE160543", "OIPN\nGSE126773", "NC\nGSE246156", "CCI\nGSE212311")
    limit <- max(0.5, ceiling(max(abs(nes)) * 2) / 2)
    draw <- function(idx, filename, title, groups = NULL) {
      if (!length(idx)) {
        stale <- file.path(out_dir, paste0(filename, c(".pdf", ".png")))
        unlink(stale[file.exists(stale)])
        writeLines("No pathways met the declared criteria; no heatmap generated.", file.path(out_dir, paste0(filename, "_empty.txt")))
        message(filename, ": no eligible pathways")
        return(invisible(NULL))
      }
      empty_note <- file.path(out_dir, paste0(filename, "_empty.txt"))
      if (file.exists(empty_note)) unlink(empty_note)
      long <- do.call(rbind, lapply(seq_along(inputs), function(i) {
        data.frame(pathway = pathways[idx], study = labels[i], NES = nes[idx, i],
                   label = paste0(sprintf("%.2f", nes[idx, i]), ifelse(fdr[idx, i] < cutoff, "*", "")),
                   group = if (is.null(groups)) "" else groups)
      }))
      long$study <- factor(long$study, levels = labels)
      long$pathway <- factor(long$pathway, levels = rev(pathways[idx]))
      pretty_names <- function(x) gsub("_", " ", sub("^HALLMARK_", "", x))
      p <- ggplot2::ggplot(long, ggplot2::aes(study, pathway, fill = NES)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.4) +
        ggplot2::geom_text(ggplot2::aes(label = label), size = 3.5) +
        ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0, limits = c(-limit, limit), name = "NES") +
        ggplot2::scale_y_discrete(labels = pretty_names) +
        ggplot2::labs(title = title, subtitle = "* within-study FDR < 0.05; independent bulk DRG contrasts", x = NULL, y = NULL) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_text(color = "black"),
                       plot.title = ggplot2::element_text(face = "bold"))
      if (!is.null(groups)) p <- p + ggplot2::facet_grid(group ~ ., scales = "free_y", space = "free_y")
      height <- max(3.8, 2.4 + 0.42 * length(idx))
      ggplot2::ggsave(file.path(out_dir, paste0(filename, ".pdf")), p, width = 10, height = height)
      ggplot2::ggsave(file.path(out_dir, paste0(filename, ".png")), p, width = 10, height = height, dpi = 600, bg = "white")
    }
    draw(which(wide$shared_positive), "Fig1_shared_positive", "Concordant positive Hallmark enrichment")
    o <- wide[[paste0("oipn_", divergence_mode)]]
    ph <- wide[[paste0("physical_", divergence_mode)]]
    idx <- which(o | ph)
    groups <- if (divergence_mode == "direction") {
      ifelse(o[idx], "Positive in OIPN / negative in NC and CCI", "Positive in NC and CCI / negative in OIPN")
    } else {
      ifelse(o[idx], "OIPN-favoring support", "Physical-injury-favoring support")
    }
    draw(idx, paste0("Fig2_divergent_", divergence_mode),
         if (divergence_mode == "direction") "Divergent Hallmark enrichment directions" else "Different patterns of positive enrichment support", groups)
  }
  message("Manuscript Hallmark extraction complete: ", out_dir)
  invisible(wide)
}

# Set the environment variable HALLMARK_TABLES_ONLY=true to skip drawing.
# Set HALLMARK_DIVERGENCE_MODE=support only for the separately defined FDR pattern.
hallmark_manuscript_result <- run_manuscript_hallmarks(
  draw_plots = tolower(Sys.getenv("HALLMARK_TABLES_ONLY", "false")) != "true",
  divergence_mode = Sys.getenv("HALLMARK_DIVERGENCE_MODE", "direction")
)
