#!/usr/bin/env Rscript
# Descriptive sample-level expression of existing priorities, not cell fractions.
run_sample_expression_context <- function(draw_plots = TRUE) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) stop("Install DESeq2 for NC size factors.")
  if (draw_plots && !requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2.")
  bases <- c("results/shared_Hallmark_leading_edge_three_dataset",
             "results/shared_negative_Hallmark_leading_edge_three_dataset",
             "results/divergent_Hallmark_leading_edge_three_dataset")
  priority_files <- c(file.path(bases[1:2], "gene_prioritization", "priority_OIPN_plus_physical.csv"),
                      file.path(bases[3], "priority_OIPN_plus_physical.csv"))
  evidence_files <- file.path(bases, "gene_evidence_long.csv")
  files <- c(priority_files, evidence_files)
  if (any(!file.exists(files))) stop("Run scripts 05/06/07/08 first; missing priority/evidence tables.")
  priorities <- lapply(priority_files, read.csv, stringsAsFactors = FALSE)
  genes <- sort(unique(unlist(lapply(priorities, function(d) d$symbol))))
  evidence <- do.call(rbind, lapply(evidence_files, read.csv, stringsAsFactors = FALSE))
  mapping <- unique(evidence[evidence$symbol %in% genes & evidence$rank_available,
                            c("symbol", "study", "feature_id")])
  if (anyDuplicated(mapping[c("symbol", "study")])) stop("Inconsistent original GSEA feature mapping.")
  studies <- c("OIPN_GSE160543", "NC_GSE246156", "CCI_GSE212311")
  samples <- list(c(paste0("GSM487500", 3:6), paste0("GSM487501", 1:4)),
                  c(paste0("GSM786379", 4:6), paste0("GSM786377", 0:2)),
                  c(paste0("GSM652375", 1:3), paste0("GSM652374", 8:9), "GSM6523750"))
  count_files <- c("results/GSE160543_Oxaliplatin_vs_Vehicle/normalized_counts_complete_case.csv",
                   "results/GSE246156_NC_L5_day7/raw_count_matrix.csv",
                   "results/GSE212311_CCI_L4L6_day11/DESeq2/normalized_counts_filtered.csv")
  if (any(!file.exists(count_files))) stop("Missing canonical count matrices.")
  result <- diagnostics <- list()
  for (i in seq_along(studies)) {
    tab <- read.csv(count_files[i], check.names = FALSE, stringsAsFactors = FALSE)
    ids <- as.character(tab[[1L]])
    if (anyDuplicated(ids) || !all(samples[[i]] %in% names(tab))) stop("Invalid count feature/sample keys.")
    counts <- as.matrix(tab[, samples[[i]], drop = FALSE]); storage.mode(counts) <- "numeric"
    rownames(counts) <- ids
    if (any(!is.finite(counts)) || any(counts < 0)) stop("Invalid count values.")
    condition <- rep(c("Control", "Neuropathy"), each = length(samples[[i]]) / 2L)
    if (i == 2L) {
      # Match the original NC prefilter; estimating size factors does not refit DE.
      keep <- rowSums(counts >= 10) >= 3L
      if (any(counts != round(counts))) stop("NC raw counts must be integers.")
      dds <- DESeq2::DESeqDataSetFromMatrix(counts[keep, , drop = FALSE],
        data.frame(condition = condition, row.names = samples[[i]]), design = ~ condition)
      dds <- DESeq2::estimateSizeFactors(dds)
      counts <- DESeq2::counts(dds, normalized = TRUE)
    }
    for (gene in genes) {
      feature <- mapping$feature_id[mapping$study == studies[i] & mapping$symbol == gene]
      if (length(feature) != 1L || !feature %in% rownames(counts)) stop("Priority feature missing: ", studies[i], " / ", gene)
      x <- as.numeric(counts[feature, ]); logx <- log2(x + 1)
      z <- if (stats::sd(logx) > 0) (logx - mean(logx)) / stats::sd(logx) else rep(0, length(logx))
      origin <- paste(c("shared_positive", "shared_negative", "divergent")[vapply(priorities,
        function(d) gene %in% d$symbol, logical(1))], collapse = ";")
      result[[length(result)+1L]] <- data.frame(study = studies[i], symbol = gene, feature_id = feature,
        sample = samples[[i]], condition = condition, origin = origin,
        normalized_count = x, log2_count_plus_1 = logx, within_study_z = z)
      for (group in c("Control", "Neuropathy")) {
        v <- x[condition == group]; lv <- logx[condition == group]
        diagnostics[[length(diagnostics)+1L]] <- data.frame(study = studies[i], symbol = gene, condition = group,
          mean_log2_count_plus_1 = mean(lv), min_count = min(v), max_count = max(v),
          max_sample_fraction_of_group_counts = if (sum(v) > 0) max(v) / sum(v) else NA_real_)
      }
    }
  }
  d <- do.call(rbind, result); out <- "results/three_dataset_sample_expression_context"
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.csv(d, file.path(out, "priority_sample_expression.csv"), row.names = FALSE)
  write.csv(do.call(rbind, diagnostics), file.path(out, "group_expression_diagnostics.csv"), row.names = FALSE)
  inputs <- c(files, count_files)
  write.csv(data.frame(input = inputs, md5 = unname(tools::md5sum(inputs))), file.path(out, "input_checksums.csv"), row.names = FALSE)
  if (draw_plots) {
    d$symbol <- factor(d$symbol, levels = rev(genes))
    d$study <- factor(d$study, levels = studies)
    p <- ggplot2::ggplot(d, ggplot2::aes(sample, symbol, fill = within_study_z)) +
      ggplot2::geom_tile() + ggplot2::facet_grid(. ~ study, scales = "free_x", space = "free_x") +
      ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
      ggplot2::labs(title = "Sample-level expression of existing priority genes", x = NULL, y = NULL,
        fill = "Within-study z", caption = "log2(normalized count + 1); each gene standardized within its study. No cell fractions or cell-type assignment.") +
      ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
    ggplot2::ggsave(file.path(out, "priority_sample_expression.pdf"), p, width = 12, height = max(6, .2*length(genes)+2))
    ggplot2::ggsave(file.path(out, "priority_sample_expression.png"), p, width = 12, height = max(6, .2*length(genes)+2), dpi = 300)
  }
  writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
  writeLines(c("Exploratory display of already selected genes; no new significance selection.",
    "OIPN and CCI use canonical normalized counts; NC size factors use original raw counts/filter.",
    "Compare sample patterns within a study, not calibrated count magnitudes across studies.",
    "A dominant sample is descriptive, not an automatic outlier-removal criterion.",
    "These candidate panels cannot estimate cell composition or prove neuronal localization."
  ), file.path(out, "analysis_notes.txt"))
  message("Sample-level expression context complete: ", out)
  invisible(d)
}
sample_expression_context <- run_sample_expression_context(
  draw_plots = tolower(Sys.getenv("CONTEXT_TABLES_ONLY", "false")) != "true")
