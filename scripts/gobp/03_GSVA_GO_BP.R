#!/usr/bin/env Rscript
source("scripts/gobp/helpers.R", local = TRUE)
gobp_fit_gsva <- function(scores, meta) {
  if (!is.matrix(scores) || !is.numeric(scores) || !nrow(scores) || any(!is.finite(scores)) ||
      anyDuplicated(rownames(scores)) || !identical(colnames(scores), meta$sample)) stop("Invalid GSVA score matrix.")
  md <- meta; md$condition <- factor(md$condition, levels = c("Control", "Neuropathy"))
  design <- model.matrix(~ condition, md)
  fit <- limma::eBayes(limma::lmFit(scores, design), trend = FALSE, robust = FALSE)
  tt <- limma::topTable(fit, coef = "conditionNeuropathy", number = Inf, sort.by = "none", adjust.method = "BH", confint = .95)
  tt <- tt[match(rownames(scores), rownames(tt)), , drop = FALSE]
  ctrl <- md$condition == "Control"; inj <- !ctrl
  delta <- rowMeans(scores[, inj, drop = FALSE]) - rowMeans(scores[, ctrl, drop = FALSE])
  if (!isTRUE(all.equal(unname(delta), unname(tt$logFC), tolerance = 1e-8)) || any(!is.finite(tt$adj.P.Val))) stop("GSVA contrast/fit validation failed.")
  result <- data.frame(pathway = rownames(scores), control_mean = rowMeans(scores[, ctrl, drop = FALSE]),
    neuropathy_mean = rowMeans(scores[, inj, drop = FALSE]), delta_GSVA = unname(delta),
    CI95_low = tt$CI.L, CI95_high = tt$CI.R, moderated_t = tt$t, p_value = tt$P.Value, FDR_BH = tt$adj.P.Val,
    n_control = sum(ctrl), n_neuropathy = sum(inj), tested_family_size = nrow(scores), row.names = NULL)
  result$direction <- ifelse(delta > 0, "higher_in_neuropathy", ifelse(delta < 0, "lower_in_neuropathy", "no_difference"))
  result$complete_separation <- vapply(seq_len(nrow(scores)), function(i) {
    a <- scores[i, inj]; b <- scores[i, ctrl]
    if (min(a) > max(b)) "neuropathy_above_control" else if (max(a) < min(b)) "neuropathy_below_control" else "overlap"
  }, character(1))
  result
}
run_gobp_gsva <- function(draw_plots = TRUE) {
  gobp_require(c("GSVA", "limma", "BiocParallel", "DESeq2"))
  if (!"gsvaParam" %in% getNamespaceExports("GSVA")) stop("Update GSVA: gsvaParam API is required.")
  sets <- gobp_load_sets()
  selection_file <- file.path(gobp_root, "comparison/all_GO_BP_pathways_classified.csv")
  selection <- gobp_read(selection_file, c("pathway", "shared_positive", "shared_negative", "oipn_direction", "physical_direction"))
  selected <- selection$pathway[selection$shared_positive | selection$shared_negative | selection$oipn_direction | selection$physical_direction]
  all_stats <- selected_samples <- list()
  for (study in gobp_studies) {
    message("GO:BP GSVA: ", study)
    s <- gobp_spec(study); a <- gobp_expression(study)
    # A uniform transform avoids requiring the parked OIPN rerun / unarchived VST.
    variable <- apply(a$expr, 1L, function(z) max(z) > min(z))
    expr <- a$expr[variable, , drop = FALSE]
    sizes <- vapply(sets, function(z) length(intersect(z, rownames(expr))), integer(1))
    eligible <- sizes >= gobp_min_size & sizes <= gobp_max_size
    if (!any(eligible)) stop("No eligible GSVA GO:BP sets: ", study)
    param <- GSVA::gsvaParam(exprData = expr, geneSets = sets[eligible], minSize = gobp_min_size,
      maxSize = gobp_max_size, kcdf = "Gaussian", tau = 1, maxDiff = TRUE, absRanking = FALSE)
    set.seed(s$seed)
    scores <- as.matrix(GSVA::gsva(param, verbose = FALSE, BPPARAM = BiocParallel::SerialParam()))
    if (!setequal(rownames(scores), names(sets)[eligible]) || !identical(colnames(scores), a$meta$sample) || any(!is.finite(scores))) stop("Unexpected GSVA output: ", study)
    scores <- scores[sort(rownames(scores)), a$meta$sample, drop = FALSE]
    fit <- gobp_fit_gsva(scores, a$meta); fit$study <- study
    gp <- file.path(gobp_root, "GSEA", study, "GSEA_GO_BP_results.csv")
    g <- gobp_read(gp); idx <- match(fit$pathway, g$pathway)
    fit$GSEA_NES <- g$NES[idx]; fit$GSEA_FDR <- g$padj[idx]
    fit$direction_agrees_with_GSEA <- ifelse(is.finite(fit$GSEA_NES), sign(fit$delta_GSVA) == sign(fit$GSEA_NES), NA)
    fit$selected_for_display <- fit$pathway %in% selected
    out <- file.path(gobp_root, "GSVA", study)
    gobp_write(data.frame(pathway = rownames(scores), scores, check.names = FALSE), file.path(out, "GSVA_GO_BP_scores.csv"))
    gobp_write(fit, file.path(out, "GSVA_limma_all_GO_BP.csv"))
    gobp_write(data.frame(pathway = names(sets), mapped_rat_genes = lengths(sets), variable_expression_overlap = sizes, eligible = eligible), file.path(out, "pathway_coverage.csv"))
    gobp_write(a$genes, file.path(out, "feature_to_symbol_audit.csv")); gobp_write(a$meta, file.path(out, "sample_metadata.csv"))
    gobp_audit(c(s$counts, s$de, s$mapping, s$ranks, s$meta, gp, selection_file, gobp_lock_files(), "scripts/gobp/03_GSVA_GO_BP.R"), out)
    all_stats[[study]] <- fit
    display_ids <- intersect(selected, rownames(scores))
    if (length(display_ids)) {
      m <- scores[display_ids, , drop = FALSE]
      sd <- apply(m, 1L, stats::sd); z <- sweep(sweep(m, 1L, rowMeans(m), "-"), 1L, ifelse(sd > 0, sd, 1), "/")
      selected_samples[[study]] <- do.call(rbind, lapply(seq_len(ncol(m)), function(j) data.frame(pathway = rownames(m), study = study,
        sample = colnames(m)[j], condition = a$meta$condition[j], source_group = a$meta$source_group[j],
        sample_label = paste0(ifelse(j <= ncol(m) / 2, "Ctrl", "Inj"), ifelse(j <= ncol(m) / 2, j, j - ncol(m) / 2)),
        GSVA_score = unname(m[, j]), within_study_z = unname(z[, j]), mark = "", row.names = NULL)))
    }
  }
  out <- file.path(gobp_root, "GSVA")
  combined <- do.call(rbind, all_stats); rownames(combined) <- NULL
  gobp_write(combined, file.path(out, "GSVA_limma_all_GO_BP_three_studies.csv"))
  gobp_write(combined[combined$selected_for_display, , drop = FALSE], file.path(out, "GSVA_GSEA_selected_comparison.csv"))
  long <- if (length(selected_samples)) do.call(rbind, selected_samples) else data.frame(pathway = character(), study = character(), sample = character(), condition = character(), source_group = character(), sample_label = character(), GSVA_score = numeric(), within_study_z = numeric(), mark = character())
  gobp_write(long, file.path(out, "GSVA_selected_sample_scores.csv"))
  writeLines(c("All cohorts rescored using log2(normalized count +1), Gaussian GSVA; no pooled samples.",
    "OIPN/CCI use archived normalized counts; NC size factors estimated on original >=10 in >=3 sample filter; no DE refit.",
    "Original GSEA representative features are retained for symbol alignment; selection can depend on original DE statistic.",
    "This differs from historical VST GSVA: GO-vs-Hallmark differences cannot be attributed solely to gene-set choice.",
    "Constant genes removed before calculating eligibility; BH across all scored GO:BP sets per study.",
    "limma eBayes trend=FALSE robust=FALSE, Neuropathy minus Control; model-based exploratory estimates.",
    "GSEA-selected displays are from the same data, not independent validation; raw scores are within-study only.",
    "Plot z-scores standardized within each study and pathway; missing rows are not assigned zero scores."), file.path(out, "analysis_notes.txt"))
  if (draw_plots) {
    limit <- if (nrow(long)) max(1, ceiling(max(abs(long$within_study_z)) * 2) / 2) else 1
    gobp_plot_pages(long, intersect(selected, unique(long$pathway)), file.path(out, "figures"), "GSVA_selected_samples", "GO:BP sample-level GSVA; within-study z scores", "within_study_z", limit, x = "sample_label", facet_studies = TRUE)
    fig <- file.path(out, "figures")
    unlink(list.files(fig, pattern = "^GSVA_sample_dotplots_page_", full.names = TRUE))
    ids <- intersect(selected, unique(long$pathway))
    pages <- split(ids, ceiling(seq_along(ids) / 6L))
    for (j in seq_along(pages)) {
      d <- long[long$pathway %in% pages[[j]], , drop = FALSE]
      d$condition <- factor(d$condition, levels = c("Control", "Neuropathy"))
      labels <- function(z) vapply(gobp_pretty(z), function(a) paste(strwrap(a, width = 45), collapse = "\n"), character(1))
      panels <- as.vector(t(outer(labels(pages[[j]]), gobp_studies, paste, sep = "\n")))
      d$panel <- factor(paste(labels(d$pathway), d$study, sep = "\n"), levels = panels)
      p <- ggplot2::ggplot(d, ggplot2::aes(condition, GSVA_score, color = condition)) +
        ggplot2::geom_point(position = ggplot2::position_jitter(width = .08, height = 0, seed = 160543), size = 1.5) +
        ggplot2::stat_summary(fun = mean, geom = "point", shape = 95, size = 6, color = "black") +
        ggplot2::facet_wrap(~ panel, ncol = 3L, scales = "free_y") +
        ggplot2::scale_color_manual(values = c(Control = "#607D8B", Neuropathy = "#B2182B")) +
        ggplot2::labs(title = "GO:BP sample-level GSVA", subtitle = "Each point is a sample; black bar = mean; separate y scales", x = NULL, y = "GSVA score") +
        ggplot2::theme_bw(base_size = 8) + ggplot2::theme(legend.position = "none")
      stem <- file.path(fig, sprintf("GSVA_sample_dotplots_page_%03d", j)); ht <- max(5, 3 * length(pages[[j]]))
      ggplot2::ggsave(paste0(stem, ".pdf"), p, width = 14, height = ht, limitsize = FALSE)
      ggplot2::ggsave(paste0(stem, ".png"), p, width = 14, height = ht, dpi = 200, bg = "white", limitsize = FALSE)
    }
  }
  message("GO:BP GSVA complete: ", out)
  invisible(combined)
}
if (!exists("GOBP_AUTORUN", inherits = FALSE) || isTRUE(GOBP_AUTORUN)) run_gobp_gsva(tolower(Sys.getenv("GOBP_TABLES_ONLY", "false")) != "true")
