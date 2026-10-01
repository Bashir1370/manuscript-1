#!/usr/bin/env Rscript
source("scripts/gobp/helpers.R", local = TRUE)
run_gobp_comparison <- function(draw_plots = TRUE) {
  sets <- gobp_load_sets(); ids <- names(sets)
  nes <- fdr <- sizes <- matrix(NA_real_, length(ids), 3L, dimnames = list(ids, gobp_studies))
  inputs <- character()
  for (i in seq_along(gobp_studies)) {
    study <- gobp_studies[i]; file <- file.path(gobp_root, "GSEA", study, "GSEA_GO_BP_results.csv")
    g <- gobp_read(file, c("pathway", "NES", "padj", "size", "test_status")); inputs <- c(inputs, file)
    if (anyDuplicated(g$pathway) || !all(g$pathway %in% ids) || any(!is.na(g$padj) & (g$padj < 0 | g$padj > 1))) stop("Invalid GO:BP results: ", study)
    pos <- match(g$pathway, ids); nes[pos, i] <- g$NES; fdr[pos, i] <- g$padj; sizes[pos, i] <- g$size
  }
  wide <- data.frame(pathway = ids)
  for (i in seq_along(gobp_studies)) {
    wide[[paste0(gobp_studies[i], "_NES")]] <- nes[, i]
    wide[[paste0(gobp_studies[i], "_FDR")]] <- fdr[, i]
    wide[[paste0(gobp_studies[i], "_set_size")]] <- sizes[, i]
  }
  wide <- cbind(wide, gobp_classify(nes, fdr)); out <- file.path(gobp_root, "comparison")
  gobp_write(wide, file.path(out, "all_GO_BP_pathways_classified.csv"))
  flags <- c("shared_positive", "shared_negative", "oipn_direction", "physical_direction", "oipn_support", "physical_support")
  for (flag in flags) gobp_write(wide[wide[[flag]], , drop = FALSE], file.path(out, paste0(flag, ".csv")))
  gobp_write(wide[wide$oipn_direction | wide$physical_direction, , drop = FALSE], file.path(out, "Fig2_direction_selection.csv"))
  gobp_write(data.frame(metric = c("locked_pathways", "available_all3", flags),
    count = c(length(ids), sum(wide$available_all3), vapply(wide[flags], sum, integer(1)))), file.path(out, "selection_summary.csv"))
  writeLines(c("Shared +/-: strict same NES sign and pathway BH FDR <0.05 in all three studies.",
    "Opposite direction: OIPN NES has the opposite sign to both NC and CCI; no FDR selection filter.",
    "Support: positive significant enrichment in the favored study/studies and no positive significant enrichment in comparators.",
    "Missing, undersized, oversized or numerically failed pathways are unavailable, not zero or nonsignificant evidence.",
    "Each study's FDR is over its full eligible GO:BP family, not selected plots.",
    "Cross-study recurrence is descriptive; no combined FDR, formal interaction test or OIPN-specificity claim.",
    "GO terms overlap and nest; counts of significant terms are not independent biological programs."), file.path(out, "analysis_notes.txt"))
  gobp_audit(c(inputs, gobp_lock_files(), "scripts/gobp/02_compare_GO_BP.R"), out)
  if (draw_plots) {
    long <- do.call(rbind, lapply(seq_along(gobp_studies), function(i) data.frame(pathway = ids,
      study = gobp_studies[i], NES = nes[, i], mark = ifelse(is.finite(fdr[, i]) & fdr[, i] < gobp_fdr, "*", ""))))
    finite <- nes[is.finite(nes)]; limit <- if (length(finite)) max(.5, ceiling(max(abs(finite)) * 2) / 2) else 1
    fig <- file.path(out, "figures")
    gobp_plot_pages(long, ids, fig, "GO_BP_all", "All GO:BP pathways", "NES", limit)
    for (flag in c("shared_positive", "shared_negative")) gobp_plot_pages(long, ids[wide[[flag]]], fig, paste0("Fig1_", flag), paste("GO:BP", flag), "NES", limit)
    gobp_plot_pages(long, ids[wide$oipn_direction | wide$physical_direction], fig, "Fig2_divergent_direction", "Opposite GO:BP NES directions", "NES", limit)
  }
  message("GO:BP comparison complete: ", out)
  invisible(wide)
}
if (!exists("GOBP_AUTORUN", inherits = FALSE) || isTRUE(GOBP_AUTORUN)) run_gobp_comparison(tolower(Sys.getenv("GOBP_TABLES_ONLY", "false")) != "true")
