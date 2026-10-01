#!/usr/bin/env Rscript
source("scripts/gobp/helpers.R", local = TRUE)
run_gobp_sample_context <- function(draw_plots = TRUE) {
  categories <- c("shared_positive", "shared_negative", "divergent")
  files <- file.path(gobp_root, "leading_edge", categories, "priority_OIPN_plus_physical.csv")
  priorities <- lapply(files, function(f) gobp_read(f, "symbol")); names(priorities) <- categories
  symbols <- sort(unique(unlist(lapply(priorities, function(d) d$symbol))), method = "radix")
  out <- file.path(gobp_root, "sample_expression")
  rows <- diagnostics <- list()
  if (length(symbols)) for (study in gobp_studies) {
    a <- gobp_expression(study)
    if (!all(symbols %in% rownames(a$expr))) stop("Priority genes unavailable in expression: ", study)
    for (symbol in symbols) {
      x <- a$counts[symbol, ]; lx <- a$expr[symbol, ]
      z <- if (sd(lx) > 0) (lx - mean(lx)) / sd(lx) else rep(0, length(lx))
      origin <- paste(categories[vapply(priorities, function(d) symbol %in% d$symbol, logical(1))], collapse = ";")
      rows[[length(rows) + 1L]] <- data.frame(study = study, symbol = symbol,
        feature_id = a$genes$feature_id[match(symbol, a$genes$symbol)], sample = a$meta$sample,
        condition = a$meta$condition, source_group = a$meta$source_group, origin = origin,
        normalized_count = unname(x), log2_count_plus_1 = unname(lx), within_study_z = unname(z),
        sample_label = paste0(rep(c("Ctrl", "Inj"), each = length(x) / 2), rep(seq_len(length(x) / 2), 2L)))
      for (group in c("Control", "Neuropathy")) {
        v <- x[a$meta$condition == group]
        diagnostics[[length(diagnostics) + 1L]] <- data.frame(study = study, symbol = symbol, condition = group,
          mean_log2_count_plus_1 = mean(log2(v + 1)), min_count = min(v), max_count = max(v),
          max_sample_fraction_of_group_counts = if (sum(v) > 0) max(v) / sum(v) else NA_real_)
      }
    }
  }
  d <- if (length(rows)) do.call(rbind, rows) else data.frame(study = character(), symbol = character(), feature_id = character(), sample = character(), condition = character(), source_group = character(), origin = character(), normalized_count = numeric(), log2_count_plus_1 = numeric(), within_study_z = numeric(), sample_label = character())
  diag <- if (length(diagnostics)) do.call(rbind, diagnostics) else data.frame(study = character(), symbol = character(), condition = character(), mean_log2_count_plus_1 = numeric(), min_count = numeric(), max_count = numeric(), max_sample_fraction_of_group_counts = numeric())
  gobp_write(d, file.path(out, "priority_sample_expression.csv"))
  gobp_write(diag, file.path(out, "group_expression_diagnostics.csv"))
  gobp_audit(c(files, unlist(lapply(gobp_studies, function(study) {
    s <- gobp_spec(study); c(s$counts, s$de, s$mapping, s$ranks, s$meta)
  })), "scripts/gobp/05_sample_expression_GO_BP.R"), out)
  writeLines(c("Descriptive expression context for already prioritized GO:BP genes; no new significance filtering.",
    "log2(normalized count +1), z-standardized separately for each gene and study for display only.",
    "No cell fraction estimates, automatic sample exclusion, neuronal localization or mechanistic inference."), file.path(out, "analysis_notes.txt"))
  if (draw_plots) {
    d$pathway <- d$symbol; d$mark <- ""
    limit <- if (nrow(d)) max(1, ceiling(max(abs(d$within_study_z)) * 2) / 2) else 1
    gobp_plot_pages(d, symbols, file.path(out, "figures"), "priority_sample_expression", "GO:BP priority genes: sample-level expression", "within_study_z", limit, x = "sample_label", facet_studies = TRUE)
  }
  message("GO:BP sample expression complete: ", length(symbols), " genes")
  invisible(d)
}
if (!exists("GOBP_AUTORUN", inherits = FALSE) || isTRUE(GOBP_AUTORUN)) run_gobp_sample_context(tolower(Sys.getenv("GOBP_TABLES_ONLY", "false")) != "true")
