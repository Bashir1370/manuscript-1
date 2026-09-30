# Shared helpers for the manuscript GSVA workflow; sourced by the two entry scripts.
# The selected pathways are post hoc GSEA-derived displays, not independent tests.
gsva_selected_pathways <- character()

gsva_check_scores <- function(scores, expected_samples) {
  if (!is.matrix(scores) || !is.numeric(scores) || nrow(scores) != 50L ||
      anyNA(scores) || any(!is.finite(scores)) ||
      is.null(rownames(scores)) || anyDuplicated(rownames(scores)) ||
      !all(startsWith(rownames(scores), "HALLMARK_")) ||
      anyDuplicated(colnames(scores)) ||
      !setequal(colnames(scores), expected_samples) ||
      ncol(scores) != length(expected_samples) ||
      !all(gsva_selected_pathways %in% rownames(scores))) {
    stop("Invalid complete 50-pathway GSVA matrix or sample identities.")
  }
  scores[sort(rownames(scores)), expected_samples, drop = FALSE]
}

gsva_fit_scores <- function(scores, metadata) {
  if (!requireNamespace("limma", quietly = TRUE)) stop("Install limma.")
  if (!all(c("sample", "condition") %in% names(metadata)) ||
      anyNA(metadata[, c("sample", "condition")]) || anyDuplicated(metadata$sample) ||
      !all(metadata$condition %in% c("Control", "Neuropathy")) ||
      any(table(factor(metadata$condition, levels = c("Control", "Neuropathy"))) < 3L)) {
    stop("Metadata must identify >=3 independent samples per group.")
  }
  scores <- gsva_check_scores(scores, metadata$sample)
  metadata$condition <- factor(metadata$condition, levels = c("Control", "Neuropathy"))
  design <- stats::model.matrix(~ condition, data = metadata)
  if (qr(design)$rank != ncol(design)) stop("GSVA design matrix is not full rank.")
  # Apply the same moderated linear model to all 50 pathways in each study.
  # No row scaling, pooling, DEG filtering, or FDR-based pathway selection.
  fit <- limma::eBayes(limma::lmFit(scores, design), trend = FALSE, robust = FALSE)
  tab <- limma::topTable(fit, coef = "conditionNeuropathy", number = Inf,
                         sort.by = "none", adjust.method = "BH", confint = 0.95)
  tab <- tab[match(rownames(scores), rownames(tab)), , drop = FALSE]
  control <- metadata$condition == "Control"
  injury <- metadata$condition == "Neuropathy"
  delta <- rowMeans(scores[, injury, drop = FALSE]) - rowMeans(scores[, control, drop = FALSE])
  if (!isTRUE(all.equal(unname(tab$logFC), unname(delta), tolerance = 1e-8))) {
    stop("limma coefficient does not match Neuropathy minus Control.")
  }
  if (any(!is.finite(tab$P.Value)) || any(!is.finite(tab$adj.P.Val))) {
    stop("Non-finite GSVA statistical results; inspect score variances.")
  }
  result <- data.frame(
    pathway = rownames(scores), control_mean = rowMeans(scores[, control, drop = FALSE]),
    neuropathy_mean = rowMeans(scores[, injury, drop = FALSE]), delta_GSVA = delta,
    CI95_low = tab$CI.L, CI95_high = tab$CI.R, moderated_t = tab$t,
    p_value = tab$P.Value, FDR_BH_50 = tab$adj.P.Val,
    n_control = sum(control), n_neuropathy = sum(injury),
    stringsAsFactors = FALSE, row.names = NULL
  )
  result$direction <- ifelse(result$delta_GSVA > 0, "higher_in_neuropathy",
                             ifelse(result$delta_GSVA < 0, "lower_in_neuropathy", "no_difference"))
  result$complete_separation <- vapply(seq_len(nrow(scores)), function(i) {
    a <- scores[i, injury]; b <- scores[i, control]
    if (min(a) > max(b)) "neuropathy_above_control"
    else if (max(a) < min(b)) "neuropathy_below_control" else "overlap"
  }, character(1))
  result
}
