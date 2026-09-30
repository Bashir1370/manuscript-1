#!/usr/bin/env Rscript
# Post hoc gene prioritization; retains original within-study DE FDR.
run_three_dataset_gene_prioritization <- function() {
  base <- "results/shared_Hallmark_leading_edge_three_dataset"
  input <- file.path(base, "shared_all3.csv")
  if (!file.exists(input)) stop("Missing input: ", input,
    "; run scripts/cross_model/05_three_dataset_leading_edge.R from repository root.")
  d <- read.csv(input, check.names = FALSE, stringsAsFactors = FALSE)
  studies <- c("OIPN_GSE160543", "NC_GSE246156", "CCI_GSE212311")
  effect_cols <- paste0(studies, "_log2FC")
  fdr_cols <- paste0(studies, "_gene_FDR")
  le_cols <- paste0(studies, "_LE")
  required <- c("pathway", "symbol", "shared_all3", effect_cols, fdr_cols, le_cols)
  if (!all(required %in% names(d))) stop("Invalid shared_all3 schema.")
  if (!nrow(d) || anyNA(d$symbol) || any(!nzchar(d$symbol)) ||
      anyNA(d$pathway) || any(!nzchar(d$pathway)) ||
      anyDuplicated(d[c("pathway", "symbol")])) stop("Invalid pathway-gene keys.")
  if (!is.logical(d$shared_all3) || anyNA(d$shared_all3) ||
      !all(d$shared_all3)) stop("Input contains non-shared genes.")
  for (col in le_cols) {
    if (!is.logical(d[[col]]) || anyNA(d[[col]]) || !all(d[[col]]))
      stop("All-three leading-edge membership required: ", col)
  }
  for (col in effect_cols) {
    if (!is.numeric(d[[col]]) || any(!is.finite(d[[col]])))
      stop("Finite original log2FC required: ", col)
  }
  for (col in fdr_cols) {
    x <- d[[col]]
    if (!is.numeric(x) && !all(is.na(x))) stop("Numeric gene FDR required: ", col)
    if (any(!is.na(x) & (!is.finite(x) | x < 0 | x > 1)))
      stop("Invalid gene FDR: ", col)
  }
  # Gene evidence must be identical when the same symbol occurs in several pathways.
  same <- function(x) {
    if (all(is.na(x))) return(TRUE)
    if (anyNA(x)) return(FALSE)
    all(x == x[1L])
  }
  symbols <- sort(unique(d$symbol), method = "radix")
  for (symbol in symbols) {
    rows <- d[d$symbol == symbol, , drop = FALSE]
    if (!all(vapply(rows[c(effect_cols, fdr_cols)], same, logical(1))))
      stop("Inconsistent repeated gene evidence: ", symbol)
  }
  genes <- d[match(symbols, d$symbol), c("symbol", effect_cols, fdr_cols), drop = FALSE]
  rownames(genes) <- NULL
  genes$n_shared_pathways <- vapply(symbols,
    function(s) as.integer(sum(d$symbol == s)), integer(1))
  genes$shared_pathways <- vapply(symbols, function(s)
    paste(sort(unique(d$pathway[d$symbol == s]), method = "radix"), collapse = ";"),
    character(1))
  effects <- as.matrix(genes[effect_cols])
  fdr <- as.matrix(genes[fdr_cols])
  significant <- !is.na(fdr) & fdr < 0.05
  genes$positive_all3 <- rowSums(effects > 0) == 3L
  genes$n_gene_FDR_lt_0_05 <- rowSums(significant)
  genes$OIPN_significant <- significant[, 1L]
  genes$NC_significant <- significant[, 2L]
  genes$CCI_significant <- significant[, 3L]
  genes$selected_OIPN_plus_physical <- genes$positive_all3 &
    genes$OIPN_significant & (genes$NC_significant | genes$CCI_significant)
  genes$strict_significant_all3 <- genes$positive_all3 &
    genes$n_gene_FDR_lt_0_05 == 3L
  genes$physical_only_significant <- genes$positive_all3 &
    !genes$OIPN_significant & genes$NC_significant & genes$CCI_significant
  genes$support_pattern <- ifelse(genes$strict_significant_all3,
    "OIPN_NC_CCI", ifelse(genes$OIPN_significant & genes$NC_significant,
    "OIPN_NC", ifelse(genes$OIPN_significant & genes$CCI_significant,
    "OIPN_CCI", ifelse(genes$NC_significant & genes$CCI_significant,
    "NC_CCI", "fewer_than_two_significant"))))
  selected <- genes[genes$selected_OIPN_plus_physical, , drop = FALSE]
  strict <- genes[genes$strict_significant_all3, , drop = FALSE]
  two <- selected[selected$n_gene_FDR_lt_0_05 == 2L, , drop = FALSE]
  physical <- genes[genes$physical_only_significant, , drop = FALSE]
  memberships <- d[d$symbol %in% selected$symbol, , drop = FALSE]
  memberships <- memberships[order(memberships$symbol, memberships$pathway,
    method = "radix"), , drop = FALSE]
  stopifnot(all(strict$symbol %in% selected$symbol),
    nrow(selected) == nrow(strict) + nrow(two),
    !any(selected$symbol %in% physical$symbol))
  out <- file.path(base, "gene_prioritization")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write <- function(x, name) write.csv(x, file.path(out, name), row.names = FALSE, na = "")
  write(genes, "all_shared_genes_with_priority_flags.csv")
  write(selected, "priority_OIPN_plus_physical.csv")
  write(strict, "strict_significant_all3.csv")
  write(two, "priority_significant_exactly2.csv")
  write(physical, "excluded_physical_only_significant.csv")
  write(memberships, "priority_pathway_memberships.csv")
  writeLines(selected$symbol, file.path(out, "priority_genes_STRING.txt"))
  writeLines(strict$symbol, file.path(out, "strict_genes_STRING.txt"))
  counts <- data.frame(metric = c("shared_unique_genes", "positive_all3",
    "FDR_ge2_any_pair_positive_all3", "priority_OIPN_plus_physical",
    "strict_significant_all3", "priority_significant_exactly2",
    "excluded_physical_only", "priority_pathway_memberships"),
    count = c(nrow(genes), sum(genes$positive_all3),
    sum(genes$positive_all3 & genes$n_gene_FDR_lt_0_05 >= 2L),
    nrow(selected), nrow(strict), nrow(two), nrow(physical), nrow(memberships)))
  write(counts, "selection_summary.csv")
  write(data.frame(input = input, md5 = unname(tools::md5sum(input))),
    "input_checksums.csv")
  writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
  message("OIPN plus physical-injury prioritization complete: ", nrow(selected),
    " unique genes (", nrow(strict), " significant in all three); ", out)
  invisible(genes)
}
three_dataset_gene_prioritization <- run_three_dataset_gene_prioritization()
