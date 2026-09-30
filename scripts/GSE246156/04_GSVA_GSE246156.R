#!/usr/bin/env Rscript

# GSE246156 NC L5 DRG, day 7: sample-level Hallmark GSVA.
# From the repository root, after steps 02 and 03:
# source("scripts/GSE246156/04_GSVA_GSE246156.R")
# Positive score difference means Compression > Sham within this study.

options(stringsAsFactors = FALSE)

required <- c("GSVA", "msigdbr", "pheatmap")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install missing packages: ", paste(missing, collapse = ", "))
}

base_dir <- file.path("results", "GSE246156_NC_L5_day7")
vst_file <- file.path(base_dir, "DESeq2", "GSE246156_vst_expression_matrix.csv")
meta_file <- file.path(base_dir, "locked_primary_samples.csv")
mapping_file <- file.path(base_dir, "pathway_analysis",
                          "Ensembl_annotation_all_mappings.csv")
outdir <- file.path(base_dir, "pathway_analysis")
for (input in c(vst_file, meta_file, mapping_file)) {
  if (!file.exists(input)) stop("Missing input: ", input,
                                ". Run the NC preparation, DESeq2 and GSEA/GO steps first.")
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

meta <- read.csv(meta_file, check.names = FALSE)
expected_sham <- paste0("GSM786379", 4:6)
expected_compression <- paste0("GSM786377", 0:2)
expected <- c(expected_sham, expected_compression)
if (!all(c("gsm", "condition", "tissue", "drg_level", "day") %in% names(meta)) ||
    nrow(meta) != 6L || anyDuplicated(meta$gsm) ||
    !setequal(meta$gsm, expected)) {
  stop("Locked metadata must contain exactly the six NC L5 day-7 GSMs.")
}
meta <- meta[match(expected, meta$gsm), , drop = FALSE]
if (!identical(as.character(meta$condition), rep(c("Sham", "Compression"), each = 3L)) ||
    !all(meta$tissue == "DRG") || !all(meta$drg_level == "L5") ||
    !all(as.character(meta$day) == "7")) {
  stop("Locked metadata condition, tissue, level or day do not match the primary contrast.")
}

vst_data <- read.csv(vst_file, row.names = 1L, check.names = FALSE)
if (nrow(vst_data) < 5000L || ncol(vst_data) != 6L ||
    anyDuplicated(rownames(vst_data)) || !setequal(colnames(vst_data), expected)) {
  stop("VST matrix has unexpected gene IDs or sample columns.")
}
vst <- as.matrix(vst_data[, expected, drop = FALSE])
storage.mode(vst) <- "numeric"
if (anyNA(vst) || any(!is.finite(vst)) ||
    !all(grepl("^ENSRNOG[0-9]+(\\.[0-9]+)?$", rownames(vst)))) {
  stop("VST expression must be finite and indexed by rat Ensembl gene IDs.")
}
gene_id <- sub("\\.[0-9]+$", "", rownames(vst))
if (anyDuplicated(gene_id)) stop("Duplicate Ensembl IDs after version removal.")

# Use the complete annotation audit produced by step 03. Exclude Ensembl
# IDs with multiple symbols; for duplicate symbols choose the Ensembl row
# with the greatest mean VST expression, without using treatment effects.
mapping <- read.csv(mapping_file, check.names = FALSE, na.strings = c("", "NA"))
if (!all(c("ENSEMBL", "SYMBOL") %in% names(mapping))) {
  stop("Annotation audit lacks ENSEMBL and SYMBOL columns.")
}
pairs <- unique(mapping[!is.na(mapping$SYMBOL) & nzchar(mapping$SYMBOL) &
                          !is.na(mapping$ENSEMBL), c("ENSEMBL", "SYMBOL")])
pairs <- pairs[pairs$ENSEMBL %in% gene_id, , drop = FALSE]
ambiguous <- names(which(table(pairs$ENSEMBL) > 1L))
pairs <- pairs[!pairs$ENSEMBL %in% ambiguous, , drop = FALSE]
pos <- match(pairs$ENSEMBL, gene_id)
if (anyNA(pos)) stop("Annotation and expression gene IDs failed to align.")
pairs$mean_vst <- rowMeans(vst[pos, , drop = FALSE])
pairs <- pairs[order(-pairs$mean_vst, pairs$ENSEMBL), , drop = FALSE]
duplicate_symbols <- sum(duplicated(pairs$SYMBOL))
pairs <- pairs[!duplicated(pairs$SYMBOL), , drop = FALSE]
if (nrow(pairs) < 5000L || anyDuplicated(pairs$SYMBOL)) {
  stop("Insufficient unambiguous rat symbols for GSVA: ", nrow(pairs))
}
pos <- match(pairs$ENSEMBL, gene_id)
expr_symbol <- vst[pos, , drop = FALSE]
rownames(expr_symbol) <- pairs$SYMBOL
if (anyDuplicated(rownames(expr_symbol)) || !identical(colnames(expr_symbol), expected)) {
  stop("Symbol expression matrix failed uniqueness or sample alignment checks.")
}

write.csv(pairs, file.path(outdir, "GSVA_symbol_mapping.csv"), row.names = FALSE)

# Match the previous GSE160543 analysis: rat Hallmark sets and GSVA's
# gsvaParam defaults applied to VST-scale (continuous) expression.
hallmark_table <- msigdbr::msigdbr(species = "Rattus norvegicus",
                                   collection = "H")
if (!all(c("gs_name", "gene_symbol") %in% names(hallmark_table))) {
  stop("Unexpected rat Hallmark annotation columns from msigdbr.")
}
hallmark <- lapply(split(hallmark_table$gene_symbol, hallmark_table$gs_name), unique)
overlap <- vapply(hallmark, function(x) length(intersect(x, rownames(expr_symbol))),
                  integer(1))
if (length(hallmark) != 50L || any(overlap < 15L)) {
  stop("Expected 50 rat Hallmark pathways with at least 15 mapped genes each.")
}
write.csv(data.frame(pathway = names(hallmark),
                     hallmark_genes = lengths(hallmark), mapped_genes = overlap),
          file.path(outdir, "GSVA_pathway_overlap.csv"), row.names = FALSE)

param <- GSVA::gsvaParam(exprData = expr_symbol, geneSets = hallmark)
scores <- as.matrix(GSVA::gsva(param, verbose = FALSE))
if (!setequal(rownames(scores), names(hallmark)) ||
    !identical(colnames(scores), expected) || any(!is.finite(scores))) {
  stop("GSVA returned unexpected pathways, sample order or non-finite scores.")
}
scores <- scores[sort(rownames(scores)), expected, drop = FALSE]
write.csv(data.frame(pathway = rownames(scores), scores, check.names = FALSE),
          file.path(outdir, "GSVA_Hallmark_scores.csv"), row.names = FALSE)

# The same eight pathways tested by the OIPN GSVA analysis in main.
selected <- c(
  "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MITOTIC_SPINDLE", "HALLMARK_P53_PATHWAY",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_IL6_JAK_STAT3_SIGNALING"
)
if (!all(selected %in% rownames(scores))) {
  stop("Missing a preselected OIPN comparison pathway from GSVA output.")
}
selected_scores <- scores[selected, , drop = FALSE]
write.csv(data.frame(pathway = selected, selected_scores, check.names = FALSE),
          file.path(outdir, "selected_GSVA_program_scores.csv"), row.names = FALSE)

stats <- do.call(rbind, lapply(selected, function(pathway) {
  sham <- as.numeric(scores[pathway, expected_sham])
  compression <- as.numeric(scores[pathway, expected_compression])
  # For 3 versus 3, the smallest attainable two-sided exact p is 0.1.
  # An asymptotic test is identified explicitly if exact ranks have ties.
  tied <- anyDuplicated(c(sham, compression)) > 0L
  test <- stats::wilcox.test(compression, sham, alternative = "two.sided",
                            exact = !tied, correct = FALSE)
  data.frame(
    pathway = pathway,
    Sham_mean = mean(sham), Compression_mean = mean(compression),
    mean_difference_Compression_minus_Sham = mean(compression) - mean(sham),
    Sham_median = median(sham), Compression_median = median(compression),
    wilcox_W = unname(test$statistic),
    p_value = unname(test$p.value),
    p_method = if (tied) "asymptotic_tied_ranks" else "exact_two_sided",
    complete_separation = if (min(compression) > max(sham)) "Compression_above_Sham"
                          else if (max(compression) < min(sham)) "Compression_below_Sham"
                          else "overlap",
    stringsAsFactors = FALSE
  )
}))
stats$FDR_BH_8 <- p.adjust(stats$p_value, method = "BH")
write.csv(stats, file.path(outdir, "GSVA_program_statistics.csv"), row.names = FALSE)

annotation <- data.frame(condition = factor(meta$condition,
                                             levels = c("Sham", "Compression")),
                         row.names = meta$gsm)
pheatmap::pheatmap(selected_scores, scale = "row", cluster_cols = FALSE,
                   annotation_col = annotation,
                   main = "GSE246156 L5 day 7: selected Hallmark GSVA",
                   filename = file.path(outdir, "GSVA_selected_programs_heatmap.png"),
                   width = 9, height = 6)

writeLines(c(
  "GSE246156: L5 day 7 DRG, three Compression versus three Sham samples.",
  "Input: DESeq2 VST matrix; rat Hallmark msigdbr; GSVA gsvaParam defaults.",
  "Symbol mapping: ambiguous Ensembl-symbol mappings excluded; duplicate symbols resolved by highest mean VST.",
  paste("VST genes:", nrow(vst)),
  paste("Ambiguous Ensembl IDs excluded:", length(ambiguous)),
  paste("Duplicate symbol rows discarded:", duplicate_symbols),
  paste("Unique symbols scored:", nrow(expr_symbol)),
  paste("Hallmark sets scored:", nrow(scores)),
  "Eight pathways selected before this NC GSVA test, matching the OIPN GSVA script.",
  "Selected pathways: two-sided Wilcoxon; exact when untied; BH adjustment across eight tests.",
  "With 3 versus 3 untied scores, the smallest two-sided exact p is 0.1.",
  "GSVA scores are relative within this cohort; do not compare raw score magnitudes across studies.",
  "Bulk DRG scores do not identify the cell of origin, causality or neuronal cell-cycle re-entry."
), file.path(outdir, "GSVA_analysis_summary.txt"))
writeLines(capture.output(sessionInfo()),
           file.path(outdir, "GSVA_sessionInfo.txt"))
message("GSE246156 NC Hallmark GSVA completed: ", outdir)
