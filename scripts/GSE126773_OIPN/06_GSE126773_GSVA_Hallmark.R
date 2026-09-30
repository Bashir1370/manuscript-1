#!/usr/bin/env Rscript
# Run from repository root; uses the existing RMA matrix and probe annotation.
# Produces scores/tables only. The four-study script below handles figures.
run_GSE126773_GSVA <- function() {
  required <- c("GSVA", "msigdbr", "limma", "BiocParallel")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
  helper <- "scripts/cross_model/GSVA_helpers.R"
  if (!file.exists(helper)) stop("Run from the repository root after pulling the branch.")
  source(helper, local = TRUE)
  root <- "results/GSE126773_OIPN"
  inputs <- c(
    expression = file.path(root, "GSE126773_RMA_expression_matrix.csv"),
    metadata = file.path(root, "GSE126773_OIPN_analysis_samples.csv"),
    annotation = file.path(root, "annotation", "GSE126773_Oxaliplatin_vs_Vehicle_limma_annotated.csv")
  )
  if (any(!file.exists(inputs))) stop("Missing inputs: ", paste(inputs[!file.exists(inputs)], collapse = "; "))
  expected <- paste0("GSM36126", 35:40)
  meta <- read.csv(inputs["metadata"], check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("GSM", "group") %in% names(meta)) || nrow(meta) != 6L ||
      anyDuplicated(meta$GSM) || !setequal(meta$GSM, expected)) stop("Invalid six-sample manifest.")
  meta <- meta[match(expected, meta$GSM), , drop = FALSE]
  if (!identical(meta$group, rep(c("Vehicle", "Oxaliplatin"), each = 3L))) {
    stop("Manifest differs from the locked Control+OX minus Control contrast.")
  }
  expr <- as.matrix(read.csv(inputs["expression"], row.names = 1L, check.names = FALSE))
  if (!is.numeric(expr) || nrow(expr) < 5000L || anyDuplicated(rownames(expr)) ||
      ncol(expr) != 6L || !setequal(colnames(expr), expected) || any(!is.finite(expr))) {
    stop("Invalid RMA matrix: expected finite probe expression and the six locked GSMs.")
  }
  expr <- expr[, expected, drop = FALSE]
  annotation <- read.csv(inputs["annotation"], check.names = FALSE,
                          stringsAsFactors = FALSE, na.strings = c("", "NA"))
  if (!all(c("probe_id", "SYMBOL") %in% names(annotation))) stop("Missing probe_id/SYMBOL mapping.")
  # Read only the mapping columns: do not select probes by DE significance or |t|.
  pairs <- unique(annotation[, c("probe_id", "SYMBOL")])
  pairs$SYMBOL <- trimws(pairs$SYMBOL)
  pairs <- unique(pairs)
  pairs <- pairs[!is.na(pairs$probe_id) & !is.na(pairs$SYMBOL) &
                   nzchar(pairs$SYMBOL) & pairs$probe_id %in% rownames(expr), , drop = FALSE]
  ambiguous <- names(which(table(pairs$probe_id) > 1L))
  pairs <- pairs[!pairs$probe_id %in% ambiguous &
                   !grepl("[;,/|[:space:]]", pairs$SYMBOL), , drop = FALSE]
  if (anyDuplicated(pairs$probe_id)) stop("Ambiguous probe annotation remains.")
  pairs <- pairs[order(pairs$SYMBOL, pairs$probe_id), , drop = FALSE]
  mapped <- expr[match(pairs$probe_id, rownames(expr)), , drop = FALSE]
  # Average log2 RMA expression over probes per gene, independent of group effects.
  sums <- rowsum(mapped, group = pairs$SYMBOL, reorder = TRUE)
  probe_counts <- table(pairs$SYMBOL)
  expression <- sweep(sums, 1L, as.numeric(probe_counts[rownames(sums)]), "/")
  variable <- apply(expression, 1L, stats::sd) > 0
  expression <- expression[variable, , drop = FALSE]
  if (nrow(expression) < 5000L || anyDuplicated(rownames(expression)) || any(!is.finite(expression))) {
    stop("Too few finite, variable unique symbols after annotation.")
  }
  hallmark_table <- msigdbr::msigdbr(db_species = "HS", species = "Rattus norvegicus", collection = "H")
  if (!all(c("gs_name", "gene_symbol") %in% names(hallmark_table))) stop("Unexpected msigdbr columns.")
  sets <- lapply(split(hallmark_table$gene_symbol, hallmark_table$gs_name), unique)
  overlap <- vapply(sets, function(x) length(intersect(x, rownames(expression))), integer(1))
  if (length(sets) != 50L || any(overlap < 15L)) stop("Expected all 50 Hallmarks with >=15 measured genes.")
  # Gaussian is explicit for continuous log2 RMA; remaining score parameters
  # match the defaults used by the existing RNA-seq GSVA workflows.
  param <- GSVA::gsvaParam(exprData = expression, geneSets = sets, kcdf = "Gaussian",
                           minSize = 1, maxSize = Inf, tau = 1, maxDiff = TRUE, absRanking = FALSE)
  scores <- as.matrix(GSVA::gsva(param, verbose = TRUE, BPPARAM = BiocParallel::SerialParam()))
  scores <- gsva_check_scores(scores, expected)
  metadata <- data.frame(sample = expected, source_group = meta$group,
                          condition = rep(c("Control", "Neuropathy"), each = 3L))
  statistics <- gsva_fit_scores(scores, metadata)
  out <- file.path(root, "GSVA_Hallmark")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.csv(data.frame(pathway = rownames(scores), scores, check.names = FALSE),
            file.path(out, "GSVA_Hallmark_scores.csv"), row.names = FALSE)
  write.csv(metadata, file.path(out, "GSVA_sample_metadata.csv"), row.names = FALSE)
  write.csv(statistics, file.path(out, "GSVA_limma_all_50.csv"), row.names = FALSE)
  write.csv(statistics[match(gsva_selected_pathways, statistics$pathway), ],
            file.path(out, "GSVA_limma_selected_8.csv"), row.names = FALSE)
  pairs$gene_retained <- pairs$SYMBOL %in% rownames(expression)
  write.csv(pairs, file.path(out, "GSVA_probe_mapping.csv"), row.names = FALSE)
  write.csv(data.frame(pathway = names(sets), set_genes = lengths(sets), measured_genes = overlap),
            file.path(out, "GSVA_pathway_overlap.csv"), row.names = FALSE)
  provenance_columns <- intersect(c("gs_name", "gene_symbol", "db_gene_symbol", "db_version",
                                    "ortholog_sources", "num_ortholog_sources"), names(hallmark_table))
  write.csv(unique(hallmark_table[, provenance_columns, drop = FALSE]),
            file.path(out, "GSVA_rat_Hallmark_gene_sets.csv"), row.names = FALSE)
  write.csv(data.frame(input = unname(inputs), md5 = unname(tools::md5sum(inputs))),
            file.path(out, "GSVA_input_checksums.csv"), row.names = FALSE)
  saveRDS(scores, file.path(out, "GSVA_Hallmark_scores.rds"))
  writeLines(capture.output(sessionInfo()), file.path(out, "GSVA_sessionInfo.txt"))
  writeLines(c(
    "Rat bulk DRG, GSE126773: Control+OX minus Control; 3 versus 3.",
    "Existing RMA expression; no CEL download or differential-expression filtering.",
    "Single-symbol probes; duplicates collapsed by mean log2 RMA, independent of treatment effect.",
    paste("RMA probes:", nrow(expr)), paste("Ambiguous probe IDs excluded:", length(ambiguous)),
    paste("Retained variable gene symbols:", nrow(expression)),
    "Human MSigDB H projected to rat by msigdbr; exported gene-set membership/version fields.",
    "GSVA: Gaussian; minSize=1, maxSize=Inf; tau=1; maxDiff=TRUE; absRanking=FALSE; serial execution.",
    "limma: Neuropathy-Control, eBayes trend=FALSE/robust=FALSE; BH over all 50 sets.",
    "delta_GSVA is a score difference, not a log2 fold change of expression.",
    "Selected 8 sets come from GSEA; these are complementary analyses of the same data.",
    "Compare directions within studies; raw score magnitudes are not calibrated between studies."
  ), file.path(out, "GSVA_analysis_summary.txt"))
  message("GSE126773 GSVA complete: ", out)
  invisible(scores)
}
GSE126773_GSVA_result <- run_GSE126773_GSVA()
