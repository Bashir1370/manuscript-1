#!/usr/bin/env Rscript
# GSE212311: annotate significant CCI versus Sham DESeq2 features.
# Run from the manuscript-1 repository root after step 02.
# source("scripts/GSE212311/03_GSE212311_annotation.R")

required_packages <- c("AnnotationDbi", "org.Rn.eg.db")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                             logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install these R packages before annotation: ",
       paste(missing_packages, collapse = ", "))
}

base_dir <- file.path("results", "GSE212311_CCI_L4L6_day11")
input_file <- file.path(base_dir, "DESeq2",
                        "GSE212311_significant_DEGs_FDR_lt_0.05.csv")
out_dir <- file.path(base_dir, "annotation")
if (!file.exists(input_file)) {
  stop("DESeq2 DEG file not found. Run step 02 from the repository root: ",
       input_file)
}

deg <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c("gene_id", "gene_class", "log2FoldChange", "padj")
if (!all(required_columns %in% names(deg))) {
  stop("DESeq2 output lacks: ",
       paste(setdiff(required_columns, names(deg)), collapse = ", "))
}
if (anyNA(deg$gene_id) || any(!nzchar(deg$gene_id)) ||
    anyDuplicated(deg$gene_id)) {
  stop("DESeq2 output contains missing, empty, or duplicated gene_id values.")
}
if (any(!grepl("^(ENSRNOG[0-9]+(\\.[0-9]+)?|MSTRG\\.[0-9]+)$",
               deg$gene_id))) {
  stop("Unexpected gene_id format in DESeq2 output.")
}

is_ensembl <- startsWith(deg$gene_id, "ENSRNOG")
if (any(is_ensembl != (deg$gene_class == "annotated_Ensembl"))) {
  stop("gene_class disagrees with the gene_id prefix in DESeq2 output.")
}
# A version suffix is not part of the org.Rn.eg.db ENSEMBL key.
ensembl <- sub("\\.[0-9]+$", "", deg$gene_id[is_ensembl])
if (anyDuplicated(ensembl)) {
  stop("Multiple DEG rows collapse to the same Ensembl ID after version removal.")
}
keys <- unique(ensembl)

if (length(keys)) {
  mapping <- suppressMessages(AnnotationDbi::select(
    org.Rn.eg.db::org.Rn.eg.db,
    keys = keys, keytype = "ENSEMBL", columns = c("SYMBOL", "ENTREZID")
  ))
  mapping <- unique(mapping[, c("ENSEMBL", "SYMBOL", "ENTREZID"),
                            drop = FALSE])
  if (!all(keys %in% mapping$ENSEMBL)) {
    stop("Annotation database did not return every requested Ensembl key.")
  }
} else {
  mapping <- data.frame(ENSEMBL = character(), SYMBOL = character(),
                        ENTREZID = character())
}

# AnnotationDbi can return several rows per ID. Retain every original DEG
# exactly once. Assign a symbol or Entrez ID only when that field is unique.
by_id <- split(mapping, mapping$ENSEMBL)
single_value <- function(id, field) {
  values <- unique(as.character(by_id[[id]][[field]]))
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 1L) values else NA_character_
}
symbol <- vapply(keys, single_value, character(1), field = "SYMBOL")
entrez <- vapply(keys, single_value, character(1), field = "ENTREZID")
symbol_count <- vapply(keys, function(id) {
  values <- unique(as.character(by_id[[id]]$SYMBOL))
  sum(!is.na(values) & nzchar(values))
}, integer(1))

annotated <- deg
annotated$SYMBOL <- NA_character_
annotated$ENTREZID <- NA_character_
annotated$annotation_status <- "MSTRG_novel_locus"
if (length(keys)) {
  idx <- which(is_ensembl)
  pos <- match(ensembl, keys)
  annotated$SYMBOL[idx] <- symbol[pos]
  annotated$ENTREZID[idx] <- entrez[pos]
  annotated$annotation_status[idx] <- ifelse(
    symbol_count[pos] > 1L, "ambiguous_symbol",
    ifelse(is.na(symbol[pos]), "no_symbol", "mapped_symbol")
  )
}
if (nrow(annotated) != nrow(deg) ||
    !identical(annotated$gene_id, deg$gene_id)) {
  stop("Annotation changed the number or order of DEG features.")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(annotated, file.path(out_dir, "GSE212311_DEG_annotated.csv"),
          row.names = FALSE)
write.csv(mapping, file.path(out_dir, "annotation_mapping_table.csv"),
          row.names = FALSE)
write.csv(annotated[is.na(annotated$SYMBOL), , drop = FALSE],
          file.path(out_dir, "unmapped_genes.csv"), row.names = FALSE)

mapped <- sum(!is.na(annotated$SYMBOL))
summary_text <- c(
  paste("Total DEG features:", nrow(deg)),
  paste("Ensembl DEG features:", length(keys)),
  paste("MSTRG novel loci (no Ensembl lookup):", sum(!is_ensembl)),
  paste("Ensembl features with an unambiguous symbol:", mapped),
  paste("Ensembl features with ambiguous symbols:",
        sum(annotated$annotation_status == "ambiguous_symbol")),
  paste("Ensembl features without a symbol:",
        sum(annotated$annotation_status == "no_symbol")),
  paste("Features without an unambiguous symbol:",
        sum(is.na(annotated$SYMBOL))),
  paste("Ensembl symbol mapping percentage:",
        if (length(keys)) round(100 * mapped / length(keys), 2) else NA_real_)
)
writeLines(summary_text, file.path(out_dir, "annotation_summary.txt"))
message("GSE212311 annotation completed: ", mapped, " of ", length(keys),
        " Ensembl DEGs have unambiguous symbols. Output: ", out_dir)
