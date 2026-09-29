#!/usr/bin/env Rscript
# Audit original GEO feature annotations before interpreting GSE212311 pathways.
# Run from the repository root after steps 01 and 02:
# source("scripts/GSE212311/05_GSE212311_feature_audit.R")
# This script does not reinterpret MSTRG gene_name values as validated symbols.

root <- file.path("results", "GSE212311_CCI_L4L6_day11")
source_file <- file.path("data", "source-cache", "GSE212311",
                         "GSE212311_genes_fpkm_expression.txt.gz")
de_file <- file.path(root, "DESeq2",
                     "GSE212311_CCI_vs_Sham_DE_all_genes.csv")
counts_file <- file.path(root, "raw_count_matrix.csv")
out <- file.path(root, "feature_audit")
for (p in c(source_file, de_file, counts_file)) {
  if (!file.exists(p)) stop("Required input missing; run steps 01 and 02: ", p)
}
expected_md5 <- "71033155b073b1af600c51a20cc9b33a"
if (!identical(tolower(unname(tools::md5sum(source_file))), expected_md5)) {
  stop("GEO source MD5 differs from the source verified in step 01.")
}
tab <- read.delim(gzfile(source_file), check.names = FALSE, quote = "",
                  comment.char = "", stringsAsFactors = FALSE,
                  na.strings = character())
needed <- c("gene_id", "gene_name", "trans_type", "Description",
            paste0("count.Sham_", 1:3), paste0("count.CCI_", 1:3))
if (nrow(tab) != 43350L || !all(needed %in% names(tab)) ||
    anyDuplicated(tab$gene_id)) {
  stop("Source feature layout or IDs differ from step 01.")
}
counts <- read.csv(counts_file, row.names = 1L, check.names = FALSE)
de <- read.csv(de_file, check.names = FALSE)
if (!identical(rownames(counts), tab$gene_id) ||
    !all(c("gene_id", "gene_class", "padj", "log2FoldChange") %in% names(de)) ||
    anyDuplicated(de$gene_id)) {
  stop("Step 01/02 outputs do not match the verified GEO feature IDs.")
}
source_columns <- c(paste0("count.Sham_", 1:3),
                    paste0("count.CCI_", 1:3))
source_counts <- as.matrix(tab[, source_columns, drop = FALSE])
if (!isTRUE(all.equal(unname(as.matrix(counts)), unname(source_counts),
                      check.attributes = FALSE, tolerance = 0))) {
  stop("Locked count matrix differs from GEO integer count columns.")
}
kept <- rowSums(source_counts >= 10L) >= 3L
if (!setequal(de$gene_id, tab$gene_id[kept])) {
  stop("DESeq2 tested features disagree with the step 01 count filter.")
}

is_ens <- startsWith(tab$gene_id, "ENSRNOG")
is_mstrg <- grepl("^MSTRG\\.[0-9]+$", tab$gene_id)
if (any(!is_ens & !is_mstrg)) stop("Unexpected source feature ID class.")
name <- trimws(as.character(tab$gene_name))
name_present <- !is.na(name) & nzchar(name) &
  !tolower(name) %in% c("na", "n/a", "none", "null", "-", ".", "--") &
  name != tab$gene_id
classes <- list(ENSRNOG = is_ens, MSTRG = is_mstrg)
inventory <- do.call(rbind, lapply(names(classes), function(label) {
  idx <- classes[[label]]
  data.frame(
    feature_class = label,
    source_features = sum(idx),
    zero_count_features = sum(idx & rowSums(source_counts) == 0),
    tested_features = sum(idx & kept),
    significant_features = sum(idx & tab$gene_id %in%
      de$gene_id[!is.na(de$padj) & de$padj < 0.05]),
    source_features_with_nonempty_gene_name = sum(idx & name_present),
    tested_features_with_nonempty_gene_name = sum(idx & kept & name_present),
    total_counts_across_six_samples = sum(source_counts[idx, , drop = FALSE]),
    stringsAsFactors = FALSE
  )
}))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(inventory, file.path(out, "feature_class_inventory.csv"),
          row.names = FALSE)

# These are the original source labels, not a validated MSTRG-to-gene map.
sig_mstrg <- de[de$gene_class == "MSTRG_novel_locus" &
                  !is.na(de$padj) & de$padj < 0.05, , drop = FALSE]
at <- match(sig_mstrg$gene_id, tab$gene_id)
if (anyNA(at)) stop("Significant MSTRG IDs missing from source.")
sig_mstrg$source_gene_name <- tab$gene_name[at]
sig_mstrg$source_trans_type <- tab$trans_type[at]
sig_mstrg$source_description <- tab$Description[at]
sig_mstrg$source_name_present <- name_present[at]
write.csv(sig_mstrg, file.path(out, "significant_MSTRG_source_labels.csv"),
          row.names = FALSE)

named <- tab[is_mstrg & name_present, c("gene_id", "gene_name",
                                       "trans_type", "Description")]
write.csv(named, file.path(out, "MSTRG_source_names.csv"),
          row.names = FALSE)
writeLines(c(
  "GSE212311 source feature audit; no MSTRG names were accepted as official gene symbols.",
  paste("Verified source MD5:", expected_md5),
  paste("Tested Ensembl/MSTRG:", inventory$tested_features[1],
        "/", inventory$tested_features[2]),
  paste("Ensembl fraction of six-sample raw counts:",
        round(100 * inventory$total_counts_across_six_samples[1] /
              sum(inventory$total_counts_across_six_samples), 3), "%"),
  paste("Significant MSTRG with a nonempty source gene_name:",
        sum(sig_mstrg$source_name_present), "of", nrow(sig_mstrg)),
  "Inspect source labels and duplicate names before changing gene-level or pathway analysis."
), file.path(out, "feature_audit_summary.txt"))
message("GSE212311 source annotation audit complete: ", out)
