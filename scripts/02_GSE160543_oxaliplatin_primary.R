#!/usr/bin/env Rscript

# Primary bulk-DRG analysis: GSE160543 Oxaliplatin versus Vehicle (4 versus 4).
# Run from the repository root: Rscript scripts/02_GSE160543_oxaliplatin_primary.R
# The GEO files contain fractional expected_count estimates and unequal gene
# inventories. This script rounds expected counts and analyzes complete-case
# gene IDs only; it audits the genes excluded for missingness. No missing value
# is interpreted as a zero count. Paclitaxel samples are not analyzed.

options(stringsAsFactors = FALSE, timeout = max(300L, getOption("timeout")))

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install("DESeq2", ask = FALSE, update = FALSE)
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("DESeq2 could not be loaded. Install it with BiocManager::install('DESeq2').")
}

samples <- data.frame(
  gsm = c("GSM4875003", "GSM4875004", "GSM4875005", "GSM4875006",
          "GSM4875011", "GSM4875012", "GSM4875013", "GSM4875014"),
  group = factor(rep(c("Vehicle", "Oxaliplatin"), each = 4L),
                 levels = c("Vehicle", "Oxaliplatin")),
  stringsAsFactors = FALSE
)
stopifnot(nrow(samples) == 8L, !anyDuplicated(samples$gsm))
rownames(samples) <- samples$gsm

cache_dir <- file.path("data", "source-cache")
output_dir <- file.path("results", "GSE160543_Oxaliplatin_vs_Vehicle")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
tar_path <- file.path(cache_dir, "GSE160543_RAW.tar")
geo_url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/",
                  "GSE160nnn/GSE160543/suppl/GSE160543_RAW.tar")

if (!file.exists(tar_path) || file.info(tar_path)$size == 0) {
  message("Downloading official GEO supplementary archive...")
  tmp <- paste0(tar_path, ".partial")
  if (file.exists(tmp)) unlink(tmp)
  tryCatch(
    download.file(geo_url, destfile = tmp, mode = "wb", method = "libcurl",
                  quiet = FALSE),
    error = function(e) {
      if (file.exists(tmp)) unlink(tmp)
      stop("GEO download failed: ", conditionMessage(e))
    }
  )
  if (!file.exists(tmp) || file.info(tmp)$size < 1000000) {
    stop("Downloaded archive appears incomplete; remove the .partial file and retry.")
  }
  if (!file.rename(tmp, tar_path)) stop("Could not move downloaded archive into cache.")
}

archive_members <- tryCatch(untar(tar_path, list = TRUE),
                            error = function(e) stop("Invalid GEO tar archive: ", conditionMessage(e)))
expected_names <- c(paste0(samples$gsm[1:4], "_Vehicle_", 1:4, ".txt.gz"),
                    paste0(samples$gsm[5:8], "_oxaliplatin_", 1:4, ".txt.gz"))
if (!all(expected_names %in% archive_members)) {
  stop("The tar archive does not contain all eight locked GSM files. Missing: ",
       paste(setdiff(expected_names, archive_members), collapse = ", "))
}
raw_dir <- file.path(cache_dir, "GSE160543_raw_files")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
missing_files <- expected_names[!file.exists(file.path(raw_dir, expected_names))]
if (length(missing_files)) {
  untar(tar_path, files = missing_files, exdir = raw_dir)
}
if (!all(file.exists(file.path(raw_dir, expected_names)))) {
  stop("Extraction failed for one or more selected samples.")
}

read_sample <- function(path) {
  # header=FALSE avoids a malformed header in one GEO file with extra tabs.
  tab <- read.delim(gzfile(path), header = FALSE, skip = 1L, sep = "\t",
                    quote = "", fill = TRUE, check.names = FALSE,
                    comment.char = "", stringsAsFactors = FALSE)
  if (ncol(tab) < 5L) stop("Expected at least five columns in ", path)
  ids <- trimws(as.character(tab[[2L]]))
  symbols <- trimws(as.character(tab[[1L]]))
  values <- suppressWarnings(as.numeric(tab[[5L]]))
  if (!length(ids) || anyNA(ids) || any(ids == "") || anyDuplicated(ids)) {
    stop("Missing or duplicated gene IDs in ", path)
  }
  if (anyNA(values) || any(!is.finite(values)) || any(values < 0)) {
    stop("Missing, non-finite or negative expected_count in ", path)
  }
  data.frame(gene_id = ids, symbol = symbols, expected_count = values,
             stringsAsFactors = FALSE)
}

message("Reading selected GEO files...")
tables <- lapply(file.path(raw_dir, expected_names), read_sample)
names(tables) <- samples$gsm
all_ids <- sort(unique(unlist(lapply(tables, `[[`, "gene_id"), use.names = FALSE)))
complete_ids <- sort(Reduce(intersect, lapply(tables, `[[`, "gene_id")))
if (length(complete_ids) < 10000L) {
  stop("Too few genes are present in all eight samples: ", length(complete_ids))
}

sample_audit <- data.frame(
  gsm = samples$gsm, group = as.character(samples$group),
  file = expected_names,
  n_gene_ids = vapply(tables, nrow, integer(1)),
  n_fractional = vapply(tables, function(x)
    sum(abs(x$expected_count - round(x$expected_count)) > 1e-8), integer(1)),
  estimated_library_sum = vapply(tables, function(x)
    sum(x$expected_count), numeric(1)),
  n_missing_from_union = vapply(tables, function(x)
    length(setdiff(all_ids, x$gene_id)), integer(1)),
  stringsAsFactors = FALSE
)
write.csv(sample_audit, file.path(output_dir, "sample_input_audit.csv"), row.names = FALSE)

# Explicitly measure the scale of unequal gene inventories before dropping IDs.
max_when_present <- setNames(numeric(length(all_ids)), all_ids)
presence <- setNames(integer(length(all_ids)), all_ids)
for (tab in tables) {
  idx <- match(tab$gene_id, all_ids)
  presence[idx] <- presence[idx] + 1L
  max_when_present[idx] <- pmax(max_when_present[idx], tab$expected_count)
}
missing_ids <- all_ids[presence < nrow(samples)]
missing_audit <- data.frame(
  gene_id = missing_ids,
  n_samples_present = unname(presence[missing_ids]),
  max_expected_count_when_present = unname(max_when_present[missing_ids]),
  stringsAsFactors = FALSE
)
write.csv(missing_audit, file.path(output_dir, "excluded_gene_missingness.csv"),
          row.names = FALSE)

raw_counts <- matrix(0, nrow = length(complete_ids), ncol = nrow(samples),
                     dimnames = list(complete_ids, samples$gsm))
symbol_matrix <- matrix("", nrow = length(complete_ids), ncol = nrow(samples))
for (j in seq_along(tables)) {
  tab <- tables[[j]]
  idx <- match(complete_ids, tab$gene_id)
  if (anyNA(idx)) stop("Internal alignment failure for ", names(tables)[j])
  raw_counts[, j] <- tab$expected_count[idx]
  symbol_matrix[, j] <- tab$symbol[idx]
}
discordant_symbols <- rowSums(symbol_matrix != symbol_matrix[, 1L]) > 0L
if (any(discordant_symbols)) {
  stop("Gene IDs have discordant symbols across samples: ", sum(discordant_symbols))
}
symbols <- symbol_matrix[, 1L]
symbols[is.na(symbols)] <- ""

# GEO supplies RSEM-style fractional expected counts. Rounding is documented;
# FPKM is never used as a DESeq2 input. Round before applying the count filter.
counts <- round(raw_counts)
storage.mode(counts) <- "integer"
keep <- rowSums(counts) >= 10L
counts <- counts[keep, , drop = FALSE]
symbols <- symbols[keep]
if (nrow(counts) < 5000L || anyNA(counts) || any(counts < 0L)) {
  stop("Unexpected filtered count matrix.")
}

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts, colData = samples, design = ~ group
)
dds <- DESeq2::DESeq(dds, quiet = TRUE)
res <- DESeq2::results(dds,
                       contrast = c("group", "Oxaliplatin", "Vehicle"),
                       alpha = 0.05)
out <- data.frame(gene_id = rownames(res), symbol = symbols,
                  as.data.frame(res), check.names = FALSE)
out <- out[order(is.na(out$padj), out$padj, -abs(out$log2FoldChange)), , drop = FALSE]
rownames(out) <- NULL
sig <- !is.na(out$padj) & out$padj < 0.05
sig_fc <- sig & !is.na(out$log2FoldChange) & abs(out$log2FoldChange) > 1

write.csv(out, file.path(output_dir, "DE_all_genes.csv"), row.names = FALSE, na = "NA")
write.csv(out[sig, , drop = FALSE], file.path(output_dir, "DE_FDR_lt_0.05.csv"),
          row.names = FALSE, na = "NA")
write.csv(out[sig_fc, , drop = FALSE],
          file.path(output_dir, "DE_FDR_lt_0.05_abs_log2FC_gt_1.csv"),
          row.names = FALSE, na = "NA")
write.csv(out[tolower(out$symbol) %in% c("cdk1", "cdc20", "cdkn1a"), , drop = FALSE],
          file.path(output_dir, "three_original_candidate_genes.csv"),
          row.names = FALSE, na = "NA")
norm <- DESeq2::counts(dds, normalized = TRUE)
write.csv(data.frame(gene_id = rownames(norm), symbol = symbols, norm,
                     check.names = FALSE),
          file.path(output_dir, "normalized_counts_complete_case.csv"), row.names = FALSE)

vst <- DESeq2::vst(dds, blind = TRUE)
vmat <- SummarizedExperiment::assay(vst)

write.csv(
  vmat,
  file.path(
    output_dir,
    "vst_expression_matrix.csv"
  ),
  row.names = TRUE
)
                                 
pcs <- prcomp(t(vmat), center = TRUE, scale. = FALSE)
variance <- 100 * pcs$sdev^2 / sum(pcs$sdev^2)
pca_table <- data.frame(gsm = samples$gsm, group = as.character(samples$group),
                        PC1 = pcs$x[, 1L], PC2 = pcs$x[, 2L])
write.csv(pca_table, file.path(output_dir, "PCA_coordinates.csv"), row.names = FALSE)
png(file.path(output_dir, "PCA.png"), width = 1400, height = 1100, res = 160)
plot(pcs$x[, 1L], pcs$x[, 2L],
     col = ifelse(samples$group == "Vehicle", "#2474A6", "#D35400"),
     pch = 19, cex = 1.5,
     xlab = sprintf("PC1 (%.1f%%)", variance[1L]),
     ylab = sprintf("PC2 (%.1f%%)", variance[2L]),
     main = "GSE160543: Oxaliplatin versus Vehicle")
text(pcs$x[, 1L], pcs$x[, 2L], labels = samples$gsm, pos = 3, cex = 0.65)
legend("topright", legend = c("Vehicle", "Oxaliplatin"),
       col = c("#2474A6", "#D35400"), pch = 19, bty = "n")
dev.off()

distance <- as.matrix(dist(t(vmat)))
dimnames(distance) <- list(samples$gsm, samples$gsm)
write.csv(distance, file.path(output_dir, "sample_distances.csv"))
png(file.path(output_dir, "sample_distance_heatmap.png"),
    width = 1300, height = 1100, res = 160)
heatmap(distance, Rowv = NA, Colv = NA, scale = "none",
        col = grDevices::colorRampPalette(c("#FFFFFF", "#22577A"))(100),
        margins = c(10, 10), main = "VST sample distances")
dev.off()

summary_lines <- c(
  "GSE160543 primary differential-expression analysis",
  "Contrast: Oxaliplatin versus Vehicle; four independent GEO samples per group.",
  "Paclitaxel samples: excluded from the model.",
  paste("GEO archive:", geo_url),
  paste("Downloaded archive MD5:", unname(tools::md5sum(tar_path))),
  "Input column: expected_count (fractional estimates), rounded to nearest integer.",
  "Gene IDs: complete-case intersection of the eight sample files; missing is never set to zero.",
  paste("Gene IDs in any selected file:", length(all_ids)),
  paste("Gene IDs present in all eight:", length(complete_ids)),
  paste("Excluded for unequal gene inventories:", length(missing_ids)),
  paste("Excluded genes with max count >=100 in an available sample:",
        sum(missing_audit$max_expected_count_when_present >= 100)),
  paste("Genes retained after total rounded count >=10:", nrow(counts)),
  paste("FDR < 0.05:", sum(sig)),
  paste("FDR < 0.05 and absolute log2FC > 1:", sum(sig_fc)),
  paste("  up:", sum(sig_fc & out$log2FoldChange > 0, na.rm = TRUE)),
  paste("  down:", sum(sig_fc & out$log2FoldChange < 0, na.rm = TRUE)),
  "Important: Sample-specific gene absence is substantial and includes high-count genes.",
  "The complete-case result is an auditable sensitivity-limited analysis, not proof that excluded genes did not change.",
  "Bulk DRG mRNA does not establish cell type, protein activity, senescence or causality."
)
writeLines(summary_lines, file.path(output_dir, "analysis_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message(paste(summary_lines, collapse = "\n"))
