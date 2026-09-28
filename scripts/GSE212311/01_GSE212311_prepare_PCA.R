#!/usr/bin/env Rscript
# Run from the manuscript-1 repository root:
# Rscript scripts/GSE212311/01_GSE212311_prepare_PCA.R
# GEO count-table import + blind exploratory QC. No DE model is fitted here.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})
options(timeout = max(600L, getOption("timeout")))
if (!file.exists("README.md") || !dir.exists("scripts"))
  stop("Run this script from the manuscript-1 repository root.")

out <- file.path("results", "GSE212311_CCI_L4L6_day11")
qc <- file.path(out, "PCA_6samples")
cache <- file.path("data", "source-cache", "GSE212311")
dir.create(qc, recursive = TRUE, showWarnings = FALSE)
dir.create(cache, recursive = TRUE, showWarnings = FALSE)
url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE212nnn/",
              "GSE212311/suppl/GSE212311_genes_fpkm_expression.txt.gz")
path <- file.path(cache, "GSE212311_genes_fpkm_expression.txt.gz")
expected_md5 <- "71033155b073b1af600c51a20cc9b33a"

download_verified <- function() {
  if (file.exists(path)) {
    current <- unname(tools::md5sum(path))
    if (identical(tolower(current), expected_md5)) return(invisible(NULL))
    stop("Existing GEO cache differs from the verified source (MD5 ", current,
         "). Investigate before replacing it: ", path)
  }
  tmp <- paste0(path, ".partial")
  for (attempt in seq_len(5L)) {
    if (file.exists(tmp)) unlink(tmp)
    message("Downloading GSE212311 GEO file, attempt ", attempt, "/5")
    status <- tryCatch(suppressWarnings(utils::download.file(
      url, tmp, mode = "wb", method = "libcurl", quiet = TRUE)),
      error = function(e) e)
    if (is.numeric(status) && length(status) == 1L && status == 0L &&
        file.exists(tmp) && file.info(tmp)$size > 1000000L &&
        identical(tolower(unname(tools::md5sum(tmp))), expected_md5)) {
      if (!file.rename(tmp, path)) stop("Cannot move verified GEO download to cache.")
      return(invisible(NULL))
    }
    if (attempt < 5L) Sys.sleep(attempt)
  }
  stop("GEO download could not be verified. Source: ", url,
       "; expected MD5: ", expected_md5)
}
download_verified()

# The source is ordered CCI then Sham; analysis columns are locked Sham then CCI.
samples <- data.frame(
  gsm = c("GSM6523751", "GSM6523752", "GSM6523753",
          "GSM6523748", "GSM6523749", "GSM6523750"),
  title = c("Sham_rep1", "Sham_rep2", "Sham_rep3",
            "CCI_rep1", "CCI_rep2", "CCI_rep3"),
  condition = c(rep("Sham", 3L), rep("CCI", 3L)),
  source_count_column = c(paste0("count.Sham_", 1:3),
                          paste0("count.CCI_", 1:3)),
  source_fpkm_column = c(paste0("FPKM.Sham_", 1:3),
                         paste0("FPKM.CCI_", 1:3)),
  tissue = "ipsilateral_L4-L6_DRG", day = 11L,
  strain = "Sprague-Dawley", sex = "Male",
  stringsAsFactors = FALSE
)
if (nrow(samples) != 6L || anyDuplicated(samples$gsm) ||
    any(table(samples$condition) != 3L)) stop("Locked design changed.")

tab <- read.delim(gzfile(path), check.names = FALSE, quote = "",
                  comment.char = "", stringsAsFactors = FALSE,
                  na.strings = character())
expected_columns <- c("gene_id", "chr", "start", "end", "strand",
                      "gene_name", "transcript_id", "GO", "KEGG", "KO_ENTRY",
                      "EC", "Description", "trans_type",
                      paste0("FPKM.CCI_", 1:3), paste0("FPKM.Sham_", 1:3),
                      paste0("count.CCI_", 1:3), paste0("count.Sham_", 1:3))
if (nrow(tab) != 43350L || !identical(names(tab), expected_columns))
  stop("Unexpected GEO feature inventory or column layout; audit the source.")
genes <- as.character(tab$gene_id)
if (anyNA(genes) || anyDuplicated(genes) ||
    any(!grepl("^(ENSRNOG[0-9]+|MSTRG\\.[0-9]+)$", genes)))
  stop("Invalid or duplicated gene identifiers.")

count_input <- tab[, samples$source_count_column, drop = FALSE]
counts <- as.matrix(count_input)
if (!is.numeric(counts) || anyNA(counts) || any(!is.finite(counts)) ||
    any(counts < 0) || any(counts != round(counts)) ||
    any(counts > .Machine$integer.max))
  stop("GEO count columns do not contain valid nonnegative integer counts.")
storage.mode(counts) <- "integer"
dimnames(counts) <- list(genes, samples$gsm)
if (any(colSums(counts) == 0)) stop("At least one sample has zero counts.")
fpkm <- as.matrix(tab[, samples$source_fpkm_column, drop = FALSE])
if (!is.numeric(fpkm) || anyNA(fpkm) ||
    any(!is.finite(fpkm)) || any(fpkm < 0))
  stop("GEO FPKM columns are invalid.")

samples$source_url <- url
samples$source_md5 <- expected_md5
samples$assigned_count_sum <- as.numeric(colSums(counts))
samples$genes_nonzero <- as.integer(colSums(counts > 0L))
write.csv(samples, file.path(out, "locked_primary_samples.csv"), row.names = FALSE)
write.csv(counts, file.path(out, "raw_count_matrix.csv"))
write.csv(samples[, c("gsm", "title", "condition", "source_count_column",
                      "source_fpkm_column", "assigned_count_sum", "genes_nonzero",
                      "source_md5")],
          file.path(out, "sample_input_audit.csv"), row.names = FALSE)

# Blind transform: no group information enters the PCA.
keep <- rowSums(counts >= 10L) >= 3L
if (sum(keep) < 1000L) stop("Too few genes remain after count filtering.")
meta <- samples
rownames(meta) <- meta$gsm
meta$condition <- factor(meta$condition, levels = c("Sham", "CCI"))
dds <- DESeqDataSetFromMatrix(countData = counts[keep, , drop = FALSE],
                              colData = meta, design = ~ 1)
vsd <- vst(dds, blind = TRUE)
mat <- SummarizedExperiment::assay(vsd)
pc <- prcomp(t(mat), center = TRUE, scale. = FALSE)
variance <- 100 * pc$sdev^2 / sum(pc$sdev^2)
coords <- data.frame(gsm = meta$gsm, condition = meta$condition,
                     PC1 = unname(pc$x[meta$gsm, 1L]),
                     PC2 = unname(pc$x[meta$gsm, 2L]))
write.csv(coords, file.path(qc, "PCA_coordinates.csv"), row.names = FALSE)
write.csv(data.frame(PC = paste0("PC", seq_along(variance)),
                     variance_percent = variance),
          file.path(qc, "PCA_variance.csv"), row.names = FALSE)
write.csv(as.matrix(dist(t(mat))), file.path(qc, "sample_distances.csv"))

p <- ggplot(coords, aes(PC1, PC2, colour = condition, label = gsm)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8, show.legend = FALSE, check_overlap = TRUE) +
  scale_colour_manual(values = c(Sham = "#2171B5", CCI = "#D95F02")) +
  labs(title = "GSE212311: rat ipsilateral L4-L6 DRG, day 11",
       subtitle = "Blind VST PCA: 3 Sham versus 3 CCI",
       x = sprintf("PC1 (%.1f%%)", variance[1L]),
       y = sprintf("PC2 (%.1f%%)", variance[2L])) +
  theme_bw(base_size = 12)
ggsave(file.path(qc, "PCA_6samples.png"), p,
       width = 8, height = 6, dpi = 300)

writeLines(c(
  "GSE212311 independent CCI vs Sham analysis; 3 biological replicates per group.",
  "Adult male Sprague-Dawley rat ipsilateral L4-L6 DRG; day 11.",
  paste("Source URL:", url), paste("Source MD5:", expected_md5),
  paste("Input genes:", nrow(counts), "including",
        sum(startsWith(genes, "ENSRNOG")), "annotated ENSRNOG and",
        sum(startsWith(genes, "MSTRG.")), "MSTRG loci."),
  paste("PCA retained >=10 raw counts in >=3 samples:", sum(keep), "genes."),
  "GEO file contains FPKM and integer count columns; only counts enter DESeq2/VST.",
  "Sham nerve exposure without ligation is described in the source publication.",
  "No sample was removed and no DE model was fitted in this step."
), file.path(qc, "QC_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(qc, "sessionInfo.txt"))
message("GSE212311 input and PCA complete: ", nrow(counts), " genes x 6 samples; ", qc)
