#!/usr/bin/env Rscript

# Run from the repository root:
# Rscript scripts/05_GSE246156_prepare_NC_L5_day7.R
# This script locks and audits the primary NC contrast. It does not run DESeq2.

options(stringsAsFactors = FALSE, timeout = max(600L, getOption("timeout")))

manifest_path <- file.path("data", "metadata", "geo_sample_manifest.csv")
if (!file.exists(manifest_path)) stop("Run this script from the repository root.")
manifest <- read.csv(manifest_path, check.names = FALSE, na.strings = c("", "NA"))
required <- c("accession", "gsm", "condition", "tissue", "drg_level", "day",
              "replicate_label", "sra_accession")
if (!all(required %in% names(manifest))) {
  stop("Missing manifest columns: ", paste(setdiff(required, names(manifest)), collapse = ", "))
}

nc <- manifest[manifest$accession == "GSE246156", , drop = FALSE]
drg <- nc[nc$tissue == "DRG", , drop = FALSE]
if (nrow(nc) != 48L || nrow(drg) != 36L || anyDuplicated(nc$gsm) ||
    anyNA(drg[, c("gsm", "condition", "drg_level", "day", "replicate_label")])) {
  stop("GSE246156 sample inventory differs from the audited 48/36 design.")
}

cells <- expand.grid(condition = c("Sham", "Compression"),
                     drg_level = c("L4", "L5", "L6"),
                     day = c(3L, 7L), stringsAsFactors = FALSE)
cells$n_samples <- vapply(seq_len(nrow(cells)), function(i) {
  sum(drg$condition == cells$condition[i] &
        drg$drg_level == cells$drg_level[i] &
        drg$day == cells$day[i])
}, integer(1))
if (any(cells$n_samples != 3L) ||
    any(!drg$condition %in% cells$condition) ||
    any(!drg$drg_level %in% cells$drg_level) ||
    any(!drg$day %in% cells$day)) {
  stop("The DRG condition x level x day design is not 12 cells of 3 samples.")
}
for (i in seq_len(nrow(cells))) {
  rows <- drg$condition == cells$condition[i] &
    drg$drg_level == cells$drg_level[i] & drg$day == cells$day[i]
  if (!identical(sort(as.integer(drg$replicate_label[rows])), 1:3)) {
    stop("Unexpected replicate labels in ", cells$condition[i], " ",
         cells$drg_level[i], " day ", cells$day[i])
  }
}

# Prespecified primary comparison: DRG at the compressed L5 level on day 7.
# This is 3 versus 3, not a pooled 9 versus 9 across spinal levels.
primary_ids <- c("GSM7863794", "GSM7863795", "GSM7863796",
                 "GSM7863770", "GSM7863771", "GSM7863772")
primary <- drg[match(primary_ids, drg$gsm), , drop = FALSE]
expected_group <- rep(c("Sham", "Compression"), each = 3L)
if (anyNA(primary$gsm) || !identical(primary$gsm, primary_ids) ||
    !identical(primary$condition, expected_group) ||
    any(primary$drg_level != "L5") || any(primary$day != 7L)) {
  stop("The six locked L5 day-7 GSMs disagree with the audited manifest.")
}
primary$analysis_role <- "primary_L5_day7"
primary <- primary[, c("gsm", "condition", "tissue", "drg_level", "day",
                       "replicate_label", "sra_accession", "analysis_role")]

output_dir <- file.path("results", "GSE246156_NC_L5_day7")
cache_dir <- file.path("data", "source-cache", "GSE246156")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(cells, file.path(output_dir, "DRG_strata_inventory.csv"), row.names = FALSE)
write.csv(primary, file.path(output_dir, "locked_primary_samples.csv"), row.names = FALSE)

# File names are taken from the GSE246156 GEO series-matrix supplementary_file_1
# field. Never pair count files to groups by alphabetical file order.
source_names <- c(
  GSM7863794 = "GSM7863794_featurecounts_D7dui1L5.txt.gz",
  GSM7863795 = "GSM7863795_featurecounts_D7dui2L5.txt.gz",
  GSM7863796 = "GSM7863796_featurecounts_D7dui3L5.txt.gz",
  GSM7863770 = "GSM7863770_featurecounts_D7shi1L5.txt.gz",
  GSM7863771 = "GSM7863771_featurecounts_D7shi2L5.txt.gz",
  GSM7863772 = "GSM7863772_featurecounts_D7shi3L5.txt.gz"
)
stopifnot(identical(names(source_names), primary$gsm))

get_source_url <- function(gsm, filename) {
  paste0("https://ftp.ncbi.nlm.nih.gov/geo/samples/",
         substr(gsm, 1L, 7L), "nnn/", gsm, "/suppl/", filename)
}
source_urls <- mapply(get_source_url, primary$gsm, unname(source_names),
                      USE.NAMES = FALSE)
paths <- file.path(cache_dir, unname(source_names))

for (i in seq_along(paths)) {
  if (file.exists(paths[i]) && file.info(paths[i])$size > 100L) next
  tmp <- paste0(paths[i], ".partial")
  if (file.exists(tmp)) unlink(tmp)
  message("Downloading ", primary$gsm[i], " from GEO...")
  tryCatch(download.file(source_urls[i], tmp, mode = "wb", method = "libcurl"),
           error = function(e) {
             if (file.exists(tmp)) unlink(tmp)
             stop("Download failed for ", primary$gsm[i], ": ", conditionMessage(e))
           })
  if (!file.exists(tmp) || file.info(tmp)$size <= 100L ||
      !file.rename(tmp, paths[i])) {
    stop("Incomplete download for ", primary$gsm[i], ". Remove the partial file and retry.")
  }
}

read_featurecounts <- function(path) {
  tab <- tryCatch(read.delim(gzfile(path), header = TRUE, sep = "\t",
                             comment.char = "#", quote = "", check.names = FALSE),
                  error = function(e) stop("Cannot read ", path, ": ", conditionMessage(e)))
  if (ncol(tab) < 2L || nrow(tab) < 1000L) stop("Unexpected featureCounts layout: ", path)
  ids <- as.character(tab[[1L]])
  values <- suppressWarnings(as.numeric(as.character(tab[[ncol(tab)]])))
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) ||
      anyNA(values) || any(!is.finite(values)) || any(values < 0) ||
      any(abs(values - round(values)) > 1e-8) || any(values > .Machine$integer.max)) {
    stop("Invalid gene IDs or non-integer raw counts in ", path)
  }
  data.frame(gene_id = ids, count = as.integer(values))
}
tables <- lapply(paths, read_featurecounts)
names(tables) <- primary$gsm
all_ids <- sort(unique(unlist(lapply(tables, `[[`, "gene_id"), use.names = FALSE)))
shared_ids <- Reduce(intersect, lapply(tables, `[[`, "gene_id"))
audit <- data.frame(gsm = primary$gsm, condition = primary$condition,
                    source_url = source_urls, file = paths,
                    md5 = unname(tools::md5sum(paths)),
                    n_features = vapply(tables, nrow, integer(1)),
                    library_count_sum = vapply(tables, function(x) sum(x$count), numeric(1)),
                    n_missing_from_union = vapply(tables, function(x)
                      length(setdiff(all_ids, x$gene_id)), integer(1)))
write.csv(audit, file.path(output_dir, "sample_input_audit.csv"), row.names = FALSE)

if (length(shared_ids) != length(all_ids)) {
  stop("Gene inventories differ among samples; audit saved. Do not silently fill missing genes with zero.")
}
counts <- vapply(tables, function(x) x$count[match(all_ids, x$gene_id)],
                 integer(length(all_ids)))
dimnames(counts) <- list(all_ids, primary$gsm)
write.csv(counts, file.path(output_dir, "raw_count_matrix.csv"))
writeLines(c("GSE246156 NC primary: L5 DRG, day 7, Compression versus Sham (3 versus 3).",
             "GSM order in count matrix: three Sham, then three Compression.",
             "Data source: per-sample GEO featureCounts supplementary files.",
             "No DESeq2, GSEA, GSVA or cross-model analysis has been run.",
             "Animal IDs and independence across DRG levels remain unresolved."),
           file.path(output_dir, "selection_summary.txt"))
message("Input preparation complete: ", nrow(counts), " genes x ", ncol(counts), " samples.")
print(primary, row.names = FALSE)
