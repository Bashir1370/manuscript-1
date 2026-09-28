#!/usr/bin/env Rscript
# Run at repository root: Rscript scripts/GSE246156/05_PCA_all_samples_GSE246156.R
# QC only. Do not infer animal independence or change the locked primary contrast from PCA.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})
options(timeout = max(600L, getOption("timeout")))

manifest_file <- file.path("data", "metadata", "geo_sample_manifest.csv")
if (!file.exists(manifest_file)) stop("Run from the repository root.")
meta <- read.csv(manifest_file, stringsAsFactors = FALSE, check.names = FALSE)
need <- c("accession", "gsm", "title", "condition", "tissue", "drg_level", "day")
if (!all(need %in% names(meta))) stop("Incomplete GEO manifest.")
meta <- meta[meta$accession == "GSE246156", , drop = FALSE]
meta <- meta[order(meta$gsm), , drop = FALSE]
if (nrow(meta) != 48L || anyDuplicated(meta$gsm) ||
    sum(meta$tissue == "DRG") != 36L ||
    sum(meta$tissue == "Sciatic nerve") != 12L ||
    anyNA(meta[, c("gsm", "title", "condition", "tissue", "day")])) {
  stop("GSE246156 metadata does not match the audited 48/36 design.")
}
meta$day <- as.integer(meta$day)
drg <- meta[meta$tissue == "DRG", , drop = FALSE]
if (anyNA(meta$day) || any(!meta$day %in% c(3L, 7L)) ||
    any(!meta$condition %in% c("Sham", "Compression")) ||
    any(!drg$drg_level %in% c("L4", "L5", "L6"))) {
  stop("Unexpected day, condition or DRG level.")
}
cells <- table(drg$condition, drg$drg_level, drg$day)
if (any(dim(cells) != c(2L, 3L, 2L)) || any(cells != 3L)) {
  stop("Unexpected DRG condition x level x day sample counts.")
}

cache <- file.path("data", "source-cache", "GSE246156")
out <- file.path("results", "GSE246156_all_samples_PCA")
dir.create(cache, recursive = TRUE, showWarnings = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
matrix_url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/",
                     "GSE246156/matrix/GSE246156_series_matrix.txt.gz")
matrix_file <- file.path(cache, "GSE246156_series_matrix.txt.gz")

download_checked <- function(url, path) {
  if (file.exists(path) && file.info(path)$size > 100L) return(invisible(path))
  if (file.exists(path)) unlink(path)
  tmp <- paste0(path, ".partial")
  curl_exe <- Sys.which(if (.Platform$OS.type == "windows") "curl.exe" else "curl")
  last_error <- "unknown transfer error"
  for (attempt in seq_len(8L)) {
    message("Downloading ", basename(path), " (attempt ", attempt, "/8)")
    if (nzchar(curl_exe)) {
      # -C - resumes a previously interrupted .partial download when GEO supports ranges.
      # curl exits nonzero if the server closes the connection before the full response.
      args <- c("--fail", "--location", "--silent", "--show-error",
                "--connect-timeout", "30", "--max-time", "600",
                "--continue-at", "-", "--output", shQuote(tmp), shQuote(url))
      output <- suppressWarnings(system2(curl_exe, args = args,
                                         stdout = TRUE, stderr = TRUE))
      status <- attr(output, "status")
      if (is.null(status)) status <- 0L
      ok <- identical(as.integer(status), 0L)
      if (!ok) {
        last_error <- paste(tail(output, 3L), collapse = " ")
        # Exit 33 means the remote server cannot resume; retry from byte zero.
        if (identical(as.integer(status), 33L) && file.exists(tmp)) unlink(tmp)
      }
    } else {
      # Base R fallback: it cannot resume, so retry complete downloads.
      if (file.exists(tmp)) unlink(tmp)
      result <- tryCatch(utils::download.file(url, tmp, mode = "wb",
                                               method = "libcurl", quiet = TRUE),
                         error = function(e) e)
      ok <- is.numeric(result) && length(result) == 1L && result == 0L
      if (!ok) last_error <- if (inherits(result, "error"))
        conditionMessage(result) else paste("download.file returned", result)
    }
    if (ok && file.exists(tmp) && file.info(tmp)$size > 100L) {
      if (!file.rename(tmp, path)) stop("Cannot move downloaded file to ", path)
      return(invisible(path))
    }
    message("Transfer interrupted: ", last_error)
    if (attempt < 8L) Sys.sleep(min(2L * attempt, 10L))
  }
  stop("Download failed after 8 attempts: ", url, "; ", last_error,
       ". Partial data remains at ", tmp)
}
download_checked(matrix_url, matrix_file)
geo_lines <- readLines(gzfile(matrix_file), warn = FALSE)
read_geo_field <- function(key) {
  line <- geo_lines[startsWith(geo_lines, paste0(key, "\t"))]
  if (length(line) != 1L) stop("Missing/duplicated GEO field: ", key)
  values <- as.character(read.table(text = line, sep = "\t", quote = "\"",
                                    comment.char = "", header = FALSE,
                                    check.names = FALSE, stringsAsFactors = FALSE)[1L, -1L])
  if (length(values) != 48L) stop("GEO field has unexpected length: ", key)
  values
}
geo_gsm <- read_geo_field("!Sample_geo_accession")
geo_title <- read_geo_field("!Sample_title")
geo_url <- read_geo_field("!Sample_supplementary_file_1")
if (anyDuplicated(geo_gsm) || !setequal(geo_gsm, meta$gsm)) {
  stop("GEO GSMs differ from the tracked manifest.")
}
ix <- match(meta$gsm, geo_gsm)
if (!identical(meta$title, geo_title[ix])) {
  stop("GEO titles differ from the tracked sample metadata.")
}
urls <- sub("^ftp://", "https://", geo_url[ix])
if (any(!grepl("^https://ftp\\.ncbi\\.nlm\\.nih\\.gov/geo/samples/", urls)) ||
    any(!grepl("_featurecounts_.*\\.txt\\.gz$", urls)) ||
    any(!startsWith(basename(urls), paste0(meta$gsm, "_")))) {
  stop("Unexpected supplementary file URL or GSM-to-file mapping.")
}
files <- file.path(cache, basename(urls))
for (i in seq_along(files)) download_checked(urls[i], files[i])

read_counts <- function(path) {
  x <- tryCatch(read.delim(gzfile(path), sep = "\t", header = TRUE,
                           comment.char = "#", quote = "", check.names = FALSE),
                error = function(e) stop("Cannot read ", path, ": ", conditionMessage(e)))
  if (nrow(x) < 1000L || ncol(x) < 2L) stop("Bad featureCounts layout: ", path)
  id <- as.character(x[[1L]])
  value <- suppressWarnings(as.numeric(as.character(x[[ncol(x)]])))
  if (anyNA(id) || any(!nzchar(id)) || anyDuplicated(id) ||
      anyNA(value) || any(!is.finite(value)) || any(value < 0) ||
      any(value != round(value)) || any(value > .Machine$integer.max)) {
    stop("Invalid gene IDs or non-integer counts: ", path)
  }
  list(id = id, value = as.integer(value))
}
tables <- lapply(files, read_counts)
names(tables) <- meta$gsm
reference <- sort(tables[[1L]]$id)
missing <- vapply(tables, function(x) length(setdiff(reference, x$id)), integer(1))
extra <- vapply(tables, function(x) length(setdiff(x$id, reference)), integer(1))
audit <- data.frame(gsm = meta$gsm, condition = meta$condition,
                    tissue = meta$tissue, drg_level = meta$drg_level,
                    day = meta$day, source_url = urls, file = files,
                    md5 = unname(tools::md5sum(files)),
                    n_features = vapply(tables, function(x) length(x$id), integer(1)),
                    library_count_sum = vapply(tables, function(x) sum(x$value), numeric(1)),
                    missing_genes = missing, extra_genes = extra)
write.csv(audit, file.path(out, "sample_input_audit.csv"), row.names = FALSE)
if (any(missing != 0L) || any(extra != 0L)) {
  stop("Gene sets differ across files. See sample_input_audit.csv; no zero imputation.")
}
counts <- vapply(tables, function(x) x$value[match(reference, x$id)],
                 integer(length(reference)))
dimnames(counts) <- list(reference, meta$gsm)
if (anyNA(counts)) stop("Missing aligned count value.")

run_pca <- function(selected, label) {
  m <- counts[, selected$gsm, drop = FALSE]
  keep <- rowSums(m >= 10L) >= 3L
  if (sum(keep) < 1000L) stop("Too few genes after filtering in ", label)
  m <- m[keep, , drop = FALSE]
  sample_data <- selected
  rownames(sample_data) <- sample_data$gsm
  dds <- DESeqDataSetFromMatrix(countData = m, colData = sample_data, design = ~ 1)
  transformed <- vst(dds, blind = TRUE)
  p <- prcomp(t(assay(transformed)), center = TRUE, scale. = FALSE)
  variance <- 100 * p$sdev^2 / sum(p$sdev^2)
  coord <- cbind(selected, as.data.frame(p$x[, seq_len(min(5L, ncol(p$x))), drop = FALSE]))
  write.csv(coord, file.path(out, paste0("PCA_", label, "_coordinates.csv")),
            row.names = FALSE)
  write.csv(data.frame(PC = paste0("PC", seq_along(variance)),
                       variance_percent = variance),
            file.path(out, paste0("PCA_", label, "_variance.csv")), row.names = FALSE)
  list(data = coord, xlab = sprintf("PC1 (%.1f%%)", variance[1L]),
       ylab = sprintf("PC2 (%.1f%%)", variance[2L]), genes = nrow(m))
}
plot_pca <- function(obj, label, title, map, facet = NULL, point_size = 3) {
  p <- ggplot(obj$data, map) + geom_point(size = point_size, alpha = 0.85) +
       labs(title = title, x = obj$xlab, y = obj$ylab) + theme_bw(base_size = 12) +
       theme(plot.title = element_text(face = "bold"))
  if (!is.null(facet)) p <- p + facet
  ggsave(file.path(out, paste0("PCA_", label, ".png")), plot = p,
         width = 10, height = 7, dpi = 300)
}

all48 <- run_pca(meta, "all48")
plot_pca(all48, "all48", "GSE246156: all 48 samples (tissue overview)",
         aes(PC1, PC2, colour = tissue, shape = factor(day)))
all36 <- run_pca(drg, "DRG36")
plot_pca(all36, "DRG36", "GSE246156: 36 DRG samples",
         aes(PC1, PC2, colour = condition, shape = drg_level),
         facet_grid(. ~ day))
day7 <- drg[drg$day == 7L, , drop = FALSE]
day7_pca <- run_pca(day7, "DRG_day7_18")
plot_pca(day7_pca, "DRG_day7_18", "GSE246156: DRG day 7 (9 Sham, 9 Compression)",
         aes(PC1, PC2, colour = condition, shape = drg_level))
writeLines(c("PCA is descriptive QC; it does not establish animal independence.",
             "Input: GEO per-sample integer featureCounts; source URLs and MD5 in sample_input_audit.csv.",
             "Each PCA uses size-factor-normalized blind VST, design ~ 1, genes with >=10 counts in >=3 samples.",
             paste0("All 48: ", all48$genes, " genes; DRG 36: ", all36$genes,
                    " genes; DRG day 7: ", day7_pca$genes, " genes."),
             "Interpret day and DRG level before any pooling; animal IDs are missing from GEO.",
             "Existing locked L5 day-7 3 vs 3 contrast remains unchanged."),
           file.path(out, "PCA_readme.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
message("Done: ", out)
