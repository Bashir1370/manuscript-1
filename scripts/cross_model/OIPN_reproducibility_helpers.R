# Freeze the currently acquired membership for future reruns. Historical releases
# were not archived: this lock is NOT proof of the historical membership/version.
oipn_locked_hallmark <- function() {
  lock <- "data/gene_sets/Hallmark_rat_locked.csv"
  meta <- "data/gene_sets/Hallmark_rat_lock_metadata.csv"
  if (!file.exists(lock)) {
    if (!requireNamespace("msigdbr", quietly = TRUE)) stop("Install msigdbr to acquire the initial lock.")
    f <- names(formals(msigdbr::msigdbr))
    args <- list(species = "Rattus norvegicus")
    args[[if ("collection" %in% f) "collection" else "category"]] <- "H"
    x <- as.data.frame(do.call(msigdbr::msigdbr, args))
    if (!all(c("gs_name", "gene_symbol") %in% names(x))) stop("Unexpected msigdbr schema.")
    dir.create(dirname(lock), recursive = TRUE, showWarnings = FALSE)
    write.csv(x, lock, row.names = FALSE, na = "")
    version <- if ("db_version" %in% names(x)) paste(unique(x$db_version), collapse = ";") else "not_exposed"
    write.csv(data.frame(acquired_UTC = format(Sys.time(), tz = "UTC", usetz = TRUE),
      msigdbr_version = as.character(utils::packageVersion("msigdbr")), db_version = version,
      species = "Rattus norvegicus", collection = "H", md5 = unname(tools::md5sum(lock)),
      historical_release_recovered = FALSE), meta, row.names = FALSE)
  }
  if (!file.exists(meta)) stop("Missing gene-set lock metadata; do not silently recreate it.")
  m <- read.csv(meta, stringsAsFactors = FALSE)
  if (nrow(m) != 1L || !identical(m$md5, unname(tools::md5sum(lock)))) stop("Gene-set lock checksum mismatch.")
  x <- read.csv(lock, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("gs_name", "gene_symbol") %in% names(x)) || anyNA(x$gs_name) || anyNA(x$gene_symbol) ||
      any(!nzchar(x$gs_name)) || any(!nzchar(x$gene_symbol))) stop("Invalid locked gene sets.")
  sets <- lapply(split(x$gene_symbol, x$gs_name), unique)
  if (length(sets) != 50L || any(!startsWith(names(sets), "HALLMARK_"))) stop("Expected exactly 50 Hallmarks.")
  sets
}

oipn_rerun_paths <- function() {
  root <- "results/GSE160543_Oxaliplatin_vs_Vehicle"
  primary <- Sys.getenv("OIPN_PRIMARY_OUTDIR", root)
  out <- Sys.getenv("OIPN_PATHWAY_OUTDIR", file.path(root, "reproducibility_rerun", "Pathway_analysis"))
  if (normalizePath(out, mustWork = FALSE) == normalizePath(file.path(root, "Pathway_analysis"), mustWork = FALSE) ||
      normalizePath(out, mustWork = FALSE) == normalizePath(file.path(root, "pathway_analysis"), mustWork = FALSE))
    stop("Use a separate OIPN_PATHWAY_OUTDIR; canonical historical outputs must be preserved.")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  list(root = root, primary = primary, out = out)
}

oipn_write_provenance <- function(out, inputs, parameters, stem) {
  inputs <- unique(c(inputs, "data/gene_sets/Hallmark_rat_locked.csv", "data/gene_sets/Hallmark_rat_lock_metadata.csv"))
  write.csv(data.frame(input = inputs, md5 = unname(tools::md5sum(inputs))),
    file.path(out, paste0(stem, "_input_checksums.csv")), row.names = FALSE)
  write.csv(parameters, file.path(out, paste0(stem, "_parameters.csv")), row.names = FALSE)
  packages <- sort(loadedNamespaces())
  write.csv(data.frame(package = packages,
    version = vapply(packages, function(p) as.character(utils::packageVersion(p)), character(1)),
    R_version = as.character(getRversion())),
    file.path(out, paste0(stem, "_package_versions.csv")), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(out, paste0(stem, "_sessionInfo.txt")))
}

oipn_compare_gsea <- function(current, baseline, out) {
  if (!file.exists(baseline)) stop("Missing archived GSEA baseline: ", baseline)
  old <- read.csv(baseline, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("pathway", "NES", "padj", "leadingEdge")
  if (!all(required %in% names(old)) || anyDuplicated(old$pathway) || anyDuplicated(current$pathway))
    stop("Invalid GSEA comparison keys/schema.")
  keys <- sort(union(old$pathway, current$pathway))
  a <- match(keys, old$pathway); b <- match(keys, current$pathway)
  overlap <- vapply(seq_along(keys), function(i) {
    if (is.na(a[i]) || is.na(b[i])) return(NA_real_)
    x <- strsplit(old$leadingEdge[a[i]], ";", fixed = TRUE)[[1L]]
    y <- strsplit(current$leadingEdge[b[i]], ";", fixed = TRUE)[[1L]]
    length(intersect(x, y)) / length(union(x, y))
  }, numeric(1))
  d <- data.frame(pathway = keys, archived_present = !is.na(a), rerun_present = !is.na(b),
    archived_NES = old$NES[a], rerun_NES = current$NES[b],
    archived_FDR = old$padj[a], rerun_FDR = current$padj[b], leading_edge_Jaccard = overlap)
  d$delta_NES <- d$rerun_NES - d$archived_NES
  d$direction_agrees <- sign(d$archived_NES) == sign(d$rerun_NES)
  d$FDR_0_05_class_agrees <- (d$archived_FDR < .05) == (d$rerun_FDR < .05)
  write.csv(d, file.path(out, "GSEA_archived_comparison.csv"), row.names = FALSE, na = "")
  d
}

oipn_compare_gsva <- function(current, baseline, out) {
  if (!file.exists(baseline)) stop("Missing archived GSVA baseline: ", baseline)
  old <- read.csv(baseline, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"pathway" %in% names(old) || anyDuplicated(old$pathway) ||
      !setequal(setdiff(names(old), "pathway"), colnames(current))) stop("GSVA baseline sample/key mismatch.")
  if (!setequal(old$pathway, rownames(current))) stop("GSVA pathway identities changed.")
  rows <- lapply(sort(rownames(current)), function(p) {
    previous <- as.numeric(unlist(old[match(p, old$pathway), colnames(current)], use.names = FALSE))
    data.frame(pathway = p, sample = colnames(current), archived_score = previous,
      rerun_score = as.numeric(current[p, ]), delta = as.numeric(current[p, ]) - previous)
  })
  d <- do.call(rbind, rows)
  write.csv(d, file.path(out, "GSVA_archived_comparison.csv"), row.names = FALSE)
  d
}
