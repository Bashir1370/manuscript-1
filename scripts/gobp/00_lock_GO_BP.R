#!/usr/bin/env Rscript
source("scripts/gobp/helpers.R", local = TRUE)
run_gobp_lock <- function() {
  manifest_file <- file.path(gobp_lock_dir, "manifest.csv")
  if (file.exists(manifest_file)) {
    sets <- gobp_load_sets()
    message("Reusing checksum-validated GO:BP lock: ", length(sets), " sets")
    return(invisible(sets))
  }
  if (dir.exists(gobp_lock_dir) && length(list.files(gobp_lock_dir, all.files = FALSE))) stop("Incomplete lock exists; inspect it before rebuilding.")
  gobp_require("msigdbr")
  if (!all(c("db_species", "collection", "subcollection") %in% names(formals(msigdbr::msigdbr)))) stop("Update msigdbr; the current collection API is required.")
  message("Downloading Human MSigDB C5:GO:BP, mapped to rat orthologs...")
  audit <- as.data.frame(msigdbr::msigdbr(db_species = "HS", species = "Rattus norvegicus", collection = "C5", subcollection = "GO:BP"))
  needed <- c("gs_name", "gene_symbol", "gs_collection", "gs_subcollection", "db_version")
  if (!all(needed %in% names(audit)) || !nrow(audit) || any(audit$gs_collection != "C5") ||
      any(audit$gs_subcollection != "GO:BP") || any(!startsWith(audit$gs_name, "GOBP_"))) stop("Unexpected MSigDB GO:BP table.")
  version <- unique(as.character(audit$db_version))
  if (length(version) != 1L || anyNA(version) || !nzchar(version)) stop("Missing or ambiguous MSigDB release.")
  valid <- !is.na(audit$gene_symbol) & nzchar(audit$gene_symbol)
  members <- unique(audit[valid, c("gs_name", "gene_symbol")])
  members <- members[order(members$gs_name, members$gene_symbol, method = "radix"), , drop = FALSE]
  if (!nrow(members)) stop("No mapped rat membership.")
  columns <- intersect(c("gs_name", "gs_id", "gs_description", "gs_exact_source", "gs_url", "db_version"), names(audit))
  meta <- unique(audit[, columns, drop = FALSE])
  if (anyDuplicated(meta$gs_name)) stop("Inconsistent pathway metadata.")
  meta$mapped_rat_genes <- as.integer(table(members$gs_name)[meta$gs_name])
  meta$mapped_rat_genes[is.na(meta$mapped_rat_genes)] <- 0L
  dir.create(gobp_lock_dir, recursive = TRUE, showWarnings = FALSE)
  f <- gobp_lock_files()
  con <- gzfile(f[1L], "wt"); tryCatch(write.csv(members, con, row.names = FALSE), finally = close(con))
  con <- gzfile(f[2L], "wt"); tryCatch(write.csv(audit, con, row.names = FALSE, na = ""), finally = close(con))
  gobp_write(meta, f[3L])
  manifest <- data.frame(file = basename(f), md5 = unname(tools::md5sum(f)),
    db_species = "HS", target_species = "Rattus norvegicus", collection = "C5", subcollection = "GO:BP",
    db_version = version, msigdbr_version = as.character(utils::packageVersion("msigdbr")),
    created_UTC = format(Sys.time(), tz = "UTC", usetz = TRUE),
    min_size = gobp_min_size, max_size = gobp_max_size)
  gobp_write(manifest, manifest_file)
  writeLines(capture.output(sessionInfo()), file.path(gobp_lock_dir, "R_sessionInfo.txt"))
  message("Locked GO:BP: ", nrow(meta), " pathways; MSigDB ", version)
  invisible(gobp_load_sets())
}
if (!exists("GOBP_AUTORUN", inherits = FALSE) || isTRUE(GOBP_AUTORUN)) run_gobp_lock()
