# Shared implementation for the isolated three-study GO:BP experiment.
# Run entry scripts from repository root. No historical files are overwritten.
gobp_root <- "results/GO_BP_three_dataset"
gobp_lock_dir <- "data/gene_sets/GO_BP_rat_locked"
gobp_min_size <- 15L
gobp_max_size <- 500L
gobp_fdr <- 0.05
gobp_studies <- c("OIPN_GSE160543", "NC_GSE246156", "CCI_GSE212311")

gobp_require <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install missing R packages: ", paste(missing, collapse = ", "))
}
gobp_read <- function(path, columns = character()) {
  if (!file.exists(path)) stop("Missing input: ", path, "; run from repository root.")
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("NA", ""))
  if (!all(columns %in% names(x))) stop("Missing columns in ", path, ": ", paste(setdiff(columns, names(x)), collapse = ", "))
  x
}
gobp_write <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE, na = "")
}
gobp_audit <- function(inputs, out) {
  inputs <- unique(c(inputs, "scripts/gobp/helpers.R"))
  if (any(!file.exists(inputs))) stop("Cannot checksum missing inputs.")
  gobp_write(data.frame(input = inputs, md5 = unname(tools::md5sum(inputs))), file.path(out, "input_checksums.csv"))
  writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
}
gobp_spec <- function(study) {
  roots <- c(OIPN_GSE160543 = "results/GSE160543_Oxaliplatin_vs_Vehicle",
             NC_GSE246156 = "results/GSE246156_NC_L5_day7", CCI_GSE212311 = "results/GSE212311_CCI_L4L6_day11")
  root <- roots[[study]]
  if (is.null(root)) stop("Unknown study: ", study)
  if (study == gobp_studies[1L]) return(list(root = root,
    samples = c(paste0("GSM487500", 3:6), paste0("GSM487501", 1:4)),
    groups = rep(c("Vehicle", "Oxaliplatin"), each = 4L),
    meta = "data/metadata/gse160543_analysis_selection.csv", group_col = "group",
    ranks = file.path(root, "Pathway_analysis/ranked_statistics.csv"),
    de = file.path(root, "DE_all_genes.csv"), mapping = NULL,
    counts = file.path(root, "normalized_counts_complete_case.csv"), seed = 160543L))
  if (study == gobp_studies[2L]) return(list(root = root,
    samples = c(paste0("GSM786379", 4:6), paste0("GSM786377", 0:2)),
    groups = rep(c("Sham", "Compression"), each = 3L),
    meta = file.path(root, "locked_primary_samples.csv"), group_col = "condition",
    ranks = file.path(root, "pathway_analysis/ranked_statistics.csv"),
    de = file.path(root, "DESeq2/GSE246156_Compression_vs_Sham_DE_all_genes.csv"),
    mapping = file.path(root, "pathway_analysis/Ensembl_annotation_all_mappings.csv"),
    counts = file.path(root, "raw_count_matrix.csv"), seed = 246156L))
  list(root = root, samples = c(paste0("GSM652375", 1:3), paste0("GSM652374", 8:9), "GSM6523750"),
    groups = rep(c("Sham", "CCI"), each = 3L),
    meta = file.path(root, "locked_primary_samples.csv"), group_col = "condition",
    ranks = file.path(root, "pathway_analysis_source_aware/ranked_statistics.csv"),
    de = file.path(root, "pathway_analysis_source_aware/feature_to_symbol_audit.csv"), mapping = NULL,
    counts = file.path(root, "DESeq2/normalized_counts_filtered.csv"), seed = 212311L)
}
gobp_metadata <- function(study) {
  s <- gobp_spec(study); m <- gobp_read(s$meta, c("gsm", s$group_col))
  if (anyDuplicated(m$gsm) || !all(s$samples %in% m$gsm)) stop("Invalid sample metadata: ", study)
  m <- m[match(s$samples, m$gsm), , drop = FALSE]
  if (!identical(as.character(m[[s$group_col]]), s$groups)) stop("Locked groups differ: ", study)
  if (study == gobp_studies[2L] && (!all(m$day == 7L) || !all(m$drg_level == "L5"))) stop("NC stratum differs.")
  if (study == gobp_studies[3L] && (!all(m$day == 11L) || !all(m$tissue == "ipsilateral_L4-L6_DRG"))) stop("CCI stratum differs.")
  data.frame(sample = s$samples, source_group = s$groups,
    condition = rep(c("Control", "Neuropathy"), each = length(s$samples) / 2L))
}
gobp_gene_table <- function(study) {
  s <- gobp_spec(study)
  d <- gobp_read(s$de)
  if (study == gobp_studies[1L]) {
    needed <- c("symbol", "gene_id", "stat", "log2FoldChange", "padj")
    if (!all(needed %in% names(d))) stop("Invalid OIPN DE schema.")
    d <- d[!is.na(d$symbol) & nzchar(d$symbol) & is.finite(d$stat), , drop = FALSE]
    d <- d[order(-d$stat, method = "radix"), , drop = FALSE]
    d <- d[!duplicated(d$symbol), , drop = FALSE]
    g <- data.frame(symbol = d$symbol, feature_id = as.character(d$gene_id), statistic = d$stat,
      log2FC = d$log2FoldChange, gene_FDR = d$padj, mapping_source = "original_DE_symbol")
  } else if (study == gobp_studies[2L]) {
    if (!all(c("gene", "stat", "log2FoldChange", "padj") %in% names(d))) stop("Invalid NC DE schema.")
    m <- gobp_read(s$mapping, c("ENSEMBL", "SYMBOL"))
    d$ensembl <- sub("\\.[0-9]+$", "", d$gene)
    pairs <- unique(m[!is.na(m$SYMBOL) & nzchar(m$SYMBOL), c("ENSEMBL", "SYMBOL")])
    ambiguous <- names(which(table(pairs$ENSEMBL) > 1L))
    pairs <- pairs[!pairs$ENSEMBL %in% ambiguous, , drop = FALSE]
    d <- merge(d, pairs, by.x = "ensembl", by.y = "ENSEMBL", sort = FALSE)
    d <- d[is.finite(d$stat), , drop = FALSE]
    d <- d[order(-abs(d$stat), d$ensembl), , drop = FALSE]
    d <- d[!duplicated(d$SYMBOL), , drop = FALSE]
    g <- data.frame(symbol = d$SYMBOL, feature_id = d$gene, statistic = d$stat,
      log2FC = d$log2FoldChange, gene_FDR = d$padj, mapping_source = "original_unambiguous_Ensembl_symbol")
  } else {
    if (!all(c("symbol", "gene_id", "stat", "log2FoldChange", "padj", "source") %in% names(d))) stop("Invalid CCI audit schema.")
    g <- data.frame(symbol = d$symbol, feature_id = as.character(d$gene_id), statistic = d$stat,
      log2FC = d$log2FoldChange, gene_FDR = d$padj, mapping_source = d$source)
  }
  if (!nrow(g) || anyNA(g[c("symbol", "feature_id")]) || anyDuplicated(g$symbol) ||
      anyDuplicated(g$feature_id) || any(!is.finite(g$statistic)) || any(!is.finite(g$log2FC)) ||
      any(!is.na(g$gene_FDR) & (!is.finite(g$gene_FDR) | g$gene_FDR < 0 | g$gene_FDR > 1))) stop("Invalid gene evidence: ", study)
  ranks <- gobp_read(s$ranks, c("symbol", "statistic"))
  if (anyDuplicated(ranks$symbol) || !setequal(ranks$symbol, g$symbol) ||
      !isTRUE(all.equal(ranks$statistic, g$statistic[match(ranks$symbol, g$symbol)], tolerance = 1e-8))) stop("Gene representatives differ from archived ranks: ", study)
  g[match(ranks$symbol, g$symbol), , drop = FALSE]
}
gobp_lock_files <- function() file.path(gobp_lock_dir, c("membership.csv.gz", "mapping_audit.csv.gz", "pathway_metadata.csv"))
gobp_load_sets <- function() {
  f <- gobp_lock_files(); manifest <- gobp_read(file.path(gobp_lock_dir, "manifest.csv"),
    c("file", "md5", "db_version", "db_species", "target_species", "collection", "subcollection", "min_size", "max_size"))
  if (anyNA(manifest) || anyDuplicated(manifest$file) || length(unique(manifest$db_version)) != 1L ||
      any(manifest$db_species != "HS") || any(manifest$target_species != "Rattus norvegicus") ||
      any(manifest$collection != "C5") || any(manifest$subcollection != "GO:BP") ||
      any(manifest$min_size != gobp_min_size) || any(manifest$max_size != gobp_max_size)) stop("Lock release/species/parameters differ.")
  if (!setequal(manifest$file, basename(f)) || any(!file.exists(f)) ||
      !identical(unname(tools::md5sum(f)), manifest$md5[match(basename(f), manifest$file)])) stop("GO:BP lock checksum mismatch; do not silently rebuild.")
  x <- gobp_read(f[1L], c("gs_name", "gene_symbol"))
  if (!nrow(x) || anyNA(x[c("gs_name", "gene_symbol")]) || anyDuplicated(x[c("gs_name", "gene_symbol")]) ||
      any(!startsWith(x$gs_name, "GOBP_"))) stop("Invalid locked GO:BP membership.")
  sets <- lapply(split(x$gene_symbol, x$gs_name), function(z) sort(unique(z), method = "radix"))
  meta <- gobp_read(f[3L], "gs_name")
  if (anyDuplicated(meta$gs_name) || !all(names(sets) %in% meta$gs_name)) stop("Lock pathway identities differ.")
  ids <- sort(meta$gs_name, method = "radix")
  setNames(lapply(ids, function(id) if (id %in% names(sets)) sets[[id]] else character()), ids)
}
gobp_expression <- function(study) {
  s <- gobp_spec(study); meta <- gobp_metadata(study)
  tab <- gobp_read(s$counts)
  ids <- as.character(tab[[1L]])
  if (anyNA(ids) || anyDuplicated(ids) || !all(s$samples %in% names(tab))) stop("Invalid count feature/sample keys: ", study)
  counts <- as.matrix(tab[, s$samples, drop = FALSE])
  if (!is.numeric(counts) || any(!is.finite(counts)) || any(counts < 0)) stop("Invalid counts: ", study)
  rownames(counts) <- ids
  if (study == gobp_studies[2L]) {
    gobp_require("DESeq2")
    if (any(counts != round(counts)) || any(counts > .Machine$integer.max)) stop("NC raw counts must be integers.")
    counts <- counts[rowSums(counts >= 10) >= 3L, , drop = FALSE]
    storage.mode(counts) <- "integer"
    md <- data.frame(condition = factor(meta$condition, levels = c("Control", "Neuropathy")), row.names = meta$sample)
    dds <- DESeq2::DESeqDataSetFromMatrix(counts, md, ~ condition)
    dds <- DESeq2::estimateSizeFactors(dds)
    counts <- DESeq2::counts(dds, normalized = TRUE)
  }
  g <- gobp_gene_table(study)
  if (!all(g$feature_id %in% rownames(counts))) stop("Rank representatives missing from expression matrix: ", study)
  counts <- counts[match(g$feature_id, rownames(counts)), s$samples, drop = FALSE]
  rownames(counts) <- g$symbol
  list(counts = counts, expr = log2(counts + 1), genes = g, meta = meta)
}
gobp_classify <- function(nes, fdr) {
  stopifnot(is.matrix(nes), is.matrix(fdr), identical(dim(nes), dim(fdr)), ncol(nes) == 3L)
  valid <- rowSums(is.finite(nes) & is.finite(fdr)) == 3L
  pos_sig <- is.finite(nes) & is.finite(fdr) & nes > 0 & fdr < gobp_fdr
  neg_sig <- is.finite(nes) & is.finite(fdr) & nes < 0 & fdr < gobp_fdr
  data.frame(available_all3 = valid,
    shared_positive = valid & rowSums(nes > 0, na.rm = TRUE) == 3L & rowSums(pos_sig) == 3L,
    shared_negative = valid & rowSums(nes < 0, na.rm = TRUE) == 3L & rowSums(neg_sig) == 3L,
    oipn_direction = valid & !is.na(nes[, 1L]) & nes[, 1L] > 0 & rowSums(nes[, 2:3, drop = FALSE] < 0, na.rm = TRUE) == 2L,
    physical_direction = valid & !is.na(nes[, 1L]) & nes[, 1L] < 0 & rowSums(nes[, 2:3, drop = FALSE] > 0, na.rm = TRUE) == 2L,
    oipn_support = valid & pos_sig[, 1L] & rowSums(pos_sig[, 2:3, drop = FALSE]) == 0L,
    physical_support = valid & !pos_sig[, 1L] & rowSums(pos_sig[, 2:3, drop = FALSE]) == 2L,
    n_positive_significant = rowSums(pos_sig), n_negative_significant = rowSums(neg_sig))
}
gobp_pretty <- function(x) gsub("_", " ", sub("^GOBP_", "", x))
gobp_plot_pages <- function(long, ids, out, stem, title, value, limit, x = "study", page_size = 40L, facet_studies = FALSE) {
  # Recreate this plot family's files to prevent stale pages after selection changes.
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  old <- list.files(out, full.names = TRUE)
  old <- old[startsWith(basename(old), paste0(stem, "_page_")) | basename(old) == paste0(stem, "_empty.txt")]
  unlink(old)
  if (!length(ids)) {
    writeLines("No rows met the declared criteria.", file.path(out, paste0(stem, "_empty.txt")))
    return(invisible(NULL))
  }
  gobp_require("ggplot2")
  if (!is.finite(limit) || limit <= 0) stop("Plot limit must be finite and positive.")
  pages <- split(ids, ceiling(seq_along(ids) / page_size))
  for (j in seq_along(pages)) {
    d <- long[long$pathway %in% pages[[j]], , drop = FALSE]
    d$pathway <- factor(d$pathway, levels = rev(pages[[j]]))
    d$study <- factor(d$study, levels = gobp_studies)
    if (x == "sample_label") d$sample_label <- factor(d$sample_label, levels = c(paste0("Ctrl", 1:4), paste0("Inj", 1:4)))
    d$display_value <- pmax(-limit, pmin(limit, d[[value]]))
    if (!"mark" %in% names(d)) d$mark <- ""
    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[x]], y = pathway, fill = display_value)) +
      ggplot2::geom_tile(color = "white", linewidth = .25) + ggplot2::geom_text(ggplot2::aes(label = mark), size = 2.5) +
      ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
        limits = c(-limit, limit), na.value = "#BDBDBD", name = value) +
      ggplot2::scale_y_discrete(labels = function(z) vapply(gobp_pretty(z), function(a) paste(strwrap(a, width = 60), collapse = "\n"), character(1))) +
      ggplot2::labs(title = title, subtitle = "Independent DRG cohorts; * within-study FDR <0.05; gray = unavailable", x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = 9) + ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
    if (facet_studies) p <- p + ggplot2::facet_grid(. ~ study, scales = "free_x", space = "free_x")
    fn <- file.path(out, sprintf("%s_page_%03d", stem, j)); ht <- max(4, 2 + .42 * length(pages[[j]]))
    ggplot2::ggsave(paste0(fn, ".pdf"), p, width = 14, height = ht, limitsize = FALSE)
    ggplot2::ggsave(paste0(fn, ".png"), p, width = 14, height = ht, dpi = 200, bg = "white", limitsize = FALSE)
  }
  invisible(NULL)
}
