#!/usr/bin/env Rscript
source("scripts/gobp/helpers.R", local = TRUE)
run_gobp_gsea <- function() {
  gobp_require("fgsea")
  sets <- gobp_load_sets()
  for (study in gobp_studies) {
    message("GO:BP GSEA: ", study)
    s <- gobp_spec(study); meta <- gobp_metadata(study); g <- gobp_gene_table(study)
    ranks <- setNames(g$statistic, g$symbol)
    # Stable secondary ordering preserves each original representative and statistic.
    ranks <- ranks[order(-ranks, names(ranks), method = "radix")]
    sizes <- vapply(sets, function(z) length(intersect(z, names(ranks))), integer(1))
    eligible <- sizes >= gobp_min_size & sizes <= gobp_max_size
    if (!any(eligible)) stop("No eligible GO:BP sets: ", study)
    warnings <- character()
    set.seed(s$seed)
    fg <- withCallingHandlers(as.data.frame(fgsea::fgseaMultilevel(pathways = sets[eligible], stats = ranks,
      minSize = gobp_min_size, maxSize = gobp_max_size, eps = 1e-10, nproc = 1L,
      nPermSimple = 10000L, sampleSize = 101L, gseaParam = 1, scoreType = "std")),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
    if (!all(c("pathway", "pval", "NES", "ES", "size", "leadingEdge") %in% names(fg)) ||
        anyDuplicated(fg$pathway) || !setequal(fg$pathway, names(sets)[eligible])) stop("Unexpected fgsea output: ", study)
    fg <- fg[match(names(sets)[eligible], fg$pathway), , drop = FALSE]
    # Full tested family, not the displayed or significant subset.
    fg$padj <- p.adjust(fg$pval, method = "BH", n = sum(eligible))
    fg$leadingEdge <- vapply(fg$leadingEdge, function(z) paste(sort(unique(z)), collapse = ";"), character(1))
    fg$study <- study
    fg$test_status <- ifelse(is.finite(fg$pval) & is.finite(fg$NES) & is.finite(fg$padj), "tested", "numerical_failure")
    out <- file.path(gobp_root, "GSEA", study)
    gobp_write(fg, file.path(out, "GSEA_GO_BP_results.csv"))
    coverage <- data.frame(pathway = names(sets), mapped_rat_genes = lengths(sets), rank_overlap = sizes,
      eligible = eligible, status = ifelse(sizes < gobp_min_size, "below_min_size", ifelse(sizes > gobp_max_size, "above_max_size", "eligible")))
    gobp_write(coverage, file.path(out, "pathway_coverage.csv"))
    gobp_write(g, file.path(out, "gene_feature_evidence.csv"))
    gobp_write(data.frame(symbol = names(ranks), statistic = unname(ranks)), file.path(out, "ranked_statistics.csv"))
    gobp_write(meta, file.path(out, "sample_metadata.csv"))
    gobp_write(data.frame(parameter = c("seed", "engine", "minSize", "maxSize", "eps", "nPermSimple", "gseaParam", "tested_sets", "tied_statistics"),
      value = c(s$seed, "fgseaMultilevel", gobp_min_size, gobp_max_size, 1e-10, 10000, 1, sum(eligible), sum(duplicated(ranks)))), file.path(out, "parameters.csv"))
    writeLines(unique(warnings), file.path(out, "warnings.txt"))
    gobp_audit(c(s$ranks, s$de, s$mapping, s$meta, gobp_lock_files(), file.path(gobp_lock_dir, "manifest.csv"), "scripts/gobp/01_GSEA_GO_BP.R"), out)
  }
  invisible(TRUE)
}
if (!exists("GOBP_AUTORUN", inherits = FALSE) || isTRUE(GOBP_AUTORUN)) run_gobp_gsea()
