#!/usr/bin/env Rscript
# Meaningful regression checks on archived inputs and known historical outputs.
# No GO result claims; no plots; temporary fixture outputs are removed on exit.
validate_GO_BP_workflow <- function() {
  env <- new.env(parent = globalenv()); env$GOBP_AUTORUN <- FALSE
  for (f in list.files("scripts/gobp", pattern = "[.]R$", full.names = TRUE)) parse(f)
  for (f in c("02_compare_GO_BP.R", "03_GSVA_GO_BP.R", "04_leading_edge_GO_BP.R", "05_sample_expression_GO_BP.R")) source(file.path("scripts/gobp", f), local = env)
  tmp <- tempfile("gobp-validation-"); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  env$gobp_root <- tmp
  source_files <- c(
    "results/GSE160543_Oxaliplatin_vs_Vehicle/Pathway_analysis/GSEA_Hallmark_results.csv",
    "results/GSE246156_NC_L5_day7/pathway_analysis/GSEA_Hallmark_results.csv",
    "results/GSE212311_CCI_L4L6_day11/pathway_analysis_source_aware/GSEA_Hallmark_results.csv")
  tabs <- lapply(source_files, env$gobp_read); ids <- sort(tabs[[1L]]$pathway)
  env$gobp_load_sets <- function() setNames(rep(list(character()), length(ids)), ids)
  env$gobp_lock_files <- function() character()
  for (i in seq_along(env$gobp_studies)) {
    study <- env$gobp_studies[i]; m <- env$gobp_metadata(study); g <- env$gobp_gene_table(study)
    stopifnot(nrow(m) == c(8L, 6L, 6L)[i], nrow(g) == c(16614L, 12424L, 13799L)[i])
    tabs[[i]]$test_status <- "tested"
    env$gobp_write(tabs[[i]], file.path(tmp, "GSEA", study, "GSEA_GO_BP_results.csv"))
    env$gobp_write(g, file.path(tmp, "GSEA", study, "gene_feature_evidence.csv"))
    s <- env$gobp_spec(study); c <- env$gobp_read(s$counts)
    stopifnot(all(s$samples %in% names(c)), !anyDuplicated(as.character(c[[1L]])), all(g$feature_id %in% as.character(c[[1L]])))
    cat(study, "sample, rank, feature and count alignment PASS\n")
  }
  wide <- env$run_gobp_comparison(draw_plots = FALSE)
  stopifnot(sum(wide$shared_positive) == 10L, sum(wide$shared_negative) == 2L,
    sum(wide$oipn_direction | wide$physical_direction) == 7L)
  env$run_gobp_leading_edge(draw_plots = FALSE)
  expected <- list(shared_positive = c(memberships = 109L, unique = 78L, priority = 20L, strict = 8L),
    shared_negative = c(memberships = 44L, unique = 39L, priority = 7L, strict = 1L),
    divergent = c(opposite_unique = 126L, priority = 12L))
  for (category in names(expected)) {
    base <- file.path(tmp, "leading_edge", category)
    all <- env$gobp_read(file.path(base, "gene_membership_summary.csv"))
    shared <- all[all$shared_all3, , drop = FALSE]
    priority <- env$gobp_read(file.path(base, "priority_OIPN_plus_physical.csv"))
    strict <- env$gobp_read(file.path(base, "strict_significant_all3.csv"))
    values <- c(memberships = nrow(shared), unique = length(unique(shared$symbol)), priority = nrow(priority),
      strict = nrow(strict), opposite_unique = length(unique(all$symbol[all$opposite_both_physical])))
    stopifnot(identical(as.integer(values[names(expected[[category]])]), as.integer(expected[[category]])))
    stopifnot(!anyDuplicated(priority$symbol))
    cat(category, "known historical membership/prioritization counts PASS\n")
  }
  # Missing values and zero NES must not become evidence for sharing/divergence.
  nes <- rbind(c(1, 1, 1), c(-1, -1, -1), c(1, -1, -1), c(1, NA, -1), c(0, -1, -1))
  fdr <- matrix(.01, nrow(nes), 3L)
  z <- env$gobp_classify(nes, fdr)
  stopifnot(identical(z$shared_positive, c(TRUE, FALSE, FALSE, FALSE, FALSE)),
    identical(z$shared_negative, c(FALSE, TRUE, FALSE, FALSE, FALSE)),
    identical(z$oipn_direction, c(FALSE, FALSE, TRUE, FALSE, FALSE)), !anyNA(z))
  fdr[1, 2] <- NA; z <- env$gobp_classify(nes, fdr); stopifnot(!z$shared_positive[1])
  cat("Missing/zero NES and missing FDR handling PASS\n")
  # Exercise a real compressed lock, including a term with no mapped genes.
  lockenv <- new.env(parent = globalenv()); source("scripts/gobp/helpers.R", local = lockenv)
  lockenv$gobp_lock_dir <- file.path(tmp, "synthetic_lock"); dir.create(lockenv$gobp_lock_dir)
  lf <- lockenv$gobp_lock_files()
  con <- gzfile(lf[1L], "wt"); write.csv(data.frame(gs_name = rep("GOBP_TEST_A", 2L), gene_symbol = c("Ahsp", "Bard1")), con, row.names = FALSE); close(con)
  con <- gzfile(lf[2L], "wt"); write.csv(data.frame(test = "mapping_fixture"), con, row.names = FALSE); close(con)
  lockenv$gobp_write(data.frame(gs_name = c("GOBP_TEST_A", "GOBP_TEST_B")), lf[3L])
  mf <- data.frame(file = basename(lf), md5 = unname(tools::md5sum(lf)), db_version = "test",
    db_species = "HS", target_species = "Rattus norvegicus", collection = "C5", subcollection = "GO:BP", min_size = 15L, max_size = 500L)
  lockenv$gobp_write(mf, file.path(lockenv$gobp_lock_dir, "manifest.csv"))
  sets <- lockenv$gobp_load_sets(); stopifnot(identical(lengths(sets), c(GOBP_TEST_A = 2L, GOBP_TEST_B = 0L)))
  cat("Compressed lock roundtrip and zero-mapped term handling PASS\n")
  cat("\n", file = lf[3L], append = TRUE)
  stopifnot(inherits(try(lockenv$gobp_load_sets(), silent = TRUE), "try-error"))
  cat("Corrupt lock correctly rejected PASS\n")
  # All-empty selections are expected and produce valid empty tables, not errors.
  for (flag in c("shared_positive", "shared_negative", "oipn_direction", "physical_direction")) wide[[flag]] <- FALSE
  env$gobp_write(wide, file.path(tmp, "comparison/all_GO_BP_pathways_classified.csv"))
  env$run_gobp_leading_edge(draw_plots = FALSE)
  empty <- env$run_gobp_sample_context(draw_plots = FALSE)
  stopifnot(nrow(empty) == 0L)
  cat("Empty pathway/priority/sample-context handling PASS\n")
  # Changing the collection size must change the BH family, never assume 50.
  if (requireNamespace("limma", quietly = TRUE)) {
    set.seed(13); scores <- matrix(rnorm(137L * 8L), 137L, 8L)
    rownames(scores) <- paste0("GOBP_TEST_", seq_len(137L)); meta <- env$gobp_metadata(env$gobp_studies[1L]); colnames(scores) <- meta$sample
    fit <- env$gobp_fit_gsva(scores, meta)
    stopifnot(nrow(fit) == 137L, all(fit$tested_family_size == 137L),
      isTRUE(all.equal(fit$FDR_BH, p.adjust(fit$p_value, "BH"))),
      isTRUE(all.equal(fit$delta_GSVA, unname(rowMeans(scores[, 5:8]) - rowMeans(scores[, 1:4])))))
    cat("GSVA limma full-family BH and contrast PASS\n")
  } else cat("GSVA limma functional test SKIPPED: limma unavailable\n")
  cat("Validation complete. This is not a full GO:BP GSEA/GSVA execution.\n")
  invisible(TRUE)
}
validate_GO_BP_workflow()
