#!/usr/bin/env Rscript
# One complete entry point; run from the repository root.
# Sys.setenv(GOBP_TABLES_ONLY="true") skips figures but computes all tables.
run_GO_BP_three_dataset <- function() {
  if (!file.exists("scripts/gobp/helpers.R")) stop("Set working directory to the repository root.")
  env <- new.env(parent = globalenv()); env$GOBP_AUTORUN <- FALSE
  for (file in c("00_lock_GO_BP.R", "01_GSEA_GO_BP.R", "02_compare_GO_BP.R", "03_GSVA_GO_BP.R", "04_leading_edge_GO_BP.R", "05_sample_expression_GO_BP.R")) {
    source(file.path("scripts/gobp", file), local = env)
  }
  plots <- tolower(Sys.getenv("GOBP_TABLES_ONLY", "false")) != "true"
  env$gobp_require(c("msigdbr", "fgsea", "GSVA", "limma", "DESeq2", "BiocParallel", if (plots) "ggplot2"))
  # Validate all archived inputs before expensive scoring or creating a lock.
  for (study in env$gobp_studies) {
    env$gobp_metadata(study); env$gobp_gene_table(study)
    if (!file.exists(env$gobp_spec(study)$counts)) stop("Missing counts: ", study)
  }
  env$run_gobp_lock()
  env$run_gobp_gsea()
  env$run_gobp_comparison(draw_plots = plots)
  env$run_gobp_gsva(draw_plots = plots)
  env$run_gobp_leading_edge(draw_plots = plots)
  env$run_gobp_sample_context(draw_plots = plots)
  message("GO:BP three-study experiment complete: results/GO_BP_three_dataset")
  invisible(TRUE)
}
run_GO_BP_three_dataset()
