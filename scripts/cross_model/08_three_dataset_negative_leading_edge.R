#!/usr/bin/env Rscript
# Original negative-NES leading edges; never reverse or recompute the ranks.
local({
  env <- new.env(parent = globalenv())
  env$LE_AUTORUN <- FALSE
  env$GENE_PRIORITY_AUTORUN <- FALSE
  sys.source("scripts/cross_model/05_three_dataset_leading_edge.R", envir = env)
  sys.source("scripts/cross_model/06_three_dataset_gene_prioritization.R", envir = env)
  env$run_shared_leading_edge(direction = "negative",
    draw_plots = tolower(Sys.getenv("LE_TABLES_ONLY", "false")) != "true")
  env$run_three_dataset_gene_prioritization(direction = "negative")
})
