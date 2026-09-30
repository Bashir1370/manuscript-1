#!/usr/bin/env Rscript
# Regenerate missing VST and compare DE/GSEA/GSVA without promoting new results.
run_oipn_reproducibility <- function() {
  root <- "results/GSE160543_Oxaliplatin_vs_Vehicle"
  primary <- file.path(root, "reproducibility_rerun", "primary")
  out <- file.path(root, "reproducibility_rerun", "Pathway_analysis")
  vars <- c("OIPN_PRIMARY_OUTDIR", "OIPN_PATHWAY_OUTDIR")
  previous <- Sys.getenv(vars, unset = NA_character_)
  on.exit({
    for (i in seq_along(vars)) {
      if (is.na(previous[i])) Sys.unsetenv(vars[i]) else
        do.call(Sys.setenv, stats::setNames(list(previous[i]), vars[i]))
    }
  }, add = TRUE)
  Sys.setenv(OIPN_PRIMARY_OUTDIR = primary, OIPN_PATHWAY_OUTDIR = out)
  env <- new.env(parent = globalenv())
  set.seed(160543L)
  sys.source("scripts/02_GSE160543_oxaliplatin_primary.R", envir = env)
  old <- read.csv(file.path(root, "DE_all_genes.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  new <- read.csv(file.path(primary, "DE_all_genes.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  if (anyDuplicated(old$gene_id) || anyDuplicated(new$gene_id)) stop("Duplicate DE comparison keys.")
  keys <- sort(union(old$gene_id, new$gene_id))
  a <- match(keys, old$gene_id); b <- match(keys, new$gene_id)
  compare <- data.frame(gene_id = keys, archived_present = !is.na(a), rerun_present = !is.na(b),
    archived_symbol = old$symbol[a], rerun_symbol = new$symbol[b],
    archived_log2FC = old$log2FoldChange[a], rerun_log2FC = new$log2FoldChange[b],
    archived_stat = old$stat[a], rerun_stat = new$stat[b],
    archived_FDR = old$padj[a], rerun_FDR = new$padj[b])
  compare$delta_log2FC <- compare$rerun_log2FC - compare$archived_log2FC
  compare$delta_stat <- compare$rerun_stat - compare$archived_stat
  compare$FDR_0_05_class_agrees <- (compare$archived_FDR < .05) == (compare$rerun_FDR < .05)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.csv(compare, file.path(out, "DE_archived_comparison.csv"), row.names = FALSE, na = "")
  sys.source("scripts/03_GSE160543_pathway_analysis.R", envir = env)
  sys.source("scripts/04_GSE160543_GO_GSVA_analysis.R", envir = env)
  writeLines(c(
    "Historical canonical DE/GSEA/GSVA tables were not replaced.",
    "Current acquired Hallmark membership is frozen; historical release was not recovered.",
    "A fixed seed cannot eliminate differences caused by software/gene-set/annotation versions.",
    "Review DE, GSEA and GSVA archived-comparison tables before choosing a sensitivity analysis.",
    "Cross-study scripts still use original canonical tables, not this rerun."
  ), file.path(out, "rerun_notes.txt"))
  message("OIPN reproducibility rerun complete: ", out)
}
run_oipn_reproducibility()
