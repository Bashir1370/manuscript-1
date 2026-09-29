#!/usr/bin/env Rscript

# Three-model Hallmark pathway specificity analysis
#
# Purpose:
# Compare already generated Hallmark GSEA outputs from:
#   1) OIPN  (chemical neuropathy; GSE160543)
#   2) NC    (physical injury neuropathy; GSE246156)
#   3) CCI   (physical nerve injury; GSE212311)
#
# No raw counts are merged. Comparison is performed only at pathway level.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
})

root <- getwd()
out_dir <- file.path(root, "results", "cross_model_hallmark_specificity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

find_hallmark <- function(patterns) {
  candidates <- list.files(
    file.path(root, "results"),
    pattern = "csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  hit <- candidates[sapply(candidates, function(x) {
    any(grepl(paste(patterns, collapse = "|"), x, ignore.case = TRUE))
  })]

  if (length(hit) == 0) {
    stop("Could not identify Hallmark GSEA result file for: ", paste(patterns, collapse = ", "))
  }

  hit[1]
}

read_hallmark <- function(path, model) {
  x <- read_csv(path, show_col_types = FALSE)

  pathway_col <- names(x)[grepl("pathway|term|name|geneset", names(x), ignore.case = TRUE)][1]
  nes_col <- names(x)[grepl("NES", names(x), ignore.case = TRUE)][1]
  fdr_col <- names(x)[grepl("FDR|padj|adj", names(x), ignore.case = TRUE)][1]

  if (is.na(pathway_col) || is.na(nes_col) || is.na(fdr_col)) {
    stop("Required columns not found in ", path)
  }

  x %>%
    transmute(
      pathway = .data[[pathway_col]],
      NES = as.numeric(.data[[nes_col]]),
      FDR = as.numeric(.data[[fdr_col]]),
      model = model
    )
}

# Adjust these patterns if filenames change in future
op <- read_hallmark(find_hallmark(c("GSE160543", "hallmark", "gsea")), "OIPN")
nc <- read_hallmark(find_hallmark(c("GSE246156", "hallmark", "gsea")), "NC")
cci <- read_hallmark(find_hallmark(c("GSE212311", "hallmark", "gsea")), "CCI")

all <- bind_rows(op, nc, cci) %>%
  mutate(
    significant = FDR < 0.05,
    direction = case_when(
      NES > 0 ~ "UP",
      NES < 0 ~ "DOWN",
      TRUE ~ "NA"
    )
  )

write_csv(all, file.path(out_dir, "hallmark_three_model_matrix.csv"))

wide <- all %>%
  select(pathway, model, significant, direction, NES, FDR) %>%
  pivot_wider(names_from = model, values_from = c(significant, direction, NES, FDR))

write_csv(wide, file.path(out_dir, "hallmark_three_model_wide.csv"))

shared <- wide %>%
  filter(significant_OIPN == TRUE,
         significant_NC == TRUE,
         significant_CCI == TRUE)

write_csv(shared, file.path(out_dir, "hallmark_shared_neuropathy_programs.csv"))

opipn_specific <- wide %>%
  filter(significant_OIPN == TRUE,
         significant_NC != TRUE,
         significant_CCI != TRUE)

write_csv(opipn_specific, file.path(out_dir, "hallmark_OIPN_specific.csv"))

physical_specific <- wide %>%
  filter(significant_OIPN != TRUE,
         significant_NC == TRUE,
         significant_CCI == TRUE)

write_csv(physical_specific, file.path(out_dir, "hallmark_physical_injury_specific.csv"))

heat <- all %>%
  mutate(score = ifelse(FDR < 0.05, NES, 0)) %>%
  select(pathway, model, score) %>%
  pivot_wider(names_from = model, values_from = score) %>%
  tibble::column_to_rownames("pathway")

pheatmap(as.matrix(heat),
         filename = file.path(out_dir, "hallmark_three_model_heatmap.png"))

writeLines(
  c(
    "Three-model Hallmark specificity analysis completed.",
    "Comparison performed at pathway level only.",
    "Raw expression matrices were not merged.",
    "Models: OIPN, NC, CCI."
  ),
  file.path(out_dir, "analysis_summary.txt")
)

message("Completed: ", out_dir)
