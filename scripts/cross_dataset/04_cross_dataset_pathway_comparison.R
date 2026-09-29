# Cross-dataset Hallmark pathway comparison
# Purpose: compare pathway-level NES across independent CIPN datasets

message("Starting cross-dataset Hallmark pathway comparison")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
})

inputs <- list(
  GSE160543 = "results/GSE160543_Paclitaxel/pathway_analysis/GSE160543_Hallmark_ranked_GSEA_results.csv",
  GSE126773 = "results/GSE126773_OIPN/pathway_analysis/GSE126773_Hallmark_ranked_GSEA_results.csv",
  GSE212311 = "results/GSE212311_CCI_L4L6_day11/pathway_analysis/GSE212311_Hallmark_ranked_GSEA_results.csv"
)

available <- inputs[file.exists(unlist(inputs))]

if(length(available) == 0){
  stop("No GSEA result files found")
}

read_gsea <- function(path, dataset){
  x <- read.csv(path, stringsAsFactors = FALSE)
  
  nes_col <- intersect(c("NES", "nes"), colnames(x))
  pathway_col <- intersect(c("Description", "ID", "pathway", "gs_name"), colnames(x))
  
  if(length(nes_col) == 0 || length(pathway_col) == 0){
    stop("Could not identify pathway/NES columns in ", path)
  }
  
  data.frame(
    pathway = x[[pathway_col[1]]],
    NES = x[[nes_col[1]]],
    dataset = dataset,
    stringsAsFactors = FALSE
  )
}

all_results <- bind_rows(
  lapply(names(available), function(nm){
    read_gsea(available[[nm]], nm)
  })
)

out_dir <- "results/cross_dataset_pathway_comparison"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  all_results,
  file.path(out_dir, "Hallmark_NES_all_datasets.csv"),
  row.names = FALSE
)

matrix_df <- dplyr::select(
  all_results,
  pathway,
  dataset,
  NES
) %>%
  tidyr::pivot_wider(
    names_from = dataset,
    values_from = NES
  )

write.csv(
  matrix_df,
  file.path(out_dir, "Hallmark_NES_comparison_matrix.csv"),
  row.names = FALSE
)

heatmap_data <- as.data.frame(matrix_df)
rownames(heatmap_data) <- heatmap_data$pathway
heatmap_data$pathway <- NULL

pdf(file.path(out_dir, "Hallmark_NES_heatmap.pdf"), width = 8, height = 12)
pheatmap(
  as.matrix(heatmap_data),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Hallmark pathway NES comparison"
)
dev.off()

message("Cross-dataset pathway comparison completed")
