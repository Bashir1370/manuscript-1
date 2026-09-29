############################################################
# GSE126773 OIPN
# Ranked Hallmark GSEA
#
# Full limma moderated t-statistic ranking.
# No DEG filtering is performed.
############################################################

message("Starting GSE126773 ranked Hallmark GSEA")

suppressPackageStartupMessages({
  library(dplyr)
  library(clusterProfiler)
  library(msigdbr)
})

rank_file <- "results/GSE126773_OIPN/annotation/GSE126773_ranked_gene_statistics.csv"

if(!file.exists(rank_file)){
  stop("Missing ranked statistics file: ", rank_file)
}

rank_df <- read.csv(
  rank_file,
  stringsAsFactors = FALSE
)

required_cols <- c("SYMBOL", "stat")
missing_cols <- setdiff(required_cols, colnames(rank_df))

if(length(missing_cols) > 0){
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

rank_df <- rank_df %>%
  filter(
    !is.na(SYMBOL),
    SYMBOL != "",
    !is.na(stat)
  ) %>%
  group_by(SYMBOL) %>%
  slice_max(
    order_by = abs(stat),
    n = 1
  ) %>%
  ungroup()

geneList <- rank_df$stat
names(geneList) <- rank_df$SYMBOL

geneList <- sort(
  geneList,
  decreasing = TRUE
)

message("Ranked genes: ", length(geneList))

hallmark <- msigdbr(
  species = "Rattus norvegicus",
  collection = "H"
)

hallmark_sets <- dplyr::select(
  hallmark,
  gs_name,
  gene_symbol
)

out_dir <- "results/GSE126773_OIPN/pathway_analysis"

if(!dir.exists(out_dir)){
  dir.create(out_dir, recursive = TRUE)
}

set.seed(123)

gsea_res <- GSEA(
  geneList = geneList,
  TERM2GENE = hallmark_sets,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  verbose = FALSE
)

write.csv(
  as.data.frame(gsea_res),
  file.path(out_dir, "GSE126773_Hallmark_ranked_GSEA_results.csv"),
  row.names = FALSE
)

saveRDS(
  gsea_res,
  file.path(out_dir, "GSE126773_Hallmark_GSEA_object.rds")
)

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "sessionInfo.txt")
)

message("GSE126773 Hallmark GSEA completed")
