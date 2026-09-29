############################################################
# GSE126773 OIPN
# Probe annotation and ranked statistics preparation
#
# Input:
#   GSE126773_Oxaliplatin_vs_Vehicle_limma_all_results.csv
#
# Output:
#   Annotated limma results
#   Ranked statistics for GSEA
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(annotate)
  library(rat2302.db)
})

message("Starting GSE126773 probe annotation")

input_file <-
  "results/GSE126773_OIPN/GSE126773_Oxaliplatin_vs_Vehicle_limma_all_results.csv"

if(!file.exists(input_file)){
  stop("Missing limma results file: ", input_file)
}

res <- read.csv(
  input_file,
  stringsAsFactors = FALSE
)

if(!"ID" %in% colnames(res)){
  stop("Expected probe identifier column 'ID' not found")
}

probe_ids <- as.character(res$ID)

annotation <- data.frame(
  ID = probe_ids,
  SYMBOL = getSYMBOL(
    probe_ids,
    "rat2302.db"
  ),
  stringsAsFactors = FALSE
)

annotated <- res %>%
  left_join(
    annotation,
    by = "ID"
  )

out_dir <-
  "results/GSE126773_OIPN/annotation"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  annotated,
  file.path(
    out_dir,
    "GSE126773_Oxaliplatin_vs_Vehicle_limma_annotated.csv"
  ),
  row.names = FALSE
)

rank_df <- annotated %>%
  filter(
    !is.na(SYMBOL),
    SYMBOL != "",
    !is.na(t)
  ) %>%
  group_by(SYMBOL) %>%
  slice_max(
    order_by = abs(t),
    n = 1
  ) %>%
  ungroup()

geneList <- rank_df$t
names(geneList) <- rank_df$SYMBOL

geneList <- sort(
  geneList,
  decreasing = TRUE
)

saveRDS(
  geneList,
  file.path(
    out_dir,
    "GSE126773_ranked_t_statistics.rds"
  )
)

write.csv(
  rank_df,
  file.path(
    out_dir,
    "GSE126773_ranked_gene_statistics.csv"
  ),
  row.names = FALSE
)

message("Annotation completed")
message("Ranked genes prepared: ", length(geneList))
