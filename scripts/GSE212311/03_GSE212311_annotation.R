# GSE212311 DEG annotation
# Map DESeq2 results to rat gene symbols

suppressPackageStartupMessages({
  library(tidyverse)
  library(org.Rn.eg.db)
  library(AnnotationDbi)
})

base_dir <- "results/GSE212311_CCI_L4L6_day11"
out_dir <- file.path(base_dir, "annotation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(
  base_dir,
  "DESeq2",
  "GSE212311_significant_DEGs_FDR_lt_0.05.csv"
)

if (!file.exists(input_file)) {
  stop("DESeq2 DEG file not found: ", input_file)
}

deg <- read.csv(input_file, stringsAsFactors = FALSE)

if (!"gene" %in% colnames(deg)) {
  stop("Expected gene column was not found.")
}

mapping <- AnnotationDbi::select(
  org.Rn.eg.db,
  keys = unique(deg$gene),
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "ENTREZID")
) %>%
  distinct()

annotated <- deg %>%
  left_join(mapping, by = c("gene" = "ENSEMBL"))

write.csv(
  annotated,
  file.path(out_dir, "GSE212311_DEG_annotated.csv"),
  row.names = FALSE
)

write.csv(
  mapping,
  file.path(out_dir, "annotation_mapping_table.csv"),
  row.names = FALSE
)

unmapped <- annotated %>%
  filter(is.na(SYMBOL))

write.csv(
  unmapped,
  file.path(out_dir, "unmapped_genes.csv"),
  row.names = FALSE
)

summary_text <- c(
  paste0("Total DEG input: ", nrow(deg)),
  paste0("Mapped symbols: ", sum(!is.na(annotated$SYMBOL))),
  paste0("Unmapped genes: ", sum(is.na(annotated$SYMBOL))),
  paste0("Mapping percentage: ", round(mean(!is.na(annotated$SYMBOL))*100, 2), "%")
)

writeLines(
  summary_text,
  file.path(out_dir, "annotation_summary.txt")
)

message("GSE212311 annotation completed")
