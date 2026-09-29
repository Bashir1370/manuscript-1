############################################################
# GSE126773 OIPN
# RMA normalization and limma differential expression
# Contrast: Control + Oxaliplatin vs Control
############################################################

suppressPackageStartupMessages({
  library(affy)
  library(limma)
  library(dplyr)
})

message("Starting GSE126773 RMA + limma analysis")

manifest_file <- "results/GSE126773_OIPN/GSE126773_OIPN_analysis_samples.csv"

if(!file.exists(manifest_file)){
  stop("Missing analysis manifest: ", manifest_file)
}

manifest <- read.csv(
  manifest_file,
  stringsAsFactors = FALSE
)

if(!all(file.exists(manifest$file))){
  stop("One or more CEL files are missing")
}

message("Reading CEL files")

raw_data <- ReadAffy(
  filenames = manifest$file
)

message("Running RMA normalization")

eset <- rma(raw_data)

expr_matrix <- exprs(eset)

colnames(expr_matrix) <- manifest$GSM

out_dir <- "results/GSE126773_OIPN"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  expr_matrix,
  file.path(out_dir,"GSE126773_RMA_expression_matrix.csv"),
  row.names = TRUE
)

message("Building limma design matrix")

condition <- factor(
  manifest$group,
  levels = c("Vehicle","Oxaliplatin")
)

design <- model.matrix(~0 + condition)
colnames(design) <- levels(condition)

contrast_matrix <- makeContrasts(
  Oxaliplatin - Vehicle,
  levels = design
)

message("Running limma")

fit <- lmFit(expr_matrix, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(
  fit2,
  number = Inf,
  adjust.method = "BH"
)

results$probe_id <- rownames(results)

write.csv(
  results,
  file.path(out_dir,"GSE126773_Oxaliplatin_vs_Vehicle_limma_all_results.csv"),
  row.names = FALSE
)

sig <- results %>%
  filter(adj.P.Val < 0.05)

write.csv(
  sig,
  file.path(out_dir,"GSE126773_Oxaliplatin_vs_Vehicle_significant_DEGs_FDR_lt_0.05.csv"),
  row.names = FALSE
)

saveRDS(
  fit2,
  file.path(out_dir,"GSE126773_limma_object.rds")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir,"sessionInfo.txt")
)

message("GSE126773 limma analysis completed")
message("Significant probes: ", nrow(sig))
