############################################################
# GSE126773 OIPN
# CEL file audit and sample manifest construction
#
# Contrast:
# Control + Oxaliplatin vs Control
############################################################

suppressPackageStartupMessages({
  library(dplyr)
})

message("Starting GSE126773 CEL audit")

cel_dir <- "data/source-cache/GSE126773"

cel_files <- list.files(
  cel_dir,
  pattern = "\\.CEL\\.gz$",
  full.names = TRUE
)

if(length(cel_files) == 0){
  stop("No CEL files found in: ", cel_dir)
}

sample_manifest <- data.frame(
  file = cel_files,
  sample = basename(cel_files),
  stringsAsFactors = FALSE
)

sample_manifest <- sample_manifest %>%
  mutate(
    GSM = sub("_.*", "", sample),
    group = case_when(
      grepl("Control\\+OX", sample) ~ "Oxaliplatin",
      grepl("Control_", sample) ~ "Vehicle",
      grepl("Cancer\\+OX", sample) ~ "Cancer_Oxaliplatin",
      grepl("Cancer_", sample) ~ "Cancer",
      TRUE ~ "Unknown"
    )
  )

out_dir <- "results/GSE126773_OIPN"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  sample_manifest,
  file.path(
    out_dir,
    "GSE126773_CEL_sample_manifest.csv"
  ),
  row.names = FALSE
)

analysis_manifest <- sample_manifest %>%
  filter(
    group %in% c("Vehicle", "Oxaliplatin")
  )

write.csv(
  analysis_manifest,
  file.path(
    out_dir,
    "GSE126773_OIPN_analysis_samples.csv"
  ),
  row.names = FALSE
)

message("CEL audit completed")
message("Samples selected for OIPN analysis: ", nrow(analysis_manifest))

print(table(analysis_manifest$group))
