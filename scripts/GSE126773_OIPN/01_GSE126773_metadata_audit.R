# GSE126773 OIPN metadata audit
# Purpose: verify sample groups before microarray processing

message('Starting GSE126773 OIPN metadata audit')

out_dir <- 'results/GSE126773_OIPN'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- data.frame(
  sample = character(),
  group = character(),
  notes = character(),
  stringsAsFactors = FALSE
)

write.csv(metadata,
          file.path(out_dir,'GSE126773_sample_metadata_template.csv'),
          row.names = FALSE)

message('Metadata template created. Populate after GEO sample verification.')
