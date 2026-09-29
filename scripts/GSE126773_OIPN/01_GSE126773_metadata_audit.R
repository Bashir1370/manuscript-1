# GSE126773 OIPN metadata audit
# Rat DRG oxaliplatin validation study

library(GEOquery)
library(dplyr)

message('Starting GSE126773 metadata audit')

outdir <- 'results/GSE126773_OIPN/metadata'
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

gse <- getGEO('GSE126773', GSEMatrix = TRUE)
eset <- gse[[1]]

pheno <- Biobase::pData(eset)
write.csv(pheno, file.path(outdir,'GSE126773_raw_phenoData.csv'), row.names = FALSE)

sample_groups <- data.frame(
  sample = rownames(pheno),
  title = pheno$title,
  source = pheno$source_name_ch1,
  stringsAsFactors = FALSE
)

sample_groups$group <- ifelse(grepl('oxaliplatin', sample_groups$title, ignore.case=TRUE),
                              'Oxaliplatin',
                              ifelse(grepl('control|vehicle|naive', sample_groups$title, ignore.case=TRUE),
                                     'Control','Review'))

write.csv(sample_groups,
          file.path(outdir,'GSE126773_sample_groups.csv'),
          row.names=FALSE)

message('Metadata audit completed')
