# GSE126773 CEL download preparation

message('Preparing GSE126773 CEL download')

cache_dir <- 'data/source-cache/GSE126773'
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

message('Register GEO supplementary CEL files before download.')
