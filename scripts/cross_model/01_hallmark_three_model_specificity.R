############################################################
# Hallmark three-model specificity analysis
#
# Goal:
# Compare Hallmark pathway enrichment across three neuropathy models:
#
# 1) OIPN  : Oxaliplatin-induced peripheral neuropathy
# 2) NC    : Physical injury neuropathy model
# 3) CCI   : Chronic constriction injury model
#
# Important:
# - Raw counts are NOT merged.
# - Comparison is performed only at Hallmark pathway level.
# - Each dataset remains an independent analysis.
#
############################################################


suppressPackageStartupMessages({
  
  library(tidyverse)
  library(pheatmap)
  
})


############################################################
# 1. INPUT FILES
############################################################

# CHANGE ONLY THIS SECTION IF FILE NAMES DIFFER

OIPN_file <- 
  "results/OIPN/Hallmark_GSEA_results.csv"


NC_file <- 
  "results/NC/Hallmark_GSEA_results.csv"


CCI_file <- 
  "results/GSE212311_CCI_L4L6_day11/pathway_analysis_source_aware/Hallmark_GSEA_results.csv"



input_files <- c(
  OIPN = OIPN_file,
  NC   = NC_file,
  CCI  = CCI_file
)



############################################################
# 2. CHECK INPUT FILES
############################################################


missing_files <- input_files[!file.exists(input_files)]


if(length(missing_files) > 0){
  
  stop(
    paste(
      "Missing input files:",
      paste(missing_files, collapse="\n")
    )
  )
  
}



############################################################
# 3. READ FUNCTION
############################################################


read_hallmark <- function(file, model){
  
  
  message("Reading: ", model)
  
  
  x <- read.csv(
    file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  
  required_columns <- c(
    "pathway",
    "NES",
    "padj"
  )
  
  
  missing_columns <- setdiff(
    required_columns,
    colnames(x)
  )
  
  
  if(length(missing_columns)>0){
    
    stop(
      paste(
        model,
        "missing columns:",
        paste(missing_columns, collapse=", ")
      )
    )
    
  }
  
  
  
  x %>%
    
    select(
      pathway,
      NES,
      padj
    ) %>%
    
    mutate(
      model = model
    )
  
  
}



############################################################
# 4. LOAD THREE DATASETS
############################################################


OIPN <- read_hallmark(
  OIPN_file,
  "OIPN"
)


NC <- read_hallmark(
  NC_file,
  "NC"
)


CCI <- read_hallmark(
  CCI_file,
  "CCI"
)



############################################################
# 5. INPUT QUALITY CONTROL
############################################################


message("Running input QC...")


if(
  identical(OIPN$NES, NC$NES) ||
  identical(OIPN$NES, CCI$NES) ||
  identical(NC$NES, CCI$NES)
){
  
  stop(
    paste(
      "QC FAILED:",
      "Two or more models have identical NES vectors.",
      "Input files are probably duplicated."
    )
  )
  
}



############################################################
# 6. COMBINE RESULTS
############################################################


combined <- bind_rows(
  OIPN,
  NC,
  CCI
)



############################################################
# 7. CREATE WIDE MATRIX
############################################################


wide <- combined %>%
  
  select(
    pathway,
    model,
    NES,
    padj
  ) %>%
  
  pivot_wider(
    
    names_from = model,
    
    values_from = c(
      NES,
      padj
    )
    
  )



############################################################
# 8. CLASSIFICATION
############################################################


FDR_cutoff <- 0.05



classified <- wide %>%
  
  rowwise() %>%
  
  mutate(
    
    
    OIPN_sig =
      !is.na(padj_OIPN) &
      padj_OIPN < FDR_cutoff,
    
    
    NC_sig =
      !is.na(padj_NC) &
      padj_NC < FDR_cutoff,
    
    
    CCI_sig =
      !is.na(padj_CCI) &
      padj_CCI < FDR_cutoff,
    
    
    
    classification = case_when(
      
      
      OIPN_sig &
        NC_sig &
        CCI_sig ~
        
        "Shared neuropathy program",
      
      
      
      OIPN_sig &
        !NC_sig &
        !CCI_sig ~
        
        "OIPN-specific",
      
      
      
      !OIPN_sig &
        NC_sig &
        CCI_sig ~
        
        "Physical injury-specific",
      
      
      
      TRUE ~
        
        "Other"
      
    )
    
  ) %>%
  
  ungroup()



############################################################
# 9. OUTPUT DIRECTORY
############################################################


output_dir <-
  "results/cross_model_hallmark_specificity"



dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)



############################################################
# 10. SAVE TABLES
############################################################



write.csv(
  
  wide,
  
  file.path(
    output_dir,
    "hallmark_three_model_matrix.csv"
  ),
  
  row.names = FALSE
  
)



write.csv(
  
  classified %>%
    filter(
      classification ==
        "Shared neuropathy program"
    ),
  
  file.path(
    output_dir,
    "hallmark_shared_neuropathy_programs.csv"
  ),
  
  row.names = FALSE
  
)



write.csv(
  
  classified %>%
    filter(
      classification ==
        "OIPN-specific"
    ),
  
  file.path(
    output_dir,
    "hallmark_OIPN_specific.csv"
  ),
  
  row.names = FALSE
  
)



write.csv(
  
  classified %>%
    filter(
      classification ==
        "Physical injury-specific"
    ),
  
  file.path(
    output_dir,
    "hallmark_physical_injury_specific.csv"
  ),
  
  row.names = FALSE
  
)




############################################################
# 11. HEATMAP
############################################################


heatmap_matrix <- wide %>%
  
  select(
    pathway,
    NES_OIPN,
    NES_NC,
    NES_CCI
  ) %>%
  
  column_to_rownames(
    "pathway"
  )


png(
  
  filename =
    file.path(
      output_dir,
      "hallmark_three_model_heatmap.png"
    ),
  
  width = 1800,
  
  height = 2500,
  
  res = 250
  
)



pheatmap(
  
  as.matrix(
    heatmap_matrix
  ),
  
  cluster_rows = TRUE,
  
  cluster_cols = FALSE,
  
  fontsize_row = 8
  
)


dev.off()



############################################################
# 12. AUDIT REPORT
############################################################


sink(
  
  file.path(
    output_dir,
    "input_audit.txt"
  )
  
)


cat(
  "Hallmark three model specificity analysis\n\n"
)



for(i in names(input_files)){
  
  cat(
    "\nMODEL:",
    i,
    "\nFILE:",
    input_files[i],
    "\n"
  )
  
}



cat("\n\nNumber of pathways:\n")

print(
  
  combined %>%
    
    count(model)
  
)



cat("\n\nClassification summary:\n")


print(
  
  classified %>%
    
    count(classification)
  
)



sink()



message(
  "Hallmark three-model specificity analysis completed."
)