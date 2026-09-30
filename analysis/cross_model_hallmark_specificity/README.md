# Three-model Hallmark pathway specificity analysis

## Objective

Compare Hallmark pathway enrichment across three independent neuropathy models:

- OIPN: oxaliplatin-induced peripheral neuropathy (chemical injury)
- NC: neuropathy model (physical injury)
- CCI: chronic constriction injury (physical nerve injury)

## Principle

Raw RNA-seq matrices are not merged.

Each dataset is analyzed independently. Comparison is performed only at the pathway level using existing Hallmark GSEA outputs.

This avoids direct integration of heterogeneous transcriptomic datasets while allowing comparison of biological programs.

## Input datasets

Hallmark GSEA outputs used:

- OIPN:
  `results/GSE160543_Oxaliplatin_vs_Vehicle/Pathway_analysis/GSEA_Hallmark_results.csv`

- NC:
  `results/GSE246156_NC_L5_day7/pathway_analysis/GSEA_Hallmark_results.csv`

- CCI:
  `results/GSE212311_CCI_L4L6_day11/pathway_analysis_source_aware/GSEA_Hallmark_results.csv`

## Classification

The analysis separates:

1. Shared neuropathy programs
   - significant in OIPN + NC + CCI

2. OIPN-enriched programs
   - significant in OIPN but not NC or CCI

3. Physical injury-associated programs
   - significant in NC + CCI but not OIPN

## Outputs

Generated files are stored in:

`results/cross_model_hallmark_specificity/`

Including:

- three-model Hallmark matrix
- shared pathway list
- OIPN-enriched pathway list
- physical injury-associated pathway list
- heatmap visualization
- input audit report

## Script

`scripts/cross_model/01_hallmark_three_model_specificity.R`
