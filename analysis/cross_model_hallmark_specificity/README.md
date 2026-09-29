# Three-model Hallmark pathway specificity analysis

## Objective

Compare Hallmark pathway enrichment across three independent neuropathy models:

- OIPN: oxaliplatin-induced peripheral neuropathy (chemical injury)
- NC: neuropathy model (physical injury)
- CCI: chronic constriction injury (physical nerve injury)

## Principle

Raw RNA-seq matrices are not merged.

Each dataset is analyzed independently. Comparison is performed only at the pathway level using existing Hallmark GSEA outputs.

## Outputs

The analysis separates:

1. Shared neuropathy programs
   - significant in OIPN + NC + CCI

2. OIPN-specific programs
   - significant only in OIPN

3. Physical injury-specific programs
   - significant in NC + CCI but not OIPN

## Script

`scripts/cross_model/01_hallmark_three_model_specificity.R`
