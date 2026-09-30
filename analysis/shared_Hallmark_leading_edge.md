# Current shared Hallmark leading-edge evidence

See [three_dataset_protocol.md](three_dataset_protocol.md) for the current three-study method, source provenance and results. The 10 shared-positive pathways contain 109 all-three leading-edge memberships representing 78 distinct genes. Outputs are in `results/shared_Hallmark_leading_edge_three_dataset/`; use `scripts/cross_model/05_three_dataset_leading_edge.R` for local reproduction and figures.

## OIPN-required follow-up prioritization (2026-09-30)
See [three_dataset_gene_prioritization.md](three_dataset_gene_prioritization.md). Keeping the 78 all-three leading-edge genes and positive effects in all three, requiring gene FDR <0.05 in OIPN and at least one of NC/CCI selects 20 unique genes (8 strict all-three and 12 OIPN+NC). Five NC+CCI-only significant genes are excluded from this priority list. This post hoc descriptive selection retains the full 78-gene evidence and original within-study gene FDR.
