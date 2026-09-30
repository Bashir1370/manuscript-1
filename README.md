# OIPN manuscript revision

Working repository for the computational redesign of the oxaliplatin-induced peripheral neuropathy (OIPN) project. This repository is public. Do not upload unpublished manuscript drafts or confidential data here.

## Current files

- [Reviewer comments](review/reviewer-comments-original.md) — reviewer text supplied by the author.
- [Structured reviewer response map](review/reviewer-feedback-summary.md) — working summary of the major issues.
- [Computational redesign plan](analysis/redesign-plan.md) — staged analyses and decision gates.
- [Stage 0 metadata audit](analysis/stage0-metadata-audit.md) — verified GEO sample design and unresolved assumptions.
- [GEO sample manifest](data/metadata/geo_sample_manifest.csv), [group counts](data/metadata/sample_groups.csv), and [source checksums](data/metadata/source_provenance.csv).
- [Metadata extraction script](scripts/01_audit_geo_metadata.py) — reproduces the audit tables from GEO series matrices.
- [Primary OIPN DESeq2 script](scripts/02_GSE160543_oxaliplatin_primary.R) — Oxaliplatin versus Vehicle only; download, input audit, complete-case gene alignment, DE and QC plots. Requires R and DESeq2.
- [Primary result tables and QC figures](results/GSE160543_Oxaliplatin_vs_Vehicle/) — complete output from the author's run of the primary script.
- [Primary result audit](analysis/gse160543-primary-results.md) — observed results, missing-gene caveat and interpretation limits.

## Current status

The submitted manuscript was removed and the default branch history was rebuilt without it. GitHub may retain old commits or cached views accessible by their direct commit IDs until purged by GitHub.

Official GEO sample metadata have been audited. The author's executed GSE160543 differential-expression tables and QC figures are available with source accession, GSM identifiers, script and archive MD5. The original submitted DEG and cross-model overlap numbers have not been independently reproduced. Animal-level behavioral/qPCR tables have not been added.

## Proposed organization as work proceeds

- `data/metadata/`: sample inventory, provenance and analysis-ready metadata.
- `data/processed/`: compact derived tables with documented inputs.
- `scripts/`: runnable analysis code and dependency lock information.
- `results/`: DEG, pathway, sensitivity and cell-signature outputs.
- `figures/`: generated figures.
- `review/`: response planning without unpublished manuscript text.

Public input datasets should be documented with accession, source URL, download date and checksum. Large raw archives can be referenced by accession rather than duplicated in Git.


## Current three-study branch

This branch is `three-dataset-neuropathy-analysis`, copied from the complete [four-study archive](https://github.com/Bashir1370/manuscript-1/tree/CCI-GSE212311-analysis). The retained studies are OIPN GSE160543, NC GSE246156 and CCI GSE212311. Read the [protocol and exclusion rationale](analysis/three_dataset_protocol.md) before interpreting the restricted analysis.

Run these scripts from repository root, in order:

1. [GSEA selection and manuscript heatmaps](scripts/cross_model/03_manuscript_Hallmark_heatmaps.R).
2. [Three-study GSVA comparison](scripts/cross_model/04_three_dataset_GSVA_comparison.R), using [helpers](scripts/cross_model/GSVA_helpers.R).
3. [Three-study leading-edge extraction and gene heatmaps](scripts/cross_model/05_three_dataset_leading_edge.R).
4. [OIPN-required gene prioritization](scripts/cross_model/06_three_dataset_gene_prioritization.R): original all-three membership and positive effects, gene FDR <0.05 in OIPN and NC or CCI; [method and results](analysis/three_dataset_gene_prioritization.md).

Audited results:

- [GSEA](results/manuscript_hallmark_three_dataset/): 10 shared positive, two shared negative, seven opposite-direction sets. The original >=3 significant-study threshold is retained and now requires all three.
- [GSVA](results/GSVA_three_dataset/): 150 original within-study limma/BH50 fits, 57 selected pathway-study rows and 380 sample-pathway rows. Models remain independent per study.
- [Leading edges](results/shared_Hallmark_leading_edge_three_dataset/): 109 all-three pathway-gene memberships representing 78 distinct genes from the 10 shared-positive programs.
- [Source provenance](analysis/three_dataset_input_provenance.csv). Table validation passed; R execution and figure review remain local. No new figures have been generated here.

Previous reports for retained studies are preserved as historical evidence; canonical selection rules and source-aware CCI outputs for this branch are specified in the protocol. Results of this post hoc restricted analysis do not establish agreement across all OIPN studies.

The post hoc OIPN-required priority list contains 20 unique genes, including the strict eight significant in all three studies. See [STRING list](results/shared_Hallmark_leading_edge_three_dataset/gene_prioritization/priority_genes_STRING.txt) and [full gene evidence](results/shared_Hallmark_leading_edge_three_dataset/gene_prioritization/priority_OIPN_plus_physical.csv).

### Opposite-direction Hallmark follow-up (three studies)

Run `source("scripts/cross_model/07_three_dataset_divergent_leading_edge.R")` from the repository root. This follows the seven Fig2 opposite-NES pathways using original leading-edge unions, evaluates gene direction, and prioritizes gene FDR < 0.05 in OIPN plus NC or CCI. [Methods and archived results](analysis/three_dataset_divergent_leading_edge.md); outputs: `results/divergent_Hallmark_leading_edge_three_dataset/`. The script generates plots locally; archived tables were independently reconstructed without R execution.

Gene log2FC heatmaps from scripts 05 and 07 use one common linear color range of -6 to +6 by default. `LE_COLOR_LIMIT` can set another shared range; each output folder records `plot_settings.csv`. See the [display audit](analysis/leading_edge_color_scale.md). This affects rendering only.
