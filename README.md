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


## Manuscript Hallmark evidence

- [Evidence archive, selection rules and interpretation limits](analysis/manuscript_hallmark_evidence.md).
- [Complete R script for both manuscript heatmaps](scripts/cross_model/03_manuscript_Hallmark_heatmaps.R): run from repository root; default direction rule and separate support-pattern mode.
- [Audited classification tables](results/manuscript_hallmark/): six shared positive sets, no shared negative set, two opposite-direction sets (HEME METABOLISM and COMPLEMENT; FDR annotated rather than used for direction selection), and separate support selections. No new figures have been generated; execute the script locally in R.

## Four-study GSVA

- [Protocol, input audit and interpretation limits](analysis/GSVA_four_dataset_protocol.md).
- Run [GSE126773 GSVA scoring](scripts/GSE126773_OIPN/06_GSE126773_GSVA_Hallmark.R), then [four-study comparison and figure script](scripts/cross_model/04_four_dataset_GSVA_comparison.R), from the repository root. Both use [GSVA helpers](scripts/cross_model/GSVA_helpers.R).
- Tests use all 50 Hallmarks per study with limma/BH; eight GSEA-selected pathways are displayed. New outputs go to `results/GSE126773_OIPN/GSVA_Hallmark/` and `results/GSVA_four_dataset/`; historical Wilcoxon files are retained. Input alignment was audited; R execution and figure review remain to be performed locally.

## Shared Hallmark leading-edge genes

- [Selection, representative matching, observed counts and limitations](analysis/shared_Hallmark_leading_edge.md).
- [Archived gene evidence and overlaps](results/shared_Hallmark_leading_edge/): 116 distinct genes shared in at least three leading edges of the same pathway; nine distinct genes shared in all four, with gene and pathway FDR recorded separately.
- [Complete R extraction and heatmap script](scripts/cross_model/05_shared_Hallmark_leading_edge.R): run from repository root; `LE_TABLES_ONLY=true` suppresses figures. Source CSV extraction was independently validated; R runtime and figure review remain local.
