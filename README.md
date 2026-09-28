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
