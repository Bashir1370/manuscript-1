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

## Current status

The submitted manuscript was removed and the default branch history was rebuilt without it. GitHub may retain old commits or cached views accessible by their direct commit IDs until purged by GitHub.

Official GEO sample metadata have been audited. No raw count matrices, animal-level behavioral/qPCR tables, executed differential-expression outputs, or new figures have yet been added. The reported DEG and overlap numbers have not been independently reproduced. Future results must be traceable to input accession/sample metadata and code.

## Proposed organization as work proceeds

- `data/metadata/`: sample inventory, provenance and analysis-ready metadata.
- `data/processed/`: compact derived tables with documented inputs.
- `scripts/`: runnable analysis code and dependency lock information.
- `results/`: DEG, pathway, sensitivity and cell-signature outputs.
- `figures/`: generated figures.
- `review/`: response planning without unpublished manuscript text.

Public input datasets should be documented with accession, source URL, download date and checksum. Large raw archives can be referenced by accession rather than duplicated in Git.
