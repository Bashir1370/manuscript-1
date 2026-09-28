# OIPN manuscript revision

Working repository for the computational redesign of the oxaliplatin-induced peripheral neuropathy (OIPN) project. This repository is public. Do not upload unpublished manuscript drafts or confidential data here.

## Current files

- [Reviewer comments](review/reviewer-comments-original.md) — reviewer text supplied by the author.
- [Structured reviewer response map](review/reviewer-feedback-summary.md) — working summary of the major issues.
- [Computational redesign plan](analysis/redesign-plan.md) — staged analyses and decision gates; **a plan, not completed results**.

## Current status

The submitted manuscript has been removed from the current branch. Earlier Git commits may still contain it until repository history is rewritten or the repository is replaced.

No raw count matrices, animal-level behavioral/qPCR tables, analysis scripts, or new figures have yet been added. The reported numbers have not been independently reproduced. Future results must be traceable to input accession/sample metadata and code.

## Proposed organization as work proceeds

- `data/metadata/`: sample inventory, provenance and analysis-ready metadata.
- `data/processed/`: compact derived tables with documented inputs.
- `scripts/`: runnable analysis code and dependency lock information.
- `results/`: DEG, pathway, sensitivity and cell-signature outputs.
- `figures/`: generated figures.
- `review/`: response planning without unpublished manuscript text.

Public input datasets should be documented with accession, source URL, download date and checksum. Large raw archives can be referenced by accession rather than duplicated in Git.
