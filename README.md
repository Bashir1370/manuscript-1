# OIPN manuscript revision

## Current GO:BP exploratory branch

This branch is `three-dataset-gobp-exploration`, based on the three-study archive
at `99d5a8e31e25126ff99e4392de1657ce71dc16cc`. The current experimental workflow
uses **Human MSigDB C5:GO:BP mapped to rat** for OIPN GSE160543, NC GSE246156
and CCI GSE212311. Inherited Hallmark files below remain historical archives.

Read the [GO:BP protocol and local commands](analysis/GO_BP_three_dataset_protocol.md).
From repository root, run the complete entry script:

```r
source("scripts/gobp/run_GO_BP_three_dataset.R")
```

It locks one GO release, performs GSEA and classifications, rescores GSVA,
extracts shared positive/negative and divergent leading edges, prioritizes genes
with OIPN-required evidence, and displays sample expression. Outputs are isolated
under `results/GO_BP_three_dataset/`; the new gene-set archive is
`data/gene_sets/GO_BP_rat_locked/`. Figures are drawn by the author's local R run.
Set `GOBP_TABLES_ONLY=true` to compute the complete tables without plots.

For the author-selected 11 shared GO:BP representatives with adjacent gene
evidence, run `source("scripts/gobp/06_representative_GO_BP_heatmap.R")` after
stages 02 and 04. See the [selection and plotting notes](analysis/GO_BP_representative_11.md).
This optional stage reads completed results; it does not repeat the experiment.

Important: full-family GO BH replaces BH50. Shared pathway significance still
requires all three studies. Opposite NES selection still has no FDR filter.
GSVA uses a uniform log2(normalized count +1) transform in this experiment;
this differs from historical VST GSVA and does not require the parked script 09.
GO terms overlap and positive NES does not establish functional activation.

[Functional regression validation](analysis/GO_BP_workflow_validation.txt) and
[repeatable validation script](scripts/gobp/validate_GO_BP_workflow.R) cover
archived ranks/samples/features, known Hallmark classification and gene counts,
missing and empty selections, and full-family limma BH. **They are not a completed
GO:BP experiment:** the gene-set lock, new GO results and figures require local R.

---


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


## Inherited three-study Hallmark archive

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

Latest [integrated repository/results review](analysis/three_dataset_integrated_review.md) audits uploaded results at `6d199c8`, distinguishes shared and exploratory evidence, records upstream reproduction gaps, and prioritizes the next analyses. [Numeric checks](analysis/three_dataset_integrated_review_checks.json).

### Shared-negative completion and reproducibility preparation

[Methods, verified results and local execution](analysis/negative_leading_edge_completion.md):
the two shared-negative Hallmarks contribute 44 all-three leading-edge memberships,
39 unique genes with negative effects in all three, seven OIPN-required priorities,
and one strict all-three gene (Hsph1). [Archived tables](results/shared_negative_Hallmark_leading_edge_three_dataset/)
were reconstructed from original inputs without R execution or plot generation.

- Run `scripts/cross_model/09_OIPN_reproducibility_rerun.R` to regenerate primary/VST,
  GSEA and GSVA in separate rerun directories and compare them with the archive.
  It repairs old pathway/CSV handling, freezes the currently acquired [Hallmark membership](data/gene_sets/README.md),
  and records seed, parameters, package versions and checksums. Historical membership
  is not recovered; canonical cross-study inputs are not replaced.
- Run `scripts/cross_model/08_three_dataset_negative_leading_edge.R` for the negative
  tables, gene priorities and two local heatmaps. Positive script defaults remain intact.
- Run `scripts/cross_model/10_three_dataset_sample_expression_context.R` for descriptive
  sample expression of existing priorities. Cell composition and neuronal localization
  are not inferred. Context numerical results and all new figures remain pending local R execution.

The older integrated review remains a dated record; its negative-leading-edge and
upstream-code gaps are addressed by this addition, with runtime verification still pending.
