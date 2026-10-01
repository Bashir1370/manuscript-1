# Selected GO BP evidence and sample sensitivity

This stage completes sample-level evidence for the manuscript's locked 37 shared
and five opposite-direction genes, and compares GSEA with GSVA for the 11 shared
and nine opposite-direction representative pathways. It addresses the central
question of conserved DRG injury responses and OIPN's position relative to NC and
CCI. The five opposite genes are supported primarily by OIPN and NC.

## Run

From the repository root after stages 06 and 07:

```r
Sys.unsetenv("GOBP_TABLES_ONLY")
source("scripts/gobp/08_selected_evidence_and_sensitivity.R")
```

The normal entry point uses the existing `gobp_expression()` helper, including
DESeq2 size-factor estimation for NC. It does not refit differential expression,
GSEA, GSVA scoring or limma. Its default output is
`results/GO_BP_three_dataset/selected_evidence_42_genes_20_pathways/`.
The historical stage-05 expression archive remains intact. Use the new output
for the revised 42-gene selection, since stage 05 covers the previous priorities
and omits 11 of the currently selected shared genes.

## Scope and outputs

- `selected_genes.csv`: 42 genes, category and original log2FC/gene FDR.
- `selected_pathways.csv`: 20 ordered pathways and original NES/pathway FDR.
- `selected_gene_sample_expression.csv`: 840 gene-sample records, normalized
  counts, log2(count+1), display z-scores, original feature IDs and study labels.
- `gene_group_diagnostics.csv`: 252 gene-study-group records, means/ranges and
  largest sample's fraction of group counts.
- `GSEA_GSVA_selected_comparison.csv`: 60 pathway-study comparisons with original
  NES, GSVA delta, confidence intervals and FDR; no selected-subset adjustment.
- `GSEA_GSVA_concordance_summary.csv`: directional/significance summaries.
- `gene_leave_one_out.csv`, `gene_sensitivity_summary.csv`: 840 omissions and
  126 gene-study summaries.
- `pathway_leave_one_out.csv`, `pathway_sensitivity_summary.csv`: 400 omissions
  and 60 pathway-study summaries using the frozen GSVA scores.
- Two PDF/300-dpi PNG figures, input checksums, analysis notes and R session.

Selections are rederived from same-path Leading-edge membership and original
whole-study gene FDR. Missing genes, duplicate keys, altered group assignments,
inconsistent feature representatives/statistics and divergent score contrasts
stop execution before scientific tables are written. The present 37/5 and 11/9
counts are locked deliberately; changed selections require review.

## Definition of sensitivity

Each sample is omitted once, within its own study. All remaining control and
neuropathy samples are retained. For genes the effect is the difference between
group means of log2(normalized count+1). This descriptive effect is **not** a
refitted DESeq2 log2FC. Normalization is fixed. For pathways the effect is the
difference between group means of the archived per-sample GSVA scores. GSVA
scoring and normalization are fixed, so this is not full-pipeline resampling.

Direction stability means every single-sample omission retains the strict sign
of the descriptive baseline. A zero baseline is not classified as stable.
Agreement with the original effect is exported separately. No additional
p-values, FDR, automatic exclusions or new gene selections are produced.
Stability of direction does not establish significance, cellular origin or
causality. GSEA and GSVA on the same samples are complementary descriptions,
not independent validation. Differences between models remain confounded with
study, sampling time and design. The prior post hoc study restriction applies.

## Verified results

All 126 gene-study descriptive baselines agree in sign with the original DE
estimates. The single-sample omission results are:

| Evidence | OIPN | NC | CCI |
|---|---:|---:|---:|
| Shared genes retaining direction in every omission | 37/37 | 36/37 | 33/37 |
| Opposite genes retaining direction in every omission | 5/5 | 5/5 | 1/5 |
| Shared pathway GSVA deltas retaining direction | 11/11 | 10/11 | 11/11 |
| Opposite pathway GSVA deltas retaining direction | 9/9 | 4/9 | 1/9 |

Cdk1 and Cdkn1a retain their positive descriptive direction in all three
studies for every omission. The unstable shared genes are Tln1 in NC, and
Cd14, Parp14, Dhcr7 and Syt11 in CCI. In CCI the opposite genes Cav1, Cdh5,
Tns2 and Wnt6 can change sign after omission; Col4a2 remains negative.
The shared pathway tissue remodeling is direction-sensitive in NC.

GSEA versus GSVA direction agreement is 11/11 in every study for shared terms,
and 9/9, 6/9 and 3/9 in OIPN, NC and CCI for opposite terms. GSVA significance
for the shared representatives is 0/11, 7/11 and 0/11; for opposite terms it is
0/9, 2/9 and 0/9. The pathway stability rows above refer to each GSVA baseline,
not necessarily to the GSEA direction.

## Local production execution and verification

The default entry point was executed locally with R 4.5.2 on Windows 11 and
DESeq2 1.50.2. The scientific tables and four binary figures were archived at
commit `b08296108e276225ce168f36c4bd53ba903f86e1`. Numeric values match the
previous independently reproduced normalization results within 1e-10 tolerance.
All 33 recorded input checksums match repository inputs allowing for LF versus
Windows CRLF serialization. No original DE/GSEA/GSVA source was changed by
this upload. All 840 gene and 400 pathway omission contrasts, log transforms
and within-study sample z-scores were independently recomputed and verified.

The earlier validation used archived OIPN/CCI normalized matrices and an
independently reproduced NC median-ratio normalization because DESeq2 was not
available in that runtime. The subsequent local production execution now
confirms those scientific values with the default DESeq2 helper. Synthetic
stable, sign-changing and zero effects; a one-row matrix; missing data;
reordered samples; negative count rejection; and loss of mandatory OIPN
Leading-edge support were also checked in the earlier validation.

## Manuscript presentation update

The stage-08 plotting function now draws category boundaries in both figures,
control/neuropathy boundaries in the gene figure, and explicit group headers.
The revised figures are drawn directly from the uploaded production tables;
no normalization, scoring, selection or scientific table is rerun for this
presentation update. `figure_inputs_checksums.csv` identifies those tables and
the updated script, separately from the original production `input_checksums.csv`
and `R_sessionInfo.txt`, which remain the provenance of the local table run.
`figure_render_sessionInfo.txt` records the rendering runtime.

Use [GO_BP_manuscript_results_methods_legends.md](GO_BP_manuscript_results_methods_legends.md)
for the manuscript-ready computational section and legends. Supplementary
Table S3 contains only the compact selected-evidence sensitivity and direction
agreement summary; the full opposite-direction pathway list is not appended.
The stable opposite-category GSVA term in CCI is Schwann cell differentiation;
its positive GSVA effect opposes its negative GSEA NES. Stability must not be
reported as agreement with GSEA or renewed statistical significance.
