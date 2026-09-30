# Shared-positive Hallmark leading-edge evidence

## Question and selection

Do the six shared-positive Hallmark programs draw on the same genes across neuropathy studies? Extract existing leading-edge/core-enrichment lists for E2F TARGETS, G2M CHECKPOINT, IL6 JAK STAT3 SIGNALING, INTERFERON ALPHA RESPONSE, INTERFERON GAMMA RESPONSE and TNFA SIGNALING VIA NFKB. These pathways have positive NES in all four studies and pathway FDR <0.05 in at least three. No GSEA is rerun and no new gene-set release is introduced.

This is exploratory follow-up of pathways selected using the same data. Leading-edge membership, pathway significance, gene differential-expression significance and direction are separate evidence dimensions. A shared leading edge is not independent validation or a causal/cell-type-specific mechanism.

## Inputs and representative matching

Inputs are archived repository results at commit `6b18bfe716f7b6d2079aa875fa355bdc67a58779`. `results/shared_Hallmark_leading_edge/input_blob_shas.csv` records the exact 12 source CSV paths and Git blob identifiers. The newly uploaded GSVA outputs are preserved and not used for leading-edge selection.

| Study | Leading-edge field | Gene evidence and original representative |
|---|---|---|
| OIPN GSE160543 | `leadingEdge`, semicolon separated | Full DE table, highest signed statistic per symbol as in original GSEA; reconstructed rows verified against archived ranks |
| OIPN GSE126773 | `core_enrichment`, slash separated | Existing ranked probe table, selected by largest absolute moderated t in the original analysis; 14,658 unique symbols, no duplicate/tied representatives |
| NC GSE246156 | `leadingEdge`, semicolon separated | Archived Ensembl annotation plus full DE table; exclude ambiguous Ensembl-to-symbol mappings, select largest absolute statistic per symbol, Ensembl-ID tie break; verified against archived ranks |
| CCI GSE212311 | `leadingEdge`, semicolon separated | Existing source-aware feature-to-symbol audit, including original feature IDs and mapping source; verified against archived ranks |

The original ranked symbol universes contain 16,614 / 14,658 / 12,424 / 13,799 symbols, respectively. Gene log2FC and FDR are taken from the same feature/probe representative as the rank. They retain the original gene-level FDR family from DESeq2 or limma; no BH recalculation over the selected leading-edge genes is performed. DESeq2 and microarray estimates arise from different assays/designs; compare direction and within-study support rather than calibrated cross-study magnitudes. The original GSEA gene-set releases, mappings and ranked universes may differ; no harmonized rescore is implied.

## Membership rules

For each pathway separately, take the union of its four leading edges and enumerate each gene in all four studies:

- `leading_edge`: present in that study's original leading-edge list.
- `ranked_not_leading_edge`: present in the original GSEA rank but absent from that leading edge.
- `unavailable_to_GSEA`: absent from the original rank. Membership, effect and FDR are exported as missing, never zero. This can reflect assay coverage, annotation or statistical filtering; it does **not** prove the gene was unmeasured or unexpressed.

`shared_ge3` requires leading-edge membership in at least three studies for the **same pathway**. `shared_all4` requires all four. These are descriptive thresholds and do not supply an overlap p value. Membership counts do not filter on gene FDR or individual log2FC sign. Those are annotated separately. `all_four_ranked` identifies complete rank coverage; `positive_all4` additionally requires positive log2FC in all four. Counts of membership in significant pathways are provided without replacing the primary rule.

## Observed results

| Hallmark | Union leading-edge genes | Shared in ≥3 | Shared in ≥3 and ranked in all 4 | Shared in all 4 | Shared in ≥3 and positive log2FC in all 4 |
|---|---:|---:|---:|---:|---:|
| E2F TARGETS | 176 | 49 | 46 | 6 | 21 |
| G2M CHECKPOINT | 152 | 33 | 29 | 3 | 11 |
| IL6 JAK STAT3 SIGNALING | 61 | 11 | 11 | 1 | 3 |
| INTERFERON ALPHA RESPONSE | 83 | 16 | 15 | 0 | 5 |
| INTERFERON GAMMA RESPONSE | 151 | 29 | 29 | 2 | 14 |
| TNFA SIGNALING VIA NFKB | 145 | 19 | 18 | 0 | 5 |

The ≥3 selection contains 157 pathway–gene memberships representing **116 distinct genes**, not 157 independent genes. Eleven of those pathway–gene rows lack rank coverage in the fourth study. The all-four selection contains 12 pathway–gene memberships representing **nine distinct genes**:

- E2F TARGETS: Cdk1, Cdkn3, Ddx39a, E2f8, Kif22, Tcf19.
- G2M CHECKPOINT: Cdkn3, Kif11, Kif22.
- IL6 JAK STAT3 SIGNALING: Ccl7.
- INTERFERON GAMMA RESPONSE: Ccl7, Tnfaip6.
- No all-four leading-edge intersection for INTERFERON ALPHA RESPONSE or TNFA SIGNALING VIA NFKB.

All nine all-four genes have positive log2FC in every study. This does not mean each is individually significant in every study: inspect gene FDR and pathway FDR in the exported evidence table. Repeated membership across overlapping Hallmarks is not independent replication. A zero all-four intersection does not negate shared pathway enrichment; different leading-edge genes can contribute across studies, and rank coverage differs.

The resulting evidence supports partial convergence at the gene level within the shared programs, with substantial variation in contributing genes. It does not establish that these genes operate in the same DRG cell population or cause neuropathic pain.

## Outputs

In `results/shared_Hallmark_leading_edge/`:

- `gene_evidence_long.csv`: 3,072 pathway–gene–study records with availability, membership, original feature, rank statistic, log2FC, gene FDR and pathway NES/FDR.
- `gene_membership_summary.csv`: all 768 union pathway–gene rows with study-specific membership/effect/FDR and counts.
- `shared_ge3.csv` and `shared_all4.csv`: descriptive selections, 157 and 12 rows.
- `pathway_overlap_summary.csv` and `pathway_evidence.csv`: overlap counts and the 24 pathway–study evidence rows.
- `shared_gene_pathway_counts.csv`: 116 distinct ≥3 genes with the number of selected pathways in which each meets that rule.
- `input_blob_shas.csv`: exact source snapshot. Local R execution additionally exports input MD5 checksums and R session information.

## Reproduce and draw figures locally

```r
setwd("D:/manuscript_1/manuscript-1-four-dataset")
Sys.unsetenv("LE_TABLES_ONLY")
source("scripts/cross_model/05_shared_Hallmark_leading_edge.R")
```

Tables use base R; figures require ggplot2. Set `Sys.setenv(LE_TABLES_ONLY="true")` to regenerate tables without figures. The script validates all leading-edge identifiers against original ranks and stops if representative reconstruction differs. It leaves GSEA, DE, GSVA and prior manuscript results unchanged, but overwrites its own derived output files.

The script creates one PDF and PNG per Hallmark for the ≥3 genes, with original log2FC (no row z scaling), asterisk for original gene FDR <0.05, and grey/NA for unavailable-to-GSEA cells. The fill is clipped at ±2 for legibility only; numerical estimates remain unchanged in CSVs. All genes satisfying the membership rule are shown, including those whose effects are inconsistent in the fourth study. Figures are intended for local review before manuscript use.

## Validation and remaining limits

The archived tables were calculated from the original CSV blobs using JavaScript, then independently checked in Python: representative/rank matching in the three RNA-seq analyses, membership/count identities, unique keys, exact selected subsets, missing-data handling, feature effect/FDR joins and unique-gene totals. All 3,072 evidence rows and 768 union rows passed these checks; 197 study–pathway–gene cells are unavailable to GSEA. Script strings/delimiters and matching/selection logic were reviewed. R is unavailable in this environment, so the R script and figure rendering have not been executed here. No figures have been created. Local R execution is the remaining runtime check.

Leading-edge interpretation follows the [GSEA user guide](https://docs.gsea-msigdb.org/GSEA/GSEA_User_Guide/). Keep this analysis framed as data-derived evidence for shared transcriptional programs, consistent with the manuscript's revised question.
