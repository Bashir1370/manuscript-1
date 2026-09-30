# Three-study results review
Reviewed 2026-09-30. Results snapshot: `49421d3006fc3d9550aea9932ba941f541861bde`.
Studies: OIPN GSE160543 (4 control, 4 neuropathy), NC GSE246156 (3+3), CCI GSE212311 (3+3).

## Validation of uploaded R outputs
Downloaded the committed CSVs and checked them independently against the original study evidence used to prepare this branch. Checks passed:
- 50 GSEA pathways and classification rules: 10 shared positive, 2 shared negative, 7 divergent by NES sign.
- 150 GSVA study-pathway fits, 57 selected comparison rows, 380 selected sample-pathway rows. Verified BH correction over all 50 pathways separately per study, original-score group contrasts, and within-study display z scores.
- Leading-edge evidence: 3,075 rows, 1,025 pathway-gene union memberships, 109 all-three leading-edge memberships, 78 unique genes. Verified original gene/rank availability, leading-edge membership, rank statistic, log2FC and original gene FDR. Shared_ge3 and shared_all3 are identical for three studies.
Validation concerns numeric tables; the uploaded PDF/PNG figures have not been visually inspected in this review. No new figures were generated.

## Shared pathway results
Every pathway below meets the GSEA rule: same NES sign and within-study FDR <0.05 in all three studies. GSVA FDR values below use the 50-pathway testing family in each study.

| Pathway | NES OIPN / NC / CCI | GSVA FDR OIPN | GSVA FDR NC | GSVA FDR CCI | Shared leading-edge genes |
|---|---|---:|---:|---:|---:|
| ALLOGRAFT_REJECTION | 1.751 / 2.405 / 2.246 | 0.6038 | 0.0025 | 0.0066 | 20 |
| APOPTOSIS | 1.736 / 2.055 / 1.652 | 0.2566 | 0.0055 | 0.2147 | 6 |
| COAGULATION | 1.761 / 1.533 / 1.843 | 0.1654 | 0.1659 | 0.2147 | 4 |
| EPITHELIAL_MESENCHYMAL_TRANSITION | 1.996 / 2.075 / 1.832 | 0.0618 | 0.4311 | 0.3554 | 8 |
| G2M_CHECKPOINT | 2.902 / 2.646 / 1.483 | 0.0043 | 0.0005 | 0.3554 | 19 |
| IL6_JAK_STAT3_SIGNALING | 1.615 / 2.317 / 2.450 | 0.2740 | 0.0047 | 0.0018 | 8 |
| INTERFERON_ALPHA_RESPONSE | 2.217 / 1.705 / 1.809 | 0.0106 | 0.1530 | 0.4576 | 5 |
| INTERFERON_GAMMA_RESPONSE | 2.189 / 2.187 / 2.353 | 0.0415 | 0.0122 | 0.0292 | 18 |
| P53_PATHWAY | 2.041 / 1.678 / 1.551 | 0.1345 | 0.1530 | 0.2988 | 6 |
| TNFA_SIGNALING_VIA_NFKB | 1.690 / 2.571 / 2.470 | 0.1393 | 0.0021 | 0.0035 | 15 |
| FATTY_ACID_METABOLISM | -2.092 / -2.096 / -1.457 | 0.0618 | 0.0122 | 0.3554 | Not analyzed |
| OXIDATIVE_PHOSPHORYLATION | -2.774 / -2.115 / -2.149 | 0.0069 | 0.0277 | 0.0538 | Not analyzed |

All 12 shared pathways have concordant GSVA group-contrast direction in all studies. Among the ten shared positive pathways, INTERFERON_GAMMA_RESPONSE alone has GSVA FDR <0.05 in all three; G2M_CHECKPOINT, ALLOGRAFT_REJECTION, IL6_JAK_STAT3_SIGNALING and TNFA_SIGNALING_VIA_NFKB meet that threshold in two. OXIDATIVE_PHOSPHORYLATION meets it in OIPN and NC; CCI FDR 0.0538 does not meet the prespecified threshold.
E2F_TARGETS is positive in all three GSEA studies but CCI FDR approximately 0.3954 excludes it from the shared-significant set.

GSVA is complementary analysis of the same samples, with pathways selected using GSEA; it is not independent replication. A nonsignificant GSVA result does not establish absence of a change. The reused score workflows have different assay universes/mappings; raw deltas are descriptive within studies, not calibrated comparative effect sizes.

## Divergent directions
There are four pathways positive in OIPN and negative in both NC and CCI:
APICAL_JUNCTION, HEME_METABOLISM, MYOGENESIS, NOTCH_SIGNALING.
There are three with the reverse pattern:
COMPLEMENT, MYC_TARGETS_V1, REACTIVE_OXYGEN_SPECIES_PATHWAY.

**None of these seven is GSEA-significant in all three studies.** They support a descriptive opposite-sign panel with FDR annotations, not a statistically established between-model interaction or OIPN-specific mechanism.
NOTCH_SIGNALING is nonsignificant in every study. MYC_TARGETS_V1 has CCI NES 0.632 and FDR 0.9992, so its positive sign is weak evidence.
MYOGENESIS and APICAL_JUNCTION have opposite signs between OIPN and NC with GSEA FDR <0.05 in those two studies; CCI does not meet that threshold. HEME_METABOLISM has the same two-study pattern.
Among all 19 display pathways, GSVA direction agrees with GSEA for 19/19 in OIPN, 19/19 in NC, and 17/19 in CCI. The two CCI disagreements are APICAL_JUNCTION and MYOGENESIS; neither method is significant for either pathway in CCI.

## Shared genes
The 109 shared pathway-gene memberships comprise 78 unique genes. All have positive log2FC in all three studies. Original gene-level FDR is <0.05 in all three studies for **8** unique genes:
Bard1, Ccnf, Cdkn1a, Csf1, Il4r, Itgal, Kif22, Mcm5.
Among the remaining unique genes, 17 meet gene-level FDR <0.05 in two studies and 53 in one.
Leading-edge recurrence is therefore not equivalent to three-study gene-level significance.

Tgfb1 recurs in five selected pathways; Cdkn1a, Csf1 and Il4r recur in four each. Tgfb1 is gene-significant in only one study, while the latter three are significant in all three. Pathway counts reflect overlapping gene sets and must not be interpreted as independent replication, hub centrality or causal prioritization.

## Manuscript interpretation
The results support recurring injury-associated DRG transcriptional programs across these three cohorts, encompassing immune/inflammatory gene sets, cell-cycle/stress gene sets and reduced metabolic gene sets. These are interpretations of enrichment labels and expression, not proof of biochemical activation.
Bulk DRG cannot resolve cell origin or separate altered cell proportions from within-cell regulation. ALLOGRAFT_REJECTION does not indicate transplantation; EPITHELIAL_MESENCHYMAL_TRANSITION does not demonstrate literal epithelial transition in DRG; G2M does not demonstrate neuronal proliferation.

GSE126773 was excluded after observing divergent behavior, without an established sample-label or technical defect. State this post hoc decision explicitly; preserve the four-study analysis as a sensitivity comparison. With only one retained OIPN cohort, conclusions cannot be generalized as universal OIPN signatures. Cohort/model/time differences are confounded, and no direct between-model interaction test has been performed.

## Focused next analysis
1. Extend leading-edge extraction to the two shared negative pathways using the same source-aware rank mapping. Their genes were not included in the current ten-positive-pathway script.
2. For candidate follow-up, keep the 8 shared leading-edge genes with gene FDR <0.05 in all three separate from the wider 78-gene recurrence set; report full evidence rather than selecting by pathway recurrence alone.
3. Inspect the uploaded main figures visually before finalizing manuscript panel choice. A defensible primary narrative emphasizes the shared 10-positive/2-negative GSEA result, with transparent GSVA support and a separately labeled exploratory divergent panel.

Source directories:
- ../results/manuscript_hallmark_three_dataset/
- ../results/GSVA_three_dataset/
- ../results/shared_Hallmark_leading_edge_three_dataset/
