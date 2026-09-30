# Integrated three-study repository review

Review date: 2026-09-30. Reviewed uploaded-results commit: `6d199c8c0207267e2c60d9c4812e798b066c7bde`.
Scope: full recursive tree (283 files, 25 scripts, 64 PDF/PNG artifacts), 101 script/protocol/reviewer/summary/session text files, central GSEA/GSVA/leading-edge tables, original study GSEA/rank files and original gene evidence. This is a repository and numeric-evidence audit, not a new DESeq2/GSVA run. Figures are present but have not been visually inspected here.

## What has been completed

Three independent rat DRG contrasts: OIPN GSE160543 4+4, NC GSE246156 L5 day 7 3+3, CCI GSE212311 ipsilateral L4-L6 day 11 3+3. GSE126773 was excluded after observing discordance; its four-study archive is preserved. This is a post hoc cohort restriction, not a verified sample-label correction.

Canonical cross-study scripts are 03 (GSEA), 04 with GSVA_helpers (GSVA), 05 (shared-positive leading edges), 06 (gene prioritization) and 07 (opposite-direction leading-edge unions). Source-aware CCI outputs are canonical; older Ensembl-only and old specificity reports are historical.

Uploaded results include GSEA 50-pathway supplementary and selected panels, GSVA sample displays, ten shared-positive and seven divergent gene heatmaps, runtime session files and input checksums. Both gene-plot settings files report the common linear -6 to +6 color scale. A clean remote upload is confirmed by the current commit message and file inventory.

## Revalidated numeric results

| Analysis | Finding |
|---|---|
| GSEA shared positive | 10 pathways, positive NES and FDR <0.05 in every retained study |
| GSEA shared negative | 2 pathways, negative NES and FDR <0.05 in every retained study |
| Opposite-NES Fig2 panel | 7 pathways, direction selection without a pathway FDR filter |
| GSVA | 150 independent study-pathway fits, BH over 50 per study |
| Shared-positive leading edges | 109 same-pathway/all-three memberships, 78 unique genes |
| Shared gene priorities | 20 unique genes significant in OIPN plus NC/CCI; strict all-three subset 8 |
| Divergent leading edges | 616 union memberships; 133 opposite-direction memberships, 126 unique genes |
| Divergent gene priorities | 12 unique genes significant in OIPN plus NC/CCI; strict all-three subset 0 |

All 50 GSEA NES/FDR cells were checked against original source outputs. BH correction was independently recalculated for all 150 GSVA fits; selected group-score contrasts and within-study display z scores were checked using 380 sample records. Across shared and divergent long tables, 4593 rank-available evidence rows matched original feature-specific log2FC and adjusted P values, treating blank and NA as missing. Rank availability, signed statistics and original leading-edge membership also matched. Original NA adjusted P values remain missing, not significant. Detailed checks are in the accompanying JSON.

## Scientific interpretation

The shared-positive programs span immune/inflammatory labels (ALLOGRAFT_REJECTION, IL6_JAK_STAT3_SIGNALING, INTERFERON_ALPHA_RESPONSE, INTERFERON_GAMMA_RESPONSE, TNFA_SIGNALING_VIA_NFKB), cell-cycle/stress labels (G2M_CHECKPOINT, P53_PATHWAY, APOPTOSIS), and tissue-remodeling labels (COAGULATION, EPITHELIAL_MESENCHYMAL_TRANSITION). These are overlapping expression signatures, not independent mechanisms or functional proof of their names.

The two negative programs are FATTY_ACID_METABOLISM and OXIDATIVE_PHOSPHORYLATION. All twelve shared programs also have a concordant GSVA group-contrast direction. INTERFERON_GAMMA_RESPONSE alone among the shared-positive programs has GSVA FDR <0.05 in every study (OIPN 0.0415, NC 0.0122, CCI 0.0292). G2M_CHECKPOINT is GSVA-significant in OIPN and NC, not CCI. OXIDATIVE_PHOSPHORYLATION is GSVA-significant in OIPN and NC; CCI FDR 0.0538 remains above the threshold. GSVA complements GSEA on the same samples; it is not independent replication.

The 20 shared priorities comprise the strict eight (Bard1, Ccnf, Cdkn1a, Csf1, Il4r, Itgal, Kif22, Mcm5) and twelve additional OIPN+NC genes with positive but nonsignificant CCI effects. Their composition supports discussing cell-cycle-associated and immune-associated expression together, without neuronal localization or cell-cycle re-entry claims. The same symbol recurring in several Hallmarks does not provide independent replication or establish network centrality.

Four divergent pathways are positive in OIPN and negative in NC/CCI (APICAL_JUNCTION, HEME_METABOLISM, MYOGENESIS, NOTCH_SIGNALING); three have the reverse pattern (COMPLEMENT, MYC_TARGETS_V1, REACTIVE_OXYGEN_SPECIES_PATHWAY). None is GSEA-significant in every study. NOTCH_SIGNALING is nonsignificant throughout. CCI GSVA direction disagrees with CCI GSEA for APICAL_JUNCTION and MYOGENESIS; both methods are nonsignificant in those CCI comparisons.

The twelve divergent priorities are Plat/Plaur (COMPLEMENT), Ahsp/Alas2/Epb42/Fbxo9/Kel/Slc4a1 (HEME_METABOLISM), Psmd8 (MYC_TARGETS_V1), Col4a2/Smtn (MYOGENESIS), and Glrx (ROS). Six come from the heme program. The official Hallmark description includes erythroblast differentiation, and AHSP is an erythroid-associated protein. Therefore, a possible erythroid/tissue-composition contribution is an interpretation to investigate, not an established diagnosis of contamination. These data cannot establish OIPN specificity, ferroptosis or neuronal oxidative injury.
Sources: [MSigDB HEME_METABOLISM](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_HEME_METABOLISM.html), [original AHSP study](https://pubmed.ncbi.nlm.nih.gov/12066189/).

## Sensitivity and study limitations

The four-study archive has six shared-positive pathways and no shared-negative pathway under its original >=3 significance/all-four direction rule. Five positives persist in the current ten: G2M_CHECKPOINT, IL6_JAK_STAT3_SIGNALING, INTERFERON_ALPHA_RESPONSE, INTERFERON_GAMMA_RESPONSE and TNFA_SIGNALING_VIA_NFKB. E2F_TARGETS is excluded here because CCI FDR does not meet the all-three threshold; five additional positives are gained because the excluded OIPN study has discordant signs. Report both the cohort restriction and threshold consequence. Three-study results cannot establish universal OIPN reproducibility.

OIPN source files exclude 3213 IDs under the complete-case intersection, including 66 IDs with an observed count >=100 in at least one available file. NC has much larger effects for some genes; the original estimates match the plotted evidence, and display saturation is not proof of a technical error. CCI's eligible ranking uses explicit source-aware MSTRG annotation, not Ensembl features alone. Different rank universes/representative rules, tissue levels, time points, models and sample sizes limit calibrated cross-study magnitude comparisons. Bulk tissue cannot separate within-cell regulation from composition changes.

## Reproduction gaps found in historical OIPN scripts

These are static findings; no R execution was available to demonstrate runtime behavior and no source scripts were changed during this review.

- scripts/03_GSE160543_pathway_analysis.R selects/plots gs_name after fgsea, although its result key is pathway. It also passes the list-valued leadingEdge column straight to write.csv rather than serializing it, and has no local random seed.
- OIPN downstream scripts use lowercase pathway_analysis while archived canonical outputs use Pathway_analysis. Windows may resolve this case difference; case-sensitive systems will not.
- The OIPN GSVA script reads vst_expression_matrix.csv, which is absent from the uploaded tree. The primary DE script contains its export, so a full primary rerun can regenerate it; a downstream-only rerun from the current checkout cannot.
- Gene sets are obtained from live msigdbr calls; exact gene membership/version manifests and a dependency lock are not archived. Session files show the same msigdbr package version but do not by themselves lock the retrieved gene-set release.
- Older study-design and specificity reports describe earlier goals/mappings. Follow project-goal.md and the canonical three-study protocol; do not mix historical OIPN-specific labels with current exploratory divergence criteria.

These gaps limit end-to-end reproducibility of the historical upstream workflow. They do not invalidate the checked source-to-current-table joins, but must be repaired before claiming a fully reproducible pipeline.

## Recommended next work, in order

1. Repair historical OIPN upstream scripts and paths; lock exact Hallmark memberships, rat mappings, software versions and random seeds. Keep archived results intact and compare any rerun transparently. Define a harmonized sensitivity run rather than silently replacing historical estimates or choosing mappings to improve significance.
2. Complete the shared-response analysis by extracting original leading edges for both negative pathways. Preserve negative NES handling, missingness and original rank representatives. Evaluate concordant negative gene log2FC and separately report strict all-three versus OIPN-plus-physical original gene FDR support. Do not assume the number of genes or direction before examining the data.
3. Examine sample-level counts for the heme/erythroid, muscle/remodeling and immune signals, using documented reference markers as exploratory context. Avoid inferring actual cell fractions or contamination without suitable evidence.
4. Then use separately defined shared and divergent candidate lists for STRING/PPI, retain rat identifiers and document evidence channels/background. With small gene lists, sparse networks or degree are descriptive; a hub is not a causal target.
5. Build a compact manuscript evidence table: program, GSEA direction/FDR, GSVA direction/FDR, leading-edge support, candidate gene evidence and caveats. Put the shared injury programs in the main narrative, divergent pathways as exploratory observations, and the four-study comparison as a sensitivity result. Existing bulk-tissue qPCR validates measured transcripts in that experiment, not all GEO signatures or neuronal mechanism.

The immediate scientific gap is the two negative leading-edge programs; the immediate engineering gap is historical OIPN reproducibility. More datasets or more network methods should follow these focused steps rather than obscure them.
