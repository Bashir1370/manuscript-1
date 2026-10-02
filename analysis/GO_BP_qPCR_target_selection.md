# Selected qPCR targets for the experimental rat OIPN comparison

Initial proposal: 2026-10-01. Final six-target selection: 2026-10-02. Branch: `three-dataset-gobp-exploration`.

## Selected target panel

The final primary target panel is **Cdkn1a, Cdk1, Csf1, Hmgcs1, Cav1 and Col4a2**.
**Cdh5** remains an unselected computational candidate. Reference genes used for normalization
are additional to this target panel and require assessment of stability in the
experimental samples.

This panel connects the comparative GO:BP findings to the actual experimental
contrast: rat OIPN versus control in bulk DRG tissue. It samples shared-positive,
shared-negative and selected opposite-direction gene evidence. It is a
purposeful follow-up panel, not the result of a new statistical ranking,
a network-hub selection, or an independent gene-set analysis. Assay specificity
and amplification performance remain to be verified before final primer use.
No experimental qPCR result is asserted in this report.

## Quantitative evidence and selection rationale

Fold changes below are calculated as 2 raised to the original OIPN log2FC.
They describe the public RNA-seq cohort and are hypotheses for the experimental
comparison, not expected guaranteed qPCR outcomes. FDR retains the original
whole-study gene adjustment.

| Target | Evidence category and pathway link | OIPN log2FC | OIPN fold change versus control | OIPN gene FDR | Reason for inclusion |
|---|---|---:|---:|---:|---|
| Cdkn1a | Shared positive; response to wounding and signaling in response to DNA damage | 2.366 | 5.15 | 1.35e-34 | Significant positive expression and same-pathway leading-edge support in all three studies; positive direction stable after every sample omission |
| Cdk1 | Shared positive; the same two representative terms | 1.121 | 2.17 | 7.34e-5 | Significant positive expression and paired leading-edge support in OIPN and NC; matching positive CCI direction; stable descriptive direction in all three studies |
| Csf1 | Shared positive; cytokine-mediated signaling and leukocyte migration | 0.503 | 1.42 | 2.66e-3 | Significant increase and same-pathway leading-edge membership in all three studies; stable positive descriptive direction after every omission; adds direct coverage of shared immune evidence |
| Hmgcs1 | Shared negative; sterol biosynthesis | -0.416 | 0.75 | 1.97e-2 | Adds the shared-negative component; significant reduction in OIPN and NC, leading-edge membership in all three and stable negative descriptive direction in all three |
| Cav1 | Opposite-direction category; selected WNT signaling leading edge | 0.407 | 1.33 | 2.53e-2 | Significant, directionally stable increase in OIPN and decrease in NC; adds evidence from a different selected term than Col4a2 |
| Col4a2 | Opposite-direction category; collagen fibril organization | 0.332 | 1.26 | 2.21e-2 | Significant, directionally stable increase in OIPN and decrease in NC; the only one of the five opposite-category genes with stable negative descriptive direction in CCI |
| Cdh5, not selected | Opposite-direction category; tight junction organization | 0.529 | 1.44 | 3.25e-2 | Alternative candidate for tight junction annotation; significant and directionally stable OIPN increase and NC decrease |

Cdkn1a has stronger three-study gene-level recurrence than Cdk1. Cdk1 has CCI
gene FDR 0.894 and is absent from the CCI leading edges of the two representative
terms supporting its selection. Hmgcs1 has CCI gene FDR 0.545 despite a
directionally stable reduction and leading-edge membership. Directional
recurrence must therefore be distinguished from significance in every study.

All five original opposite-category genes, including Cav1, Col4a2 and Cdh5,
have nonsignificant CCI gene changes. Col4a2's CCI stability does not make its
CCI change statistically significant. Cav1 and Col4a2 have relatively small
OIPN effect sizes; their expression assessment should be interpreted with
effect estimates and uncertainty rather than an expectation of large changes.

## Why this panel differs from the highlighted candidate list

The manuscript highlights Cdk1 and Cdkn1a and introduces the five selected
opposite-category genes Cav1, Cdh5, Col4a2, Tns2 and Wnt6. That seven-gene list
describes computational candidates, not a requirement to assay every candidate.

Hmgcs1 already belongs to the locked 37 shared genes. Adding it to the qPCR
panel covers negative enrichment and avoids a follow-up panel containing only
genes increased in the public OIPN cohort. This does not change the 37-gene,
five-gene or 11/9-pathway computational selections.

Hmgcs1 also has relevant external biological context: Wang and colleagues
reported HMGCS1 localization in satellite glial cells of rat DRG and reduced
expression following spinal nerve ligation. This supports its relevance as a
tissue-response candidate; it does not establish the same response in OIPN or
cellular origin in the present bulk samples.

The earlier five-target proposal was expanded by adding Csf1 to cover the shared
immune response, supported by significant expression and leading-edge evidence
in all three studies. The final panel has six targets. Cdh5, Wnt6 and Tns2 remain
documented computational candidates outside this primary panel.
Their omission is a follow-up prioritization choice, not a claim of biological
irrelevance or failure of the selection rule. No new significance threshold or
automatic composite gene score was introduced.

## Network context and decision chronology

The rat STRING network at score >= 0.4 contains 55 edges among 30 of the 37
shared genes, with seven isolates retained in the complete representation.
Cdk1, Cdkn1a and Csf1 have degree 9, 6 and 4, respectively; Hmgcs1 has degree 3
in the separate four-gene sterol-associated component. The Cdk1-Cdkn1a edge
score is 0.853. Threshold filtering at 0.7/0.9 retains 21/12 edges; centrality
is threshold-dependent. Degree is contextual evidence, not a gene-selection
cutoff or a causal-importance score. Cav1 and Col4a2 are opposite-category
candidates, so their absence from this shared-only network is expected.

The initial five targets were proposed before the new STRING analysis. Csf1
was subsequently added to improve pathway coverage, rather than to optimize
a network rank or any experimental result. No qPCR outcome was used in this
selection. The network does not retroactively turn the panel into a set of hubs.

## Experimental and manuscript interpretation

The experimental analysis assesses the OIPN expression component of the
comparative framework. Results should report the actual change, uncertainty
and statistical analysis for each measured target, including discordant or
nonsignificant findings. No target should be removed solely because its qPCR
result differs from the RNA-seq prediction.

Concordant qPCR results would support tissue-level expression changes in the
experimental rat OIPN setting. OIPN-versus-control qPCR does not independently
confirm the NC/CCI directions, demonstrate OIPN specificity, assign expression
to sensory neurons, or establish a cell-cycle, senescence or causal mechanism.
Only targets actually measured are described as experimentally assessed.

Reference-gene suitability, primer specificity, amplification performance and
the experimental sample design must be documented with the eventual qPCR
methods and results. RNA-seq normalized counts are not a direct prediction of
qPCR Cq values. This report records target selection, not a validated primer
set or a completed qPCR analysis.

## Evidence sources

- [Shared gene evidence](../results/GO_BP_three_dataset/representative_11/displayed_gene_evidence.csv).
- [Shared pathway-gene memberships](../results/GO_BP_three_dataset/representative_11/selected_gene_membership_evidence.csv).
- [Five opposite-category genes](../results/GO_BP_three_dataset/divergent_representative_9/priority_final_genes.csv).
- [Opposite pathway-gene audit](../results/GO_BP_three_dataset/divergent_representative_9/priority_membership_audit.csv).
- [Selected gene sensitivity](../results/GO_BP_three_dataset/selected_evidence_42_genes_20_pathways/gene_sensitivity_summary.csv).
- [Sample expression](../results/GO_BP_three_dataset/selected_evidence_42_genes_20_pathways/selected_gene_sample_expression.csv).
- [Six-target evidence table](../results/String/qpcr_selected_panel_6.csv).
- [STRING network audit and provenance](GO_BP_STRING_network_and_qPCR_closure.md).
- [Manuscript Results, Methods and Discussion](GO_BP_manuscript_results_methods_legends.md).
- [Study design](study-design-lock.md).

External source: Wang F et al. HMG-CoA synthase isoenzymes 1 and 2 localize to
satellite glial cells in dorsal root ganglia and are differentially regulated
by peripheral nerve injury. Brain Research. 2016;1652:62-70.
[PMID 27671501](https://pubmed.ncbi.nlm.nih.gov/27671501/);
[full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC5441544/).

The numerical source tables and sample-sensitivity results were archived before
target selection. Existing differential-expression estimates, pathway selections
and sensitivity analyses are retained. The new network audit and publication
figure provide context for the final six-target follow-up decision.

