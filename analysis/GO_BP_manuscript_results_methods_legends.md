# Shared and divergent DRG transcriptional responses across neuropathy models

## Results

### Shared GO Biological Process enrichment across the three studies

We compared rat dorsal root ganglion (DRG) transcriptomes from oxaliplatin-induced peripheral neuropathy (OIPN; GSE160543), nerve compression (NC; GSE246156) and chronic constriction injury (CCI; GSE212311). Each study was analyzed independently. Of 7,535 locked GO Biological Process (GO:BP) gene sets, 3,406 were available for enrichment comparison in all three studies. We identified 102 terms with positive normalized enrichment scores (NES) and 18 with negative NES at a false discovery rate (FDR) below 0.05 in every study (Supplementary Figures S1 and S2). These recurrent enrichments support shared transcriptional responses within the three included cohorts.

Eleven representative terms were selected for display: six shared-positive terms covering response to wounding, cytokine-mediated signaling, leukocyte migration, phagocytosis, signaling in response to DNA damage and tissue remodeling, and five shared-negative terms covering oxidative phosphorylation, the tricarboxylic acid cycle, regulation of trans-synaptic signaling, potassium ion transport and sterol biosynthesis. These terms summarize overlapping GO annotations and should not be counted as independent mechanisms. Leading-edge membership and original gene-level differential-expression statistics identified 37 genes satisfying the selected evidence rule (Supplementary Table S1). All 37 had concordant log2 fold-change directions across the three studies; selection required same-pathway leading-edge membership and gene FDR below 0.05 in OIPN and at least one physical-injury study, rather than gene significance in all three studies.

### Selected opposite-direction responses have strongest gene support in OIPN and NC

The initial directional comparison identified 457 terms whose OIPN NES opposed both NC and CCI: 246 were positive in OIPN and negative in the two physical-injury studies, and 211 showed the reverse pattern. This classification imposed no significance requirement. Of these terms, 24 had OIPN pathway FDR below 0.05, comprising 20 positive and four negative terms. Nine author-selected representatives were displayed, covering collagen fibril organization, extracellular matrix assembly, WNT signaling, tight junction organization, Schwann cell differentiation, peripheral nervous system axon ensheathment, renal system development, exocytosis and cellular response to unfolded protein. The renal system development label denotes a GO annotation and does not imply a renal phenotype in DRG.

The revised same-pathway leading-edge and gene-significance rule identified five genes within the nine representatives: Cav1, Cdh5, Col4a2, Tns2 and Wnt6. All five were increased in OIPN and decreased in NC, with original gene FDR below 0.05 in both studies. Their CCI estimates were also negative, but small and nonsignificant (gene FDR 0.888-0.988). Thus, the strongest gene-level evidence for this selected divergence was between OIPN and NC. The original 39-gene selection across all 457 opposite-direction terms is retained separately, with the final five genes marked (Supplementary Table S2). Directional divergence across separate studies does not establish OIPN specificity or a formal interaction between injury models.

### GSEA and GSVA agree more consistently for the shared representatives

GSEA and GSVA effect directions agreed for all 11 shared representatives in each study. However, GSVA reached full-family FDR below 0.05 for seven of these terms in NC and for none in OIPN or CCI. The shared GSEA findings therefore had complete directional agreement with GSVA but incomplete statistical agreement. Among the nine opposite-direction representatives, GSEA-GSVA direction agreement was observed for nine terms in OIPN, six in NC and three in CCI. GSVA significance was observed for zero, two and zero terms, respectively. Because both analyses used the same samples, their agreement provides complementary evidence rather than independent validation.

### Leading-edge gene candidates connect shared and divergent pathway findings

Among the 37 shared-category genes, Cdk1 and Cdkn1a were highlighted as candidates for follow-up because both contributed leading-edge evidence to response to wounding and signaling in response to DNA damage in OIPN and NC. Cdkn1a showed positive log2 fold changes of 2.366, 2.356 and 1.253 in OIPN, NC and CCI, respectively, with gene FDR below 0.05 and leading-edge membership in both terms in every study. Cdk1 showed corresponding positive log2 fold changes of 1.121, 3.025 and 0.455. Its changes were significant in OIPN and NC (gene FDR 7.34e-5 and 8.77e-4), whereas its CCI change was nonsignificant (gene FDR 0.894) and it was absent from the CCI leading edges of these two terms. Thus, Cdkn1a provided stronger three-study gene-level recurrence, while Cdk1 met the declared shared-gene selection rule through OIPN and NC support and a matching CCI direction. Both retained positive descriptive effects after every sample omission in all three studies.

The five opposite-category candidates linked the selected enrichment patterns to specific leading-edge genes: Cav1 and Wnt6 to WNT signaling, Cdh5 to tight junction organization, Col4a2 to collagen fibril organization and Tns2 to renal system development. These links met the same-pathway leading-edge and gene-FDR rule in OIPN and NC. All five genes showed significant positive OIPN and negative NC changes and retained those descriptive directions after every omission. Their CCI changes were nonsignificant; only Col4a2 retained a negative direction after every omission. These five genes represent gene-level evidence within the displayed opposite-direction pathways, rather than validated representatives of all nine terms.

These findings nominate Cdk1 and Cdkn1a as shared-response candidates and Cav1, Cdh5, Col4a2, Tns2 and Wnt6 as candidates associated with the selected OIPN-NC divergence. This nomination derives from pathway membership, original differential-expression evidence and descriptive direction stability. The experimental rat OIPN-versus-control comparison addresses the OIPN expression component of this framework.

### Sample omission supports the shared response and qualifies the opposite-direction findings

Descriptive leave-one-sample-out analyses evaluated the 42 selected genes and 20 pathways using fixed normalization and archived GSVA scores (Supplementary Table S3). Among the 37 shared genes, the direction of the difference in mean log2(normalized count + 1) was retained after every omission for 37 genes in OIPN, 36 in NC and 33 in CCI. Tln1 was direction-sensitive in NC; Cd14, Parp14, Dhcr7 and Syt11 were direction-sensitive in CCI. Cdk1 and Cdkn1a retained positive differences after every omission in all three studies.

All five opposite-direction genes retained their descriptive directions in OIPN and NC. In CCI, Cav1, Cdh5, Tns2 and Wnt6 could change sign after sample omission; only Col4a2 consistently remained negative. The shared pathway GSVA differences retained their baseline directions for 11 of 11 terms in OIPN, 10 of 11 in NC and 11 of 11 in CCI; tissue remodeling was direction-sensitive in NC. The corresponding counts for opposite-direction terms were nine of nine, four of nine and one of nine. The sole stable opposite-category GSVA term in CCI was Schwann cell differentiation, whose positive GSVA difference disagreed with its negative GSEA NES. Stability therefore refers to each method's descriptive baseline, not necessarily agreement with GSEA.

Together, these results support recurrent bulk DRG transcriptional responses across the included neuropathy cohorts, while the selected differences involving OIPN are less consistently supported in CCI. Direction stability does not demonstrate renewed significance after sample omission, cellular origin or a causal mechanism.

## Methods

### Study contrasts and analysis scope

The study combined comparative analysis of public bulk DRG transcriptomes with an experimental rat OIPN-versus-control component. The computational comparison addressed recurrence and directional differences across neuropathy models. The experimental component addressed DRG responses within OIPN; it was not designed as an experimental comparison of OIPN, NC and CCI. Experimental procedures and measured outcomes are reported in the corresponding experimental sections.

The analysis included three independent rat bulk DRG RNA-sequencing studies: OIPN GSE160543 (four oxaliplatin and four vehicle samples), NC GSE246156 (three compression and three sham samples, L5 DRG at day 7) and CCI GSE212311 (three CCI and three sham samples, ipsilateral L4-L6 DRG at day 11). Contrasts were neuropathy minus control. Original differential-expression results, signed gene ranks and feature representatives were retained. Samples were not pooled across studies, and no cross-study batch correction was performed. The restriction to these three cohorts followed inspection of the earlier four-study analysis and exclusion of GSE126773 because of divergent profiles. This was a post hoc scope restriction, without established technical failure or sample mislabeling. Accordingly, recurrence is interpreted within the included studies rather than as replication across all OIPN datasets.

### Locked GO gene sets and enrichment analysis

Human MSigDB 2026.1.Hs C5:GO:BP gene sets were obtained through msigdbr 26.1.1 and mapped to Rattus norvegicus orthologs. The 7,535-term collection, distinct rat memberships, mapping audit, pathway metadata and checksums were locked before GO scoring. Each study's existing ranked signed statistics were analyzed with fgseaMultilevel using two-sided enrichment, gseaParam = 1, eps = 1e-10, nPermSimple = 10,000, sampleSize = 101 and one worker. Study-specific seeds were 160543, 246156 and 212311. Eligible gene sets contained 15-500 genes overlapping the study's rank. Benjamini-Hochberg adjustment covered the full eligible GO:BP family separately within each study; FDR was not recalculated over displayed terms. Coverage and unavailable terms were retained in the archive.

Shared-positive and shared-negative terms required the same strict NES sign and pathway FDR below 0.05 in all three studies. Opposite-direction terms required the OIPN NES sign to oppose both NC and CCI, without an initial significance filter. OIPN pathway FDR below 0.05 was subsequently used to define the 24-term display candidate pool. Eleven shared and nine opposite-direction representatives were selected by the authors to summarize relevant annotations and overlapping gene evidence. This was an exploratory presentation choice, not automated redundancy clustering or a new inferential test.

### Selected gene evidence

For the shared representatives, selected genes required available, nonzero log2 fold changes with the same sign in all three studies and a direction consistent with the pathway. In addition, the gene had to belong to the same pathway's leading edge and have original whole-study gene FDR below 0.05 in OIPN and at least one of NC or CCI. Leading-edge membership and significance had to coincide in the same supporting physical-injury study. The other physical-injury study was required to retain the common expression direction, but not gene significance or leading-edge membership. This yielded 37 unique genes.

For the nine opposite-direction representatives, the gene log2 fold-change sign had to match its pathway NES sign in each study, with OIPN opposing both physical-injury studies. The same leading-edge and original gene-FDR requirement applied in OIPN and one physical-injury study, yielding five unique genes. The prior 39-gene list used the union of leading edges across all 457 opposite-direction terms and a less restrictive membership rule. It was retained as an archive and distinguished from the revised five-gene selection. No FDR adjustment was performed over either selected gene list.

### GSVA and sample-level expression

GSVA used log2(normalized count + 1) expression with the archived GSEA feature representatives. OIPN and CCI used archived normalized counts. NC size factors were estimated with DESeq2 from the original raw counts after retaining features with at least 10 counts in at least three samples; differential expression was not refitted. Constant expression rows were removed before assessing size eligibility. GSVA used Gaussian scoring, tau = 1, maxDiff = TRUE, absRanking = FALSE and serial execution. Scores were analyzed separately in each study using limma, with control as the reference and eBayes(trend = FALSE, robust = FALSE). Original score differences, confidence intervals, P values and Benjamini-Hochberg FDR across all scored GO:BP sets were retained.

The GSEA-GSVA comparison used effect directions and original full-family FDR. NES and GSVA score differences were not treated as interchangeable effect magnitudes, and absolute GSVA magnitudes were not calibrated across studies. Sample expression was displayed as a gene-specific z-score of log2(normalized count + 1) within each study. Display standardization did not alter the expression values used for descriptive sensitivity analyses.

### Descriptive leave-one-sample-out analysis

Each sample was omitted once within its own study, retaining all remaining control and neuropathy samples. For genes, the descriptive effect was the neuropathy minus control difference in group means of log2(normalized count + 1). For pathways, it was the corresponding difference in archived GSVA scores. Normalization and GSVA scoring remained fixed; neither DESeq2 nor limma was refitted, and no new P values or FDR estimates were produced. Direction was classified as stable only when every omission retained the strict sign of the descriptive full-sample baseline. A zero baseline was not stable. Agreement with the original DE or GSVA estimate was recorded separately. All 840 gene and 400 pathway omission contrasts were checked by independent recomputation.

## Discussion

### Shared DRG responses and the position of OIPN

The comparative analysis places OIPN within a recurrent bulk DRG injury-response pattern encompassing immune and wound-associated enrichment together with negative enrichment of selected metabolic and neuronal-function annotations. The shared representative terms showed concordant GSEA and GSVA directions, and most selected gene and pathway effects retained their direction after single-sample omission. Cdkn1a provided significant gene-level recurrence across all three studies, whereas Cdk1 showed stronger evidence in OIPN and NC. This distinction supports interpretation at the level of tissue transcriptional remodeling without assigning the changes to a particular DRG cell population.

The opposite-direction analysis identified candidate differences, particularly between OIPN and NC. The five selected genes had significant and directionally stable changes in those two studies, while their CCI evidence was weak. These results generate hypotheses about model-associated responses rather than establishing oxaliplatin specificity. OIPN-versus-control measurements in the experimental DRG samples address the OIPN component of those hypotheses. Agreement for a measured transcript would support its tissue-level expression change in that experimental setting; it would not independently verify the opposite direction in physical-injury models or validate every member of the computational candidate list.

### Scope of inference

The public studies differ in injury model, sampling time and design, which limits attribution of between-study differences to injury type alone. The three-study restriction followed inspection of the earlier four-study results, so the recurrent findings apply to the included cohorts. GSEA and GSVA use the same samples, and the sample-omission analysis assesses descriptive direction stability with normalization and scoring fixed. The study therefore supports candidate transcriptional responses in bulk DRG; cell-specific mechanisms and causal involvement require evidence beyond expression measurements.

## Figure legends

### Shared pathway and gene evidence

**Eleven representative shared GO:BP pathways and their selected gene evidence.** Panel A shows original GSEA NES for six shared-positive and five shared-negative terms in OIPN, NC and CCI. NES uses a common symmetric scale without row standardization. All terms have pathway FDR below 0.05 in every study. Panel B shows original log2 fold changes for 37 unique genes meeting the same-pathway leading-edge and gene-significance rule in OIPN and at least one physical-injury study. The remaining study is required to have a matching expression direction, but not gene significance. Asterisks denote original within-study FDR below 0.05 for the corresponding pathway or gene. Gene colors use a separate symmetric scale clipped at +/-6 for display; raw values are retained in the tables. The terms are author-selected representatives of overlapping annotations.

### Opposite pathway and gene evidence

**Nine representative GO:BP pathways with opposite enrichment directions in OIPN and physical-injury studies.** Panel A shows seven OIPN-positive and two OIPN-negative terms selected from 24 opposite-direction terms with OIPN pathway FDR below 0.05. NC and CCI significance was not required for pathway selection. Panel B shows Cav1, Cdh5, Col4a2, Tns2 and Wnt6, which match the pathway directions and meet the same-pathway leading-edge and gene-FDR rule in OIPN and one physical-injury study. The five genes have significant opposite changes in OIPN and NC; their CCI changes are nonsignificant. Asterisks denote original within-study pathway or gene FDR below 0.05. NES and log2 fold changes use separate symmetric scales; gene colors are clipped at +/-6 for display. Renal system development is a GO annotation label, not evidence of a renal phenotype in DRG. This descriptive comparison does not establish OIPN specificity.

### Selected gene sample expression

**Sample expression of the 37 shared and five opposite-direction genes.** Columns show individual controls followed by neuropathy samples within each study. Cells display gene-specific within-study z-scores of log2(normalized count + 1), with colors clipped at -3 and +3. Horizontal lines separate the shared and opposite categories; vertical lines separate control and neuropathy samples. Scaling is performed independently within each study and does not permit comparison of absolute expression across studies. No statistical significance is encoded by the cell colors. OIPN: four vehicle and four oxaliplatin samples; NC and CCI: three sham and three injury samples each.

### GSEA and GSVA direction comparison

**Effect directions and original FDR for 20 selected GO:BP representatives.** Paired columns show GSEA NES direction and GSVA neuropathy minus control score direction for OIPN, NC and CCI. Red denotes a positive effect and blue a negative effect; color intensity does not encode magnitude. Asterisks indicate original within-study FDR below 0.05, adjusted across the full eligible GSEA or scored GSVA family, respectively. The horizontal line separates 11 shared from nine opposite-direction representatives. Shared terms have complete directional agreement, whereas the opposite-direction terms show weaker agreement in NC and CCI. Both methods use the same samples and do not constitute independent validation.

## Repository evidence and figure mapping

This GO:BP section is supported by the local execution archived at commit `b08296108e276225ce168f36c4bd53ba903f86e1`. Final journal figure numbering can be assigned when this section is integrated with the existing Hallmark, network and experimental sections.

| Legend | Repository figure |
|---|---|
| Shared pathway and gene evidence | `results/GO_BP_three_dataset/representative_11/Fig_GO_BP_11_pathways_and_genes` |
| Opposite pathway and gene evidence | `results/GO_BP_three_dataset/divergent_representative_9/Fig_GO_BP_9_opposite_pathways_and_genes` |
| Selected gene sample expression | `results/GO_BP_three_dataset/selected_evidence_42_genes_20_pathways/Fig_selected_42_gene_sample_expression` |
| GSEA and GSVA direction comparison | `results/GO_BP_three_dataset/selected_evidence_42_genes_20_pathways/Fig_selected_20_GSEA_GSVA_directions` |

The numeric sources are the comparison, representative selection and stage-08 tables under `results/GO_BP_three_dataset/`. Supplementary Tables S1 and S2 retain the 37-gene and archived 39-gene evidence. Table S3 summarizes selected direction stability and GSEA-GSVA concordance. The complete opposite-direction pathway list remains in the repository rather than being appended to the Word supplementary material.
