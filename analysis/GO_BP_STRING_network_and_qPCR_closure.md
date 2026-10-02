# Systems-biology closure and qPCR transition

Date: 2026-10-02. Branch: `three-dataset-gobp-exploration`.
Reviewed source commit: `fa278673ab02ec4eb023326352d802e80ae14783`.

## Evidence audit

The uploaded STRING TSV, Cytoscape node table and original network PNG were reviewed against the locked 37 shared genes and the existing GO:BP evidence. The TSV contains 55 unique undirected associations, with no duplicate pairs or self-loops; all STRING identifiers are rat taxonomy 10116 and all endpoints belong to the locked query. The original Cytoscape table contains 30 connected genes. Every archived degree agrees with independent recomputation. Seven query genes absent from that table are retained as isolates: Amigo1, Ephb2, Fxyd7, Kcna2, Myof, Nefh and Syt11.

The primary display uses combined score >= 0.4, the full association network and no added partners. This is a medium-confidence threshold; it must not be described as STRING's high-confidence threshold. Association confidence is not biochemical interaction strength. Functional associations and transferred evidence do not establish direct interactions in DRG.

| Minimum score | Query genes retained | Edges | Genes with an edge | Isolates |
|---:|---:|---:|---:|---:|
| 0.4 | 37 | 55 | 30 | 7 |
| 0.7 | 37 | 21 | 21 | 16 |
| 0.9 | 37 | 12 | 14 | 23 |

At 0.4, the components contain 26 genes, four sterol-associated genes and seven singletons. Cdk1 has degree 9, Icam1 8, Foxm1 7, Cdkn1a 6, Csf1 4 and Hmgcs1 3. Cdk1 degree falls to 5 and 2 at 0.7 and 0.9; Cdkn1a falls to 1 and 0. The separate sterol component contains Hmgcs1, Dhcr7, Insig1 and Sqle. Connectivity is descriptive context, not an inferential hub test or independent validation. No interaction-enrichment P value was supplied or reconstructed.

## Final primary qPCR panel

| Target | Evidence covered | Original OIPN RNA-seq fold change |
|---|---|---:|
| Cdkn1a | Shared-positive injury response; strongest three-study recurrence of the two highlighted injury candidates | 5.15 |
| Cdk1 | Shared-positive injury response with OIPN and NC statistical support | 2.17 |
| Csf1 | Shared immune evidence with significant increase and leading-edge support in all three studies | 1.42 |
| Hmgcs1 | Shared-negative sterol biosynthesis evidence | 0.75 |
| Cav1 | Selected WNT annotation with opposite OIPN-NC gene direction | 1.33 |
| Col4a2 | Collagen annotation with opposite OIPN-NC gene direction | 1.26 |

All six satisfy original OIPN gene FDR < 0.05 and preserve their OIPN and NC descriptive direction after every sample omission. Cdkn1a and Csf1 are significant in all three studies. Cdk1 and Hmgcs1 have nonsignificant CCI changes; Cav1 and Col4a2 also have nonsignificant CCI changes, and Cav1 is direction-sensitive in CCI. These distinctions are retained in the manuscript and machine-readable panel table. Cav1 and Col4a2 do not belong to the shared-network query.

The earlier five-target proposal was expanded by Csf1 to cover shared immune evidence, before assessment of experimental qPCR outcomes. Selection covers evidence axes rather than maximizing degree or exhausting the selected genes. Cdh5, Wnt6 and Tns2 remain documented computational candidates. The primary panel contains six targets; normalization genes are additional and require validation in the actual samples. No laboratory constraints are used as an explanation in the manuscript.

The next section measures these transcripts in rat OIPN versus control DRG. It can assess the OIPN expression hypotheses, including negative or discordant results. It cannot establish NC/CCI responses, oxaliplatin specificity, cellular origin or causality. The RNA-seq fold changes above are not promised qPCR effect sizes.

## Manuscript and figures

The [canonical manuscript](GO_BP_manuscript_results_methods_legends.md) now contains the leading-edge gene discussion after the GSEA-GSVA comparison, followed by sample sensitivity, STRING context and the final qPCR transition. Materials and methods cover network processing, threshold sensitivity, figure rendering and purposeful target selection. Discussion retains the limitations of independent public cohorts, shared input samples and bulk tissue measurements.

Working figure labels are Figure 1 (11 shared representatives and 37 genes), Figure 2 (nine opposite representatives and five genes), Figure 3 (shared-gene network), Supplementary Figures S1-S2 (all significant shared terms), S3 (42-gene sample expression) and S4 (GSEA-GSVA directions). Full captions and exact file mappings are in the manuscript. Global numbering may be remapped when this section is integrated with the complete paper.

The new network is available as PDF, SVG and 600-dpi PNG under `results/String/Fig_shared_37_STRING_network`. It was redrawn from the original edge list using coordinates transcribed from the author's figure. It retains 37 labels and 55 edges, rather than enlarging the original raster. Node colors denote components/isolates, node size is constant, and edge width represents confidence. PDF and SVG preserve vector output. The original uploads are retained unchanged. Existing GO figures, differential-expression estimates, representative selections and supplementary Word files are retained. No complete opposite-pathway list is added to the Word supplement.

## Reproduction and checks

From the repository root, run:

```bash
python scripts/gobp/09_shared_STRING_network.py
```

Python and matplotlib are required. The script validates the query, pairs, scores, species, archived degrees and 37 layout coordinates, then writes:

- `results/String/network_node_audit_37.csv`: 111 rows, one per gene and threshold.
- `results/String/network_threshold_sensitivity.csv`: three thresholds with all query genes retained.
- `results/String/qpcr_selected_panel_6.csv`: six targets joined to original expression statistics and descriptive sensitivity.
- `results/String/network_provenance.json`: source commit, input SHA256 checksums, rendering versions and missing-metadata record.
- The PDF, SVG and PNG network figures.

The script was executed successfully. Source-file bytes were checked against the uploaded Git blob hashes before generating provenance. All 30 original degrees matched; the six-target statistics agree with the existing source tables. The network preview was inspected for label readability and clipping. This stage does not refit expression or pathway models.

## Remaining publication metadata

The actual STRING release and Cytoscape version are absent from the supplied export and must be recorded before submission. The author-described query settings are distinguished from values verified directly in the edge table; the original settings screenshot is not part of the three newly uploaded files. No Cytoscape session file or original enrichment output was supplied. No enrichment claim relies on those missing files.

The computational section and target prioritization are complete. Actual animal numbers, tissue collection time, primers, efficiency, reference-gene validation, qPCR analysis and measured outcomes belong to the experimental section and cannot be filled from the computational archive.

See also the [target-selection report](GO_BP_qPCR_target_selection.md) and [study-design lock](study-design-lock.md).
