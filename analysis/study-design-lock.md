# Study design lock for shared DRG injury responses

Status: current three-study GO:BP manuscript revision, following the local results archived at `b08296108e276225ce168f36c4bd53ba903f86e1`.

## Central question

Does the DRG show a common transcriptional injury response across neuropathy conditions, and where does OIPN fit within that spectrum?

## Included contrasts

| Study | Neuropathy minus control | Samples |
|---|---|---:|
| OIPN GSE160543 | Oxaliplatin minus Vehicle | 4 + 4 |
| NC GSE246156 | Compression minus Sham, L5 day 7 | 3 + 3 |
| CCI GSE212311 | CCI minus Sham, ipsilateral L4-L6 day 11 | 3 + 3 |

Studies are analyzed independently. Original DE results, gene representatives and signed ranks are retained; samples are not pooled across studies. The earlier post hoc exclusion of GSE126773 after inspecting discordant profiles must be reported. No technical failure or mislabeling has been established. The [four-study archive](https://github.com/Bashir1370/manuscript-1/tree/CCI-GSE212311-analysis) retains that evidence.

## Interpretation and manuscript priorities

The unit of interpretation is bulk DRG transcriptional remodeling. Shared enrichment supports recurrence within the three included cohorts. Opposite directions are descriptive, and the selected five-gene evidence is strongest in OIPN versus NC; CCI does not provide significant replication of those genes. Model differences remain confounded with study, time and sampling design. GSEA and GSVA on the same samples are complementary analyses, not independent validation.

No inference of neuronal cell-cycle re-entry, neuronal senescence, causality or therapeutic targets follows from bulk expression alone. The experimental arm comprises rat OIPN versus control in DRG tissue. Its role is to assess OIPN tissue-level expression changes for transcripts actually measured. It does not experimentally compare neuropathy models or establish neuronal origin or causality. The seven highlighted computational candidates are not assumed to have all been measured.

## Completed analysis and current next step

Per-study DE, Hallmark comparison, locked GO:BP enrichment, GSVA, representative pathway/gene selection and descriptive sample sensitivity are archived. The selected evidence comprises 37 shared and five opposite-direction genes, and 11 shared and nine opposite-direction pathways. No additional selection or sample exclusion is proposed.

Integrate the [GO:BP Results, Methods and figure legends](GO_BP_manuscript_results_methods_legends.md) into the full manuscript, assign final figure numbers, and align the Discussion and reviewer responses with these evidential limits. The computational section is ready for integration; the whole manuscript and experimental sections still require their own source files and review.

See [GO_BP_three_dataset_protocol.md](GO_BP_three_dataset_protocol.md), [three_dataset_protocol.md](three_dataset_protocol.md) and [GO_BP_selected_evidence_and_sensitivity.md](GO_BP_selected_evidence_and_sensitivity.md) for the full analysis rules and provenance.
