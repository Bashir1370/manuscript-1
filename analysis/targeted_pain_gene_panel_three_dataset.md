# Requested 30-gene pain panel: archived differential expression

Source branch: three-dataset-neuropathy-analysis.
Source commit: 820bfe2733ab8c545289ab83cc5350f665290260. Date: 2026-10-01.

This is a targeted extraction of existing gene-level DESeq2 results, independent
of Hallmark or leading-edge membership. No DE, GSEA or GSVA model was rerun.
Positive log2FC means neuropathy greater than the study's own control:
OIPN GSE160543 Oxaliplatin vs Vehicle; NC GSE246156 L5 day7 Compression vs Sham;
CCI GSE212311 L4-L6 day11 CCI vs Sham.

Use the original genome-wide gene padj < 0.05, without recomputing FDR on
the requested panel. Match case-insensitively to the archived rat symbol.
OIPN uses its original DE symbol, NC joins the unambiguous archived Ensembl
mapping to original DE, and CCI uses its source-aware mapping joined to full DE.
No alias rescue, new annotation database or selection of the lowest P value was
used. All available panel symbols had one representative per study.

| Study | Significant | Not significant | FDR unavailable | Not recoverable from archived mapping |
|---|---:|---:|---:|---:|
| OIPN | 2 | 26 | 2 | 0 |
| NC | 11 | 7 | 0 | 12 |
| CCI | 3 | 25 | 0 | 2 |

OIPN significant genes: Scn1a down and Cxcl12 up. Tnf and Tac4 have
available original effects but NA padj; Tnf raw P=0.031 does not establish
genome-wide FDR significance. CCI significant genes: Scn1a and Kcna1 down,
Trpa1 up. NC significant genes: Bace1, Cxcl12, Ephb1, Mmp24, Trpa1, Trpv1
down; Comt, Disc1, Htr2a, Ntrk1, Tmem120a up.

Scn1a is significantly reduced in OIPN and CCI; its NC result cannot be
assessed from the available symbol mapping. Cxcl12 is significant with opposing
directions in OIPN and NC. Trpa1 is significant with opposing directions in NC
and CCI; OIPN does not reach FDR significance. No gene is significant in all
three among the 16 genes with evaluable FDR in all three. Missing evidence limits
that conclusion for other panel genes. These are within-study results, not a
formal test of differences between models or evidence of functional channel activity.

NC missing symbols: Adora1, Asic3, Cartpt, Kcna1, Ntsr1, Phf24, Prdm12,
Scn11a, Scn1a, Tac1, Tac4, Tnf. CCI missing symbols: Ntrk1, Trpv1.
The latter are also absent under the NC Ensembl IDs from full CCI DE.
Not recoverable means unavailable under the archived mapping/filter policy;
it is not evidence of zero expression, absence from DRG, or nonsignificance.
Further source annotation audit is needed to resolve these entries. Do not fill
them with zero or use annotation from another feature without verification.

Results: results/targeted_pain_gene_panel_three_dataset/gene_panel_comparison.csv
and gene_panel_long.csv. Long output retains original feature IDs, baseMean,
standard errors, Wald statistics, raw P, FDR and missing-evidence status.
source_git_blobs.json records original source identities. The four preexisting
local tables matched current Git blob SHA1. CCI values were extracted from the
current full DE blob and checked against its source-aware audit.
CCI_selected_original_DE.csv retains the exact 28 matched DE rows to reproduce
the hosted audit when the full CCI table is not downloaded. Runtime checksums
explicitly identify this subset, not the full CCI file.

Reproduce from a full checkout (Python standard library only):

```powershell
python scripts/cross_model/11_targeted_pain_gene_panel.py
```

The script was executed successfully and checked 90 gene-study entries and
CCI mapping agreement. It stops on duplicate mapped features rather than
choosing a favorable effect or P value. Existing analyses are read only.
