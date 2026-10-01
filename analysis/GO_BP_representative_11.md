# Eleven shared GO:BP representatives with gene evidence

Run this optional presentation stage after stages 02 and 04, from repository root:

```r
Sys.unsetenv("GOBP_TABLES_ONLY")
source("scripts/gobp/06_representative_GO_BP_heatmap.R")
```

The stage reads existing results and does not refit DE, GSEA, or GSVA. It is
separate from the full runner because the 11 terms are an author-selected
exploratory display, rather than an automatic part of every future experiment.
It stops if any term no longer has the declared sign and pathway FDR <0.05 in
all three studies. It also checks gene priority flags against the original
leading-edge and gene-level statistics before writing outputs.

Outputs are under `results/GO_BP_three_dataset/representative_11/`:

- `Fig_GO_BP_11_pathways_and_genes.pdf`: vector figure, two panels.
- `Fig_GO_BP_11_pathways_and_genes.png`: 300 dpi figure.
- `selected_pathway_evidence.csv`: exact ordered terms, selection rationale,
  NES, original full-family pathway FDR, and existing classification.
- `selected_gene_membership_evidence.csv`: union leading-edge evidence for the
  selected terms, including genes that do not pass priority filtering.
- `displayed_gene_evidence.csv`: one row per displayed gene, original log2FC
  and gene FDR, display role, and associated priority/shared LE pathways.
- `pathway_gene_counts.csv`, `plot_settings.csv`, `input_checksums.csv`,
  `analysis_notes.txt`, and `R_sessionInfo.txt`.

Panel A uses raw NES with one symmetric scale across all studies, without
row scaling. Six shared-positive terms are followed by five shared-negative
terms. Stars denote the original pathway FDR, never a selected-term BH rerun.

Panel B shows all unique priority genes supported within these 11 terms,
with no manually added context genes. The revised rule requires matching
log2FC direction in all three, plus BOTH same-path LE membership and original
gene FDR <0.05 in OIPN and at least one same physical study (NC or CCI).
The remaining physical study need not show gene significance or LE membership,
but its log2FC must be available and have the matching sign. Gene stars denote
original gene FDR. Cdk1 and Cdkn1a are included only when they pass these rules.
Original stage-04 flags remain intact for auditing; selected_priority records
the new rule. This stage does not change stage-04 full-result priorities.
Missing log2FC is gray. Colors use a separate symmetric log2FC scale, default +/-6; CSV values
are never clipped. To change only the gene display range:

```r
Sys.setenv(LE_COLOR_LIMIT = "4")
source("scripts/gobp/06_representative_GO_BP_heatmap.R")
```

The fixed selection is response to wounding, cytokine-mediated signaling,
leukocyte migration, phagocytosis, signaling in response to DNA damage,
tissue remodeling, oxidative phosphorylation, tricarboxylic acid cycle,
regulation of trans-synaptic signaling, potassium ion transport, and sterol
biosynthesis. Immune proliferation terms are significant shared-positive
results but were not part of this original 11-term selection; they remain in
the complete results. A later figure selection should explicitly document
any addition or replacement rather than changing this one silently.

The selection is based on the author's question and inspection of recurrence,
gene evidence and overlap. It is not an automated redundancy clustering, a
confirmatory independent gene-set test, or 11 independent mechanisms.
Same-path all-three LE Jaccard examples in the author's uploaded results:
response to wounding/wound healing: 11/12; electron transport chain/oxidative
phosphorylation: 15/17; nucleoside phosphate/purine-containing compound
biosynthesis: 24/25. These compare common leading edges, not whole GO sets.

The current supplied lock manifest records MSigDB 2026.1.Hs, msigdbr 26.1.1,
human C5:GO:BP mapped to Rattus norvegicus, and size limits 15-500. Its
checksums need the actual locked files to be independently verified.

Validation before publication: parsed in R 4.3.3, table stage executed on the
author's uploaded results, the original and revised priority rules checked, and deliberately
invalid pathway significance/gene flags rejected. The plotting code uses
base R graphics and requires no additional package. Research figures are
generated only by the author's local R execution; no plot was generated here.

Revision validation on the supplied results: 37 unique selected genes,
including Cdk1 and Cdkn1a without context labels. All previous 24 selected
priorities are retained. Tests reject absent OIPN LE, opposite third-study
log2FC, and LE/significance support split across different physical models.
The same 11 pathways and original pathway/gene statistics are retained.
