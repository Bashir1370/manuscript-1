# Completion of shared-negative programs and OIPN reproducibility preparation

Date: 2026-09-30. Original input snapshot: `af4742bec54d1c49f38c09a53f7d563dc72f5ef7`.
The original canonical DE/GSEA/GSVA tables and existing positive/divergent results
are preserved. R execution remains local; this change does not claim new GSEA or
GSVA fits, generated figures, cell fractions, or functional validation.

## Shared-negative leading edges

Use original negative-NES leading edges, without reversing ranks. A selected
pathway has NES < 0 and original within-study FDR < 0.05 in all three studies.
Gene representatives are the exact features used by each original GSEA. Their
signed statistics were checked against every original rank. Missing rank/DE
evidence remains missing, not zero or nonsignificant expression.

| Pathway | Union leading-edge genes | Membership in the same pathway in all three |
|---|---:|---:|
| FATTY_ACID_METABOLISM | 98 | 10 |
| OXIDATIVE_PHOSPHORYLATION | 161 | 34 |

There are 44 all-three pathway-gene memberships representing **39 unique genes**.
All 39 have negative original log2FC in all three studies. Prioritization retains
the original genome-wide gene FDR: FDR <0.05 in OIPN and at least one physical
model selects **7 genes**. The strict all-three significant subset contains
**Hsph1** only. The other six are supported by OIPN+NC.

Priority symbols: Gabarapl1, Hmgcs1, Hsph1, Pdhb, Slc25a4, Uqcrc2, Uqcrfs1.
No effect-size threshold was added. Membership, direction and gene significance
are separate evidence levels; the 39-gene list must not be called 39 replicated
DEGs. These are bulk-tissue expression signatures, not measurements of metabolic
flux or proof of neuronal mitochondrial dysfunction.

Tables are archived under `results/shared_negative_Hallmark_leading_edge_three_dataset/`.
They were reconstructed with Python from original inputs, not produced by an R
run. `validation.json` records the 259 union memberships, 777 study evidence rows,
rank checks and regression reconstruction of the historical positive analysis
(10 pathways, 109 memberships, 78 unique genes). Original Git blob identities are
recorded in `negative_leading_edge_source_blobs.json`; all eleven input byte streams
were verified against those blobs before the archived checksums were calculated.

`05_three_dataset_leading_edge.R` and `06_three_dataset_gene_prioritization.R` now
accept positive/negative direction, retaining positive defaults and its output
schema. `08_three_dataset_negative_leading_edge.R` calls the negative mode in an
isolated environment. It generates two gene heatmaps locally, with the existing
shared linear +/-6 scale and original gene-FDR annotations. Negative membership
tables include the positive flag for compatibility and the new negative flag.

## Reproducibility repairs and limits

The old OIPN GSEA script used `gs_name` instead of fgsea's `pathway`, exported a
list-valued leading edge directly to CSV, had no fixed seed, and wrote to a path
with a case differing from canonical uploaded outputs. These are repaired.
GSEA retains highest signed-statistic symbol representatives, minSize=15,
maxSize=500, 10,000 permutations and now explicitly uses fgseaSimple, seed=160543,
and one process. Leading edges are serialized with semicolons.

Both OIPN downstream scripts acquire and freeze the same current 50-Hallmark rat
membership on their first local run, retain gene-set metadata and checksum, and
archive executed package versions, R session, parameters and input checksums.
This does **not** establish which gene-set release produced the old tables.
Changing a seed, database release or software version can change NES/FDR/edges;
new comparisons make those changes reviewable rather than silently promoting them.

The primary script accepts `OIPN_PRIMARY_OUTDIR`. Script 09 regenerates primary
DE/VST in `reproducibility_rerun/primary` and GSEA/GSVA/GO in
`reproducibility_rerun/Pathway_analysis`. It writes DE, NES/FDR/leading-edge Jaccard,
and sample-score comparison tables against the canonical archive. Environment
overrides are restored on exit. Existing cross-study analyses continue reading
the canonical historical outputs.

GSVA preserves the historical Entrez mapping/first-feature duplicate policy and
archives that mapping. New rerun statistics use the current manuscript limma/BH50
analysis; the older eight-pathway Wilcoxon statistics are not overwritten or treated
as equivalent. GO keeps its historical package-default background. A tested-gene
background and harmonized three-study gene filters/mappings are future sensitivity
analyses, not silently incorporated in this repair. The complete-case OIPN gene
exclusion policy is retained and still needs transparent reporting.

## Sample-expression context preparation

Script 10 follows the existing 20 positive, 7 negative and 12 divergent priorities
(39 unique priority genes, a different list from the 39 negative shared genes).
It joins original GSEA representative feature IDs to canonical sample counts;
all 117 study-gene mappings were checked against the three count files. OIPN/CCI
use archived normalized counts. NC regenerates size factors from its raw counts
using the original >=10 counts in >=3 samples filter, without refitting DE.

It exports log2(normalized count+1), per-gene within-study display z scores and
group diagnostics including the largest sample's fraction of group counts.
Plots are generated only by the user's local R execution. This lets us investigate
heme, immune and remodeling candidates sample by sample, but cannot estimate cell
fractions, establish contamination, automatically remove a sample, or assign a
transcript change to neurons. Numerical context outputs remain pending R execution.

## Local execution

After `git pull --ff-only`, run from the three-study repository root:

```r
source("scripts/cross_model/09_OIPN_reproducibility_rerun.R")
Sys.unsetenv("LE_TABLES_ONLY")
source("scripts/cross_model/08_three_dataset_negative_leading_edge.R")
Sys.unsetenv("CONTEXT_TABLES_ONLY")
source("scripts/cross_model/10_three_dataset_sample_expression_context.R")
```

Script 09 requires the upstream/downstream R packages and may download the cached
GEO source archive. Script 08 uses existing CSVs, plus ggplot2 for plots. Script 10
requires DESeq2 and ggplot2. To regenerate just the negative tables, set
`LE_TABLES_ONLY=true`. Python audit reproduction:

```bash
python scripts/cross_model/audit_negative_leading_edge.py
```

Review the rerun comparisons and sample-expression diagnostics before any change
to the canonical inputs, cohort selection or manuscript claims. Upload the new
gene-set lock, rerun metadata/comparisons and locally generated figures after R
execution. No R execution or visual figure review was performed in this environment.
