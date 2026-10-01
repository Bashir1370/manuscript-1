# Three-study GO Biological Process experiment

Current status: local GO:BP results and selected-evidence sensitivity are archived. See [GO_BP_manuscript_results_methods_legends.md](GO_BP_manuscript_results_methods_legends.md) for the manuscript section and [GO_BP_selected_evidence_and_sensitivity.md](GO_BP_selected_evidence_and_sensitivity.md) for production verification. The validation-at-implementation notes below are historical; subsequent local execution is documented in those files.

Branch: `three-dataset-gobp-exploration`.
Base: `three-dataset-neuropathy-analysis`, commit `99d5a8e31e25126ff99e4392de1657ce71dc16cc`.

The purpose is to explore shared and divergent DRG transcriptional responses in
OIPN GSE160543 (4 Vehicle +4 Oxaliplatin), NC GSE246156 (3 Sham +3 Compression,
L5 day 7) and CCI GSE212311 (3 Sham +3 CCI, ipsilateral L4-L6 day 11), using
**Human MSigDB C5:GO:BP mapped to rat orthologs**, replacing the 50 Hallmarks in
the new analysis. This is exploratory follow-up after examining Hallmark results.
The prior three-study restriction and its post hoc exclusion rationale still apply.

All inherited files are retained as historical archives. New scripts live under
`scripts/gobp/`; outputs go only to `results/GO_BP_three_dataset/` and the new
gene-set lock. No original DE fit, rank, sample assignment, Hallmark result or
historical script is modified. The parked script 09 is not a dependency.

## What is repeated

1. Lock the complete C5:GO:BP collection once with `msigdbr(db_species="HS",
   species="Rattus norvegicus", collection="C5", subcollection="GO:BP")`.
   Archive distinct rat memberships, full orthology/mapping audit, pathway
   descriptions and GO IDs, database/package versions, date and checksums.
   Later executions reuse the lock and verify it; never silently download a new
   release. An incomplete/corrupt lock stops the run.
2. GSEA uses the existing signed gene statistics and exact historical feature
   representatives. Each study is analyzed separately with `fgseaMultilevel`,
   standard two-sided enrichment, gseaParam=1, eps=1e-10, nPermSimple=10000,
   sampleSize=101, one worker, and study-specific seeds 160543/246156/212311.
   All GO terms are considered. Eligibility is 15-500 genes overlapping each
   study's rank, matching the historical GSEA size bounds. Coverage/exclusions
   are exported for every locked term. This new common GSEA engine differs from
   the historical OIPN fgseaSimple engine; inherited ranks are unchanged.
3. BH adjustment covers the **full eligible GO:BP family within each study**,
   including a conservative denominator for numerically failed eligible tests.
   There is no BH-over-50, DEG-only input, correction over displayed pathways, or
   cross-study combined FDR. Terms not tested or numerically failed stay missing.
4. Shared positive/negative requires the same strict NES sign and pathway FDR
   <0.05 in **all three** studies. Opposite-direction selection requires OIPN's
   strict NES sign to oppose both NC and CCI, without a pathway FDR selection
   filter. Significance is annotated separately. Separate historical positive
   support-pattern flags are also exported.
5. GSVA is recalculated for the locked GO collection and fitted independently
   per study with limma, Neuropathy minus Control, eBayes trend=FALSE,
   robust=FALSE, and BH across **all scored GO:BP sets** in that study.
   Report score means, deltas, confidence intervals and p/FDR; compare directions,
   not calibrated raw score magnitudes across studies.
6. Leading edges are extracted for shared-positive, shared-negative and
   opposite-direction GO terms. Shared membership must occur in the **same GO
   term** in all three studies. Divergent follow-up retains the historical union
   of leading edges and requires opposite gene log2FC signs between OIPN and
   both physical-injury cohorts; membership in all three is not required.
   Preserve each original feature ID, statistic, log2FC and whole-study gene FDR.
7. Gene priorities retain the prior rule: shared genes must have the declared
   common effect direction in all three; divergent genes must have opposite
   signs. Gene FDR <0.05 must hold in OIPN and at least one of NC/CCI. Export the
   strict all-three subset separately. Keep union evidence and all pathway-gene
   memberships; deduplicate STRING gene lists.
8. Display expression of these priorities sample by sample and report group
   means/ranges and the largest sample's fraction of normalized group counts.
   No automatic sample exclusions or new gene-level significance tests.

## Explicit GSVA difference

An OIPN VST matrix is not present in the repository and regenerating it would
require the parked OIPN rerun. Therefore this experiment uses the **same
log2(normalized count +1) transform in all three studies**, Gaussian GSVA with
tau=1, maxDiff=TRUE, absRanking=FALSE and SerialParam. OIPN and CCI use their
archived normalized counts; NC size factors are estimated on the original raw
counts/filter (>=10 in >=3 samples), without refitting differential expression.
The exact GSEA feature representatives are reused for GSVA to keep gene evidence
aligned; those historical representatives can depend on the original DE
statistic, and this is documented rather than independent feature selection.
Constant expression rows are removed before size eligibility is assessed.

This is a new harmonized GSVA scoring step, **not an exact rerun of historical
VST GSVA**. Differences from historical Hallmark GSVA cannot be attributed only
to replacing the gene-set collection. GSVA and GSEA use the same samples and
are complementary descriptions, not independent validation.

## Figures and empty selections

The user draws figures locally. New scripts generate paginated PDF/PNG heatmaps
(40 rows per page), rather than an unreadable many-thousand-row GO heatmap.
NES figures use a common symmetric scale without row scaling; gene log2FC uses
the common display-only clipping range +/-6 (`LE_COLOR_LIMIT`); sample scores
are z-standardized only within each study/row. Controls precede injury samples
in every panel. Raw estimates remain in the CSV files. Short indexed leading-
edge filenames avoid Windows path-length problems; the corresponding GO names
are recorded in `figure_pathway_index.csv`.

No selected pathways or genes is a legitimate result. The workflow exports
header-only tables/empty STRING lists and notes, and removes stale figure pages
for its own generated plot families. All selected memberships remain in tables;
page splitting is solely for readability, not term filtering or GO simplification.

## Interpretation

GO annotations are not necessarily direction-specific activation signatures.
Positive NES indicates concentration toward positive gene statistics, not proof
of increased pathway activity, pain, successful regeneration or cell-type origin.
GO parent/child terms and gene memberships overlap substantially, so many terms
can represent one underlying program. Do not count them as independent mechanisms.
The larger multiple-testing family can reduce significance relative to Hallmark;
this does not invalidate the historical Hallmark findings. Report mapping and
coverage, and do not change size/FDR thresholds to obtain preferred results.

## Local execution

Create a separate worktree in PowerShell:

```powershell
cd D:\manuscript_1\manuscript-1-three-dataset
git fetch origin
git worktree add ..\manuscript-1-gobp -b three-dataset-gobp-exploration origin/three-dataset-gobp-exploration
```

If that local branch already exists, omit `-b three-dataset-gobp-exploration`
and use `three-dataset-gobp-exploration` as the final argument.

Required CRAN packages: current msigdbr and ggplot2. Required Bioconductor
packages: fgsea, GSVA (gsvaParam API), limma, DESeq2 and BiocParallel.

```r
setwd("D:/manuscript_1/manuscript-1-gobp")
# Install only if missing:
# install.packages(c("BiocManager", "msigdbr", "ggplot2"))
# BiocManager::install(c("fgsea", "GSVA", "limma", "DESeq2", "BiocParallel"), ask=FALSE)
Sys.unsetenv("GOBP_TABLES_ONLY")
source("scripts/gobp/run_GO_BP_three_dataset.R")
```

For a quicker first look without writing figures:

```r
Sys.setenv(GOBP_TABLES_ONLY="true")
source("scripts/gobp/run_GO_BP_three_dataset.R")
read.csv("results/GO_BP_three_dataset/comparison/selection_summary.csv")
```

This still performs GSEA, GSVA, leading-edge extraction and expression context.
To draw plots later without rerunning GSEA/GSVA scoring:

```r
Sys.unsetenv("GOBP_TABLES_ONLY")
source("scripts/gobp/02_compare_GO_BP.R")
# 03 rescoring is necessary to draw GSVA plots through that entry script.
source("scripts/gobp/03_GSVA_GO_BP.R")
source("scripts/gobp/04_leading_edge_GO_BP.R")
source("scripts/gobp/05_sample_expression_GO_BP.R")
```

Main outputs:

- `comparison/all_GO_BP_pathways_classified.csv`, `selection_summary.csv`;
- `GSEA/<study>/GSEA_GO_BP_results.csv`, coverage, rank and feature evidence;
- `GSVA/GSVA_limma_all_GO_BP_three_studies.csv`, scores and selected sample table;
- `leading_edge/<shared_positive|shared_negative|divergent>/` evidence/priorities;
- `sample_expression/priority_sample_expression.csv` and diagnostics;
- `data/gene_sets/GO_BP_rat_locked/` membership, mapping audit and manifest.

The new lock and scientific results are generated by local R execution; their
existence is not implied by the presence of the new scripts. Archive them on this
experimental branch after reviewing the local run.

## Validation before publishing the scripts

Executed with R 4.3.3, limma 3.58.1 and statmod 1.5.0. All nine R files parse.
The new generic comparison/leading-edge functions, exercised on the inherited
Hallmark inputs, reproduce 10 shared-positive, two shared-negative and seven
opposite-direction pathways; positive membership/priorities 109/78/20/8,
negative 44/39/7/1, and 126 divergent genes with 12 priorities. Metadata,
feature-count alignment and every representative/rank statistic agree in all
three cohorts. Missing/zero NES, missing FDR, empty selections, compressed
lock roundtrip, checksum corruption rejection and a 137-set limma/BH family
were tested. See the complete validation log.

The full MSigDB GO download, fgsea GO scoring, GSVA GO scoring and new figures
were **not** run here. Those require the author's local R execution with all
packages. No GO result counts or biological conclusions have been invented.
