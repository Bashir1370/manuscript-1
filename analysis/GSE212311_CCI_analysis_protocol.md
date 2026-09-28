# GSE212311 CCI analysis: locked input and two-stage workflow

## Biological comparison

Male Sprague-Dawley rat, ipsilateral bulk L4–L6 DRG, 11 days after CCI
or sham surgery. The sham sciatic nerve was exposed without ligation according
to the original publication (Qu et al., *Heliyon* 2024, PMID 38813203).
Primary contrast: CCI minus Sham; three independent GEO libraries per arm.

| Condition | GEO samples | Source columns |
| --- | --- | --- |
| Sham | GSM6523751–GSM6523753 | `count.Sham_1`–`count.Sham_3` |
| CCI | GSM6523748–GSM6523750 | `count.CCI_1`–`count.CCI_3` |

Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212311

The GEO file `GSE212311_genes_fpkm_expression.txt.gz` has 43,350 unique
feature IDs and 25 columns. It has six FPKM columns **and six integer count
columns**. The source MD5 observed on 2026-09-28 is
`71033155b073b1af600c51a20cc9b33a`; the scripts reject a changed file.
There are 13,855 ENSRNOG annotated gene IDs and 29,495 MSTRG loci. The
unannotated MSTRG loci are retained in gene-level testing, but only mapped
ENSRNOG IDs may be interpreted in gene-symbol pathway comparisons.

## Run on the Windows checkout

From PowerShell in `D:\manuscript_1\manuscript-1`:

```powershell
git fetch origin
git switch CCI-GSE212311-analysis
```

In RStudio, set the working directory to the repository root and run:

```r
source("scripts/GSE212311/01_GSE212311_prepare_PCA.R")
```

Requires `DESeq2` and `ggplot2`. Step 01 downloads and verifies the processed
GEO source, validates integer counts and sample mapping, writes the locked
count matrix, and computes blind VST PCA and sample distances. Inspect the
plot and QC files in `results/GSE212311_CCI_L4L6_day11/PCA_6samples/` before
interpreting differential expression.

After QC is reviewed, run:

```r
source("scripts/GSE212311/02_GSE212311_DESeq2.R")
```

Step 02 uses `design = ~ condition` and `CCI - Sham`, filters genes with at
least 10 counts in at least three samples, writes all-gene and BH FDR < 0.05
tables, plus normalized counts and annotated-gene signed statistics. A positive
log2 fold change denotes higher expression in CCI. FPKM never enters DESeq2.

Do not treat L4–L6 day-11 CCI as a matched level/time replication of the
existing L5 day-7 NC study. Analyze models within study, then compare signed
genes/pathways while noting level and time differences. Small 3-vs-3 groups
have limited power. Neither study locates a signal in a particular DRG cell
type based on bulk RNA-seq alone.

The cached GEO download in `data/source-cache/GSE212311/` is ignored by Git.
Commit compact QC and DE output only after inspection.
