# GSE246156 nerve-compression analysis protocol

Branch: `NC-GSE246156-analysis` (created from `main`).
Status: prospective analysis plan; no NC differential-expression or cross-model result has been generated.

## Study question and scope

Which bulk rat DRG transcriptional programs change after L5 nerve entrapment relative to matched sham controls, and which changes have a concordant direction with the GSE160543 oxaliplatin-versus-vehicle response?

Analyze the NC study independently. Compare gene and pathway effects between studies only after each within-study analysis passes input and quality checks. Shared responses are injury-associated tissue-level signatures, not OIPN-specific mechanisms or evidence of neuronal cell-cycle re-entry.

## Verified metadata from the repository's GEO audit

GSE246156 comprises 48 samples: 36 DRG and 12 sciatic nerve. DRG samples span compression/sham, day 3/day 7, and L4/L5/L6, with three samples in each condition–day–level cell. The series metadata do not identify individual animals or tissue side. Sciatic nerve is excluded from the DRG comparison. Replicate labels across levels do not establish biological independence.

| DRG level and day | Compression GSMs | Sham GSMs |
| --- | --- | --- |
| L4, day 3 | GSM7863755–GSM7863757 | GSM7863779–GSM7863781 |
| L5, day 3 | GSM7863758–GSM7863760 | GSM7863782–GSM7863784 |
| L6, day 3 | GSM7863761–GSM7863763 | GSM7863785–GSM7863787 |
| L4, day 7 | GSM7863767–GSM7863769 | GSM7863791–GSM7863793 |
| L5, day 7 | GSM7863770–GSM7863772 | GSM7863794–GSM7863796 |
| L6, day 7 | GSM7863773–GSM7863775 | GSM7863797–GSM7863799 |

These IDs come from `data/metadata/geo_sample_manifest.csv` in `main`. Confirm against the count-file headers and GEO source before analysis.

## Input gate before fitting a model

1. Obtain the official processed count file(s) and record source URL, download date, checksum, feature identifier type, counting method, and whether values are raw integer counts or fractional estimates. Inspect the actual file layout before writing an import routine.
2. Confirm that the six selected columns map exactly to the selected GSMs, with no duplicated, missing, or unexpected samples. Match every sample by identifier, never by file order.
3. Check animal identifiers or author methods for repeated tissue sampling across levels. Do not use a pooled 9-versus-9 contrast or an additive level model as if 18 rows were independent animals without resolving this.
4. Resolve which day/level best matches the biological question and the original study protocol. Record the primary contrast **before** inspecting differential-expression results. A single level at a single day gives a clean 3-versus-3 comparison; other level/day cells can be reported separately as sensitivity analyses.
5. Audit per-sample library size, feature counts, identifier mapping, zero/NA/duplicate features, count integer status, and PCA/sample distances. Do not silently replace missing features with zero or round normalized expression/FPKM into counts.

## Analysis once the gate is met

- Within the locked DRG level/day: DESeq2, `design = ~ condition`, `contrast = c("condition", "Compression", "Sham")`. Report all tested genes with effect size, standard error, statistic, raw and adjusted p-values, plus the filtering rule.
- QC: count distributions and library sizes, VST PCA, sample distances, and an explicit outlier audit. Do not exclude an outlier solely because it weakens the desired result.
- Ranked GSEA: use the complete signed DESeq2 statistic, an auditable gene-ID-to-symbol mapping, and the same rat Hallmark gene-set collection and parameters used for GSE160543 where applicable.
- GO Biological Process: specify up/down sets and the measurable-gene background; interpret this as complementary to ranked GSEA.
- GSVA: score sample-level pathways on appropriately normalized/transformed expression, confirm gene ID overlap, compare Compression with Sham, and report effect direction, exact p-values, and multiplicity adjustment.
- Cross-model comparison: compare signed gene statistics/effect sizes and pathway directions/NES with GSE160543. Do not pool raw count matrices across studies, treat a non-significant result as evidence of specificity, or present overlap as causal validation.
- Existing rat qPCR and behavior support tissue-level mRNA and phenotype only.

## Decision rules

- If the processed NC data are not count-scale, do not pass them to DESeq2; obtain raw counts or select a method appropriate to the available data.
- If sample identity, animal independence, or a valid matched contrast cannot be verified, stop the NC differential-expression claim and report the limitation.
- A program is described as shared only with adequately supported concordant direction in both analyses; absence of evidence in one small 3-versus-3 contrast is inconclusive.
- Cell type, functional cell-cycle activation, senescence, causal pain mechanisms, and therapeutic targeting remain outside the resolution of these bulk DRG data.

## Relationship to existing scripts

The GSE160543 workflow supplies the analysis structure, outputs, and interpretation boundaries. Its archive filenames, five-column RSEM-style `expected_count` parser, complete-case rule, and hard-coded GSMs are dataset-specific and must not be copied into NC. The NC import script will be finalized only after inspecting GSE246156's actual processed files.
