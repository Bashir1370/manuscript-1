# Stage 0: GEO sample metadata audit

Date: 2026-09-28. Source: official GEO series-matrix metadata downloaded from the URLs and SHA-256 checksums in `data/metadata/source_provenance.csv`. Generated sample-level and group-level tables are in `data/metadata/`. Reproduce with `python3 scripts/01_audit_geo_metadata.py`.

## GSE160543

- Twelve bulk rat DRG RNA-seq samples: Vehicle 4, Paclitaxel 4, Oxaliplatin 4.
- GEO titles and treatment labels agree with the sample metadata in the author's earlier `GSE160543.zip` output. That archive comes from a separate ferroptosis screening analysis and does not reproduce the submitted manuscript's DEG counts.
- GEO sample metadata do not provide individual animal IDs, DRG spinal levels, tissue side, or a processing batch field. These are blank rather than inferred in the manifest.

## GSE246156

- Forty-eight samples in all: 36 DRG and 12 sciatic-nerve samples.
- DRG sampling is a full 2 (compression/sham) × 2 (day 3/day 7) × 3 (L4/L5/L6) design, with **three samples per cell**. The sciatic nerve has three samples per condition and day.
- The manuscript's “9 sham and 9 NC” could represent all three DRG levels at one day, but the manuscript does not give the 18 GSM IDs. This is a hypothesis, not an established reconstruction.
- A pooled 9-versus-9 comparison across L4/L5/L6 would mix anatomic levels; treating all 18 rows as independent animals is unwarranted unless the source records establish animal identity. The series-matrix metadata do not identify animal IDs or tissue side. Replicate number 1 at L4 and replicate number 1 at L5 must not be assumed to be the same rat or different rats.
- The source's overall design says nerve entrapment for 3 or 7 days; it does **not** describe a 4-week RNA-seq time point. Earlier secondary search snippets giving 4 weeks should not override the primary GEO metadata.

## Decision before reanalysis

1. Obtain the exact GSM IDs used for the original 9-versus-9 NC comparison and the count source/analysis code, if available.
2. Choose a prespecified coherent contrast within a single DRG level and day (3 versus 3) or a model accounting for level/day and biological dependence once animal IDs are verified. Do not treat a mixture of levels as a simple homogeneous phenotype group.
3. Reassess whether WGCNA is defensible at the actual number of independent animals. A large module from a mixed-level dataset can reflect tissue level or time rather than compression.
4. Recompute the submitted manuscript's 384 OIPN DEGs and 63-gene overlap from its exact inputs. The author's earlier `GSE160543.zip` was created under a different analysis and therefore does not verify those numbers.

No differential expression or cell-type inference was performed in this audit.
