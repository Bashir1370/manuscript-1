# Frozen Hallmark membership for reproducibility reruns

`scripts/cross_model/OIPN_reproducibility_helpers.R` creates
`Hallmark_rat_locked.csv` and `Hallmark_rat_lock_metadata.csv` on the first local
OIPN rerun. Both GSEA and GSVA then read the same CSV without querying a live
gene-set release. The CSV retains the complete returned msigdbr metadata,
including database version when exposed, and the companion file records the
acquisition time, package version and checksum. Subsequent runs verify it.

Commit both generated files after inspecting the rerun. Do not edit or regenerate
one silently. Historical Hallmark membership was not archived and is not recovered
by this procedure. The lock applies to the new OIPN rerun; the original NC/CCI
outputs are preserved, not harmonized or rescored. Package-version CSVs and R
session files describe the executed environment, not a dependency restoration lock.
