# Backup Study Design (No additional wet-lab required)

## Status

This represents the fallback study design for manuscript revision if alternative hypotheses do not provide a stronger and more defensible direction.

## Research question

**Identification of conserved molecular programs associated with oxaliplatin-induced peripheral neuropathy through integrative transcriptomic analysis of dorsal root ganglia.**

Central question:

What transcriptional programs are altered after oxaliplatin exposure in DRG tissue, and which of these programs are conserved across neuropathic conditions?

## Main conceptual revision

The study will not claim neuronal cell-cycle re-entry, neuronal senescence, causality, or therapeutic targeting.

The interpretation will be shifted from a mechanistic neuronal model:

Oxaliplatin → neuronal cell-cycle re-entry → senescence → inflammation → pain

Toward a bulk DRG transcriptomic model:

Oxaliplatin → transcriptional remodeling of DRG tissue → conserved neuropathy-associated molecular programs

## Aim 1: Characterize oxaliplatin-associated transcriptional changes

Dataset:

GSE160543

Comparison:

Oxaliplatin vs Vehicle

Analyses:

- DESeq2
- GSEA
- Hallmark pathway analysis
- GO/KEGG enrichment
- GSVA

## Aim 2: Identify conserved neuropathy-associated programs

Dataset:

GSE246156

Strategy:

Avoid direct interpretation of shared DEG lists as disease drivers.

Perform pathway-level comparison between:

- Oxaliplatin-associated transcriptional programs
- Nerve compression-associated transcriptional programs

Focus:

- stress responses
- inflammatory programs
- DNA damage responses
- metabolic alterations
- cell-cycle-associated programs if supported by data

## Aim 3: Integrate existing experimental validation

Existing qPCR data will be retained.

Cdk1, Cdc20, and Cdkn1a will be presented as representative genes from altered transcriptional programs rather than mechanistic drivers.

Inflammatory markers will support the presence of transcriptional remodeling but will not establish a causal cascade.

## Claims to remove

- Post-mitotic sensory neuron cell-cycle re-entry
- Neuronal senescence
- Cdk1/Cdc20/Cdkn1a as causal drivers
- Therapeutic target claims

## Claims supported by current data

- Oxaliplatin induces measurable transcriptional changes in bulk DRG tissue.
- Some molecular programs may be conserved across neuropathic conditions.
- Candidate genes can be validated at the mRNA level in an independent animal model.

## Additional computational analysis planned

- Pathway-level comparison instead of gene-overlap-only strategy
- Cell-type interpretation using published marker signatures (exploratory)
- Sensitivity analysis of pathway robustness

This design avoids additional wet-lab experiments and aligns interpretation with the actual resolution of bulk DRG transcriptomic data.
