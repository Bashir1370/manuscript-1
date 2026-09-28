#!/usr/bin/env Rscript

# GSE160543 pathway analysis
# Input: DESeq2 results from scripts/02_GSE160543_oxaliplatin_primary.R
# Goal: characterize bulk DRG transcriptional programs after oxaliplatin exposure.
# This script intentionally uses ranked statistics and avoids mechanistic claims.

options(stringsAsFactors = FALSE)

packages <- c("fgsea", "msigdbr", "dplyr", "ggplot2", "data.table")
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("Install package: ", p)
  }
}

library(fgsea)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(data.table)

input <- "results/GSE160543_Oxaliplatin_vs_Vehicle/DE_all_genes.csv"
outdir <- "results/GSE160543_Oxaliplatin_vs_Vehicle/pathway_analysis"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input)) stop("Missing DESeq2 result file: ", input)

res <- fread(input)

# Remove genes without statistics and duplicate symbols
rank_table <- res %>%
  filter(!is.na(stat), !is.na(symbol), symbol != "") %>%
  arrange(desc(stat)) %>%
  distinct(symbol, .keep_all = TRUE)

ranks <- rank_table$stat
names(ranks) <- rank_table$symbol
ranks <- sort(ranks, decreasing = TRUE)

write.csv(data.frame(symbol = names(ranks), statistic = ranks),
          file.path(outdir, "ranked_statistics.csv"),
          row.names = FALSE)

hallmark <- msigdbr(species = "Rattus norvegicus", category = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)

fg <- fgsea(
  pathways = hallmark,
  stats = ranks,
  minSize = 15,
  maxSize = 500,
  nperm = 10000
) %>%
  arrange(padj)

write.csv(as.data.frame(fg),
          file.path(outdir, "GSEA_Hallmark_results.csv"),
          row.names = FALSE)

selected <- fg %>%
  filter(gs_name %in% c(
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_P53_PATHWAY",
    "HALLMARK_DNA_REPAIR",
    "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
  ))

write.csv(selected,
          file.path(outdir, "prespecified_biological_programs.csv"),
          row.names = FALSE)

p <- ggplot(selected,
            aes(x = reorder(gs_name, NES), y = NES)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(title = "GSE160543 Oxaliplatin vs Vehicle: selected programs",
       x = NULL,
       y = "Normalized enrichment score")

ggsave(file.path(outdir, "selected_programs_NES.png"),
       p, width = 8, height = 5, dpi = 300)

writeLines(capture.output(sessionInfo()),
           file.path(outdir, "sessionInfo.txt"))
