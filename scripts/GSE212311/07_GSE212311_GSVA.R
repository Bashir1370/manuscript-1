#!/usr/bin/env Rscript
# GSE212311 CCI day 11: sample-level Hallmark GSVA on the locked six samples.
# Run from the repository root after steps 01, 02 and 06:
# source("scripts/GSE212311/07_GSE212311_GSVA.R")
# Positive score difference denotes higher scores in CCI than Sham.

required <- c("DESeq2", "GSVA", "msigdbr", "pheatmap")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required R packages: ", paste(missing, collapse = ", "))
if (!file.exists("README.md") || !dir.exists("scripts"))
  stop("Run from the manuscript-1 repository root.")

root <- file.path("results", "GSE212311_CCI_L4L6_day11")
meta_file <- file.path(root, "locked_primary_samples.csv")
count_file <- file.path(root, "raw_count_matrix.csv")
de_file <- file.path(root, "DESeq2", "GSE212311_CCI_vs_Sham_DE_all_genes.csv")
mapping_file <- file.path(root, "pathway_analysis_source_aware",
                          "feature_to_symbol_audit.csv")
for (input in c(meta_file, count_file, de_file, mapping_file))
  if (!file.exists(input)) stop("Missing input: ", input,
                                 ". Complete steps 01, 02 and 06 first.")

expected <- c("GSM6523751", "GSM6523752", "GSM6523753",
              "GSM6523748", "GSM6523749", "GSM6523750")
sham_ids <- expected[1:3]
cci_ids <- expected[4:6]
meta <- read.csv(meta_file, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(meta) != 6L || !identical(meta$gsm, expected) ||
    !identical(meta$condition, rep(c("Sham", "CCI"), each = 3L)) ||
    any(meta$day != 11L) || any(meta$tissue != "ipsilateral_L4-L6_DRG"))
  stop("Metadata disagree with the locked six-sample day-11 contrast.")

count_data <- read.csv(count_file, row.names = 1L, check.names = FALSE)
if (nrow(count_data) != 43350L || !identical(colnames(count_data), expected) ||
    anyDuplicated(rownames(count_data)))
  stop("Count matrix differs from step 01 output.")
counts <- as.matrix(count_data)
if (!is.numeric(counts) || anyNA(counts) || any(!is.finite(counts)) ||
    any(counts < 0) || any(counts != round(counts)) ||
    any(counts > .Machine$integer.max))
  stop("Invalid integer count matrix.")
storage.mode(counts) <- "integer"
keep <- rowSums(counts >= 10L) >= 3L
de <- read.csv(de_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("gene_id", "stat") %in% names(de)) ||
    anyDuplicated(de$gene_id) || !setequal(de$gene_id, rownames(counts)[keep]))
  stop("Filtered feature inventory does not match step 02 DESeq2 output.")

mapping <- read.csv(mapping_file, check.names = FALSE,
                    stringsAsFactors = FALSE)
if (!all(c("gene_id", "symbol", "stat", "source") %in% names(mapping)) ||
    nrow(mapping) < 5000L || anyNA(mapping$gene_id) ||
    anyNA(mapping$symbol) || anyDuplicated(mapping$gene_id) ||
    anyDuplicated(mapping$symbol) ||
    !all(mapping$gene_id %in% de$gene_id) ||
    !all(mapping$source %in% c("GEO_MSTRG_RGD_protein_coding", "OrgDb_ENS")))
  stop("Source-aware feature-to-symbol audit is incomplete or ambiguous.")
de_idx <- match(mapping$gene_id, de$gene_id)
if (any(!is.finite(mapping$stat)) ||
    !isTRUE(all.equal(mapping$stat, de$stat[de_idx], tolerance = 1e-8)))
  stop("Source-aware audit statistics differ from step 02 output.")

# Reconstruct the same filtered DESeq2 design and a non-blind VST for
# sample-level scoring; mapping is locked by step 06, not reselected by effect.
rownames(meta) <- meta$gsm
meta$condition <- factor(meta$condition, levels = c("Sham", "CCI"))
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts[keep, , drop = FALSE],
                                      colData = meta, design = ~ condition)
dds <- DESeq2::DESeq(dds, quiet = TRUE)
vst <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = FALSE))
if (!identical(colnames(vst), expected) ||
    !identical(rownames(vst), rownames(counts)[keep]) ||
    any(!is.finite(vst)))
  stop("VST output has unexpected features, samples or values.")
expr <- vst[match(mapping$gene_id, rownames(vst)), expected, drop = FALSE]
rownames(expr) <- mapping$symbol
if (anyNA(expr) || anyDuplicated(rownames(expr)))
  stop("Symbol expression matrix failed alignment checks.")

hallmark_table <- msigdbr::msigdbr(species = "Rattus norvegicus",
                                   collection = "H")
if (!all(c("gs_name", "gene_symbol") %in% names(hallmark_table)))
  stop("Unexpected rat Hallmark annotation columns.")
hallmark <- lapply(split(hallmark_table$gene_symbol, hallmark_table$gs_name),
                   unique)
overlap <- vapply(hallmark, function(x) length(intersect(x, rownames(expr))),
                  integer(1))
if (length(hallmark) != 50L || any(overlap < 15L))
  stop("Expected 50 rat Hallmark pathways with >=15 mapped genes each.")
scores <- as.matrix(GSVA::gsva(GSVA::gsvaParam(exprData = expr,
                                                geneSets = hallmark),
                               verbose = FALSE))
if (!setequal(rownames(scores), names(hallmark)) ||
    !identical(colnames(scores), expected) || any(!is.finite(scores)))
  stop("GSVA returned unexpected pathways, samples or values.")
scores <- scores[sort(rownames(scores)), expected, drop = FALSE]

selected <- c(
  "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MITOTIC_SPINDLE", "HALLMARK_P53_PATHWAY",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_IL6_JAK_STAT3_SIGNALING"
)
if (!all(selected %in% rownames(scores))) stop("Selected pathway is absent.")
selected_scores <- scores[selected, , drop = FALSE]
stats <- do.call(rbind, lapply(selected, function(pathway) {
  sham <- as.numeric(scores[pathway, sham_ids])
  cci <- as.numeric(scores[pathway, cci_ids])
  tied <- anyDuplicated(c(sham, cci)) > 0L
  test <- stats::wilcox.test(cci, sham, alternative = "two.sided",
                            exact = !tied, correct = FALSE)
  data.frame(
    pathway = pathway, Sham_mean = mean(sham), CCI_mean = mean(cci),
    mean_difference_CCI_minus_Sham = mean(cci) - mean(sham),
    Sham_median = median(sham), CCI_median = median(cci),
    wilcox_W = unname(test$statistic), p_value = unname(test$p.value),
    p_method = if (tied) "asymptotic_tied_ranks" else "exact_two_sided",
    complete_separation = if (min(cci) > max(sham)) "CCI_above_Sham"
                          else if (max(cci) < min(sham)) "CCI_below_Sham"
                          else "overlap", stringsAsFactors = FALSE
  )
}))
stats$FDR_BH_8 <- p.adjust(stats$p_value, method = "BH")

out <- file.path(root, "GSVA_source_aware")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(mapping[, c("gene_id", "symbol", "source")],
          file.path(out, "GSVA_symbol_mapping.csv"), row.names = FALSE)
write.csv(data.frame(pathway = names(hallmark),
                     hallmark_genes = lengths(hallmark), mapped_genes = overlap),
          file.path(out, "GSVA_pathway_overlap.csv"), row.names = FALSE)
write.csv(data.frame(pathway = rownames(scores), scores, check.names = FALSE),
          file.path(out, "GSVA_Hallmark_scores.csv"), row.names = FALSE)
write.csv(data.frame(pathway = selected, selected_scores, check.names = FALSE),
          file.path(out, "selected_GSVA_program_scores.csv"), row.names = FALSE)
write.csv(stats, file.path(out, "GSVA_program_statistics.csv"), row.names = FALSE)
annotation <- data.frame(condition = meta$condition, row.names = meta$gsm)
pheatmap::pheatmap(selected_scores, scale = "row", cluster_cols = FALSE,
                   annotation_col = annotation,
                   main = "GSE212311 day 11: selected Hallmark GSVA",
                   filename = file.path(out, "GSVA_selected_programs_heatmap.png"),
                   width = 9, height = 6)
writeLines(c(
  "GSE212311: day-11 ipsilateral L4-L6 DRG, three CCI versus three Sham.",
  "Input: integer counts, DESeq2 design-based VST, 50 rat Hallmark sets.",
  "Exactly the source-aware unique symbol mapping from step 06 is used.",
  paste("Filtered features:", sum(keep)),
  paste("Unique symbols scored:", nrow(expr)),
  paste("Hallmark sets scored:", nrow(scores)),
  "Eight pathways selected to match the OIPN and NC GSVA analyses.",
  "Two-sided Wilcoxon, exact without ties; BH adjustment across eight tests.",
  "For untied 3 versus 3, the minimum two-sided exact p value is 0.1.",
  "Interpret GSVA as within-study descriptive scores; compare directions, not raw score sizes, across studies.",
  "Bulk DRG scores do not identify cell of origin or establish neuronal cell-cycle re-entry."
), file.path(out, "GSVA_analysis_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "GSVA_sessionInfo.txt"))
message("GSE212311 sample-level Hallmark GSVA complete: ", out)
