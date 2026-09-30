#!/usr/bin/env Rscript
# Follow up the existing opposite-NES Fig2 panel; no GSEA or gene-set rerun.
run_divergent_leading_edge <- function(draw_plots = TRUE) {
  # One common linear color scale for all studies and pathways; display only.
  color_limit <- suppressWarnings(as.numeric(Sys.getenv("LE_COLOR_LIMIT", "6")))
  if (length(color_limit)!=1L || !is.finite(color_limit) || color_limit<=0)
    stop("LE_COLOR_LIMIT must be one finite positive number (default: 6).")
  read <- function(path) {
    if (!file.exists(path)) stop("Missing input: ", path, "; run from repository root.")
    read.csv(path, check.names=FALSE, stringsAsFactors=FALSE)
  }
  selection_path <- "results/manuscript_hallmark_three_dataset/Fig2_direction_selection.csv"
  classification_path <- "results/manuscript_hallmark_three_dataset/all_50_pathways_classified.csv"
  selection <- read(selection_path); classification <- read(classification_path)
  studies_check <- c("OIPN_GSE160543", "NC_GSE246156", "CCI_GSE212311")
  required <- c("pathway", paste0(studies_check, "_NES"), paste0(studies_check, "_FDR"))
  if (!all(required %in% names(classification)) || nrow(classification)!=50L ||
      anyDuplicated(classification$pathway)) stop("Invalid complete GSEA classification.")
  for (col in required[-1L]) {
    if (!is.numeric(classification[[col]]) || any(!is.finite(classification[[col]])))
      stop("Invalid classification values: ", col)
  }
  n <- as.matrix(classification[paste0(studies_check, "_NES")])
  f <- as.matrix(classification[paste0(studies_check, "_FDR")])
  if (any(f<0 | f>1)) stop("Invalid pathway FDR.")
  expected <- classification$pathway[(n[,1]>0 & n[,2]<0 & n[,3]<0) |
                                    (n[,1]<0 & n[,2]>0 & n[,3]>0)]
  if (!"pathway" %in% names(selection) || anyDuplicated(selection$pathway) ||
      !length(expected) || !setequal(selection$pathway, expected))
    stop("Fig2 selection is not the opposite-direction panel; run script 03 with HALLMARK_DIVERGENCE_MODE unset.")
  selected <- sort(expected, method="radix")
  if (draw_plots && !requireNamespace("ggplot2", quietly=TRUE))
    stop("Install ggplot2, or set DIVERGENT_LE_TABLES_ONLY=true.")
  roots <- c(OIPN_GSE160543 = "results/GSE160543_Oxaliplatin_vs_Vehicle",
    NC_GSE246156 = "results/GSE246156_NC_L5_day7",
    CCI_GSE212311 = "results/GSE212311_CCI_L4L6_day11")
  studies <- names(roots)
  inputs <- c(selection_path, classification_path); gene_tables <- gsea_tables <- list()
  # Preserve the feature/probe representative actually used by each GSEA.
  for (study in studies) {
    root <- roots[[study]]
    subdir <- switch(study, OIPN_GSE160543="Pathway_analysis",
                     NC_GSE246156="pathway_analysis", CCI_GSE212311="pathway_analysis_source_aware")
    gp <- file.path(root, subdir, "GSEA_Hallmark_results.csv")
    g <- read(gp); inputs <- c(inputs, gp)
    if (!all(c("pathway", "NES", "padj", "leadingEdge") %in% names(g)) ||
        nrow(g) != 50L || anyDuplicated(g$pathway) || !all(selected %in% g$pathway)) stop("Invalid GSEA: ", study)
    g <- g[match(selected, g$pathway), c("pathway", "NES", "padj", "leadingEdge")]
    if (any(!is.finite(g$NES)) || any(!is.finite(g$padj)) || any(g$padj < 0 | g$padj > 1) ||
        anyNA(g$leadingEdge) || any(!nzchar(g$leadingEdge))) stop("Invalid selected GSEA values: ", study)
    expected <- classification[match(g$pathway, classification$pathway), ]
    if (!isTRUE(all.equal(g$NES, expected[[paste0(study, "_NES")]], tolerance=1e-8)) ||
        !isTRUE(all.equal(g$padj, expected[[paste0(study, "_FDR")]], tolerance=1e-8)))
      stop("Classification differs from original GSEA: ", study)
    g$edges <- strsplit(g$leadingEdge, ";", fixed = TRUE)
    g$edges <- lapply(g$edges, function(x) sort(unique(trimws(x))))
    if (any(vapply(g$edges, function(x) any(!nzchar(x)), logical(1)))) stop("Blank leading-edge identifier.")
    if (study == "OIPN_GSE160543") {
      dp <- file.path(root, "DE_all_genes.csv"); d <- read(dp); inputs <- c(inputs, dp)
      d <- d[!is.na(d$symbol) & nzchar(d$symbol) & is.finite(d$stat), ]
      d <- d[order(-d$stat), ]; d <- d[!duplicated(d$symbol), ]
      dt <- data.frame(symbol=d$symbol, feature_id=d$gene_id, statistic=d$stat,
        log2FC=d$log2FoldChange, gene_FDR=d$padj, mapping_source="original_DE_symbol")
    } else if (study == "NC_GSE246156") {
      dp <- file.path(root, "DESeq2", "GSE246156_Compression_vs_Sham_DE_all_genes.csv")
      mp <- file.path(root, subdir, "Ensembl_annotation_all_mappings.csv")
      d <- read(dp); map <- read(mp); inputs <- c(inputs, dp, mp)
      d$ensembl <- sub("\\.[0-9]+$", "", d$gene)
      pairs <- unique(map[!is.na(map$SYMBOL) & nzchar(map$SYMBOL), c("ENSEMBL", "SYMBOL")])
      ambiguous <- names(which(table(pairs$ENSEMBL) > 1L))
      pairs <- pairs[!pairs$ENSEMBL %in% ambiguous, ]
      d <- merge(d, pairs, by.x="ensembl", by.y="ENSEMBL", sort=FALSE)
      d <- d[is.finite(d$stat), ]; d <- d[order(-abs(d$stat), d$ensembl), ]
      d <- d[!duplicated(d$SYMBOL), ]
      dt <- data.frame(symbol=d$SYMBOL, feature_id=d$gene, statistic=d$stat,
        log2FC=d$log2FoldChange, gene_FDR=d$padj, mapping_source="original_unambiguous_Ensembl_symbol")
    } else {
      dp <- file.path(root, subdir, "feature_to_symbol_audit.csv")
      d <- read(dp); inputs <- c(inputs, dp)
      dt <- data.frame(symbol=d$symbol, feature_id=d$gene_id, statistic=d$stat,
        log2FC=d$log2FoldChange, gene_FDR=d$padj, mapping_source=d$source)
    }
    if (anyDuplicated(dt$symbol) || anyNA(dt$symbol) || any(!is.finite(dt$statistic)) ||
        any(!is.finite(dt$log2FC))) stop("Invalid gene evidence: ", study)
    if (!is.numeric(dt$gene_FDR) || any(!is.na(dt$gene_FDR) & (!is.finite(dt$gene_FDR) | dt$gene_FDR<0 | dt$gene_FDR>1)))
      stop("Invalid original gene FDR: ", study)
    {
      rp <- file.path(root, subdir, "ranked_statistics.csv"); rank <- read(rp); inputs <- c(inputs, rp)
      if (anyDuplicated(rank$symbol) || !setequal(rank$symbol, dt$symbol) ||
          !isTRUE(all.equal(unname(rank$statistic), unname(dt$statistic[match(rank$symbol, dt$symbol)]), tolerance=1e-8))) {
        stop("Representative reconstruction does not match original GSEA ranks: ", study)
      }
    }
    if (!all(unlist(g$edges) %in% dt$symbol)) stop("Leading-edge symbols missing from original rank: ", study)
    gene_tables[[study]] <- dt; gsea_tables[[study]] <- g
  }
  # Gene FDR retains its original genome-wide DE testing family.
  sig <- function(x) !is.na(x) & x < .05
  long <- rows <- list()
  for (p in selected) {
    edges <- lapply(gsea_tables, function(g) g$edges[[match(p, g$pathway)]])
    symbols <- sort(unique(unlist(edges)), method="radix")
    for (symbol in symbols) {
      lfc <- fdr <- rep(NA_real_, 3L)
      ranked <- member <- logical(3L)
      for (i in seq_along(studies)) {
        study <- studies[i]; dt <- gene_tables[[study]]; idx <- match(symbol, dt$symbol)
        ranked[i] <- !is.na(idx); member[i] <- symbol %in% edges[[i]]
        gp <- gsea_tables[[study]]; gi <- match(p, gp$pathway)
        if (ranked[i]) { lfc[i] <- dt$log2FC[idx]; fdr[i] <- dt$gene_FDR[idx] }
        long[[length(long)+1L]] <- data.frame(pathway=p, symbol=symbol, study=study,
          rank_available=ranked[i], leading_edge=if (ranked[i]) member[i] else NA,
          status=if (!ranked[i]) "unavailable_to_GSEA" else if (member[i]) "leading_edge" else "ranked_not_leading_edge",
          feature_id=if (ranked[i]) dt$feature_id[idx] else NA_character_,
          mapping_source=if (ranked[i]) dt$mapping_source[idx] else NA_character_,
          rank_statistic=if (ranked[i]) dt$statistic[idx] else NA_real_,
          log2FC=lfc[i], gene_FDR=fdr[i], pathway_NES=gp$NES[gi], pathway_FDR=gp$padj[gi])
      }
      opposite <- all(ranked) && ((lfc[1]>0 && lfc[2]<0 && lfc[3]<0) ||
                                  (lfc[1]<0 && lfc[2]>0 && lfc[3]>0))
      pathway_nes <- gsea_tables[[1]]$NES[match(p, gsea_tables[[1]]$pathway)]
      r <- data.frame(pathway=p, symbol=symbol, n_ranked=sum(ranked), n_leading_edge=sum(member),
        OIPN_matches_pathway_direction=if (all(ranked)) lfc[1]*pathway_nes>0 else FALSE,
        opposite_both_physical=opposite, OIPN_significant=sig(fdr[1]),
        NC_significant=sig(fdr[2]), CCI_significant=sig(fdr[3]),
        priority_OIPN_plus_physical=opposite && sig(fdr[1]) && (sig(fdr[2]) || sig(fdr[3])),
        strict_significant_all3=opposite && all(sig(fdr)))
      for (i in seq_along(studies)) {
        r[[paste0(studies[i], "_LE")]] <- if (ranked[i]) member[i] else NA
        r[[paste0(studies[i], "_log2FC")]] <- lfc[i]
        r[[paste0(studies[i], "_gene_FDR")]] <- fdr[i]
      }
      rows[[length(rows)+1L]] <- r
    }
  }
  long <- do.call(rbind, long); summary <- do.call(rbind, rows)
  opposite <- summary[summary$opposite_both_physical, , drop=FALSE]
  priority <- summary[summary$priority_OIPN_plus_physical, , drop=FALSE]
  totals <- do.call(rbind, lapply(selected, function(p) {
    d <- summary[summary$pathway==p, ]
    data.frame(pathway=p, union_LE_genes=nrow(d), all_three_ranked=sum(d$n_ranked==3L),
      opposite_both_physical=sum(d$opposite_both_physical),
      priority_OIPN_plus_physical=sum(d$priority_OIPN_plus_physical),
      strict_significant_all3=sum(d$strict_significant_all3))
  }))
  out <- "results/divergent_Hallmark_leading_edge_three_dataset"
  dir.create(out, recursive=TRUE, showWarnings=FALSE)
  write <- function(x, name) write.csv(x, file.path(out, name), row.names=FALSE, na="")
  write(classification[match(selected, classification$pathway), ], "selected_pathways.csv")
  write(long, "gene_evidence_long.csv"); write(summary, "gene_membership_summary.csv")
  write(opposite, "opposite_direction_genes.csv")
  write(priority, "priority_OIPN_plus_physical.csv")
  write(summary[summary$strict_significant_all3, , drop=FALSE], "strict_significant_all3.csv")
  write(totals, "pathway_overlap_summary.csv")
  unique_priority <- sort(unique(priority$symbol), method="radix")
  writeLines(unique_priority, file.path(out, "priority_genes_STRING.txt"))
  for (label in c("positive", "negative")) {
    keep <- if (label=="positive") priority$OIPN_GSE160543_log2FC>0 else priority$OIPN_GSE160543_log2FC<0
    writeLines(sort(unique(priority$symbol[keep]), method="radix"),
      file.path(out, paste0("priority_OIPN_", label, "_STRING.txt")))
  }
  counts <- data.frame(metric=c("selected_pathways", "union_pathway_gene_rows",
    "study_evidence_rows", "opposite_pathway_gene_rows", "opposite_unique_genes",
    "priority_pathway_gene_rows", "priority_unique_genes", "strict_unique_genes"),
    count=c(length(selected), nrow(summary), nrow(long), nrow(opposite),
    length(unique(opposite$symbol)), nrow(priority), length(unique_priority),
    length(unique(summary$symbol[summary$strict_significant_all3]))))
  write(counts, "selection_summary.csv")
  write(data.frame(input=inputs, md5=unname(tools::md5sum(inputs))), "input_checksums.csv")
  write(data.frame(color_scale="linear_log2FC", lower_limit=-color_limit,
    upper_limit=color_limit, shared_across_studies=TRUE, shared_across_pathways=TRUE,
    display_clipping_only=TRUE), "plot_settings.csv")
  if (draw_plots) for (p in selected) {
    genes <- opposite$symbol[opposite$pathway==p]
    if (!length(genes)) next
    d <- long[long$pathway==p & long$symbol %in% genes, ]
    d$symbol <- factor(d$symbol, levels=rev(sort(genes, method="radix")))
    d$study <- factor(d$study, levels=studies, labels=c("OIPN 160543", "NC 246156", "CCI 212311"))
    d$display_log2FC <- pmax(-color_limit, pmin(color_limit, d$log2FC))
    d$mark <- ifelse(!is.na(d$gene_FDR) & d$gene_FDR<.05, "*", "")
    gp <- ggplot2::ggplot(d, ggplot2::aes(study, symbol, fill=display_log2FC)) +
      ggplot2::geom_tile(color="white", linewidth=.25) +
      ggplot2::geom_text(ggplot2::aes(label=mark), size=3) +
      ggplot2::scale_fill_gradient2(low="#2166AC", mid="#F7F7F7", high="#B2182B",
        midpoint=0, limits=c(-color_limit,color_limit), name=paste0("log2FC\nclipped at +/-", color_limit)) +
      ggplot2::labs(title=gsub("_", " ", sub("^HALLMARK_", "", p)),
        subtitle="Gene log2FC: OIPN opposite to both NC and CCI; * gene FDR <0.05",
        x=NULL, y=NULL, caption="Union of original leading edges; membership in all three is not required.\nNeuropathy minus control; no row scaling; exploratory cohort comparison.") +
      ggplot2::theme_minimal(base_size=10) +
      ggplot2::theme(panel.grid=ggplot2::element_blank())
    stem <- paste0(sub("^HALLMARK_", "", p), "_opposite_gene_log2FC")
    height <- max(4, 2+.20*length(genes))
    ggplot2::ggsave(file.path(out, paste0(stem, ".pdf")), gp, width=10, height=height, limitsize=FALSE)
    ggplot2::ggsave(file.path(out, paste0(stem, ".png")), gp, width=10, height=height, dpi=300, bg="white", limitsize=FALSE)
  }
  writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
  message("Divergent Hallmark leading-edge extraction complete: ", length(unique_priority),
    " priority genes; ", out)
  invisible(totals)
}
divergent_leading_edge_result <- run_divergent_leading_edge(
  draw_plots=tolower(Sys.getenv("DIVERGENT_LE_TABLES_ONLY", "false"))!="true")
