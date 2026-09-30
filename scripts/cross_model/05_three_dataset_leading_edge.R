#!/usr/bin/env Rscript
# Extract existing GSEA leading edges; do not rerun GSEA or update gene sets.
run_shared_leading_edge <- function(draw_plots = TRUE, direction = c("positive", "negative")) {
  direction <- match.arg(direction)
  effect_sign <- if (direction == "positive") 1 else -1
  selection_file <- file.path("results/manuscript_hallmark_three_dataset",
                              paste0("shared_", direction, ".csv"))
  # One common linear color scale for all studies and pathways; display only.
  color_limit <- suppressWarnings(as.numeric(Sys.getenv("LE_COLOR_LIMIT", "6")))
  if (length(color_limit)!=1L || !is.finite(color_limit) || color_limit<=0)
    stop("LE_COLOR_LIMIT must be one finite positive number (default: 6).")
  selection <- read.csv(selection_file, stringsAsFactors=FALSE)
  selected <- sort(selection$pathway)
  if (!length(selected) || anyDuplicated(selected)) stop("Run three-study GSEA extraction first.")
  if (draw_plots && !requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2, or set LE_TABLES_ONLY=true.")
  read <- function(path) {
    if (!file.exists(path)) stop("Missing input: ", path, "; run from repository root.")
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  roots <- c(OIPN_GSE160543 = "results/GSE160543_Oxaliplatin_vs_Vehicle",
    NC_GSE246156 = "results/GSE246156_NC_L5_day7",
    CCI_GSE212311 = "results/GSE212311_CCI_L4L6_day11")
  studies <- names(roots)
  inputs <- selection_file; gene_tables <- gsea_tables <- list()
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
        any(effect_sign * g$NES <= 0) || anyNA(g$leadingEdge) || any(!nzchar(g$leadingEdge))) stop("Invalid selected GSEA values: ", study)
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
  if (any(vapply(selected, function(p) sum(vapply(gsea_tables, function(g) g$padj[g$pathway==p] < .05, logical(1))) < 3L,
                 logical(1)))) stop(paste0("The three-study shared-", direction, " pathway selection no longer holds."))
  long <- summary_rows <- pathway_rows <- list(); k <- j <- 0L
  for (p in selected) {
    edges <- lapply(gsea_tables, function(g) g$edges[[match(p, g$pathway)]])
    genes <- sort(unique(unlist(edges)))
    for (gene in genes) {
      k <- k + 1L
      ranked <- vapply(gene_tables, function(d) gene %in% d$symbol, logical(1))
      member <- vapply(edges, function(x) gene %in% x, logical(1))
      lfc <- fdr <- rep(NA_real_, length(studies))
      for (i in seq_along(studies)) {
        study <- studies[i]; d <- gene_tables[[study]]; idx <- match(gene, d$symbol)
        g <- gsea_tables[[study]]; gi <- match(p, g$pathway)
        if (!is.na(idx)) { lfc[i] <- d$log2FC[idx]; fdr[i] <- d$gene_FDR[idx] }
        long[[length(long)+1L]] <- data.frame(pathway=p, symbol=gene, study=study,
          rank_available=ranked[i], leading_edge=if (ranked[i]) member[i] else NA,
          status=if (!ranked[i]) "unavailable_to_GSEA" else if (member[i]) "leading_edge" else "ranked_not_leading_edge",
          feature_id=if (ranked[i]) d$feature_id[idx] else NA_character_,
          mapping_source=if (ranked[i]) d$mapping_source[idx] else NA_character_,
          rank_statistic=if (ranked[i]) d$statistic[idx] else NA_real_,
          log2FC=lfc[i], gene_FDR=fdr[i], pathway_NES=g$NES[gi], pathway_FDR=g$padj[gi])
      }
      row <- data.frame(pathway=p, symbol=gene, n_ranked=sum(ranked), n_leading_edge=sum(member),
        n_LE_in_significant_pathways=sum(member & vapply(gsea_tables, function(g) g$padj[g$pathway==p] < .05, logical(1))),
        n_positive_log2FC=sum(lfc > 0, na.rm=TRUE), n_negative_log2FC=sum(lfc < 0, na.rm=TRUE),
        n_gene_FDR_lt_0_05=sum(fdr < .05, na.rm=TRUE), all_three_ranked=all(ranked),
        shared_ge3=sum(member)>=3L, shared_all3=all(member), positive_all3=all(ranked) && all(lfc>0))
      if (direction == "negative") row$negative_all3 <- all(ranked) && all(lfc<0)
      for (i in seq_along(studies)) {
        row[[paste0(studies[i], "_LE")]] <- if (ranked[i]) member[i] else NA
        row[[paste0(studies[i], "_log2FC")]] <- lfc[i]
        row[[paste0(studies[i], "_gene_FDR")]] <- fdr[i]
      }
      summary_rows[[k]] <- row
    }
    j <- j + 1L; ps <- do.call(rbind, summary_rows)[vapply(summary_rows, function(x) x$pathway==p, logical(1)), ]
    pathway_rows[[j]] <- data.frame(pathway=p, union_LE_genes=length(genes),
      shared_ge3=sum(ps$shared_ge3), shared_ge3_all3_ranked=sum(ps$shared_ge3 & ps$all_three_ranked),
      shared_all3=sum(ps$shared_all3), shared_ge3_positive_all3=sum(ps$shared_ge3 & ps$positive_all3))
    if (direction == "negative") pathway_rows[[j]]$shared_ge3_negative_all3 <-
      sum(ps$shared_ge3 & ps$negative_all3)
  }
  long <- do.call(rbind, long); summary <- do.call(rbind, summary_rows); totals <- do.call(rbind, pathway_rows)
  out <- if (direction == "positive") "results/shared_Hallmark_leading_edge_three_dataset" else
    "results/shared_negative_Hallmark_leading_edge_three_dataset"
  dir.create(out, recursive=TRUE, showWarnings=FALSE)
  write <- function(x, name) write.csv(x, file.path(out, name), row.names=FALSE, na="")
  write(long, "gene_evidence_long.csv"); write(summary, "gene_membership_summary.csv")
  write(summary[summary$shared_ge3, ], "shared_ge3.csv")
  write(summary[summary$shared_all3, ], "shared_all3.csv")
  write(totals, "pathway_overlap_summary.csv")
  selection <- do.call(rbind, lapply(studies, function(study) {
    g <- gsea_tables[[study]]
    data.frame(study=study, pathway=g$pathway, NES=g$NES, pathway_FDR=g$padj,
      leading_edge_genes=lengths(g$edges), ranked_symbols=nrow(gene_tables[[study]]))
  }))
  write(selection, "pathway_evidence.csv")
  unique_genes <- sort(unique(summary$symbol[summary$shared_ge3]))
  write(data.frame(symbol=unique_genes, n_selected_pathways=vapply(unique_genes,
    function(s) sum(summary$symbol==s & summary$shared_ge3), integer(1))), "shared_gene_pathway_counts.csv")
  write(data.frame(input=inputs, md5=unname(tools::md5sum(inputs))), "input_checksums.csv")
  write(data.frame(color_scale="linear_log2FC", lower_limit=-color_limit,
    upper_limit=color_limit, shared_across_studies=TRUE, shared_across_pathways=TRUE,
    display_clipping_only=TRUE), "plot_settings.csv")
  if (draw_plots) {
    for (p in selected) {
      genes <- summary$symbol[summary$pathway==p & summary$shared_ge3]
      if (!length(genes)) next
      d <- long[long$pathway==p & long$symbol %in% genes, ]
      d$symbol <- factor(d$symbol, levels=rev(sort(genes)))
      d$study <- factor(d$study, levels=studies, labels=c("OIPN 160543", "NC 246156", "CCI 212311"))
      # Clipping is display only; original effect estimates remain in CSVs.
      d$display_log2FC <- pmax(-color_limit, pmin(color_limit, d$log2FC))
      d$mark <- ifelse(!d$rank_available, "NA", ifelse(!is.na(d$gene_FDR) & d$gene_FDR<.05, "*", ""))
      plot <- ggplot2::ggplot(d, ggplot2::aes(study, symbol, fill=display_log2FC)) +
        ggplot2::geom_tile(color="white", linewidth=.25) + ggplot2::geom_text(ggplot2::aes(label=mark), size=3) +
        ggplot2::scale_fill_gradient2(low="#2166AC", mid="#F7F7F7", high="#B2182B", midpoint=0,
          limits=c(-color_limit,color_limit), na.value="#BDBDBD", name=paste0("log2FC\nclipped at +/-", color_limit)) +
        ggplot2::labs(title=gsub("_", " ", sub("^HALLMARK_", "", p)),
          subtitle="Leading-edge membership in ≥3 studies; * gene FDR <0.05; NA unavailable to GSEA",
          x=NULL, y=NULL, caption="Neuropathy minus control. No row scaling. Compare directions; assays differ.") +
        ggplot2::theme_minimal(base_size=10) + ggplot2::theme(panel.grid=ggplot2::element_blank())
      height <- max(4, 2 + .18*length(genes)); stem <- paste0(sub("^HALLMARK_", "", p), "_shared_ge3_log2FC")
      ggplot2::ggsave(file.path(out, paste0(stem, ".pdf")), plot, width=9, height=height, limitsize=FALSE)
      ggplot2::ggsave(file.path(out, paste0(stem, ".png")), plot, width=9, height=height, dpi=300, bg="white", limitsize=FALSE)
    }
  }
  writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
  message("Shared Hallmark leading-edge extraction complete: ", out)
  invisible(totals)
}
if (!exists("LE_AUTORUN", inherits = FALSE) || isTRUE(LE_AUTORUN)) {
  shared_leading_edge_result <- run_shared_leading_edge(
    draw_plots=tolower(Sys.getenv("LE_TABLES_ONLY", "false")) != "true")
}
