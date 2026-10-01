#!/usr/bin/env Rscript
source("scripts/gobp/helpers.R", local = TRUE)
run_gobp_leading_edge <- function(draw_plots = TRUE) {
  color_limit <- suppressWarnings(as.numeric(Sys.getenv("LE_COLOR_LIMIT", "6")))
  if (length(color_limit) != 1L || !is.finite(color_limit) || color_limit <= 0) stop("LE_COLOR_LIMIT must be finite and positive.")
  file <- file.path(gobp_root, "comparison/all_GO_BP_pathways_classified.csv")
  selected <- gobp_read(file, c("pathway", "shared_positive", "shared_negative", "oipn_direction", "physical_direction"))
  gsea <- genes <- list(); inputs <- file
  for (study in gobp_studies) {
    gp <- file.path(gobp_root, "GSEA", study, "GSEA_GO_BP_results.csv")
    ep <- file.path(gobp_root, "GSEA", study, "gene_feature_evidence.csv")
    gsea[[study]] <- gobp_read(gp, c("pathway", "NES", "padj", "leadingEdge"))
    genes[[study]] <- gobp_read(ep, c("symbol", "feature_id", "statistic", "log2FC", "gene_FDR", "mapping_source"))
    genes[[study]]$feature_id <- as.character(genes[[study]]$feature_id)
    original <- gobp_gene_table(study)
    if (!isTRUE(all.equal(genes[[study]], original, check.attributes = FALSE, tolerance = 1e-8))) stop("Archived gene evidence differs: ", study)
    inputs <- c(inputs, gp, ep)
  }
  # Schemas exist even when no pathways/genes qualify.
  empty_summary <- data.frame(pathway = character(), symbol = character(), n_ranked = integer(), n_leading_edge = integer(),
    shared_all3 = logical(), direction_matches_all3 = logical(), opposite_both_physical = logical(),
    OIPN_matches_pathway_direction = logical(), n_gene_FDR_lt_0_05 = integer(),
    priority_OIPN_plus_physical = logical(), strict_significant_all3 = logical())
  for (study in gobp_studies) {
    empty_summary[[paste0(study, "_LE")]] <- logical()
    empty_summary[[paste0(study, "_log2FC")]] <- numeric()
    empty_summary[[paste0(study, "_gene_FDR")]] <- numeric()
  }
  empty_long <- data.frame(pathway = character(), symbol = character(), study = character(), rank_available = logical(),
    leading_edge = logical(), status = character(), feature_id = character(), mapping_source = character(),
    rank_statistic = numeric(), log2FC = numeric(), gene_FDR = numeric(), pathway_NES = numeric(), pathway_FDR = numeric())
  for (category in c("shared_positive", "shared_negative", "divergent")) {
    ids <- if (category == "divergent") selected$pathway[selected$oipn_direction | selected$physical_direction] else selected$pathway[selected[[category]]]
    rows <- long_rows <- list()
    for (p in ids) {
      gi <- lapply(gsea, function(g) match(p, g$pathway))
      if (any(vapply(gi, is.na, logical(1)))) stop("Selected pathway missing from GSEA: ", p)
      nes <- vapply(seq_along(gsea), function(i) gsea[[i]]$NES[gi[[i]]], numeric(1))
      pfdr <- vapply(seq_along(gsea), function(i) gsea[[i]]$padj[gi[[i]]], numeric(1))
      if (any(!is.finite(nes)) || any(!is.finite(pfdr))) stop("Unavailable selected pathway: ", p)
      if (category != "divergent" && (any(pfdr >= gobp_fdr) || any(if (category == "shared_positive") nes <= 0 else nes >= 0))) stop("Invalid shared selection.")
      if (category == "divergent" && !((nes[1L] > 0 && all(nes[2:3] < 0)) || (nes[1L] < 0 && all(nes[2:3] > 0)))) stop("Invalid opposite-direction selection.")
      edges <- lapply(seq_along(gsea), function(i) {
        z <- gsea[[i]]$leadingEdge[gi[[i]]]
        if (is.na(z) || !nzchar(z)) return(character())
        unique(trimws(strsplit(z, ";", fixed = TRUE)[[1L]]))
      })
      for (i in seq_along(edges)) if (!all(edges[[i]] %in% genes[[i]]$symbol)) stop("Leading-edge gene absent from rank.")
      for (symbol in sort(unique(unlist(edges)), method = "radix")) {
        pos <- vapply(genes, function(g) match(symbol, g$symbol), integer(1))
        ranked <- !is.na(pos); member <- vapply(edges, function(e) symbol %in% e, logical(1))
        lfc <- gfdr <- rep(NA_real_, 3L)
        for (i in seq_along(gobp_studies)) {
          g <- genes[[i]]
          if (ranked[i]) { lfc[i] <- g$log2FC[pos[i]]; gfdr[i] <- g$gene_FDR[pos[i]] }
          long_rows[[length(long_rows) + 1L]] <- data.frame(pathway = p, symbol = symbol, study = gobp_studies[i],
            rank_available = ranked[i], leading_edge = if (ranked[i]) member[i] else NA,
            status = if (!ranked[i]) "unavailable_to_GSEA" else if (member[i]) "leading_edge" else "ranked_not_leading_edge",
            feature_id = if (ranked[i]) g$feature_id[pos[i]] else NA_character_,
            mapping_source = if (ranked[i]) g$mapping_source[pos[i]] else NA_character_,
            rank_statistic = if (ranked[i]) g$statistic[pos[i]] else NA_real_, log2FC = lfc[i], gene_FDR = gfdr[i], pathway_NES = nes[i], pathway_FDR = pfdr[i])
        }
        sig <- !is.na(gfdr) & gfdr < gobp_fdr
        same_direction <- all(ranked) && all(if (category == "shared_negative") lfc < 0 else lfc > 0)
        opposite <- all(ranked) && ((lfc[1L] > 0 && all(lfc[2:3] < 0)) || (lfc[1L] < 0 && all(lfc[2:3] > 0)))
        candidate <- if (category == "divergent") opposite else all(member) && same_direction
        r <- data.frame(pathway = p, symbol = symbol, n_ranked = sum(ranked), n_leading_edge = sum(member),
          shared_all3 = all(member), direction_matches_all3 = same_direction, opposite_both_physical = opposite,
          OIPN_matches_pathway_direction = ranked[1L] && lfc[1L] * nes[1L] > 0,
          n_gene_FDR_lt_0_05 = sum(sig), priority_OIPN_plus_physical = candidate && sig[1L] && any(sig[2:3]),
          strict_significant_all3 = candidate && all(sig))
        for (i in seq_along(gobp_studies)) {
          r[[paste0(gobp_studies[i], "_LE")]] <- if (ranked[i]) member[i] else NA
          r[[paste0(gobp_studies[i], "_log2FC")]] <- lfc[i]
          r[[paste0(gobp_studies[i], "_gene_FDR")]] <- gfdr[i]
        }
        rows[[length(rows) + 1L]] <- r
      }
    }
    summary <- if (length(rows)) do.call(rbind, rows) else empty_summary
    long <- if (length(long_rows)) do.call(rbind, long_rows) else empty_long
    priority <- summary[summary$priority_OIPN_plus_physical, , drop = FALSE]
    strict <- summary[summary$strict_significant_all3, , drop = FALSE]
    out <- file.path(gobp_root, "leading_edge", category)
    gobp_write(selected[match(ids, selected$pathway), , drop = FALSE], file.path(out, "selected_pathways.csv"))
    gobp_write(long, file.path(out, "gene_evidence_long.csv"))
    gobp_write(summary, file.path(out, "gene_membership_summary.csv"))
    gobp_write(summary[summary$shared_all3, , drop = FALSE], file.path(out, "shared_all3.csv"))
    gobp_write(summary[summary$opposite_both_physical, , drop = FALSE], file.path(out, "opposite_direction_genes.csv"))
    gobp_write(priority, file.path(out, "priority_pathway_memberships.csv"))
    # Deduplicate for STRING and counts; retain all pathway memberships separately.
    unique_priority <- priority[!duplicated(priority$symbol), , drop = FALSE]
    gobp_write(unique_priority, file.path(out, "priority_OIPN_plus_physical.csv"))
    gobp_write(strict[!duplicated(strict$symbol), , drop = FALSE], file.path(out, "strict_significant_all3.csv"))
    writeLines(sort(unique_priority$symbol), file.path(out, "priority_genes_STRING.txt"))
    writeLines(sort(unique(strict$symbol)), file.path(out, "strict_genes_STRING.txt"))
    totals <- if (length(ids)) do.call(rbind, lapply(ids, function(p) {
      d <- summary[summary$pathway == p, , drop = FALSE]
      data.frame(pathway = p, union_LE_genes = nrow(d), shared_all3 = sum(d$shared_all3),
        opposite_both_physical = sum(d$opposite_both_physical), priority_genes = sum(d$priority_OIPN_plus_physical), strict_genes = sum(d$strict_significant_all3))
    })) else data.frame(pathway = character(), union_LE_genes = integer(), shared_all3 = integer(), opposite_both_physical = integer(), priority_genes = integer(), strict_genes = integer())
    gobp_write(totals, file.path(out, "pathway_overlap_summary.csv"))
    gobp_write(data.frame(metric = c("selected_pathways", "union_pathway_gene_rows", "study_evidence_rows", "shared_unique_genes", "opposite_unique_genes", "priority_unique_genes", "strict_unique_genes"),
      count = c(length(ids), nrow(summary), nrow(long), length(unique(summary$symbol[summary$shared_all3])),
        length(unique(summary$symbol[summary$opposite_both_physical])), nrow(unique_priority), length(unique(strict$symbol)))), file.path(out, "selection_summary.csv"))
    gobp_write(data.frame(color_scale = "linear_log2FC", lower_limit = -color_limit, upper_limit = color_limit,
      display_clipping_only = TRUE), file.path(out, "plot_settings.csv"))
    gobp_audit(c(inputs, "scripts/gobp/04_leading_edge_GO_BP.R"), out)
    writeLines(c("Shared: leading-edge membership in the same GO term in all three; preserve gene effects and FDR separately.",
      "Divergent: union of the three leading edges; membership in every study not required, exactly as the historical divergent workflow.",
      "Priority: same-sign shared or opposite-sign divergent gene effects, gene FDR <0.05 in OIPN and NC or CCI.",
      "Original whole-study DE FDR is retained; no FDR recalculation on selected genes.",
      "GO term overlap makes repeated memberships dependent, not repeated independent validation.",
      "STRING table has one row per gene; full pathway memberships remain in a separate table."), file.path(out, "analysis_notes.txt"))
    if (draw_plots) {
      fig <- file.path(out, "figures"); dir.create(fig, recursive = TRUE, showWarnings = FALSE)
      # This directory is wholly generated by this stage; clear old plot families.
      unlink(list.files(fig, pattern = "^(GO_BP_[0-9]+_page_|GO_BP_[0-9]+_empty[.]txt)", full.names = TRUE))
      mapping <- data.frame(stem = sprintf("GO_BP_%05d", match(ids, selected$pathway)), pathway = ids)
      gobp_write(mapping, file.path(out, "figure_pathway_index.csv"))
      for (j in seq_along(ids)) {
        p <- ids[j]
        plot_genes <- summary$symbol[summary$pathway == p & if (category == "divergent") summary$opposite_both_physical else summary$shared_all3]
        d <- long[long$pathway == p & long$symbol %in% plot_genes, , drop = FALSE]
        d$pathway <- d$symbol; d$mark <- ifelse(!is.na(d$gene_FDR) & d$gene_FDR < gobp_fdr, "*", "")
        gobp_plot_pages(d, sort(plot_genes), fig, mapping$stem[j], paste(gobp_pretty(p), "gene log2FC"), "log2FC", color_limit)
      }
    }
    message("GO:BP leading edge complete: ", category, "; ", nrow(unique_priority), " priority genes")
  }
  invisible(TRUE)
}
if (!exists("GOBP_AUTORUN", inherits = FALSE) || isTRUE(GOBP_AUTORUN)) run_gobp_leading_edge(tolower(Sys.getenv("GOBP_TABLES_ONLY", "false")) != "true")
