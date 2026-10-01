#!/usr/bin/env Rscript
# Optional presentation stage; run from repository root after stages 02 and 04.
# Uses archived statistics only; does not refit DE, GSEA or GSVA.
source("scripts/gobp/helpers.R", local = TRUE)

gobp_divergent_selection <- function() {
  ids <- c("COLLAGEN_FIBRIL_ORGANIZATION", "EXTRACELLULAR_MATRIX_ASSEMBLY",
    "WNT_SIGNALING_PATHWAY", "TIGHT_JUNCTION_ORGANIZATION",
    "SCHWANN_CELL_DIFFERENTIATION", "PERIPHERAL_NERVOUS_SYSTEM_AXON_ENSHEATHMENT",
    "RENAL_SYSTEM_DEVELOPMENT", "EXOCYTOSIS", "CELLULAR_RESPONSE_TO_UNFOLDED_PROTEIN")
  data.frame(display_order = seq_along(ids), pathway = paste0("GOBP_", ids),
    category = c(rep("OIPN_positive", 7), rep("OIPN_negative", 2)),
    label = c("Collagen fibril organization", "Extracellular matrix assembly",
      "WNT signaling", "Tight junction organization", "Schwann cell differentiation",
      "Peripheral nervous system axon ensheathment", "Renal system development (GO annotation)",
      "Exocytosis", "Cellular response to unfolded protein"),
    selection_basis = "Author-selected exploratory representative; not an automated redundancy cluster")
}

gobp_divergent_gene_table <- function(members, rule, scope) {
  studies <- gobp_studies; lcols <- paste0(studies, "_log2FC"); fcols <- paste0(studies, "_gene_FDR")
  symbols <- sort(unique(members$symbol), method = "radix")
  g <- members[match(symbols, members$symbol), c("symbol", lcols, fcols), drop = FALSE]
  g$OIPN_direction <- ifelse(g[[lcols[1]]] > 0, "positive", "negative")
  g$strict_gene_FDR_all3 <- if (nrow(g)) apply(as.matrix(g[fcols]), 1,
    function(x) all(is.finite(x) & x >= 0 & x < .05)) else logical()
  g$n_priority_pathways <- vapply(symbols, function(s) sum(members$symbol == s), integer(1))
  g$priority_pathways <- vapply(symbols, function(s) paste(sort(unique(members$pathway[members$symbol == s])), collapse = ";"), character(1))
  g$supporting_physical_studies <- vapply(symbols, function(s) {
    m <- members[members$symbol == s, , drop = FALSE]
    paste(studies[2:3][vapply(2:3, function(i) {
      ok <- is.finite(m[[fcols[i]]]) & m[[fcols[i]]] < .05
      if (rule == "revised_paired_LE_and_FDR") ok <- ok & !is.na(m[[paste0(studies[i], "_LE")]]) & m[[paste0(studies[i], "_LE")]]
      any(ok)
    }, logical(1))], collapse = ";")
  }, character(1))
  g$selection_rule <- rep(rule, nrow(g)); g$pathway_scope <- rep(scope, nrow(g))
  rownames(g) <- NULL; g
}

run_gobp_divergent_representative <- function(draw_plots = TRUE, input_root = gobp_root,
    out = file.path(gobp_root, "divergent_representative_9")) {
  cf <- file.path(input_root, "comparison/all_GO_BP_pathways_classified.csv")
  mf <- file.path(input_root, "leading_edge/divergent/priority_pathway_memberships.csv")
  af <- file.path(input_root, "leading_edge/divergent/priority_OIPN_plus_physical.csv")
  lcols <- paste0(gobp_studies, "_log2FC"); fcols <- paste0(gobp_studies, "_gene_FDR")
  lecols <- paste0(gobp_studies, "_LE"); ncols <- paste0(gobp_studies, "_NES"); pcols <- paste0(gobp_studies, "_FDR")
  wide <- gobp_read(cf, c("pathway", "oipn_direction", "physical_direction", ncols, pcols))
  m <- gobp_read(mf, c("pathway", "symbol", "priority_OIPN_plus_physical", "strict_significant_all3",
    "opposite_both_physical", "OIPN_matches_pathway_direction", lcols, fcols, lecols))
  archived <- gobp_read(af, c("symbol", lcols, fcols))
  if (anyDuplicated(wide$pathway) || anyDuplicated(m[c("pathway", "symbol")]) || anyDuplicated(archived$symbol))
    stop("Duplicate pathway or gene key.")
  nes <- as.matrix(wide[ncols]); pfdr <- as.matrix(wide[pcols])
  if (any(is.finite(pfdr) & (pfdr < 0 | pfdr > 1))) stop("Invalid pathway FDR.")
  available <- apply(is.finite(nes) & is.finite(pfdr), 1, all)
  pos <- available & nes[, 1] > 0 & nes[, 2] < 0 & nes[, 3] < 0
  neg <- available & nes[, 1] < 0 & nes[, 2] > 0 & nes[, 3] > 0
  if (anyNA(wide$oipn_direction) || anyNA(wide$physical_direction) ||
      any(pos != wide$oipn_direction) || any(neg != wide$physical_direction)) stop("Pathway direction flags disagree with NES.")
  mi <- match(m$pathway, wide$pathway)
  if (anyNA(mi) || any(!(pos | neg)[mi])) stop("Priority membership outside the divergent pathway selection.")
  lfc <- as.matrix(m[lcols]); gfdr <- as.matrix(m[fcols]); le <- as.matrix(m[lecols]); le[is.na(le)] <- FALSE
  if (any(is.finite(gfdr) & (gfdr < 0 | gfdr > 1))) stop("Invalid gene FDR.")
  sig <- is.finite(gfdr) & gfdr >= 0 & gfdr < .05
  opp <- apply(is.finite(lfc), 1, all) & lfc[, 1] * lfc[, 2] < 0 & lfc[, 1] * lfc[, 3] < 0
  old <- opp & sig[, 1] & (sig[, 2] | sig[, 3])
  strict_old <- opp & apply(sig, 1, all)
  aligned <- lfc[, 1] * nes[mi, 1] > 0
  if (anyNA(old) || !all(old) || anyNA(m$priority_OIPN_plus_physical) || !all(m$priority_OIPN_plus_physical) ||
      anyNA(m$strict_significant_all3) || any(strict_old != m$strict_significant_all3) ||
      anyNA(m$opposite_both_physical) || !all(m$opposite_both_physical) ||
      anyNA(m$OIPN_matches_pathway_direction) || any(aligned != m$OIPN_matches_pathway_direction))
    stop("Archived priority flags disagree with original gene evidence; review stage 04.")
  for (s in unique(m$symbol)) for (col in c(lcols, fcols)) {
    v <- m[m$symbol == s, col]; v <- unique(v[!is.na(v)])
    if (length(v) > 1L) stop("Inconsistent gene evidence across pathways: ", s)
  }
  if (!setequal(unique(m$symbol), archived$symbol)) stop("Archived unique priorities disagree with memberships.")
  gcheck <- m[match(archived$symbol, m$symbol), c(lcols, fcols), drop = FALSE]
  if (!isTRUE(all.equal(gcheck, archived[c(lcols, fcols)], check.attributes = FALSE, tolerance = 1e-12)))
    stop("Archived unique-gene statistics disagree with memberships.")
  # Original 457 pathways are retained. OIPN-significant terms are a declared follow-up subset.
  strong <- wide[(pos | neg) & pfdr[, 1] < .05, , drop = FALSE]
  strong <- strong[order(-sign(strong[[ncols[1]]]), strong[[pcols[1]]], strong$pathway), , drop = FALSE]
  selection <- gobp_divergent_selection(); si <- match(selection$pathway, wide$pathway)
  if (anyNA(si) || !all(selection$pathway %in% strong$pathway))
    stop("A representative no longer meets opposite NES and OIPN FDR <0.05; no substitute is selected.")
  selected <- wide[si, , drop = FALSE]
  if (any(sign(selected[[ncols[1]]]) != ifelse(selection$category == "OIPN_positive", 1, -1)))
    stop("Representative pathway direction changed.")
  paired <- le & sig
  final_rule <- opp & aligned & paired[, 1] & (paired[, 2] | paired[, 3])
  strong_rows <- m$pathway %in% strong$pathway
  display_rows <- m$pathway %in% selection$pathway
  # All final genes from the 24-pathway subset must be represented within the nine terms.
  if (!setequal(m$symbol[final_rule & strong_rows], m$symbol[final_rule & display_rows]))
    stop("Nine representatives omit final genes from the OIPN-significant subset; review selection.")
  all_genes <- gobp_divergent_gene_table(m, "archived_opposite_FC_and_OIPN_plus_physical_FDR", "all_opposite_NES_pathways")
  final_genes <- gobp_divergent_gene_table(m[final_rule & display_rows, , drop = FALSE],
    "revised_paired_LE_and_FDR", "nine_OIPN_significant_representatives")
  m$pathway_OIPN_significant <- strong_rows; m$pathway_selected_for_display <- display_rows
  m$revised_paired_LE_priority <- final_rule; m$final_priority <- final_rule & display_rows
  limit <- suppressWarnings(as.numeric(Sys.getenv("LE_COLOR_LIMIT", "6")))
  if (length(limit) != 1L || !is.finite(limit) || limit <= 0) stop("LE_COLOR_LIMIT must be positive and finite.")
  gobp_write(all_genes, file.path(out, "priority_all_divergent_genes.csv"))
  gobp_write(final_genes, file.path(out, "priority_final_genes.csv"))
  gobp_write(strong, file.path(out, "OIPN_significant_opposite_pathways.csv"))
  gobp_write(cbind(selection, selected[setdiff(names(selected), "pathway")]), file.path(out, "selected_pathway_evidence.csv"))
  gobp_write(m, file.path(out, "priority_membership_audit.csv"))
  counts <- data.frame(pathway = selection$pathway,
    final_priority_genes = vapply(selection$pathway, function(p) sum(m$pathway == p & m$final_priority), integer(1)))
  gobp_write(counts, file.path(out, "pathway_gene_counts.csv"))
  summary <- data.frame(metric = c("all_opposite_pathways", "OIPN_significant_opposite_pathways", "representative_pathways",
    "archived_priority_unique_genes", "OIPN_significant_original_priority_genes", "final_priority_unique_genes"),
    count = c(sum(pos | neg), nrow(strong), nrow(selection), nrow(all_genes), length(unique(m$symbol[strong_rows & aligned])), nrow(final_genes)))
  gobp_write(summary, file.path(out, "selection_summary.csv"))
  snes <- as.matrix(selected[ncols]); sfdr <- as.matrix(selected[pcols])
  nes_limit <- max(.5, ceiling(max(abs(snes)) * 2) / 2)
  gobp_write(data.frame(panel = c("pathways", "genes"), value = c("NES", "log2FC"),
    color_limit = c(nes_limit, limit), row_scaling = FALSE), file.path(out, "plot_settings.csv"))
  writeLines(c("Exploratory opposite-direction follow-up; does not establish OIPN specificity or a between-study interaction.",
    "Full pathway selection retains opposite NES signs without a significance filter. OIPN-significant follow-up requires original pathway FDR <0.05 in OIPN only.",
    "Nine author-selected representatives are not an automated GO redundancy clustering; GO terms overlap.",
    "Original gene list: all divergent pathways, opposite log2FC in OIPN versus both physical studies, gene FDR <0.05 in OIPN and NC or CCI; same-path LE support in both is not required.",
    "Final gene list: nine OIPN-significant representatives, gene direction matches each pathway in all three studies, same-path LE AND gene FDR <0.05 in OIPN and at least one SAME physical study.",
    "The remaining physical study requires available opposite log2FC to OIPN, but neither gene significance nor LE membership.",
    "Original whole-study pathway and gene FDR are retained. No selected-subset BH recalculation or statistical refit.",
    "The 39 original genes and five final genes have different pathway scopes and LE selection rules; no identity or independence is implied.",
    "Renal system development is a GO annotation. Its enrichment does not indicate kidney development in DRG.",
    "No supplementary pathway figure, causal mechanism or cell-type localization is claimed."), file.path(out, "analysis_notes.txt"))
  gobp_audit(c(cf, mf, af, "scripts/gobp/07_divergent_GO_BP_heatmap.R"), out)
  if (draw_plots) {
    # Reuse the tested stage-06 renderer without running the shared analysis.
    penv <- new.env(parent = environment()); penv$GOBP_AUTORUN <- FALSE
    source("scripts/gobp/06_representative_GO_BP_heatmap.R", local = penv)
    labels <- vapply(selection$label, function(x) paste(strwrap(x, width = 38), collapse = "\n"), character(1))
    draw <- function() {
      graphics::layout(matrix(1:2, nrow = 1), widths = c(1.35, 1))
      penv$gobp_rep_panel(snes, sfdr, labels, "A. OIPN-significant opposite GO:BP", "NES", nes_limit, 18, significance_color = "black")
      if (nrow(final_genes)) penv$gobp_rep_panel(as.matrix(final_genes[lcols]), as.matrix(final_genes[fcols]),
        final_genes$symbol, "B. Selected opposite gene evidence", "log2FC", limit, 8)
      else { graphics::plot.new(); graphics::text(.5, .5, "No genes meet the final rule.") }
      graphics::mtext("* Original within-study FDR <0.05. A: pathway FDR; B: gene FDR. Gene LE + FDR support: OIPN and the same NC or CCI study.",
        side = 1, outer = TRUE, line = 1, cex = .7)
    }
    for (format in c("pdf", "png")) {
      path <- file.path(out, paste0("Fig_GO_BP_9_opposite_pathways_and_genes.", format))
      if (format == "pdf") grDevices::pdf(path, width = 18, height = 9, useDingbats = FALSE)
      else grDevices::png(path, width = 18, height = 9, units = "in", res = 300, bg = "white")
      tryCatch({ graphics::par(oma = c(3, 0, 0, 0)); draw() }, finally = grDevices::dev.off())
    }
  }
  message("Opposite GO:BP representatives complete: ", out, "; ", nrow(all_genes), " original and ", nrow(final_genes), " final genes.")
  invisible(list(pathways = selection, all_genes = all_genes, final_genes = final_genes, summary = summary))
}
if (!exists("GOBP_AUTORUN", inherits = FALSE) || isTRUE(GOBP_AUTORUN))
  run_gobp_divergent_representative(tolower(Sys.getenv("GOBP_TABLES_ONLY", "false")) != "true")
