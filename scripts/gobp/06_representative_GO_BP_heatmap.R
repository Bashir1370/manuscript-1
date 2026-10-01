#!/usr/bin/env Rscript
# Run from repository root. Uses completed tables; never reruns GSEA/GSVA/DE.
source("scripts/gobp/helpers.R", local = TRUE)

gobp_representative_selection <- function() {
  ids <- c("RESPONSE_TO_WOUNDING", "CYTOKINE_MEDIATED_SIGNALING_PATHWAY",
    "LEUKOCYTE_MIGRATION", "PHAGOCYTOSIS",
    "SIGNAL_TRANSDUCTION_IN_RESPONSE_TO_DNA_DAMAGE", "TISSUE_REMODELING",
    "OXIDATIVE_PHOSPHORYLATION", "TRICARBOXYLIC_ACID_CYCLE",
    "REGULATION_OF_TRANS_SYNAPTIC_SIGNALING", "POTASSIUM_ION_TRANSPORT",
    "STEROL_BIOSYNTHETIC_PROCESS")
  data.frame(display_order = seq_along(ids), pathway = paste0("GOBP_", ids),
    category = c(rep("shared_positive", 6), rep("shared_negative", 5)),
    label = c("Response to wounding", "Cytokine-mediated signaling", "Leukocyte migration",
      "Phagocytosis", "Signaling in response to DNA damage", "Tissue remodeling",
      "Oxidative phosphorylation", "Tricarboxylic acid cycle",
      "Regulation of trans-synaptic signaling", "Potassium ion transport", "Sterol biosynthesis"),
    selection_basis = "Author-selected representative; not an automated redundancy cluster")
}

gobp_rep_panel <- function(values, fdr, labels, title, legend_title, limit,
                           left_margin, gap_after = integer()) {
  stopifnot(identical(dim(values), dim(fdr)), nrow(values) == length(labels),
    ncol(values) == 3L, is.finite(limit), limit > 0)
  pal <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(201)
  n <- nrow(values); ys <- rev(seq_len(n))
  graphics::par(mar = c(3.6, left_margin, 3.3, 1), xpd = NA)
  graphics::plot.new()
  graphics::plot.window(xlim = c(.5, 5), ylim = c(.5, n + .5), xaxs = "i", yaxs = "i")
  for (i in seq_len(n)) for (j in 1:3) {
    v <- values[i, j]
    col <- if (is.finite(v)) pal[1L + round(200 * (max(-limit, min(limit, v)) + limit) / (2 * limit))] else "#BDBDBD"
    graphics::rect(j - .48, ys[i] - .48, j + .48, ys[i] + .48, col = col, border = "white")
    if (is.finite(v) && is.finite(fdr[i, j]) && fdr[i, j] < .05)
      graphics::text(j, ys[i], "*", cex = 1.05, col = if (abs(v) >= .65 * limit) "white" else "black")
  }
  graphics::axis(1, at = 1:3, labels = c("OIPN", "NC", "CCI"), tick = FALSE, cex.axis = .85)
  graphics::axis(2, at = ys, labels = labels, las = 2, tick = FALSE, cex.axis = .76)
  for (i in gap_after[gap_after > 0 & gap_after < n])
    graphics::abline(h = ys[i] - .5, col = "#333333", lwd = 1.1)
  lo <- max(.8, n / 2 - 2); hi <- min(n + .2, n / 2 + 2)
  breaks <- seq(lo, hi, length.out = 202)
  for (k in 1:201) graphics::rect(3.85, breaks[k], 4.08, breaks[k + 1], col = pal[k], border = NA)
  graphics::text(4.15, c(lo, (lo + hi) / 2, hi), c(-limit, 0, limit), adj = 0, cex = .7)
  graphics::text(3.95, hi + .6, legend_title, cex = .7)
  graphics::title(main = title, cex.main = 1)
}

run_gobp_representative <- function(draw_plots = TRUE, input_root = gobp_root,
                                  out = file.path(gobp_root, "representative_11")) {
  selection <- gobp_representative_selection()
  classification_file <- file.path(input_root, "comparison/all_GO_BP_pathways_classified.csv")
  evidence_files <- file.path(input_root, "leading_edge", c("shared_positive", "shared_negative"), "gene_membership_summary.csv")
  wide <- gobp_read(classification_file, c("pathway", "shared_positive", "shared_negative",
    paste0(gobp_studies, "_NES"), paste0(gobp_studies, "_FDR")))
  if (anyDuplicated(wide$pathway)) stop("Duplicate pathway in comparison table.")
  idx <- match(selection$pathway, wide$pathway)
  if (anyNA(idx)) stop("A representative pathway is missing; no substitute is selected automatically.")
  selected <- wide[idx, , drop = FALSE]
  nes <- as.matrix(selected[paste0(gobp_studies, "_NES")])
  pfdr <- as.matrix(selected[paste0(gobp_studies, "_FDR")])
  expected <- ifelse(selection$category == "shared_positive", 1, -1)
  if (any(!is.finite(nes)) || any(!is.finite(pfdr)) || any(pfdr < 0 | pfdr >= .05) ||
      any(sign(nes) != expected) ||
      !all(vapply(seq_len(nrow(selection)), function(i) isTRUE(selected[[selection$category[i]]][i]), logical(1))))
    stop("Representative selection no longer meets shared same-direction FDR <0.05 in all three; review the selection.")
  required <- c("pathway", "symbol", "shared_all3", "direction_matches_all3",
    "priority_OIPN_plus_physical", "strict_significant_all3",
    paste0(gobp_studies, "_LE"), paste0(gobp_studies, "_log2FC"), paste0(gobp_studies, "_gene_FDR"))
  evidence <- do.call(rbind, lapply(evidence_files, gobp_read, columns = required))
  if (anyDuplicated(evidence[c("pathway", "symbol")])) stop("Duplicate pathway/gene evidence.")
  ev <- evidence[evidence$pathway %in% selection$pathway, , drop = FALSE]
  if (!all(selection$pathway %in% ev$pathway)) stop("Missing leading-edge evidence for selected pathways.")
  lfc_cols <- paste0(gobp_studies, "_log2FC"); fdr_cols <- paste0(gobp_studies, "_gene_FDR")
  lfc <- as.matrix(ev[lfc_cols]); gfdr <- as.matrix(ev[fdr_cols])
  finite <- apply(is.finite(lfc), 1, all)
  le_all <- Reduce(`&`, lapply(ev[paste0(gobp_studies, "_LE")], function(x) !is.na(x) & x))
  dir_ok <- finite & apply(sign(lfc) == ifelse(selection$category[match(ev$pathway, selection$pathway)] == "shared_positive", 1, -1), 1, all)
  sig <- is.finite(gfdr) & gfdr >= 0 & gfdr < .05
  priority <- le_all & dir_ok & sig[, 1] & (sig[, 2] | sig[, 3])
  strict <- le_all & dir_ok & apply(sig, 1, all)
  if (anyNA(ev$priority_OIPN_plus_physical) || anyNA(ev$strict_significant_all3) ||
      any(priority != ev$priority_OIPN_plus_physical) || any(strict != ev$strict_significant_all3) ||
      any(le_all != ev$shared_all3)) stop("Gene flags disagree with original evidence; rerun stage 04.")
  # Original gene statistics must agree across pathway memberships.
  for (symbol in unique(ev$symbol)) for (col in c(lfc_cols, fdr_cols)) {
    vals <- ev[ev$symbol == symbol, col]; vals <- unique(vals[!is.na(vals)])
    if (length(vals) > 1L) stop("Inconsistent gene statistics across pathways: ", symbol)
  }
  priorities <- sort(unique(ev$symbol[priority]), method = "radix")
  focus <- c("Cdk1", "Cdkn1a")
  if (!all(focus %in% ev$symbol)) stop("Missing requested Cdk1/Cdkn1a context evidence.")
  symbols <- unique(c(priorities, focus))
  genes <- ev[match(symbols, ev$symbol), c("symbol", lfc_cols, fdr_cols), drop = FALSE]
  genes$priority_in_selected_pathways <- genes$symbol %in% priorities
  genes$display_role <- ifelse(genes$priority_in_selected_pathways, "selected-pathway priority", "context only")
  genes$priority_pathways <- vapply(symbols, function(s) paste(ev$pathway[ev$symbol == s & priority], collapse = ";"), character(1))
  genes$shared_LE_pathways <- vapply(symbols, function(s) paste(ev$pathway[ev$symbol == s & le_all], collapse = ";"), character(1))
  genes <- genes[order(-sign(genes[[lfc_cols[1]]]), genes$symbol, method = "radix"), , drop = FALSE]
  nes_limit <- max(.5, ceiling(max(abs(nes)) * 2) / 2)
  gene_limit <- suppressWarnings(as.numeric(Sys.getenv("LE_COLOR_LIMIT", "6")))
  if (length(gene_limit) != 1 || !is.finite(gene_limit) || gene_limit <= 0) stop("LE_COLOR_LIMIT must be a positive finite number.")
  gobp_write(cbind(selection, selected[setdiff(names(selected), "pathway")]), file.path(out, "selected_pathway_evidence.csv"))
  gobp_write(ev, file.path(out, "selected_gene_membership_evidence.csv"))
  gobp_write(genes, file.path(out, "displayed_gene_evidence.csv"))
  counts <- do.call(rbind, lapply(selection$pathway, function(p) {
    take <- ev$pathway == p
    data.frame(pathway = p, common_LE_all3 = sum(le_all[take]), priority_genes = sum(priority[take]), strict_all3_genes = sum(strict[take]))
  }))
  gobp_write(counts, file.path(out, "pathway_gene_counts.csv"))
  gobp_write(data.frame(panel = c("pathways", "genes"), value = c("NES", "log2FC"),
    color_limit = c(nes_limit, gene_limit), row_scaling = FALSE,
    clipped_cells = c(0, sum(abs(as.matrix(genes[lfc_cols])) > gene_limit, na.rm = TRUE))), file.path(out, "plot_settings.csv"))
  writeLines(c("Author-selected 11 representatives of the shared positive/negative GO:BP results; exploratory presentation selection.",
    "Not an automated clustering result and not 11 independent mechanisms. Full results remain in comparison/.",
    "Panel A: original NES; * original full-family pathway BH FDR <0.05. No row scaling.",
    "Panel B: unique priority genes within these selected pathways plus Cdk1/Cdkn1a context.",
    "Priority: same-path leading edge in all three, matching log2FC sign in all three, gene FDR <0.05 in OIPN plus NC or CCI.",
    "Context only: displayed for the author's question, does not meet the selected-pathway priority rule.",
    "Panel B *: original gene-level FDR <0.05; no FDR recalculation. Linear color clipping is display-only; CSVs retain original values.",
    "Cdk1 can be priority in other GO terms while context-only in this selection. Pathway significance does not imply significance of every member.",
    "No GSVA significance claim, cell proliferation measurement, neuronal localization, or causal inference."), file.path(out, "analysis_notes.txt"))
  gobp_audit(c(classification_file, evidence_files, "scripts/gobp/06_representative_GO_BP_heatmap.R"), out)
  if (draw_plots) {
    gm <- as.matrix(genes[lfc_cols]); fm <- as.matrix(genes[fdr_cols])
    labels <- ifelse(genes$priority_in_selected_pathways, genes$symbol, paste0(genes$symbol, " [context]"))
    draw <- function() {
      graphics::layout(matrix(1:2, nrow = 1), widths = c(1.35, 1))
      gobp_rep_panel(nes, pfdr, selection$label, "A. Shared GO:BP representatives", "NES", nes_limit, 17, 6)
      gobp_rep_panel(gm, fm, labels, "B. Selected gene evidence", "log2FC", gene_limit, 8)
      graphics::mtext("* FDR <0.05 within study; panel A: pathway FDR; panel B: gene FDR. [context] = not a priority in these 11 pathways.",
        side = 1, outer = TRUE, line = 1, cex = .7)
    }
    save_device <- function(format) {
      path <- file.path(out, paste0("Fig_GO_BP_11_pathways_and_genes.", format))
      height <- max(8, 2.5 + .25 * nrow(genes))
      if (format == "pdf") grDevices::pdf(path, width = 17, height = height, useDingbats = FALSE)
      else grDevices::png(path, width = 17, height = height, units = "in", res = 300, bg = "white")
      on.exit(grDevices::dev.off(), add = TRUE)
      graphics::par(oma = c(3, 0, 0, 0)); draw()
    }
    save_device("pdf"); save_device("png")
  }
  message("Representative GO:BP evidence complete: ", out, "; ", length(priorities), " unique priority genes plus context.")
  invisible(list(pathways = selection, genes = genes, counts = counts))
}
if (!exists("GOBP_AUTORUN", inherits = FALSE) || isTRUE(GOBP_AUTORUN))
  run_gobp_representative(tolower(Sys.getenv("GOBP_TABLES_ONLY", "false")) != "true")
