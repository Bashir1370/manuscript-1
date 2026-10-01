#!/usr/bin/env Rscript
# Optional evidence follow-up. Run from repository root after stages 06 and 07.
# Original DE/GSEA/GSVA statistics and selections are never refitted or replaced.
source("scripts/gobp/helpers.R", local = TRUE)

gobp08_validate_matrix <- function(x, meta, ids, nonnegative = FALSE) {
  if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x)) ||
      is.null(rownames(x)) || is.null(colnames(x)) || anyDuplicated(rownames(x)) ||
      anyDuplicated(colnames(x)) || !identical(colnames(x), as.character(meta$sample)) ||
      !all(ids %in% rownames(x)) || anyDuplicated(meta$sample) ||
      anyNA(meta$condition) || !all(meta$condition %in% c("Control", "Neuropathy")) ||
      any(table(factor(meta$condition, levels = c("Control", "Neuropathy"))) < 2L) ||
      (nonnegative && any(x < 0))) stop("Invalid selected evidence matrix or sample metadata.")
  x[match(ids, rownames(x)), , drop = FALSE]
}

gobp08_sensitivity <- function(x, meta, estimates, kind, study) {
  x <- gobp08_validate_matrix(x, meta, rownames(x))
  ctrl <- meta$condition == "Control"; inj <- !ctrl
  if (length(estimates) != nrow(x) || any(!is.finite(estimates))) stop("Missing original effect estimates.")
  baseline <- rowMeans(x[, inj, drop = FALSE]) - rowMeans(x[, ctrl, drop = FALSE])
  loo <- vapply(seq_len(ncol(x)), function(i) {
    keep <- seq_len(ncol(x)) != i
    rowMeans(x[, keep & inj, drop = FALSE]) - rowMeans(x[, keep & ctrl, drop = FALSE])
  }, numeric(nrow(x)))
  # Explicit dimensions preserve the one-row matrix case.
  loo <- matrix(loo, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
  preserves <- sign(loo) == sign(baseline)
  preserves[baseline == 0, ] <- FALSE
  long <- do.call(rbind, lapply(seq_len(nrow(x)), function(j) data.frame(
    kind = kind, study = study, id = rownames(x)[j], omitted_sample = meta$sample,
    omitted_condition = meta$condition, baseline_delta = unname(baseline[j]),
    leave_one_out_delta = unname(loo[j, ]), preserves_baseline_direction = unname(preserves[j, ]),
    original_estimate = unname(estimates[j]), stringsAsFactors = FALSE)))
  summary <- data.frame(kind = kind, study = study, id = rownames(x),
    n_control = sum(ctrl), n_neuropathy = sum(inj), control_mean = rowMeans(x[, ctrl, drop = FALSE]),
    neuropathy_mean = rowMeans(x[, inj, drop = FALSE]), baseline_delta = baseline,
    original_estimate = estimates, baseline_agrees_with_original = sign(baseline) == sign(estimates),
    leave_one_out_min = apply(loo, 1, min), leave_one_out_max = apply(loo, 1, max),
    n_direction_preserved = rowSums(preserves), n_omissions = ncol(x),
    all_omissions_preserve_direction = baseline != 0 & apply(preserves, 1, all),
    max_absolute_delta_change = apply(abs(loo - baseline), 1, max), stringsAsFactors = FALSE)
  rownames(summary) <- NULL; rownames(long) <- NULL
  list(long = long, summary = summary)
}

gobp08_selection <- function(input_root) {
  sf <- file.path(input_root, "representative_11/displayed_gene_evidence.csv")
  df <- file.path(input_root, "divergent_representative_9/priority_final_genes.csv")
  sp <- file.path(input_root, "representative_11/selected_pathway_evidence.csv")
  dp <- file.path(input_root, "divergent_representative_9/selected_pathway_evidence.csv")
  sm <- file.path(input_root, "representative_11/selected_gene_membership_evidence.csv")
  dm <- file.path(input_root, "divergent_representative_9/priority_membership_audit.csv")
  lc <- paste0(gobp_studies, "_log2FC"); fc <- paste0(gobp_studies, "_gene_FDR")
  nc <- paste0(gobp_studies, "_NES"); pc <- paste0(gobp_studies, "_FDR"); le <- paste0(gobp_studies, "_LE")
  shared <- gobp_read(sf, c("symbol", lc, fc)); opposite <- gobp_read(df, c("symbol", lc, fc))
  ps <- gobp_read(sp, c("pathway", "label", nc, pc)); pd <- gobp_read(dp, c("pathway", "label", nc, pc))
  if (anyDuplicated(shared$symbol) || anyDuplicated(opposite$symbol) ||
      anyDuplicated(ps$pathway) || anyDuplicated(pd$pathway)) stop("Duplicate selected key.")
  if (nrow(shared) != 37L || nrow(opposite) != 5L || nrow(ps) != 11L || nrow(pd) != 9L ||
      length(intersect(shared$symbol, opposite$symbol)) || length(intersect(ps$pathway, pd$pathway)))
    stop("Locked 37/5 gene or 11/9 pathway selections changed; review before extending stage 08.")
  ps$evidence_category <- "shared"; pd$evidence_category <- "opposite"
  pathways <- rbind(ps[c("pathway", "label", nc, pc, "evidence_category")],
                    pd[c("pathway", "label", nc, pc, "evidence_category")])
  for (category in c("shared", "opposite")) {
    paths <- pathways[pathways$evidence_category == category, , drop = FALSE]
    nes <- as.matrix(paths[nc]); fdr <- as.matrix(paths[pc])
    valid <- all(is.finite(nes)) && all(is.finite(fdr) & fdr >= 0 & fdr <= 1)
    if (!valid || (category == "shared" && (!all(fdr < .05) || !all(apply(nes > 0, 1, all) | apply(nes < 0, 1, all)))) ||
        (category == "opposite" && (!all(fdr[, 1] < .05) || !all(nes[, 1] * nes[, 2] < 0 & nes[, 1] * nes[, 3] < 0))))
      stop("Selected pathway direction/FDR no longer satisfies the declared rule.")
    members <- gobp_read(if (category == "shared") sm else dm, c("pathway", "symbol", lc, fc, le))
    members <- members[members$pathway %in% paths$pathway, , drop = FALSE]
    if (anyDuplicated(members[c("pathway", "symbol")])) stop("Duplicate pathway-gene membership.")
    lfc <- as.matrix(members[lc]); gf <- as.matrix(members[fc]); edges <- as.matrix(members[le]); edges[is.na(edges)] <- FALSE
    sig <- is.finite(gf) & gf >= 0 & gf < .05
    direction <- if (category == "shared") apply(lfc > 0, 1, all) | apply(lfc < 0, 1, all) else
      lfc[, 1] * lfc[, 2] < 0 & lfc[, 1] * lfc[, 3] < 0
    aligned <- lfc[, 1] * paths[[nc[1]]][match(members$pathway, paths$pathway)] > 0
    eligible <- apply(is.finite(lfc), 1, all) & direction & aligned & edges[, 1] & sig[, 1] &
      ((edges[, 2] & sig[, 2]) | (edges[, 3] & sig[, 3]))
    eligible[is.na(eligible)] <- FALSE
    genes <- if (category == "shared") shared else opposite
    if (!setequal(genes$symbol, members$symbol[eligible])) stop("Displayed genes disagree with rederived Leading-edge rule.")
    for (symbol in genes$symbol) for (col in c(lc, fc)) {
      vals <- members[members$symbol == symbol, col]
      target <- genes[genes$symbol == symbol, col]
      if (any(!is.finite(vals)) || !all(abs(vals - target) <= 1e-10 * pmax(1, abs(target))))
        stop("Selected gene statistics disagree with membership evidence: ", symbol)
    }
  }
  shared$evidence_category <- "shared"; opposite$evidence_category <- "opposite"
  genes <- rbind(shared[c("symbol", lc, fc, "evidence_category")], opposite[c("symbol", lc, fc, "evidence_category")])
  rownames(genes) <- NULL; rownames(pathways) <- NULL
  list(genes = genes, pathways = pathways, files = c(sf, df, sp, dp, sm, dm))
}

gobp08_plot <- function(expression, comparison, pathways, out) {
  plot_device <- function(stem, width, height, draw) for (format in c("pdf", "png")) {
    fn <- file.path(out, paste0(stem, ".", format))
    if (format == "pdf") grDevices::pdf(fn, width = width, height = height, useDingbats = FALSE)
    else grDevices::png(fn, width = width, height = height, units = "in", res = 300, bg = "white")
    tryCatch(draw(), finally = grDevices::dev.off())
  }
  plot_device("Fig_selected_42_gene_sample_expression", 15, 15, function() {
    graphics::layout(matrix(1:3, nrow = 1), widths = c(1.2, 1, 1))
    pal <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(201)
    for (study in gobp_studies) {
      d <- expression[expression$study == study, , drop = FALSE]; ids <- unique(d$symbol); samples <- unique(d$sample)
      graphics::par(mar = c(5, 5, 5, 1))
      graphics::plot.new(); graphics::plot.window(xlim = c(.5, length(samples) + .5), ylim = c(.5, length(ids) + .5), xaxs = "i", yaxs = "i")
      for (j in seq_along(ids)) for (i in seq_along(samples)) {
        z <- d$within_study_z[d$symbol == ids[j] & d$sample == samples[i]]
        color <- pal[1 + round((max(-3, min(3, z)) + 3) / 6 * 200)]
        graphics::rect(i - .48, length(ids) - j + .52, i + .48, length(ids) - j + 1.48, col = color, border = "white")
      }
      categories <- d$evidence_category[match(ids, d$symbol)]
      boundaries <- which(categories[-length(categories)] != categories[-1])
      if (length(boundaries)) graphics::abline(h = length(ids) - boundaries + .5, col = "#444444", lwd = 1.2)
      conditions <- d$condition[match(samples, d$sample)]
      boundaries <- which(conditions[-length(conditions)] != conditions[-1])
      if (length(boundaries)) graphics::abline(v = boundaries + .5, col = "#444444", lwd = 1.2)
      group_centers <- vapply(c("Control", "Neuropathy"), function(g) mean(which(conditions == g)), numeric(1))
      graphics::mtext(c("Control", "Neuropathy"), side = 3, at = group_centers, line = .4, cex = .75)
      graphics::axis(2, at = rev(seq_along(ids)), labels = ids, las = 2, tick = FALSE, cex.axis = .78)
      labs <- d$sample_label[match(samples, d$sample)]
      graphics::axis(1, at = seq_along(samples), labels = labs, las = 2, tick = FALSE, cex.axis = .7)
      graphics::title(main = study, cex.main = .9, line = 3)
      graphics::mtext("37 shared genes / 5 opposite-direction genes", side = 3, line = 1.6, cex = .65)
      graphics::mtext("Within-study gene z-score: blue -3, white 0, red +3", side = 1, line = 4, cex = .65)
    }
  })
  plot_device("Fig_selected_20_GSEA_GSVA_directions", 12, 10, function() {
    graphics::par(mar = c(4, 24, 4, 1))
    ids <- pathways$pathway; n <- length(ids)
    graphics::plot.new(); graphics::plot.window(xlim = c(.5, 6.5), ylim = c(.5, n + .5), xaxs = "i", yaxs = "i")
    for (j in seq_len(n)) for (s in seq_along(gobp_studies)) {
      d <- comparison[comparison$pathway == ids[j] & comparison$study == gobp_studies[s], , drop = FALSE]
      vals <- c(d$GSEA_NES, d$GSVA_delta); fs <- c(d$GSEA_FDR, d$GSVA_FDR)
      for (k in 1:2) {
        x <- (s - 1) * 2 + k; color <- if (vals[k] > 0) "#D98888" else if (vals[k] < 0) "#88AED0" else "#F7F7F7"
        graphics::rect(x - .48, n - j + .52, x + .48, n - j + 1.48, col = color, border = "white")
        if (fs[k] < .05) graphics::text(x, n - j + 1, "*", col = "black")
      }
    }
    categories <- pathways$evidence_category
    boundaries <- which(categories[-length(categories)] != categories[-1])
    if (length(boundaries)) graphics::abline(h = n - boundaries + .5, col = "#444444", lwd = 1.2)
    labels <- vapply(pathways$label, function(x) paste(strwrap(x, 42), collapse = "\n"), character(1))
    graphics::axis(2, at = rev(seq_len(n)), labels = labels, las = 2, tick = FALSE, cex.axis = .65)
    graphics::axis(1, at = 1:6, labels = rep(c("GSEA", "GSVA"), 3), tick = FALSE, cex.axis = .8)
    graphics::mtext(c("OIPN", "NC", "CCI"), side = 3, at = c(1.5, 3.5, 5.5), line = .5, cex = .8)
    graphics::title(main = "Selected pathways: effect direction and original FDR", cex.main = .9, line = 2)
    graphics::mtext("Red: positive; blue: negative. * Original within-study FDR <0.05. Colors show direction only.", side = 1, line = 2.5, cex = .7)
  })
}

run_gobp_selected_evidence <- function(draw_plots = TRUE, input_root = gobp_root,
    out = file.path(gobp_root, "selected_evidence_42_genes_20_pathways"), expression_loader = gobp_expression) {
  selected <- gobp08_selection(input_root); genes <- selected$genes; pathways <- selected$pathways
  expression_rows <- comparison_rows <- gene_sens <- pathway_sens <- diagnostics <- inputs <- list()
  for (study in gobp_studies) {
    spec <- gobp_spec(study); expected_meta <- gobp_metadata(study)
    a <- expression_loader(study)
    if (!identical(a$meta, expected_meta)) stop("Expression metadata differs from the locked study design.")
    counts <- gobp08_validate_matrix(a$counts, a$meta, genes$symbol, TRUE); lx <- log2(counts + 1)
    gf <- file.path(input_root, "GSEA", study, "gene_feature_evidence.csv")
    evidence <- gobp_read(gf, c("symbol", "feature_id", "log2FC", "gene_FDR"))
    if (anyDuplicated(evidence$symbol) || !all(genes$symbol %in% evidence$symbol)) stop("Invalid archived gene evidence.")
    ge <- evidence[match(genes$symbol, evidence$symbol), , drop = FALSE]
    gi <- a$genes[match(genes$symbol, a$genes$symbol), , drop = FALSE]
    if (!identical(as.character(ge$feature_id), as.character(gi$feature_id)) ||
        !isTRUE(all.equal(ge$log2FC, genes[[paste0(study, "_log2FC")]], tolerance = 1e-10)) ||
        !isTRUE(all.equal(ge$gene_FDR, genes[[paste0(study, "_gene_FDR")]], tolerance = 1e-10)))
      stop("Sample feature representatives disagree with selected gene statistics.")
    for (j in seq_len(nrow(genes))) {
      z <- if (stats::sd(lx[j, ]) > 0) (lx[j, ] - mean(lx[j, ])) / stats::sd(lx[j, ]) else rep(0, ncol(lx))
      expression_rows[[length(expression_rows) + 1L]] <- data.frame(study = study, symbol = genes$symbol[j],
        evidence_category = genes$evidence_category[j], feature_id = ge$feature_id[j], sample = a$meta$sample,
        condition = a$meta$condition, source_group = a$meta$source_group, normalized_count = counts[j, ],
        log2_count_plus_1 = lx[j, ], within_study_z = z,
        sample_label = paste0(rep(c("Ctrl", "Inj"), each = ncol(lx) / 2), rep(seq_len(ncol(lx) / 2), 2)),
        original_gene_log2FC = ge$log2FC[j], original_gene_FDR = ge$gene_FDR[j], stringsAsFactors = FALSE)
      for (group in c("Control", "Neuropathy")) {
        vals <- counts[j, a$meta$condition == group]
        diagnostics[[length(diagnostics) + 1L]] <- data.frame(study = study, symbol = genes$symbol[j], condition = group,
          n_samples = length(vals), mean_log2_count_plus_1 = mean(log2(vals + 1)), min_count = min(vals), max_count = max(vals),
          max_sample_fraction_of_group_counts = if (sum(vals) > 0) max(vals) / sum(vals) else NA_real_)
      }
    }
    gene_sens[[study]] <- gobp08_sensitivity(lx, a$meta, ge$log2FC, "gene_mean_log2_count_plus_1", study)
    vf <- file.path(input_root, "GSVA", study, "GSVA_limma_all_GO_BP.csv")
    sf <- file.path(input_root, "GSVA", study, "GSVA_GO_BP_scores.csv")
    mf <- file.path(input_root, "GSVA", study, "sample_metadata.csv")
    v <- gobp_read(vf, c("pathway", "delta_GSVA", "CI95_low", "CI95_high", "p_value", "FDR_BH", "GSEA_NES", "GSEA_FDR"))
    if (anyDuplicated(v$pathway) || !all(pathways$pathway %in% v$pathway)) stop("Selected pathway unavailable to GSVA.")
    v <- v[match(pathways$pathway, v$pathway), , drop = FALSE]
    nc <- pathways[[paste0(study, "_NES")]]; pc <- pathways[[paste0(study, "_FDR")]]
    if (any(!is.finite(as.matrix(v[c("delta_GSVA", "CI95_low", "CI95_high", "p_value", "FDR_BH")]))) ||
        any(v$FDR_BH < 0 | v$FDR_BH > 1) || any(v$CI95_low > v$CI95_high) ||
        !isTRUE(all.equal(nc, v$GSEA_NES, tolerance = 1e-10)) || !isTRUE(all.equal(pc, v$GSEA_FDR, tolerance = 1e-10)))
      stop("GSVA/GSEA archived statistics disagree or are invalid.")
    smeta <- gobp_read(mf, c("sample", "source_group", "condition"))
    if (!identical(smeta, expected_meta)) stop("GSVA sample assignments differ.")
    scores <- gobp_read(sf, c("pathway", a$meta$sample))
    mat <- as.matrix(scores[a$meta$sample]); rownames(mat) <- scores$pathway
    mat <- gobp08_validate_matrix(mat, smeta, pathways$pathway)
    path_s <- gobp08_sensitivity(mat, smeta, v$delta_GSVA, "pathway_frozen_GSVA_score", study)
    if (!isTRUE(all.equal(path_s$summary$baseline_delta, v$delta_GSVA, tolerance = 1e-10))) stop("Frozen GSVA scores disagree with original limma delta.")
    pathway_sens[[study]] <- path_s
    comparison_rows[[study]] <- data.frame(pathway = pathways$pathway, label = pathways$label,
      evidence_category = pathways$evidence_category, study = study, GSEA_NES = nc, GSEA_FDR = pc,
      GSEA_significant = pc < .05, GSVA_delta = v$delta_GSVA, GSVA_CI95_low = v$CI95_low,
      GSVA_CI95_high = v$CI95_high, GSVA_p_value = v$p_value, GSVA_FDR = v$FDR_BH,
      GSVA_significant = v$FDR_BH < .05, direction_agrees = sign(nc) == sign(v$delta_GSVA), stringsAsFactors = FALSE)
    inputs[[study]] <- c(spec$counts, spec$de, spec$mapping, spec$ranks, spec$meta, gf, vf, sf, mf)
  }
  expr <- do.call(rbind, expression_rows); comp <- do.call(rbind, comparison_rows)
  gs <- do.call(rbind, lapply(gene_sens, `[[`, "summary")); ps <- do.call(rbind, lapply(pathway_sens, `[[`, "summary"))
  gs$evidence_category <- genes$evidence_category[match(gs$id, genes$symbol)]
  ps$evidence_category <- pathways$evidence_category[match(ps$id, pathways$pathway)]
  summaries <- do.call(rbind, lapply(c("shared", "opposite"), function(category) do.call(rbind, lapply(gobp_studies, function(study) {
    d <- comp[comp$evidence_category == category & comp$study == study, , drop = FALSE]
    data.frame(evidence_category = category, study = study, n_pathways = nrow(d),
      n_direction_agrees = sum(d$direction_agrees), n_GSEA_significant = sum(d$GSEA_significant),
      n_GSVA_significant = sum(d$GSVA_significant), n_both_significant = sum(d$GSEA_significant & d$GSVA_significant))
  }))))
  outputs <- list(selected_genes = genes, selected_pathways = pathways, selected_gene_sample_expression = expr,
    gene_group_diagnostics = do.call(rbind, diagnostics), GSEA_GSVA_selected_comparison = comp,
    GSEA_GSVA_concordance_summary = summaries, gene_leave_one_out = do.call(rbind, lapply(gene_sens, `[[`, "long")),
    gene_sensitivity_summary = gs, pathway_leave_one_out = do.call(rbind, lapply(pathway_sens, `[[`, "long")),
    pathway_sensitivity_summary = ps)
  for (name in names(outputs)) gobp_write(outputs[[name]], file.path(out, paste0(name, ".csv")))
  writeLines(c("Descriptive evidence for the locked 37 shared and five opposite genes, and 11 shared and nine opposite pathways.",
    "Gene samples use the existing gobp_expression normalization and archived GSEA feature representatives.",
    "Gene sensitivity is the difference of group means of log2(normalized count+1), not a DESeq2 log2FC refit.",
    "GSVA sensitivity uses archived per-sample scores with scoring/normalization fixed; it is not full pipeline resampling.",
    "One sample is omitted at a time within each study. Every remaining Control and Neuropathy sample is retained.",
    "Direction preservation refers to the descriptive baseline delta; agreement with original DE/GSVA effect is reported separately.",
    "No new p-values/FDR, automatic sample exclusion, gene reprioritization, independent validation or mechanistic inference.",
    "The GSEA and GSVA figures compare directions; their magnitudes and significance are not interchangeable.",
    "Source selections and historical stage-05 outputs remain unchanged."), file.path(out, "analysis_notes.txt"))
  gobp_audit(c(selected$files, unlist(inputs), "scripts/gobp/08_selected_evidence_and_sensitivity.R"), out)
  if (draw_plots) gobp08_plot(expr, comp, pathways, out)
  message("Selected evidence complete: ", nrow(genes), " genes; ", nrow(pathways), " pathways; ", nrow(expr), " gene-sample rows.")
  invisible(outputs)
}
if (!exists("GOBP_AUTORUN", inherits = FALSE) || isTRUE(GOBP_AUTORUN))
  run_gobp_selected_evidence(tolower(Sys.getenv("GOBP_TABLES_ONLY", "false")) != "true")
