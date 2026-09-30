# Leading-edge heatmap color scale audit

The tree of the three-study branch and related NC input, differential-expression, leading-edge extraction and manuscript protocol files were reviewed. Script 05 and script 07 both previously mapped raw log2FC onto a linear scale clipped at +/-2. GSEA NES heatmaps and standardized GSVA score heatmaps use different quantities and are outside this display change.

The 133 plotted pathway-gene memberships in script 07 have median absolute log2FC of 0.1882 (OIPN), 1.2747 (NC), and 0.1719 (CCI). At the old endpoints, 1/133 OIPN, 48/133 NC and 9/133 CCI memberships reached saturation. In COMPLEMENT, 10/17 NC cells reached the endpoint. C1qa log2FC is -0.13265 in OIPN, +3.97613 in NC and +1.12219 in CCI; Olr1 is -0.00549, +5.62021 and +0.30516 respectively.

All 133 plotted NC effect values match the original DESeq2 all-gene table by the archived representative feature ID. NC DE uses results(dds, contrast=c("condition","Compression","Sham")); its exported effects are unshrunk. The six locked sample identities agree with the input and PCA metadata. Per-sample count totals range from 15.8 to 19.0 million and archived size factors from 0.868 to 1.309. These checks explain neither all biological differences nor rule out all technical effects. The precise rendering cause is larger NC estimates combined with display saturation; stronger colors do not by themselves establish an analysis error or a stronger causal mechanism.

## Display change

Both three-study gene heatmap scripts now use a common fixed linear range of -6 to +6 by default. This is an aesthetic range, not a statistical threshold. It reduces saturation without changing numerical effect sizes, gene selection, missing-data handling, FDR, ranks or pathway classifications. At this range, script 07 endpoint counts are 0/133 OIPN, 18/133 NC and 0/133 CCI; none of the 17 COMPLEMENT NC cells reach the endpoint. Extreme values beyond +/-6 still saturate, and small OIPN/CCI effects become paler as well.

Use the same range across studies and panels; do not rescale only the NC column. Scripts record plot_settings.csv and label the legend with the actual clipping limit. Sys.setenv(LE_COLOR_LIMIT="2") reproduces the earlier color range; unsetting it restores the new default of 6. This environment setting controls display only.

The change does not shrink DE estimates or imply a biological interpretation of their magnitude. If biological robustness of large NC effects is questioned, examine sample-level normalized counts, dispersion and composition separately before changing inference.

## Validation

Original NC effect joins and old/new endpoint counts were checked directly against repository CSVs. Script modifications are confined to display configuration, color mapping and a new plot-settings export; extraction and prioritization code remain unchanged. R lexical delimiter checks passed. No R runtime or image rendering is available here; local execution remains necessary. No figures were generated.
