#!/usr/bin/env Rscript
# Produces: figures/figureS16_observed_nearestgwaspeak_pattern_summary_50k_100k_comparison.pdf
#           figures/figureS24_observed_nearestgene_pattern_summary_50k_100k_comparison.pdf
# Pre-computed inputs:
#   results/observed-nearest-patterns-50k-top1_2_5_band45_50/
#   results/observed-nearest-patterns-100k-top1_2_5_band45_50/
# Run: Rscript post-process-scripts/run-figureS16-and-S24.R
# Note: uses patchwork so each window-size row (50k, 100k) has its own independent axes.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

pc_cols <- c("PC1" = "#1f78b4", "PC2" = "#33a02c", "PC3" = "#e31a1c", "PC4" = "#984EA3")
pcs <- names(pc_cols)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

load_bins <- function(kind = c("gene", "gwas")) {
  kind  <- match.arg(kind)
  files <- if (kind == "gene") {
    c(
      `50k`  = file.path("results", "observed-nearest-patterns-50k-top1_2_5_band45_50",
                         "observed_nearestgene_bin_enrichment_50k_top1_2_5_band45_50.tsv"),
      `100k` = file.path("results", "observed-nearest-patterns-100k-top1_2_5_band45_50",
                         "observed_nearestgene_bin_enrichment_100k_top1_2_5_band45_50.tsv")
    )
  } else {
    c(
      `50k`  = file.path("results", "observed-nearest-patterns-50k-top1_2_5_band45_50",
                         "observed_nearestgwaspeak_bin_enrichment_50k_top1_2_5_band45_50.tsv"),
      `100k` = file.path("results", "observed-nearest-patterns-100k-top1_2_5_band45_50",
                         "observed_nearestgwaspeak_bin_enrichment_100k_top1_2_5_band45_50.tsv")
    )
  }
  rows <- lapply(names(files), function(tag) {
    f  <- files[[tag]]
    if (!file.exists(f)) stop("Missing input: ", f)
    dt <- fread(f)
    dt[, window_tag := factor(tag, levels = c("50k", "100k"), labels = c("50 kb windows", "100 kb windows"))]
    dt
  })
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

prep_dt <- function(dt, log2r) {
  dt[, band_label := factor(band_label, levels = c("Top 1%", "Top 2%", "Top 5%", "45-50%"))]
  dt[, crosses_zero := enrich_diff_lo <= 0 & enrich_diff_hi >= 0]
  if (log2r) {
    dt[, y_val  := log2(obs_frac / null_mean)]
    dt[, y_lo   := log2(obs_frac / (obs_frac - enrich_diff_lo))]
    dt[, y_hi   := log2(obs_frac / (obs_frac - enrich_diff_hi))]
    dt[, label_y := fifelse(y_val >= 0, y_hi + 0.05, y_lo - 0.05)]
  } else {
    dt[, y_val  := enrich_diff]
    dt[, y_lo   := enrich_diff_lo]
    dt[, y_hi   := enrich_diff_hi]
    dt[, label_y := fifelse(enrich_diff >= 0, enrich_diff_hi, enrich_diff_lo)]
  }
  dt[, sig_label := fifelse(
    crosses_zero, "",
    fifelse(p_emp_directional < 0.001, "***",
            fifelse(p_emp_directional < 0.01, "**",
                    fifelse(p_emp_directional < 0.05, "*", "")))
  )]
  dt
}

make_row_plot <- function(sub, y_lab, show_legend = FALSE) {
  bin_levels <- unique(sub$bin)
  sub[, bin    := factor(bin, levels = bin_levels)]
  sub[, bin_idx := match(bin, levels(bin))]
  sub[, pc_idx  := match(pc, pcs)]
  n_pcs   <- length(pcs)
  offsets <- (seq_len(n_pcs) - (n_pcs + 1) / 2) * (0.25 / n_pcs)
  sub[, x_pos := bin_idx + offsets[pc_idx]]
  row_tag <- as.character(unique(sub$window_tag))

  pd <- position_dodge(width = 0.25)
  ggplot(sub, aes(x = bin, y = y_val, color = pc, group = pc)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_linerange(aes(ymin = y_lo, ymax = y_hi), alpha = 0.35,
                   position = pd) +
    geom_point(size = 1.6, position = pd) +
    geom_text_repel(
      data = sub[nzchar(sig_label)], inherit.aes = FALSE,
      aes(x = x_pos, y = label_y, label = sig_label, color = pc),
      size = 3.2, direction = "y", box.padding = 0.15, point.padding = 0.05,
      segment.size = 0.2, segment.alpha = 0.5, min.segment.length = 0,
      force = 0.8, max.overlaps = Inf, show.legend = FALSE
    ) +
    scale_color_manual(values = pc_cols, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    facet_wrap(~ band_label, nrow = 1) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
          panel.grid.minor = element_line(color = "grey94", linewidth = 0.2),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
          axis.text.x = element_text(angle = 30, hjust = 1, size = 10.5),
          strip.text = element_text(face = "plain", size = 11),
          legend.position = if (show_legend) "right" else "none") +
    labs(x = "Nearest-distance bin", y = y_lab, color = "PC",
         tag = row_tag)
}

plot_comp <- function(dt, out_fig, title_text, log2r = FALSE, free_scales = FALSE) {
  y_lab <- if (log2r) expression(log[2](observed/null)) else "Observed - null fraction"
  dt    <- prep_dt(copy(dt), log2r)

  if (free_scales) {
    tags  <- levels(dt$window_tag)
    plots <- lapply(seq_along(tags), function(i) {
      make_row_plot(dt[window_tag == tags[i]], y_lab, show_legend = (i == length(tags)))
    })
    fig <- Reduce(`/`, plots) + plot_layout(guides = "collect") &
      theme(legend.position = "right")
    ggsave(out_fig, fig, width = 14.5, height = 7.2)
  } else {
    dt[, bin_idx := match(bin, levels(factor(bin, levels = unique(bin))))]
    dt[, pc_idx  := match(pc, pcs)]
    n_pcs   <- length(pcs)
    offsets <- (seq_len(n_pcs) - (n_pcs + 1) / 2) * (0.25 / n_pcs)
    dt[, x_pos := bin_idx + offsets[pc_idx]]
    pd <- position_dodge(width = 0.25)
    p <- ggplot(dt, aes(x = bin, y = y_val, color = pc, group = pc)) +
      geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
      geom_linerange(aes(ymin = y_lo, ymax = y_hi), alpha = 0.35,
                     position = pd) +
      geom_point(size = 1.6, position = pd) +
      geom_text_repel(
        data = dt[nzchar(sig_label)], inherit.aes = FALSE,
        aes(x = x_pos, y = label_y, label = sig_label, color = pc),
        size = 3.2, direction = "y", box.padding = 0.15, point.padding = 0.05,
        segment.size = 0.2, segment.alpha = 0.5, min.segment.length = 0,
        force = 0.8, max.overlaps = Inf, show.legend = FALSE
      ) +
      scale_color_manual(values = pc_cols, drop = FALSE) +
      scale_x_discrete(drop = FALSE) +
      facet_grid(window_tag ~ band_label) +
      theme_minimal(base_size = 10) +
      theme(panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
            panel.grid.minor = element_line(color = "grey94", linewidth = 0.2),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
            axis.text.x = element_text(angle = 30, hjust = 1, size = 10.5),
            strip.text = element_text(face = "plain", size = 11),
            legend.position = "right") +
      labs(title = title_text, x = "Nearest-distance bin", y = y_lab, color = "PC")
    ggsave(out_fig, p, width = 14.5, height = 7.2)
  }
  message("Wrote: ", out_fig)
}

gene_dt <- load_bins("gene")
gwas_dt <- load_bins("gwas")

# Replace 50k + Top 1% rows with data computed identically to figure 4a
# (replicate-based CI for >250kb, same random seed / null replicates)
fig4_bins_file <- file.path("results", "observed-nearest-patterns",
                             "observed_topfracs_nearestgene_bin_enrichment.tsv")
fig4_reps_file <- file.path("results", "observed-nearest-patterns",
                             "observed_topfracs_nearestgene_null_replicate_stats.tsv.gz")

if (file.exists(fig4_bins_file) && file.exists(fig4_reps_file)) {
  fig4_bins <- fread(fig4_bins_file)[top_frac == 0.01]
  fig4_reps <- fread(cmd = paste("gzip -cd", shQuote(fig4_reps_file)))[top_frac == 0.01]

  # Collapse bins identically to figure 4a
  fig4_bins[, bin_group := fifelse(
    bin %in% c("0", "(0,50kb]", "(50,100kb]", "(100,250kb]"), bin, ">250kb"
  )]
  fig4_core <- fig4_bins[bin_group != ">250kb", .(
    obs_frac       = sum(obs_frac,       na.rm = TRUE),
    null_mean      = sum(null_mean,      na.rm = TRUE),
    enrich_diff    = sum(enrich_diff,    na.rm = TRUE),
    enrich_diff_lo = sum(enrich_diff_lo, na.rm = TRUE),
    enrich_diff_hi = sum(enrich_diff_hi, na.rm = TRUE)
  ), by = .(pc, top_frac, bin = bin_group)]

  far_cols <- c("(250,500kb]", "(500,750kb]", "(750kb,1Mb]", ">1Mb")
  far_obs  <- fig4_bins[bin_group == ">250kb",
                        .(obs_frac = sum(obs_frac, na.rm = TRUE)), by = .(pc, top_frac)]
  far_reps <- fig4_reps[, .(null_far = rowSums(.SD, na.rm = TRUE)),
                        by = .(pc, top_frac, replicate), .SDcols = far_cols]
  far_bins <- far_reps[, .(
    null_mean      = mean(null_far, na.rm = TRUE),
    enrich_diff_lo = far_obs[.BY, on = .(pc, top_frac)]$obs_frac -
                       quantile(null_far, 0.975, na.rm = TRUE),
    enrich_diff_hi = far_obs[.BY, on = .(pc, top_frac)]$obs_frac -
                       quantile(null_far, 0.025, na.rm = TRUE)
  ), by = .(pc, top_frac)]
  far_bins <- far_bins[far_obs, on = .(pc, top_frac)]
  far_bins[, `:=`(enrich_diff = obs_frac - null_mean, bin = ">250kb")]

  fig4_bins2 <- rbindlist(list(
    fig4_core,
    far_bins[, .(pc, top_frac, bin, obs_frac, null_mean, enrich_diff, enrich_diff_lo, enrich_diff_hi)]
  ), use.names = TRUE)

  # p_emp_directional from replicates (same logic as figure 4a)
  fig4_reps_long <- rbindlist(list(
    fig4_reps[, .(pc, top_frac, replicate, bin = "0",           null_frac = `0`)],
    fig4_reps[, .(pc, top_frac, replicate, bin = "(0,50kb]",    null_frac = `(0,50kb]`)],
    fig4_reps[, .(pc, top_frac, replicate, bin = "(50,100kb]",  null_frac = `(50,100kb]`)],
    fig4_reps[, .(pc, top_frac, replicate, bin = "(100,250kb]", null_frac = `(100,250kb]`)],
    fig4_reps[, .(pc, top_frac, replicate, bin = ">250kb",
                  null_frac = `(250,500kb]` + `(500,750kb]` + `(750kb,1Mb]` + `>1Mb`)]
  ), use.names = TRUE)

  fig4_bins2 <- fig4_reps_long[fig4_bins2, on = .(pc, top_frac, bin)]
  fig4_bins2[, p_emp_directional := if (obs_frac[1] >= mean(null_frac, na.rm = TRUE)) {
    (1 + sum(null_frac >= obs_frac[1], na.rm = TRUE)) / (1 + sum(is.finite(null_frac)))
  } else {
    (1 + sum(null_frac <= obs_frac[1], na.rm = TRUE)) / (1 + sum(is.finite(null_frac)))
  }, by = .(pc, top_frac, bin)]
  fig4_bins2 <- unique(fig4_bins2[, .(pc, bin, obs_frac, null_mean, enrich_diff,
                                       enrich_diff_lo, enrich_diff_hi, p_emp_directional)])
  fig4_bins2[, `:=`(
    band_label = "Top 1%",
    window_tag = factor("50 kb windows", levels = c("50 kb windows", "100 kb windows"))
  )]

  # Swap out the old 50k + Top 1% rows
  gene_dt <- gene_dt[!(window_tag == "50 kb windows" & band_label == "Top 1%")]
  gene_dt  <- rbindlist(list(gene_dt, fig4_bins2), use.names = TRUE, fill = TRUE)
  message("Replaced 50k Top 1% rows with figure-4a-compatible data.")
} else {
  message("Warning: figure-4a source files not found; using pre-computed 50k Top 1% data.")
}

plot_comp(gene_dt, file.path("figures", "figureS24_observed_nearestgene_pattern_summary_50k_100k_comparison.pdf"),
          "Nearest gene: binned enrichment", log2r = TRUE, free_scales = TRUE)
plot_comp(gwas_dt, file.path("figures", "figureS16_observed_nearestgwaspeak_pattern_summary_50k_100k_comparison.pdf"),
          "Nearest GWAS-significant window: binned enrichment", log2r = TRUE, free_scales = TRUE)

message("FigureS16 and FigureS24 complete.")
