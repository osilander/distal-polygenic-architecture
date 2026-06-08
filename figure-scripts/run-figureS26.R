#!/usr/bin/env Rscript
# Produces: figures/figureS26_gwasremoved_nearestgene_pattern_summary_gwasremoved100_150kb.pdf
# Pre-computed inputs:
#   results/observed-nearest-patterns-gwasremoved100kb-top1_2_5_band45_50/
#     observed_nearestgene_bin_enrichment_gwasremoved100kb_top1_2_5_band45_50.tsv
#   results/observed-nearest-patterns-gwasremoved150kb-top1_2_5_band45_50/
#     observed_nearestgene_bin_enrichment_gwasremoved150kb_top1_2_5_band45_50.tsv
#     (from compute-nearestgene-pattern-bands.R gwasremoved default gene_only {100|150})
# Run: Rscript figure-scripts/run-figureS26.R
# Note: stacks 100 kb and 150 kb removal flanks as separate rows in a single patchwork figure.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

pc_cols <- c("PC1" = "#1f78b4", "PC2" = "#33a02c", "PC3" = "#e31a1c", "PC4" = "#984EA3")
pcs     <- names(pc_cols)
flanks  <- c(100L, 150L)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

load_flank <- function(flank_kb) {
  file_tag  <- sprintf("gwasremoved%dkb", flank_kb)
  bins_file <- file.path("results",
    sprintf("observed-nearest-patterns-%s-top1_2_5_band45_50", file_tag),
    sprintf("observed_nearestgene_bin_enrichment_%s_top1_2_5_band45_50.tsv", file_tag))
  if (!file.exists(bins_file))
    stop("Missing input: ", bins_file,
         "\nRun: Rscript figure-scripts/compute-nearestgene-pattern-bands.R gwasremoved default gene_only ", flank_kb)
  dt <- fread(bins_file)
  dt[, flank_label := sprintf("\u00b1%d kb flanks removed", flank_kb)]
  dt
}

make_row_plot <- function(sub, show_legend = FALSE) {
  bin_levels <- unique(sub$bin)
  sub[, bin       := factor(bin, levels = bin_levels)]
  sub[, band_label := factor(band_label, levels = c("Top 1%", "Top 2%", "Top 5%", "45-50%"))]
  sub[, pc        := factor(pc, levels = pcs)]
  sub[, bin_idx   := match(bin, levels(bin))]
  sub[, pc_idx    := match(pc, pcs)]
  offsets <- seq(-0.75 * 0.25, 0.75 * 0.25, length.out = length(pcs))
  sub[, x_pos    := bin_idx + offsets[pc_idx]]
  sub[, y_val    := log2(obs_frac / null_mean)]
  sub[, y_lo     := log2(obs_frac / (obs_frac - enrich_diff_lo))]
  sub[, y_hi     := log2(obs_frac / (obs_frac - enrich_diff_hi))]
  sub[, label_y  := fifelse(y_val >= 0, y_hi + 0.05, y_lo - 0.05)]
  sub[, crosses_zero := enrich_diff_lo <= 0 & enrich_diff_hi >= 0]
  sub[, sig_label := fifelse(
    crosses_zero, "",
    fifelse(p_emp_directional < 0.001, "***",
            fifelse(p_emp_directional < 0.01, "**",
                    fifelse(p_emp_directional < 0.05, "*", "")))
  )]
  row_tag <- unique(sub$flank_label)

  ggplot(sub, aes(x = bin, y = y_val, color = pc, group = pc)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_linerange(aes(ymin = y_lo, ymax = y_hi), alpha = 0.35,
                   position = position_dodge(width = 0.25)) +
    geom_point(size = 1.6, position = position_dodge(width = 0.25)) +
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
          panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.3),
          axis.text.x      = element_text(angle = 30, hjust = 1, size = 10.5),
          strip.text       = element_text(face = "plain", size = 11),
          legend.position  = if (show_legend) "right" else "none") +
    labs(x = "Nearest-distance bin", y = expression(log[2](observed/null)),
         color = "PC", tag = row_tag)
}

rows <- lapply(seq_along(flanks), function(i) {
  make_row_plot(load_flank(flanks[i]), show_legend = (i == length(flanks)))
})

fig <- Reduce(`/`, rows) + plot_layout(guides = "collect") &
  theme(legend.position = "right")

out_fig <- file.path("figures",
  "figureS26_gwasremoved_nearestgene_pattern_summary_gwasremoved100_150kb.pdf")
ggsave(out_fig, fig, width = 14.5, height = 7.2)
message("Wrote: ", out_fig)
