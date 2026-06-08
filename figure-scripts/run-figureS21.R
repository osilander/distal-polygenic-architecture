#!/usr/bin/env Rscript
# run-figureS21.R
#
# Plots observed vs within-block permutation null for loading concentration
# metrics (Gini, mean run length) for the GWAS-removed PCA, PC1-4.
#
# Input:  results/gwasremoved_withinperm_peaks_iter_metrics.tsv
#         results/gwasremoved_withinperm_peaks_summary.tsv
# Output: figures/figureS21_gwasremoved_withinperm_peaks_pc1_4.pdf
#
# Run from project root:
#   Rscript figure-scripts/run-figureS21.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

iter_dt <- fread(file.path(project_root, "results",
                           "gwasremoved_withinperm_peaks_iter_metrics.tsv"))
summ_dt <- fread(file.path(project_root, "results",
                           "gwasremoved_withinperm_peaks_summary.tsv"))

null_dt <- iter_dt[scheme != "observed"]
obs_dt  <- iter_dt[scheme == "observed"]

pc_cols <- c("PC1" = "#1f78b4", "PC2" = "#33a02c", "PC3" = "#e31a1c", "PC4" = "#984EA3")

scheme_labels <- c(
  chunk250kb_withinperm = "250 kb null",
  chunk500kb_withinperm = "500 kb null"
)
scheme_colours <- c(
  chunk250kb_withinperm = "#D9D9D9",
  chunk500kb_withinperm = "#636363"
)
metric_labels <- c(
  gini        = "Gini coefficient\nof |loadings|",
  mean_run_kb = "Mean run length\nof top 1%  (kb)"
)

null_dt[, scheme_lab := factor(scheme_labels[scheme], levels = scheme_labels)]
null_dt[, pc := factor(pc, levels = paste0("PC", 1:4))]
obs_dt[,  pc := factor(pc, levels = paste0("PC", 1:4))]

base_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
    strip.text         = element_text(face = "bold", size = 9),
    axis.title.x       = element_blank(),
    legend.position    = "bottom",
    legend.title       = element_blank()
  )

stars_label <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))

make_panel <- function(metric_id) {
  nd <- null_dt[, .(pc, scheme, scheme_lab, value = get(metric_id))]
  od <- obs_dt[,  .(pc, value = get(metric_id))]
  nd[, pc := factor(pc, levels = paste0("PC", 1:4))]
  od[, pc := factor(pc, levels = paste0("PC", 1:4))]

  # Significance stars: one per (pc x null scheme), positioned above the boxplot
  pv <- summ_dt[metric == metric_id & scheme %in% names(scheme_labels),
                .(pc, scheme, empirical_p)]
  pv[, pc       := factor(pc, levels = paste0("PC", 1:4))]
  pv[, pc_idx   := as.integer(pc)]
  pv[, dodge_offset := ifelse(scheme == "chunk250kb_withinperm", -0.35, 0.35)]
  pv[, x_pos    := pc_idx + dodge_offset]
  pv[, stars    := stars_label(empirical_p)]
  pv[, pc_col   := pc_cols[as.character(pc)]]
  # y: top whisker of each null distribution (97.5th pct + small margin)
  null_tops <- nd[, .(y_top = quantile(value, 0.975, na.rm = TRUE)), by = .(pc, scheme_lab)]
  null_tops[, scheme := names(scheme_labels)[match(scheme_lab, scheme_labels)]]
  pv <- merge(pv, null_tops[, .(pc, scheme, y_top)], by = c("pc", "scheme"), all.x = TRUE)
  pv <- pv[nzchar(stars)]

  ggplot(nd, aes(x = pc, y = value, fill = scheme_lab)) +
    geom_boxplot(
      outlier.size = 0.8, outlier.alpha = 0.5,
      linewidth = 0.4, width = 0.55,
      position = position_dodge(width = 0.7)
    ) +
    geom_point(
      data = od, aes(x = pc, y = value, colour = pc),
      inherit.aes = FALSE,
      shape = 18, size = 2.25
    ) +
    { if (nrow(pv) > 0)
        geom_text(data = pv,
                  aes(x = x_pos, y = y_top, label = stars, colour = pc),
                  inherit.aes = FALSE,
                  vjust = -0.3, size = 3.5)
      else NULL } +
    scale_fill_manual(values = unname(scheme_colours), labels = scheme_labels) +
    scale_colour_manual(values = pc_cols, guide = "none") +
    labs(y = metric_labels[metric_id]) +
    base_theme
}

p_gini     <- make_panel("gini")
p_mean_run <- make_panel("mean_run_kb")

fig <- p_gini / p_mean_run +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

out_path <- file.path(project_root, "figures",
                      "figureS21_gwasremoved_withinperm_peaks_pc1_4.pdf")
ggsave(out_path, fig, width = 4.8, height = 5.6)
cat("Written:", out_path, "\n")
