#!/usr/bin/env Rscript
# Produces: figures/figure4.pdf
# Pre-computed inputs required:
#   results/observed-nearest-patterns/observed_topfracs_nearestgene_bin_enrichment.tsv
#   results/observed-nearest-patterns/observed_topfracs_nearestgene_null_replicate_stats.tsv.gz
#   results/genicclass_rank_enrichment_50k/cumulative_genicclass_enrichment_curves_pc1_pc4.tsv
#   results/nearestgene_nonzero_distance_cdf_ribbon_50k_top1pct.tsv
#   results/nearestgene_nonzero_distance_cdf_ribbon_50k_top1pct_pvalues.tsv
# Run: Rscript figure-scripts/run-figure4.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

pc_cols <- c("PC1" = "#1f78b4", "PC2" = "#33a02c", "PC3" = "#e31a1c", "PC4" = "#984EA3")
pcs <- names(pc_cols)

gene_bins_file  <- file.path("results", "observed-nearest-patterns", "observed_topfracs_nearestgene_bin_enrichment.tsv")
gene_reps_file  <- file.path("results", "observed-nearest-patterns", "observed_topfracs_nearestgene_null_replicate_stats.tsv.gz")
genic_file      <- file.path("results", "genicclass_rank_enrichment_50k", "cumulative_genicclass_enrichment_curves_pc1_pc4.tsv")
ecdf_file       <- file.path("results", "nearestgene_nonzero_distance_cdf_ribbon_50k_top1pct.tsv")
ecdf_pvals_file <- file.path("results", "nearestgene_nonzero_distance_cdf_ribbon_50k_top1pct_pvalues.tsv")
out_fig         <- file.path("figures", "figure4.pdf")
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

for (f in c(gene_bins_file, genic_file, ecdf_file, ecdf_pvals_file)) {
  if (!file.exists(f)) stop("Missing pre-computed input: ", f)
}

gene_bins  <- fread(gene_bins_file)[top_frac %in% c(0.01)]
gene_reps  <- fread(cmd = paste("gzip -cd", shQuote(gene_reps_file)))[top_frac %in% c(0.01)]
genic      <- fread(genic_file)
ecdf_dt    <- fread(ecdf_file)
ecdf_pvals <- fread(ecdf_pvals_file)

gene_bins[, bin_group := fifelse(
  bin %in% c("0", "(0,50kb]", "(50,100kb]", "(100,250kb]"), bin, ">250kb"
)]
gene_bins_core <- gene_bins[bin_group != ">250kb", .(
  obs_frac = sum(obs_frac, na.rm = TRUE),
  null_mean = sum(null_mean, na.rm = TRUE),
  enrich_diff = sum(enrich_diff, na.rm = TRUE),
  enrich_diff_lo = sum(enrich_diff_lo, na.rm = TRUE),
  enrich_diff_hi = sum(enrich_diff_hi, na.rm = TRUE)
), by = .(pc, top_frac, bin = bin_group)]

far_cols <- c("(250,500kb]", "(500,750kb]", "(750kb,1Mb]", ">1Mb")
far_obs  <- gene_bins[bin_group == ">250kb", .(obs_frac = sum(obs_frac, na.rm = TRUE)), by = .(pc, top_frac)]
far_reps <- gene_reps[, .(null_far = rowSums(.SD, na.rm = TRUE)), by = .(pc, top_frac, replicate), .SDcols = far_cols]
far_bins <- far_reps[, .(
  null_mean      = mean(null_far, na.rm = TRUE),
  enrich_diff_lo = far_obs[.BY, on = .(pc, top_frac)]$obs_frac - quantile(null_far, 0.975, na.rm = TRUE),
  enrich_diff_hi = far_obs[.BY, on = .(pc, top_frac)]$obs_frac - quantile(null_far, 0.025, na.rm = TRUE)
), by = .(pc, top_frac)]
far_bins <- far_bins[far_obs, on = .(pc, top_frac)]
far_bins[, `:=`(enrich_diff = obs_frac - null_mean, bin = ">250kb")]

gene_bins <- rbindlist(list(
  gene_bins_core,
  far_bins[, .(pc, top_frac, bin, obs_frac, null_mean, enrich_diff, enrich_diff_lo, enrich_diff_hi)]
), use.names = TRUE)

gene_reps_long <- rbindlist(list(
  gene_reps[, .(pc, top_frac, replicate, bin = "0",          null_frac = `0`)],
  gene_reps[, .(pc, top_frac, replicate, bin = "(0,50kb]",   null_frac = `(0,50kb]`)],
  gene_reps[, .(pc, top_frac, replicate, bin = "(50,100kb]", null_frac = `(50,100kb]`)],
  gene_reps[, .(pc, top_frac, replicate, bin = "(100,250kb]",null_frac = `(100,250kb]`)],
  gene_reps[, .(pc, top_frac, replicate, bin = ">250kb",
                null_frac = `(250,500kb]` + `(500,750kb]` + `(750kb,1Mb]` + `>1Mb`)]
), use.names = TRUE)

gene_bins <- gene_reps_long[gene_bins, on = .(pc, top_frac, bin)]
gene_bins[, p_emp_directional := if (obs_frac[1] >= mean(null_frac, na.rm = TRUE)) {
  (1 + sum(null_frac >= obs_frac[1], na.rm = TRUE)) / (1 + sum(is.finite(null_frac)))
} else {
  (1 + sum(null_frac <= obs_frac[1], na.rm = TRUE)) / (1 + sum(is.finite(null_frac)))
}, by = .(pc, top_frac, bin)]
gene_bins <- unique(gene_bins[, .(pc, top_frac, bin, obs_frac, null_mean, enrich_diff, enrich_diff_lo, enrich_diff_hi, p_emp_directional)])
gene_bins[, log2r    := log2(obs_frac / null_mean)]
gene_bins[, log2r_lo := log2(obs_frac / (obs_frac - enrich_diff_lo))]  # obs / null_hi
gene_bins[, log2r_hi := log2(obs_frac / (obs_frac - enrich_diff_hi))]  # obs / null_lo
gene_bins[, label_y  := fifelse(log2r >= 0, log2r_hi + 0.05, log2r_lo - 0.05)]
gene_bins[, sig_label := fifelse(p_emp_directional < 0.001, "***",
                          fifelse(p_emp_directional < 0.01, "**",
                                  fifelse(p_emp_directional < 0.05, "*", "")))]
gene_bins[, bin := factor(bin, levels = c("0", "(0,50kb]", "(50,100kb]", "(100,250kb]", ">250kb"))]
gene_bins[, pc  := factor(pc, levels = pcs)]

genic[, annotation := factor(annotation, levels = c("promoter", "exonic", "intronic", "intergenic"))]
genic <- genic[annotation %in% c("promoter", "exonic", "intronic", "intergenic")]
genic[, pc := factor(pc, levels = pcs)]
ecdf_dt[, pc    := factor(pc, levels = pcs)]
ecdf_pvals[, pc := factor(pc, levels = pcs)]
ecdf_pvals[, `:=`(x = 1950, y = 0.05, label = sprintf("p = %.3g", p_emp))]

p_gene <- ggplot(gene_bins, aes(x = bin, y = log2r, color = pc, group = pc)) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
  geom_linerange(aes(ymin = log2r_lo, ymax = log2r_hi),
                 position = position_dodge(width = 0.45), linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.45), size = 1.8) +
  geom_text(data = gene_bins[nzchar(sig_label)], aes(y = label_y, label = sig_label),
            position = position_dodge(width = 0.45), size = 3.0, show.legend = FALSE) +
  scale_color_manual(values = pc_cols, drop = FALSE) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(x = "Distance to most proximal gene body (binned)", y = expression(log[2](observed/null)), color = "PC")

p_ecdf <- ggplot() +
  geom_ribbon(data = ecdf_dt, aes(x = dist_bp / 1000, ymin = sim_lo, ymax = sim_hi),
              fill = "grey80", alpha = 0.5) +
  geom_line(data = ecdf_dt, aes(x = dist_bp / 1000, y = sim_med),
            color = "grey45", linewidth = 0.55) +
  geom_step(data = ecdf_dt, aes(x = dist_bp / 1000, y = obs_cdf, color = pc),
            linewidth = 0.9, show.legend = FALSE) +
  geom_text(data = ecdf_pvals, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 1, vjust = 0, size = 2.5, color = "black") +
  scale_color_manual(values = pc_cols, drop = FALSE) +
  scale_x_continuous(limits = c(0, 2000), breaks = seq(0, 2000, by = 500)) +
  facet_wrap(~pc, ncol = 2) +
  theme_minimal(base_size = 8.5) +
  theme(panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.2),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
        strip.text = element_text(size = 9),
        axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) +
  labs(x = "Distance to most proximal gene body (non-zero)", y = "Cumulative fraction")

p_genic <- ggplot(genic, aes(x = top_pct, y = enrichment_ratio, color = pc, fill = pc)) +
  geom_hline(yintercept = 1, color = "grey55", linewidth = 0.35, linetype = "dotted") +
  geom_ribbon(aes(ymin = null_enrichment_lo, ymax = null_enrichment_hi), alpha = 0.10, linewidth = 0, color = NA) +
  geom_line(linewidth = 0.75) +
  scale_color_manual(values = pc_cols, drop = FALSE) +
  scale_fill_manual(values = pc_cols, drop = FALSE) +
  facet_wrap(~annotation, ncol = 2) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey94", linewidth = 0.2)) +
  labs(x = "Top loading windows included (%)",
       y = "Ratio of observed fraction to expected (null) fraction",
       color = "PC")

fig <- ((p_gene / p_ecdf + plot_layout(heights = c(1, 1))) | p_genic) +
  plot_layout(widths = c(1.0, 1.45)) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = rel(3.0)))

ggsave(out_fig, fig, width = 14.6, height = 9.4)
message("Figure 4 complete: ", out_fig)
