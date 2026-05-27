#!/usr/bin/env Rscript
# Produces: figures/figureS18_gwas_removed_trait_pc_category_absrho_labelperm_pc1_pc14.pdf,
#           results/gwas_removed_trait_pc_category_absrho_labelperm_pc1_pc14.tsv
# Pre-computed inputs:
#   results/all-trait-50k-mean-{effects,pvals}.tsv, data/trait_abbrevs_categorised.txt,
#   results/trait_pc_category_absrho_labelperm_pc1_pc14.tsv (for row order),
#   results/gwas_removed_distance_topload_50k_plus100kb/reduced_pca_window_loadings.tsv,
#   results/gwas_removed_distance_topload_50k_plus150kb/reduced_pca_window_loadings.tsv
# Run: Rscript post-process-scripts/run-figureS18.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(paletteer)
  library(cowplot)
})

# ---- Top panel: absolute loading distributions (GWAS peaks vs retained) ----
# Core GWAS peaks (gws_any == TRUE) vs retained windows; flanks excluded from both groups.

loadings_file <- file.path("results", "pca_loadings_50k.tsv")
removed_file  <- file.path("results",
  "gwas_removed_distance_topload_50k_plus100kb",
  "removed_windows_gwas_plus100kb.tsv")

ld_top  <- fread(loadings_file)
rem_top <- fread(removed_file)
ld_top[, window_id := paste(chrom, start, end, sep = "_")]
ld_top[, grp := fcase(
  window_id %in% rem_top[gws_any == TRUE,  window_id], "Core peaks",
  window_id %in% rem_top[gws_any == FALSE, window_id], "Flank",
  default = "Retained"
)]

pcs_top    <- paste0("PC", 1:4)
gwas_blue  <- "#1f78b4"

for (pc in pcs_top) ld_top[, paste0("abs_", pc) := abs(get(pc))]

# Shared x limit: 99.5th percentile across all four PCs combined
all_abs <- unlist(ld_top[grp %in% c("Core peaks", "Retained"),
                          lapply(paste0("abs_", pcs_top), function(col) get(col))])
shared_xmax <- quantile(all_abs, 0.995)

make_density_panel <- function(pc_name) {
  abs_col <- paste0("abs_", pc_name)
  sub <- ld_top[grp %in% c("Core peaks", "Retained"),
                .(abs_loading = get(abs_col),
                  group = factor(grp, levels = c("Core peaks", "Retained")))]
  x_core <- sub[group == "Core peaks", abs_loading]
  x_ret  <- sub[group == "Retained",   abs_loading]
  wt     <- wilcox.test(x_core, x_ret, exact = FALSE)
  rbc    <- 2 * wt$statistic / (length(x_core) * length(x_ret)) - 1
  ann    <- sprintf("r = %.3f\np < 0.001", rbc)

  ggplot(sub, aes(x = abs_loading, colour = group, fill = group)) +
    geom_density(alpha = 0.20, linewidth = 0.7) +
    annotate("text", x = Inf, y = Inf, label = ann,
             hjust = 1.05, vjust = 1.3, size = 2.8, colour = "grey20") +
    scale_colour_manual(
      values = c("Core peaks" = gwas_blue, "Retained" = "grey45"),
      labels = c("Core peaks" = "GWAS-significant loci", "Retained" = "Non-significant loci"),
      name = NULL) +
    scale_fill_manual(
      values = c("Core peaks" = gwas_blue, "Retained" = "grey70"),
      labels = c("Core peaks" = "GWAS-significant loci", "Retained" = "Non-significant loci"),
      name = NULL) +
    # Merge colour+fill into 2 legend keys; suppress the separate fill guide
    guides(colour = guide_legend(override.aes = list(fill  = c(gwas_blue, "grey70"),
                                                      alpha = 0.5)),
           fill = "none") +
    coord_cartesian(xlim = c(0, shared_xmax)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.4),
      legend.position  = "none",
      axis.title.y     = element_blank(),
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      plot.title       = element_text(size = 10, face = "bold", hjust = 0.5,
                                      margin = margin(b = 2))
    ) +
    labs(x = "|Loading|", y = NULL, title = pc_name)
}

top_panels <- lapply(pcs_top, make_density_panel)
names(top_panels) <- pcs_top

# Shared legend: build once from a single panel rendered with legend visible
leg_source <- top_panels$PC1 +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(colour = guide_legend(override.aes = list(fill  = c(gwas_blue, "grey70"),
                                                    alpha = 0.5)),
         fill = "none")
shared_leg <- cowplot::get_legend(leg_source)

# Tag 'a' on first density panel; build density row with shared y-label
top_panels$PC1 <- top_panels$PC1 + labs(tag = "a")
top_row      <- cowplot::plot_grid(plotlist = top_panels, nrow = 1, align = "h")
y_label      <- cowplot::ggdraw() +
  cowplot::draw_label("Density", angle = 90, size = 10, colour = "grey20")
top_fig_grid <- cowplot::plot_grid(y_label, top_row, ncol = 2, rel_widths = c(0.04, 1))
top_density_block <- cowplot::plot_grid(top_fig_grid, shared_leg,
                                        ncol = 1, rel_heights = c(1, 0.15))

# ---- Middle panel: top-1% composition (GWAS core / flank / retained) ----
top_frac <- 0.01
top1pct_stats <- rbindlist(lapply(pcs_top, function(pc) {
  abs_col   <- paste0("abs_", pc)
  threshold <- quantile(ld_top[[abs_col]], 1 - top_frac, na.rm = TRUE)
  top_wins  <- ld_top[get(abs_col) >= threshold]
  n_total   <- nrow(top_wins)
  data.table(
    pc           = pc,
    `GWAS core`  = sum(top_wins$grp == "Core peaks") / n_total,
    `Flank`      = sum(top_wins$grp == "Flank")       / n_total,
    `Retained`   = sum(top_wins$grp == "Retained")    / n_total,
    n_total      = n_total
  )
}))
top1pct_long <- melt(top1pct_stats, id.vars = c("pc", "n_total"),
                     measure.vars = c("GWAS core", "Flank", "Retained"),
                     variable.name = "group", value.name = "fraction")
top1pct_long[, group := factor(group, levels = c("Retained", "Flank", "GWAS core"))]
top1pct_long[, pc    := factor(pc, levels = pcs_top)]

flank_col <- "#6aaed6"   # lighter blue for flanks
bar_cols  <- c("GWAS core" = gwas_blue, "Flank" = flank_col, "Retained" = "grey70")

# Compute totals removed for annotation
removed_fracs <- top1pct_stats[, .(pc, removed = `GWAS core` + Flank)]

mid_panel <- ggplot(top1pct_long, aes(x = pc, y = fraction, fill = group)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
  geom_text(data = removed_fracs,
            aes(x = pc, y = removed + 0.015,
                label = sprintf("%.0f%%\nremoved", removed * 100)),
            inherit.aes = FALSE, size = 2.6, colour = "grey25", vjust = 0, lineheight = 0.9) +
  scale_fill_manual(values = bar_cols,
                    labels = c("GWAS core" = "GWAS-significant", "Flank" = "+/- 100 kb flank", "Retained" = "Retained"),
                    name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.18))) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor    = element_blank(),
        panel.border        = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        legend.position     = "right",
        legend.text         = element_text(size = 9),
        axis.title.x        = element_blank()) +
  labs(tag = "b", y = "Fraction of top 1% windows", title = "Top 1% most heavily loaded windows")

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

set.seed(20260415)
pcs    <- paste0("PC", 1:14)
n_perm <- 1000L

out_tsv <- file.path("results", "gwas_removed_trait_pc_category_absrho_labelperm_pc1_pc14.tsv")
out_fig <- file.path("figures", "figureS18_gwas_removed_trait_pc_category_absrho_labelperm_pc1_pc14.pdf")
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)

# ---- Rebuild traits_t_scaled ----
message("FigureS18: rebuilding scaled trait matrix (~2-3 min)...")
source(file.path("post-process-scripts", "helpers-50k-matrix.R"))
all_betas <- fread(file.path("results", "all-trait-50k-mean-effects.tsv"))
all_pvals <- fread(file.path("results", "all-trait-50k-mean-pvals.tsv"))
prep <- prepare_X_50k(all_betas, all_pvals)
traits_t_scaled <- t(prep$X)  # traits x windows

# ---- Category colors ----
trait_categories <- fread(file.path("data", "trait_abbrevs_categorised.txt"))
custom_theme_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
  "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
  "#E5C494", "#B3B3B3"
)
all_trait_categories <- sort(unique(na.omit(as.character(trait_categories$phenotype_group))))
all_trait_categories <- all_trait_categories[nzchar(all_trait_categories)]
category_colors_all  <- setNames(custom_theme_colors[seq_along(all_trait_categories)], all_trait_categories)

conditions <- list(
  `+/-100 kb` = file.path("results", "gwas_removed_distance_topload_50k_plus100kb", "reduced_pca_window_loadings.tsv"),
  `+/-150 kb` = file.path("results", "gwas_removed_distance_topload_50k_plus150kb", "reduced_pca_window_loadings.tsv")
)

row_cor_with_vec <- function(X, y) {
  n       <- ncol(X); row_mean <- rowMeans(X); y_mean <- mean(y)
  x_ss    <- rowSums((X - row_mean)^2); y_ss <- sum((y - y_mean)^2)
  cov_num <- as.numeric(X %*% y) - n * row_mean * y_mean
  out     <- cov_num / sqrt(x_ss * y_ss)
  out[!is.finite(out)] <- NA_real_
  out
}

compute_one <- function(cond_label, loadings_file) {
  if (!file.exists(loadings_file)) stop("Missing reduced loadings: ", loadings_file)
  red        <- fread(loadings_file)
  window_col <- if ("window_id" %in% names(red)) "window_id" else "region_label"
  common_ids <- intersect(colnames(traits_t_scaled), red[[window_col]])
  if (length(common_ids) == 0L) stop("No shared windows for ", cond_label)

  trait_mat    <- as.matrix(traits_t_scaled[, common_ids, drop = FALSE])
  pc_mat       <- as.matrix(red[match(common_ids, get(window_col)), ..pcs])
  row_rank_mat <- t(apply(trait_mat, 1, rank, ties.method = "average"))
  pc_rank_mat  <- apply(pc_mat, 2, rank, ties.method = "average")
  if (is.null(dim(pc_rank_mat))) pc_rank_mat <- matrix(pc_rank_mat, ncol = length(pcs))
  colnames(pc_rank_mat) <- pcs

  rho_rows    <- vector("list", length(pcs))
  trait_names <- rownames(trait_mat)
  for (i in seq_along(pcs)) {
    rho_rows[[i]] <- data.table(trait_abbrev = trait_names, pc = pcs[[i]],
                                 spearman_rho = row_cor_with_vec(row_rank_mat, pc_rank_mat[, i]))
  }
  dt <- merge(rbindlist(rho_rows, use.names = TRUE),
              as.data.table(trait_categories)[, .(trait_abbrev, phenotype_group)],
              by = "trait_abbrev", all.x = TRUE)
  dt <- dt[is.finite(spearman_rho)]
  dt[is.na(phenotype_group) | phenotype_group == "", phenotype_group := "Uncategorized"]
  dt[, abs_rho := abs(spearman_rho)]
  dt[, abs_rho_rank := frank(abs_rho, ties.method = "average"), by = pc]

  obs <- dt[, .(n_traits = .N, observed_median_abs_rho = median(abs_rho, na.rm = TRUE),
                observed_mean_rank = mean(abs_rho_rank, na.rm = TRUE)), by = .(pc, phenotype_group)]

  trait_lab  <- unique(dt[, .(trait_abbrev, phenotype_group)])
  perm_rows  <- vector("list", n_perm)
  for (b in seq_len(n_perm)) {
    perm_map <- copy(trait_lab)
    perm_map[, phenotype_group := sample(phenotype_group, .N, replace = FALSE)]
    z <- merge(dt[, .(trait_abbrev, pc, abs_rho_rank)], perm_map, by = "trait_abbrev", all.x = FALSE)
    perm_rows[[b]] <- z[, .(perm_mean_rank = mean(abs_rho_rank, na.rm = TRUE)),
                         by = .(pc, phenotype_group)][, perm := b]
  }
  perm_dt  <- rbindlist(perm_rows, use.names = TRUE)
  null_sum <- perm_dt[, .(null_mean = mean(perm_mean_rank, na.rm = TRUE),
                           null_sd   = sd(perm_mean_rank, na.rm = TRUE)), by = .(pc, phenotype_group)]
  p_hi <- perm_dt[obs, on = .(pc, phenotype_group)][
    , .(p_emp_high = (1 + sum(perm_mean_rank >= observed_mean_rank, na.rm = TRUE)) / (.N + 1)),
    by = .(pc, phenotype_group)]

  res <- Reduce(function(x, y) merge(x, y, by = c("pc", "phenotype_group"), all = TRUE), list(obs, null_sum, p_hi))
  res[, z_score := fifelse(is.finite(null_sd) & null_sd > 0, (observed_mean_rank - null_mean) / null_sd, NA_real_)]
  res[, sig := fifelse(p_emp_high < 0.001, "***",
                       fifelse(p_emp_high < 0.01, "**", fifelse(p_emp_high < 0.05, "*", "")))]
  res[, condition := cond_label]
  res
}

message("FigureS18: computing permutations (~2 per condition)...")
res <- rbindlist(Map(compute_one, names(conditions), conditions), use.names = TRUE)

fig1_order_file <- file.path("results", "trait_pc_category_absrho_labelperm_pc1_pc14.tsv")
if (!file.exists(fig1_order_file)) stop("Missing Figure 1 heatmap table for row order: ", fig1_order_file)
fig1_order_dt <- fread(fig1_order_file)
group_order   <- unique(fig1_order_dt$phenotype_group)
res[, phenotype_group := factor(phenotype_group, levels = group_order)]
res[, pc        := factor(pc, levels = pcs)]
res[, condition := factor(condition, levels = names(conditions))]

fwrite(res[order(condition, pc, phenotype_group)], out_tsv, sep = "\t")

res[, label_rho := sprintf("%.2f", observed_median_abs_rho)]
res[, label_col := fifelse(observed_median_abs_rho > 0.12, "white", "black")]

rho_cap       <- 0.2
ag_sunset_cols <- rev(colorRampPalette(as.character(paletteer::paletteer_d("rcartocolor::ag_Sunset")))(256))

make_panel <- function(cond_label) {
  dt <- res[condition == cond_label]
  dt[, pc_index := as.integer(sub("^PC", "", as.character(pc)))]
  strip_dt <- unique(dt[, .(phenotype_group)])
  strip_dt[, strip_col := unname(category_colors_all[as.character(phenotype_group)])]
  ggplot(dt, aes(x = pc_index, y = phenotype_group, fill = observed_median_abs_rho)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_tile(data = strip_dt, aes(x = 0.35, y = phenotype_group), inherit.aes = FALSE,
              fill = strip_dt$strip_col, width = 0.18, height = 0.98, linewidth = 0) +
    geom_text(aes(label = label_rho, color = label_col), size = 3.1, na.rm = TRUE) +
    geom_text(aes(label = sig, color = label_col), size = 2.6, nudge_y = 0.28, na.rm = TRUE) +
    scale_color_identity() +
    scale_fill_gradientn(colours = ag_sunset_cols, limits = c(0, rho_cap),
                         oob = scales::squish, na.value = "grey90", name = "Median |rho|") +
    scale_x_continuous(breaks = 1:14, labels = paste0("PC", 1:14),
                       expand = expansion(add = c(0.45, 0.2))) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), axis.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 15),
          axis.text.y = element_text(size = 15),
          plot.title = element_text(size = 11, face = "bold"),
          legend.position = "none", plot.margin = margin(5.5, 5.5, 5.5, 12)) +
    labs(title = cond_label)
}

p1 <- make_panel("+/-100 kb") + labs(tag = "c")
p2 <- make_panel("+/-150 kb") + labs(tag = "d") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

legend_plot <- ggplot(data.frame(x = c(0, 1), y = c(0, 1), z = c(0, rho_cap)), aes(x, y, fill = z)) +
  geom_tile() +
  scale_fill_gradientn(colours = ag_sunset_cols, limits = c(0, rho_cap), oob = scales::squish, name = "Median |rho|") +
  guides(fill = guide_colorbar(barheight = unit(35, "mm"))) + theme_void() + theme(legend.position = "right")
leg <- cowplot::get_legend(legend_plot)

# Top row: density block (4 units) + bar chart (2 units), side by side
top_fig <- cowplot::plot_grid(top_density_block, mid_panel,
                              ncol = 2, rel_widths = c(4, 2))

bottom_fig <- (p1 | p2 | cowplot::ggdraw(leg)) + plot_layout(widths = c(1, 1, 0.16))

fig <- top_fig / bottom_fig +
  plot_layout(heights = c(1, 1.6)) +
  plot_annotation(caption = "* p<0.05, ** p<0.01, *** p<0.001")

ggsave(out_fig, fig, width = 13.97, height = 10.0)
message("Wrote: ", out_tsv)
message("FigureS18 complete: ", out_fig)
