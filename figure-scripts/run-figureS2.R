#!/usr/bin/env Rscript
# Produces: figures/figureS2_robustness_checks.pdf
# Requires: results/svd_robustness/*.tsv  (from compute-svd-robustness.R)
# Run: Rscript figure-scripts/run-figureS2.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(paletteer)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

in_dir  <- file.path("results", "svd_robustness")
fig_dir <- "figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

for (f in c("trait_subsample_pc_cosine.tsv", "trait_subsample_procrustes.tsv",
            "leave_one_out_observed.tsv", "leave_one_out_null.tsv",
            "window_subsample_pc_cosine.tsv", "window_subsample_procrustes.tsv")) {
  if (!file.exists(file.path(in_dir, f)))
    stop("Missing: ", file.path(in_dir, f), "\nRun compute-svd-robustness.R first.")
}

cosine_trait <- fread(file.path(in_dir, "trait_subsample_pc_cosine.tsv"))
proc_trait   <- fread(file.path(in_dir, "trait_subsample_procrustes.tsv"))
loo_obs  <- fread(file.path(in_dir, "leave_one_out_observed.tsv"))
loo_null <- fread(file.path(in_dir, "leave_one_out_null.tsv"))

# Compatibility: reshape wide format (old) to long format (new)
if ("pc1_cosine" %in% names(loo_obs) && !"PC" %in% names(loo_obs)) {
  pc_cols <- grep("^pc[0-9]+_cosine$", names(loo_obs), value = TRUE)
  loo_obs <- melt(loo_obs, id.vars = c("category", "n_removed", "procrustes_k4", "procrustes_k8"),
                  measure.vars = pc_cols, variable.name = "PC", value.name = "cosine")
  loo_obs[, PC := toupper(gsub("_cosine$", "", PC))]
  loo_obs[, PC := gsub("PC0*([0-9]+)", "PC\\1", PC)]
}
if ("pc1_cosine" %in% names(loo_null) && !"PC" %in% names(loo_null)) {
  pc_cols <- grep("^pc[0-9]+_cosine$", names(loo_null), value = TRUE)
  loo_null <- melt(loo_null, id.vars = c("category", "n_removed", "rep", "procrustes_k4", "procrustes_k8"),
                   measure.vars = pc_cols, variable.name = "PC", value.name = "cosine")
  loo_null[, PC := toupper(gsub("_cosine$", "", PC))]
  loo_null[, PC := gsub("PC0*([0-9]+)", "PC\\1", PC)]
}
cosine_win   <- fread(file.path(in_dir, "window_subsample_pc_cosine.tsv"))
proc_win     <- fread(file.path(in_dir, "window_subsample_procrustes.tsv"))

# ---- Shared settings ----
pc_order     <- paste0("PC", 1:8)
pc_order_loo <- paste0("PC", 1:8)
n_reps_trait <- max(cosine_trait$rep, na.rm = TRUE)
n_reps_win   <- max(cosine_win$rep,   na.rm = TRUE)
n_reps_null  <- max(loo_null$rep,     na.rm = TRUE)

cosine_trait[, PC := factor(PC, levels = pc_order)]
cosine_win[,   PC := factor(PC, levels = pc_order)]
loo_obs[,      PC := factor(PC, levels = pc_order_loo)]
loo_null[,     PC := factor(PC, levels = pc_order_loo)]

col_random   <- "#377EB8"
col_strat    <- "#E41A1C"
col_win_full <- "#4DAF4A"
col_win_gd   <- "#FF7F00"

# Match fig 1d colour palette exactly (ag_Sunset reversed = low->cool, high->warm)
ag_sunset_cols <- rev(colorRampPalette(as.character(
  paletteer::paletteer_d("rcartocolor::ag_Sunset")))(256))
# For LOO panel: invert so low cosine (big drop) = warm, high cosine = cool
loo_cosine_cols <- rev(ag_sunset_cols)

base_theme <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position  = "bottom",
        legend.key.size  = unit(0.4, "cm"))

# ============================================================
# Panel A: per-PC cosine similarity, trait subsampling
# ============================================================
cos_summ_trait <- cosine_trait[PC %in% pc_order,
  .(med = median(cosine, na.rm = TRUE),
    lo  = quantile(cosine, 0.025, na.rm = TRUE),
    hi  = quantile(cosine, 0.975, na.rm = TRUE)),
  by = .(PC, analysis)]

pA <- ggplot(cos_summ_trait,
             aes(x = PC, y = med, ymin = lo, ymax = hi,
                 colour = analysis, fill = analysis, group = analysis)) +
  geom_ribbon(alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.8) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0.01)) +
  scale_colour_manual(values = c(random_trait     = col_random,
                                 stratified_trait  = col_strat),
                      labels = c(random_trait     = paste0("Random 2/3 (n=", n_reps_trait, ")"),
                                 stratified_trait  = paste0("Stratified 2/3 (n=", n_reps_trait, ")"))) +
  scale_fill_manual(values   = c(random_trait     = col_random,
                                 stratified_trait  = col_strat),
                    labels   = c(random_trait     = paste0("Random 2/3 (n=", n_reps_trait, ")"),
                                 stratified_trait  = paste0("Stratified 2/3 (n=", n_reps_trait, ")"))) +
  labs(title  = "a   Trait subsampling: per-PC cosine similarity",
       x = NULL, y = "Cosine similarity to original",
       colour = NULL, fill = NULL) +
  base_theme

# ============================================================
# Panel B: Procrustes, trait subsampling
# ============================================================
proc_summ_trait <- proc_trait[,
  .(med = median(procrustes, na.rm = TRUE),
    lo  = quantile(procrustes, 0.025, na.rm = TRUE),
    hi  = quantile(procrustes, 0.975, na.rm = TRUE)),
  by = .(k, analysis)]
proc_summ_trait[, k_label := paste0("k = ", k)]

pB <- ggplot(proc_summ_trait,
             aes(x = analysis, y = med, ymin = lo, ymax = hi, colour = analysis)) +
  geom_errorbar(width = 0.25, linewidth = 0.7) +
  geom_point(size = 2.5) +
  facet_wrap(~k_label, nrow = 1) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete(labels = c(random_trait    = "Random",
                               stratified_trait = "Stratified")) +
  scale_colour_manual(values = c(random_trait    = col_random,
                                 stratified_trait = col_strat),
                      guide = "none") +
  labs(title = "b   Trait subsampling: Procrustes similarity",
       x = NULL, y = "Procrustes similarity to original") +
  base_theme

# ============================================================
# Panel C: per-PC cosine similarity, window subsampling
# Full matrix and GWAS-depleted on same axes
# ============================================================
datasets_present <- unique(cosine_win$dataset)
cos_summ_win <- cosine_win[PC %in% pc_order,
  .(med = median(cosine, na.rm = TRUE),
    lo  = quantile(cosine, 0.025, na.rm = TRUE),
    hi  = quantile(cosine, 0.975, na.rm = TRUE)),
  by = .(PC, dataset)]

win_col_map   <- c(full = col_win_full, gwas_depleted = col_win_gd)
win_label_map <- c(full          = paste0("Full matrix (n=", n_reps_win, ")"),
                   gwas_depleted = paste0("GWAS-depleted (n=", n_reps_win, ")"))

pC <- ggplot(cos_summ_win,
             aes(x = PC, y = med, ymin = lo, ymax = hi,
                 colour = dataset, fill = dataset, group = dataset)) +
  geom_ribbon(alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.8) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0.01)) +
  scale_colour_manual(values = win_col_map, labels = win_label_map,
                      limits = intersect(names(win_col_map), datasets_present)) +
  scale_fill_manual(values   = win_col_map, labels = win_label_map,
                    limits   = intersect(names(win_col_map), datasets_present)) +
  labs(title    = "c   Window subsampling: per-PC cosine similarity",
       subtitle = "Reference for each dataset is its own full-set PCA",
       x = NULL, y = "Cosine similarity to reference",
       colour = NULL, fill = NULL) +
  base_theme

# ============================================================
# Panel D: Procrustes, window subsampling
# ============================================================
proc_summ_win <- proc_win[,
  .(med = median(procrustes, na.rm = TRUE),
    lo  = quantile(procrustes, 0.025, na.rm = TRUE),
    hi  = quantile(procrustes, 0.975, na.rm = TRUE)),
  by = .(k, dataset)]
proc_summ_win[, k_label := paste0("k = ", k)]

pD <- ggplot(proc_summ_win,
             aes(x = k_label, y = med, ymin = lo, ymax = hi, colour = dataset)) +
  geom_errorbar(width = 0.25, linewidth = 0.7,
                position = position_dodge(width = 0.4)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.4)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = win_col_map, labels = win_label_map,
                      limits = intersect(names(win_col_map), datasets_present)) +
  labs(title = "d   Window subsampling: Procrustes similarity",
       x = NULL, y = "Procrustes similarity to reference",
       colour = NULL) +
  base_theme

# ============================================================
# Panel E: Leave-one-category-out, PC1-PC8 faceted
# Sorted by PC1 observed cosine (ascending = most disruptive at top)
# Colour encodes cosine value: low (drop) = warm; same palette as fig 1d
# ============================================================
cat_order <- loo_obs[PC == "PC1"][order(cosine), category]
loo_obs[,  category_f := factor(category, levels = cat_order)]
loo_null[, category_f := factor(category, levels = cat_order)]

n_label <- unique(loo_obs[PC == "PC1", .(category_f, n_removed)])
n_label[, label := paste0("n=", n_removed)]

null_summ <- loo_null[,
  .(null_med = median(cosine, na.rm = TRUE),
    null_lo  = quantile(cosine, 0.025, na.rm = TRUE),
    null_hi  = quantile(cosine, 0.975, na.rm = TRUE)),
  by = .(category_f, PC)]

loo_plot <- merge(loo_obs,   null_summ, by = c("category_f", "PC"))
loo_plot <- merge(loo_plot,  n_label,   by = "category_f")

x_lo <- max(0, floor(min(c(loo_plot$cosine, loo_plot$null_lo), na.rm = TRUE) * 20) / 20 - 0.02)

pE <- ggplot(loo_plot, aes(y = category_f)) +
  geom_errorbar(aes(xmin = null_lo, xmax = null_hi),
                orientation = "y", width = 0.4,
                colour = "grey65", linewidth = 0.5) +
  geom_point(aes(x = null_med), colour = "grey65", size = 1.5) +
  geom_point(aes(x = cosine, colour = cosine), size = 2.2) +
  geom_text(data = loo_plot[PC == "PC1"],
            aes(x = x_lo + 0.001, label = label),
            hjust = 0, size = 2.0, colour = "grey45") +
  facet_wrap(~PC, nrow = 2) +
  scale_x_continuous(limits = c(x_lo, 1),
                     breaks = pretty(c(x_lo, 1), n = 4),
                     expand = c(0, 0.01)) +
  scale_colour_gradientn(colours  = loo_cosine_cols,
                         limits   = c(0, 1),
                         oob      = scales::squish,
                         name     = "Cosine\nsimilarity",
                         guide    = guide_colourbar(barwidth = 0.5, barheight = 3)) +
  labs(title    = "e   Leave-one-category-out: cosine similarity (PC1-PC8)",
       subtitle = paste0("Colour = cosine to original; grey = size-matched null (",
                         n_reps_null, " reps, median + 95%)"),
       x = "Cosine similarity to original", y = NULL) +
  base_theme +
  theme(legend.position = "right",
        strip.text      = element_text(face = "bold"))

# ============================================================
# Assemble and save
# ============================================================
fig <- (pA | pB) / (pC | pD) / pE +
  plot_layout(heights = c(1, 1, 1.5)) +
  plot_annotation(
    title = "SVD robustness checks",
    theme = theme(plot.title = element_text(size = 11, face = "bold"))
  )

out_path <- file.path(fig_dir, "figureS2_robustness_checks.pdf")
ggsave(out_path, fig, width = 12, height = 14)
message("Saved: ", out_path)
