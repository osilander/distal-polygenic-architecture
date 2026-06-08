#!/usr/bin/env Rscript
# Produces: figures/figureS10_within_over_between_abs_cosine_k4_k8_normalized_expanded.pdf
#           results/pca_multiscale_anchored/within_over_between_*.tsv
# Pre-computed inputs:
#   results/pca_multiscale_anchored/pca_traits_{25k,50k,100k,200k,500k,1m}_anchored.tsv
#   results/figureS_pca_absolute_50k_trait_scores.tsv
#   results/pca_abs50k_from_windows/trait_scores.tsv
#   results/pca_mean_effect/trait_scores_{50k,100k}_mean.tsv  (from compute-mean-effect-pca.R)
#   data/trait_abbrevs_categorised.txt
# Run: Rscript figure-scripts/run-figureS10.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(irlba)
  library(ggrepel)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

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
category_colors      <- setNames(custom_theme_colors[seq_along(all_trait_categories)], all_trait_categories)

trait_map <- as.data.table(trait_categories)[, .(trait_abbrev, phenotype_group)]
setkey(trait_map, trait_abbrev)

in_pca_dir      <- file.path("results", "pca_multiscale_anchored")
in_abs_signed   <- file.path("results", "figureS_pca_absolute_50k_trait_scores.tsv")
in_abs_variant  <- file.path("results", "pca_abs50k_from_windows", "trait_scores.tsv")
out_dir         <- in_pca_dir
out_summary     <- file.path(out_dir, "within_over_between_summary_abs_cosine_k4_k8_normalized_expanded.tsv")
out_trait       <- file.path(out_dir, "within_over_between_traitlevel_abs_cosine_k4_k8_expanded.tsv")
out_summary_dist <- file.path(out_dir, "within_over_between_summary_abs_cosine_distance_k4_k8_normalized_expanded.tsv")
out_plot_cache  <- file.path(out_dir, "within_over_between_plotdata_abs_cosine_k4_k8_normalized_expanded.tsv")
out_ann_cache   <- file.path(out_dir, "within_over_between_plot_annotations_abs_cosine_k4_k8_normalized_expanded.tsv")
out_fig         <- file.path("figures", "figureS10_within_over_between_abs_cosine_k4_k8_normalized_expanded.pdf")
k_values        <- c(4L)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

scale_order    <- c("25k", "50k", "50k_mean", "100k", "100k_mean",
                    "200k", "500k", "1m",
                    "abs_from_signed_matrix", "abs_from_windows_filtered")
base_scale_order <- c("25k", "50k", "100k", "200k", "500k", "1m")
scale_labels   <- c("25k" = "25Kb", "50k" = "50Kb", "50k_mean" = "50Kb mean",
                    "100k" = "100Kb", "100k_mean" = "100Kb mean",
                    "200k" = "200Kb", "500k" = "500Kb", "1m" = "1Mb",
                    "abs_from_signed_matrix"    = "Matrix abs. value",
                    "abs_from_windows_filtered" = "Variant abs value")

score_files <- c(
  "25k"  = file.path(in_pca_dir, "pca_traits_25k_anchored.tsv"),
  "50k"  = file.path(in_pca_dir, "pca_traits_50k_anchored.tsv"),
  "100k" = file.path(in_pca_dir, "pca_traits_100k_anchored.tsv"),
  "200k" = file.path(in_pca_dir, "pca_traits_200k_anchored.tsv"),
  "500k" = file.path(in_pca_dir, "pca_traits_500k_anchored.tsv"),
  "1m"   = file.path(in_pca_dir, "pca_traits_1m_anchored.tsv"),
  "abs_from_signed_matrix"    = in_abs_signed,
  "abs_from_windows_filtered" = in_abs_variant,
  "50k_mean"  = file.path("results", "pca_mean_effect", "trait_scores_50k_mean.tsv"),
  "100k_mean" = file.path("results", "pca_mean_effect", "trait_scores_100k_mean.tsv")
)

load_scores <- function(path, k = 4L) {
  if (!file.exists(path)) stop("Missing file: ", path)
  dt        <- fread(path)
  trait_col <- intersect(c("trait_abbrev", "trait", "Trait"), names(dt))[1]
  if (is.na(trait_col)) stop("No trait id column in: ", path)
  setnames(dt, trait_col, "trait_abbrev")
  pcs  <- paste0("PC", seq_len(k))
  miss <- setdiff(pcs, names(dt))
  if (length(miss) > 0L) stop("Missing PCs in ", basename(path), ": ", paste(miss, collapse = ", "))
  out  <- dt[, c("trait_abbrev", pcs), with = FALSE]
  for (pc in pcs) set(out, j = pc, value = as.numeric(out[[pc]]))
  out  <- merge(out, trait_map, by = "trait_abbrev", all.x = FALSE, all.y = FALSE)
  out[!is.na(phenotype_group)]
}

abs_cosine <- function(a, b) {
  na <- sqrt(sum(a * a)); nb <- sqrt(sum(b * b))
  if (!is.finite(na) || !is.finite(nb) || na == 0 || nb == 0) return(NA_real_)
  abs(sum(a * b) / (na * nb))
}

darken_hex <- function(cols, factor = 0.55) {
  vapply(cols, function(cc) {
    rgb_vec <- as.numeric(grDevices::col2rgb(cc, alpha = FALSE)) / 255
    darker  <- pmax(0, pmin(1, rgb_vec * factor))
    grDevices::rgb(darker[1], darker[2], darker[3])
  }, character(1))
}

compute_trait_level <- function(dt, scale_id, k = 4L) {
  pcs <- paste0("PC", seq_len(k)); n <- nrow(dt)
  rows <- vector("list", n)
  for (i in seq_len(n)) {
    tr <- dt$trait_abbrev[i]; grp <- dt$phenotype_group[i]
    a  <- as.numeric(dt[i, ..pcs])
    others <- dt[trait_abbrev != tr]
    if (nrow(others) == 0L) next
    cosv         <- apply(others[, ..pcs], 1L, function(z) abs_cosine(a, as.numeric(z)))
    within_idx   <- which(others$phenotype_group == grp)
    between_idx  <- which(others$phenotype_group != grp)
    within_med   <- if (length(within_idx)  > 0L) median(cosv[within_idx],  na.rm = TRUE) else NA_real_
    between_med  <- if (length(between_idx) > 0L) median(cosv[between_idx], na.rm = TRUE) else NA_real_
    ratio        <- if (is.finite(within_med) && is.finite(between_med) && between_med > 0) within_med / between_med else NA_real_
    rows[[i]] <- data.table(
      scale = scale_id, scale_label = scale_labels[[scale_id]], k = k,
      trait_abbrev = tr, phenotype_group = grp,
      n_within = length(within_idx), n_between = length(between_idx),
      within_median = within_med, between_median = between_med, within_over_between = ratio
    )
  }
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

cache_ok <- file.exists(out_plot_cache) && file.exists(out_ann_cache) &&
            file.exists(out_summary)    && file.exists(out_trait)     &&
            file.exists(out_summary_dist)
if (cache_ok) {
  cached_scales <- unique(fread(out_trait, select = "scale")$scale)
  if (!all(names(score_files) %in% cached_scales)) {
    message("FigureS10: cache missing scales ", paste(setdiff(names(score_files), cached_scales), collapse = ", "), " — recomputing.")
    cache_ok <- FALSE
  }
}
if (cache_ok) {
  message("FigureS10: using cached within/between plot data...")
  plot_dt    <- fread(out_plot_cache)
  summary_dt <- fread(out_summary)
  trait_level <- fread(out_trait)
} else {
  trait_rows <- list(); ii <- 1L
  for (k in k_values) {
    for (sc in names(score_files)) {
      message("FigureS10: within/between k=", k, ": ", sc, " ...")
      dt <- load_scores(score_files[[sc]], k = k)
      trait_rows[[ii]] <- compute_trait_level(dt, sc, k = k)
      ii <- ii + 1L
    }
  }
  trait_level <- rbindlist(trait_rows, use.names = TRUE, fill = TRUE)
  fwrite(trait_level, out_trait, sep = "\t")

  summary_dt <- trait_level[is.finite(within_over_between), .(
    n_traits = .N, median_ratio = median(within_over_between, na.rm = TRUE),
    mean_ratio = mean(within_over_between, na.rm = TRUE),
    q25 = quantile(within_over_between, 0.25, na.rm = TRUE),
    q75 = quantile(within_over_between, 0.75, na.rm = TRUE),
    q05 = quantile(within_over_between, 0.05, na.rm = TRUE),
    q95 = quantile(within_over_between, 0.95, na.rm = TRUE),
    median_within = median(within_median, na.rm = TRUE),
    median_between = median(between_median, na.rm = TRUE)
  ), by = .(scale, scale_label, k, phenotype_group)]
  base_means <- summary_dt[scale %in% base_scale_order,
                            .(cat_mean = mean(median_ratio, na.rm = TRUE)), by = .(k, phenotype_group)]
  summary_dt <- merge(summary_dt, base_means, by = c("k", "phenotype_group"), all.x = TRUE)
  summary_dt[, median_ratio_norm := median_ratio / cat_mean]
  fwrite(summary_dt, out_summary, sep = "\t")

  plot_dt <- copy(summary_dt)
  plot_dt[, scale       := factor(scale, levels = scale_order)]
  plot_dt[, scale_label := factor(scale_label, levels = scale_labels[scale_order])]
  plot_dt[, k_label     := factor(paste0("k=", k), levels = paste0("k=", k_values))]
  plot_dt[, `:=`(median_within_dist = 1 - median_within, median_between_dist = 1 - median_between)]
  plot_dt[, metric_plot_raw := fifelse(
    is.finite(median_between_dist) & median_between_dist > 0,
    median_within_dist / median_between_dist, NA_real_
  )]
  dist_base_plot <- plot_dt[scale %in% base_scale_order,
                             .(cat_mean_dist_plot = mean(metric_plot_raw, na.rm = TRUE)),
                             by = .(k, phenotype_group)]
  plot_dt <- merge(plot_dt, dist_base_plot, by = c("k", "phenotype_group"), all.x = TRUE)
  plot_dt[, metric_plot := metric_plot_raw / cat_mean_dist_plot]

  med_by_scale <- plot_dt[, .(mean_norm = mean(metric_plot, na.rm = TRUE)), by = .(k_label, scale, scale_label)]
  ymax <- max(plot_dt$metric_plot, na.rm = TRUE) * 1.15
  ann  <- med_by_scale[, .(k_label, scale_label, y = ymax * 0.98, label = sprintf("mean %.2f", mean_norm))]

  fwrite(plot_dt, out_plot_cache, sep = "\t")
  fwrite(ann, out_ann_cache, sep = "\t")
}

plot_dt <- plot_dt[k == 4]
plot_dt[, scale       := factor(scale, levels = scale_order)]
plot_dt[, scale_label := factor(scale_label, levels = scale_labels[scale_order])]
plot_dt[, k_label     := factor("", levels = "")]
plot_long <- rbindlist(list(
  copy(plot_dt)[, .(scale, scale_label, k_label, phenotype_group, metric = metric_plot_raw,  panel_type = "Unnormalized")],
  copy(plot_dt)[, .(scale, scale_label, k_label, phenotype_group, metric = metric_plot,      panel_type = "Category-normalized")]
), use.names = TRUE)
plot_long[, panel_type := factor(panel_type, levels = c("Unnormalized", "Category-normalized"))]

ann <- plot_long[is.finite(metric),
                  .(mean_metric = mean(metric, na.rm = TRUE)),
                  by = .(panel_type, k_label, scale, scale_label)]
ann[, `:=`(y = 0.02, label = sprintf("%.2f", mean_metric))]

label_dt <- plot_long[scale == tail(scale_order, 1L), .(scale_label, metric, phenotype_group, k_label, panel_type)]
label_colors <- setNames(
  darken_hex(unname(category_colors[names(category_colors) %in% unique(plot_dt$phenotype_group)]), factor = 0.55),
  names(category_colors)[names(category_colors) %in% unique(plot_dt$phenotype_group)]
)
label_dt[, label_colour := label_colors[phenotype_group]]

p <- ggplot(plot_long, aes(x = scale_label, y = metric, colour = phenotype_group, group = phenotype_group)) +
  geom_line(alpha = 0.9, linewidth = 0.75) +
  geom_point(aes(fill = phenotype_group), size = 1.8, alpha = 0.95, shape = 21, stroke = 0.22, colour = "#222222") +
  geom_text(data = ann, aes(x = scale_label, y = y, label = label), inherit.aes = FALSE,
            colour = "black", size = 3.6, fontface = "plain", vjust = 1) +
  geom_text_repel(data = label_dt, aes(x = scale_label, y = metric, label = phenotype_group, colour = I(label_colour)),
                  inherit.aes = FALSE, show.legend = FALSE, direction = "y", hjust = 0, nudge_x = 1.15,
                  size = 3.0, segment.color = "grey65", segment.alpha = 0.5, segment.size = 0.25,
                  box.padding = 0.12, point.padding = 0.08, min.segment.length = 0, seed = 20260331, max.overlaps = Inf) +
  facet_grid(panel_type ~ k_label, scales = "free_y") +
  scale_color_manual(values = category_colors, drop = FALSE) +
  scale_fill_manual(values = category_colors, drop = FALSE) +
  scale_x_discrete(expand = expansion(add = c(0.1, 2.3))) +
  labs(x = NULL, y = "Within/Between |cosine| distance ratio") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 14)),
        strip.text.x = element_blank(),
        strip.text.y = element_text(face = "plain", size = 14),
        legend.position = "right", legend.title = element_blank(),
        plot.margin = margin(t = 5.5, r = 92, b = 5.5, l = 22)) +
  coord_cartesian(clip = "off")

ggsave(out_fig, p, width = 13.6, height = 7.8)

dist_summary <- copy(summary_dt)
dist_summary[, `:=`(median_within_dist = 1 - median_within, median_between_dist = 1 - median_between)]
dist_summary[, median_ratio_dist := fifelse(is.finite(median_between_dist) & median_between_dist > 0,
                                             median_within_dist / median_between_dist, NA_real_)]
dist_summary[, median_ratio_dist_cohesion := fifelse(is.finite(median_within_dist) & median_within_dist > 0,
                                                      median_between_dist / median_within_dist, NA_real_)]
dist_base <- dist_summary[scale %in% base_scale_order,
                           .(cat_mean_dist = mean(median_ratio_dist, na.rm = TRUE)), by = .(k, phenotype_group)]
dist_summary <- merge(dist_summary, dist_base, by = c("k", "phenotype_group"), all.x = TRUE)
dist_summary[, median_ratio_dist_norm := median_ratio_dist / cat_mean_dist]
dist_summary[, median_ratio_dist_cohesion_norm := 1 / median_ratio_dist_norm]
fwrite(dist_summary, out_summary_dist, sep = "\t")
fwrite(summary_dt[k == 4], file.path(out_dir, "within_over_between_summary_abs_cosine_k4_normalized_expanded.tsv"), sep = "\t")
if (nrow(summary_dt[k == 8]) > 0L) {
  fwrite(summary_dt[k == 8], file.path(out_dir, "within_over_between_summary_abs_cosine_k8_normalized_expanded.tsv"), sep = "\t")
}

message("FigureS10 complete: ", out_fig)
