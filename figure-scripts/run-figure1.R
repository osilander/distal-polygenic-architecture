#!/usr/bin/env Rscript
# Produces: figures/figure1.pdf (6-panel), figures/figureS6_trait_loading_heatmap_top60_pc1_pc5.pdf,
#           figures/figureS5_umap.pdf, results/pca_loadings_50k.tsv, results/umap-coordinates.tsv
# Pre-computed inputs required:
#   results/all-trait-50k-mean-effects.tsv, results/all-trait-50k-mean-pvals.tsv,
#   data/trait_abbrevs_categorised.txt,
#   results/figure2_cache_50k.rds (†HPC), results/trait_pc_category_absrho_labelperm_pc1_pc14.tsv

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(irlba)
  library(patchwork)
  library(RColorBrewer)
  library(scales)
  library(circlize)
  library(ComplexHeatmap)
  library(dbscan)
  library(uwot)
  library(cowplot)
  library(paletteer)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

fig2_cache <- file.path("results", "figure2_cache_50k.rds")
heat_tsv   <- file.path("results", "trait_pc_category_absrho_labelperm_pc1_pc14.tsv")
if (!file.exists(fig2_cache))
  stop("Missing pre-computed scree cache: ", fig2_cache, "\nSee README.")
if (!file.exists(heat_tsv))
  stop("Missing pre-computed category heatmap table: ", heat_tsv, "\nSee README.")

# Staleness guard: heat_tsv must not predate any PCA input source.
# If pca_loadings_50k.tsv or the raw effect/p-value tables are newer than
# heat_tsv, the heatmap was derived from an older PCA state — abort.
# Design note: this script both computes a fresh PCA and uses the pre-computed
# heat_tsv for panel d.  The guard catches any pre-existing staleness before the
# run starts.  After the PCA is written (below), heat_tsv will be timestamp-older
# than the freshly written pca_loadings_50k.tsv; the NOTE message at that point
# reminds you to re-run compute-trait-pc-category-labelperm.R if the PCA results
# actually changed.  If the underlying data has not changed, the fresh PCA is
# content-identical to the previous one and the existing heat_tsv remains valid.
.pca_inputs_for_stale_check <- c(
  file.path("results", "pca_loadings_50k.tsv"),
  file.path("results", "all-trait-50k-mean-effects.tsv"),
  file.path("results", "all-trait-50k-mean-pvals.tsv")
)
.existing_pca_inputs <- .pca_inputs_for_stale_check[file.exists(.pca_inputs_for_stale_check)]
if (length(.existing_pca_inputs) > 0L) {
  .newest_input_mtime <- max(file.mtime(.existing_pca_inputs))
  if (file.mtime(heat_tsv) < .newest_input_mtime) {
    stop("'", heat_tsv, "' predates one or more PCA input files — it may be stale.\n",
         "Re-run compute-trait-pc-category-labelperm.R to regenerate it, ",
         "then re-run this script.\nSee README Step 2.")
  }
}
rm(.pca_inputs_for_stale_check, .existing_pca_inputs)

source(file.path("figure-scripts", "helpers-plotting.R"))
source(file.path("figure-scripts", "helpers-trait-labels.R"))

trait_categories   <- fread(file.path("data", "trait_abbrevs_categorised.txt"))
large_plot_palette <- colorRampPalette(brewer.pal(8, "Set1"))(16)

custom_theme_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
  "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
  "#E5C494", "#B3B3B3"
)
all_trait_categories <- sort(unique(na.omit(as.character(trait_categories$phenotype_group))))
all_trait_categories <- all_trait_categories[nzchar(all_trait_categories)]
if (length(all_trait_categories) > length(custom_theme_colors))
  stop("Not enough colors for phenotype categories: need ",
       length(all_trait_categories), ", have ", length(custom_theme_colors))
category_colors_all <- setNames(custom_theme_colors[seq_along(all_trait_categories)], all_trait_categories)
message("Figure 1 category palette: ", length(all_trait_categories), " phenotype categories.")

fixed_pca_labels <- list(
  pc12 = c(
    "Anxious", "BasMetR", "BirthWt", "BMI_03", "Body10Y",
    "FFMTrunk", "FVitC_05", "HH_Inc", "HighStrg", "HipCirc", "HlthRating",
    "ImpedBody", "ImpedLeg", "Insomnia", "NervFeel",
    "Neurotic_30", "Obesity1", "SitHgt", "SmokeCur_16",
    "StandHgt", "UniverD", "VascHDxN", "Weight"
  ),
  pc23 = c(
    "DiaBP_A59", "EduYr", "FamSatis",
    "FFM_Leg", "FFMTrunk",
    "FVitC_05", "GripStrR", "Guilt", "HH_Inc", "HighStrg", "HipCirc",
    "HlthRating", "ImpedArm", "ImpedLeg",
    "MoodSwg", "NervFeel", "Neurotic_30", "PExpFl", "SitHgt", "StandHgt",
    "SysBP_A", "UniverD", "VascHDxN", "Weight", "WorryEmb"
  ),
  pc34 = c(
    "Anxious", "DiaBP_A59", "FFMTrunk",
    "FI3_Word", "FluIntSc", "FVitC_05",
    "Guilt", "NervFeel", "Neurotic_30",
    "SitHgt", "StandHgt", "SysBP_A", "Townsend", "ProfQual",
    "RiskTak", "SexPartn", "SmokeCur_25",
    "Speeding", "UniverD", "WorryEmb"
  )
)

cluster_levels_named <- c(
  "0"  = "Unclustered",   "1"  = "Illness",
  "2"  = "Cognitive",     "3"  = "Reproductive & Growth",
  "4"  = "Smoking & Risk","5"  = "Cancer & Sepsis",
  "6"  = "Mood & Psychotic", "7" = "Anxiety & Affective",
  "8"  = "Anthropometric & Strength", "9" = "Obesity & Inflammation",
  "10" = "Hepatic Function", "11" = "Vascular & Ocular",
  "12" = "Immune"
)

# ---- Load and filter effect matrix ----
effects_file <- "results/all-trait-50k-mean-effects.tsv"
pvals_file   <- "results/all-trait-50k-mean-pvals.tsv"
stopifnot(file.exists(effects_file), file.exists(pvals_file))
all_betas <- fread(effects_file)
all_pvals <- fread(pvals_file)

betas_filtered <- filter_effect_matrix(all_betas, all_pvals, p_thresh = NULL)
selected_traits_pca <- colnames(betas_filtered)[-(1:3)]
cols_to_keep  <- c("chrom", "start", "end", selected_traits_pca)
betas_subset  <- betas_filtered[, ..cols_to_keep]

# Exclude inversion block (17q21), MHC (6p), and 12q24 pleiotropic locus
betas_subset <- betas_subset[
  !(
    (chrom == "17" & start < 45000000  & end > 43500000) |
    (chrom == "6"  & start < 34000000  & end > 25000000) |
    (chrom == "12" & start < 113400000 & end > 111200000)
  )
]

# ---- PCA ----
max_pca_dim <- 5
pc_cols     <- paste0("PC", 1:max_pca_dim)
impute_number <- 5

trait_names <- setdiff(names(betas_subset), c("chrom", "start", "end"))
na_counts   <- rowSums(is.na(betas_subset[, ..trait_names]))
betas_subset <- betas_subset[na_counts <= impute_number]
for (col in trait_names) {
  med <- median(betas_subset[[col]], na.rm = TRUE)
  betas_subset[[col]][is.na(betas_subset[[col]])] <- med
}

win_ids <- paste(
  betas_subset$chrom,
  formatC(as.integer(betas_subset$start), format = "f", digits = 0, big.mark = ""),
  formatC(as.integer(betas_subset$end),   format = "f", digits = 0, big.mark = ""),
  sep = "_"
)
imputed_traits <- as.data.frame(betas_subset[, ..trait_names])
rownames(imputed_traits) <- win_ids
traits_t <- t(imputed_traits)

set.seed(1)
colnames(traits_t) <- rownames(imputed_traits)
traits_t_scaled <- t(scale(t(traits_t), center = TRUE, scale = TRUE))

pca_traits_scaled <- prcomp_irlba(traits_t_scaled, n = 40, center = FALSE, scale. = FALSE)
rownames(pca_traits_scaled$rotation) <- colnames(traits_t_scaled)
rownames(pca_traits_scaled$x)        <- rownames(traits_t_scaled)

pca_loadings <- as.data.frame(pca_traits_scaled$rotation)
pca_loadings$region_label <- rownames(pca_loadings)
setDT(pca_loadings)
pca_loadings[, c("chrom", "start", "end") := tstrsplit(region_label, "_", fixed = TRUE)]
pca_loadings[, chrom := as.character(chrom)]
pca_loadings[, start := as.integer(start)]
pca_loadings[, end   := as.integer(end)]
fwrite(pca_loadings, "results/pca_loadings_50k.tsv", sep = "\t")
message("Saved: results/pca_loadings_50k.tsv")
message("NOTE: pca_loadings_50k.tsv has been updated. ",
        "Re-run compute-trait-pc-category-labelperm.R before running run-figure1.R again, ",
        "otherwise panel 1d will use a stale category heatmap table.")

pca_scores <- as.data.frame(pca_traits_scaled$x)
pca_scores$Trait <- rownames(pca_scores)
pc_cols <- intersect(pc_cols, colnames(pca_scores))

trait_cat_subset <- trait_categories %>% filter(trait_abbrev %in% pca_scores$Trait)
pca_plot_data <- pca_scores[, pc_cols, drop = FALSE]
pca_plot_data$trait <- pca_scores$Trait
pca_plot_data <- left_join(pca_plot_data, trait_cat_subset, by = c("trait" = "trait_abbrev"))
pca_plot_data$phenotype_group <- factor(pca_plot_data$phenotype_group, levels = all_trait_categories)
setDT(pca_plot_data)
pca_plot_data[, total_abs := rowSums(abs(.SD)), .SDcols = pc_cols]

# Full eigendecomposition for variance-explained axis labels
G_full   <- tcrossprod(traits_t_scaled) / max(1, nrow(traits_t_scaled) - 1L)
eig_full <- eigen(G_full, symmetric = TRUE, only.values = TRUE)$values
eig_full <- pmax(eig_full, 0)
pve_pct  <- eig_full / sum(eig_full) * 100
pc1_label <- paste0("PC1 (", sprintf("%.2f", pve_pct[1]), "%)")
pc2_label <- paste0("PC2 (", sprintf("%.2f", pve_pct[2]), "%)")
pc3_label <- paste0("PC3 (", sprintf("%.2f", pve_pct[3]), "%)")

# ---- FigureS6: trait loading heatmap (top 60 traits, PC1–5) ----
color_palette <- rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100))
max_load  <- pca_plot_data[, max(abs(unlist(.SD))), .SDcols = pc_cols]
corr.lims <- c(-max_load, max_load)
top_traits <- head(pca_plot_data[order(-total_abs)], 60)$trait
sub        <- subset(pca_plot_data, trait %in% top_traits)
sub        <- sub[match(top_traits, sub$trait), ]
mat        <- as.matrix(sub[, ..pc_cols])
rownames(mat) <- sub$trait

col_fun    <- colorRamp2(seq(corr.lims[1], corr.lims[2], length.out = 100), color_palette)
row_hclust <- hclust(dist(mat), method = "ward.D2")
row_annot_data <- sub[, .(trait, phenotype_group)]
row_annot_data$phenotype_group <- droplevels(row_annot_data$phenotype_group)
group_levels  <- levels(row_annot_data$phenotype_group)
group_levels  <- group_levels[group_levels %in% as.character(unique(row_annot_data$phenotype_group))]
group_colours <- category_colors_all[group_levels]
row_ha <- rowAnnotation(
  Group = row_annot_data$phenotype_group,
  col = list(Group = group_colours),
  annotation_legend_param = list(title = "Trait group")
)
pc_trait_heat <- Heatmap(
  mat, name = "Loading", col = col_fun,
  cluster_rows = as.dendrogram(row_hclust), cluster_columns = FALSE,
  right_annotation = row_ha, row_dend_side = "left", border = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.5),
  row_names_gp    = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_title    = "",
  heatmap_legend_param = list(
    title = "PC loading", at = c(corr.lims[1], 0, corr.lims[2]),
    labels = c("Negative", "Neutral", "Positive"), legend_height = unit(3, "cm")
  )
)
pc_trait_heat_plot <- wrap_elements(full = grid::grid.grabExpr(draw(pc_trait_heat)))
ggsave("figures/figureS6_trait_loading_heatmap_top60_pc1_pc5.pdf",
       pc_trait_heat_plot, height = 10.5, width = 5.95)
message("FigureS6 complete.")

# ---- FigureS5: UMAP ----
tuned_n_neighbors <- 10
tuned_min_dist    <- 0.075
tuned_eps         <- 0.45
tuned_minPts      <- 11

umap_traits_scaled <- uwot::umap(
  traits_t_scaled, n_neighbors = tuned_n_neighbors, min_dist = tuned_min_dist,
  metric = "euclidean", n_components = 2, seed = 42
)
umap_df <- data.frame(
  trait_abbrev = rownames(traits_t_scaled),
  UMAP1 = umap_traits_scaled[, 1],
  UMAP2 = umap_traits_scaled[, 2]
)
umap_df <- left_join(umap_df, trait_cat_subset, by = "trait_abbrev")
db <- dbscan::dbscan(umap_df[, c("UMAP1", "UMAP2")], eps = tuned_eps, minPts = tuned_minPts)
umap_df$cluster <- as.character(db$cluster)
message("DBSCAN clustering complete.")

setDT(umap_df)
umap_df[, phenotype_group := factor(phenotype_group, levels = all_trait_categories)]
umap_df[, label := fifelse(
  trait_abbrev %in% selected_traits_umap_fig3a[seq(1, length(selected_traits_umap_fig3a), by = 2)],
  trait_abbrev, NA_character_
)]
cluster_ids <- sort(unique(umap_df$cluster))
cluster_colours <- c("#AAAAAA", large_plot_palette)
stopifnot(length(cluster_colours) >= length(cluster_ids))
umap_df[, cluster_f := factor(cluster, levels = cluster_ids)]

umap_post_clust_plot <- ggplot(umap_df, aes(UMAP1, UMAP2)) +
  geom_point(aes(fill = phenotype_group),
             shape = 21, colour = "black", stroke = 0.25, size = 2.5, alpha = 0.9) +
  geom_text_repel(
    data = umap_df[trait_abbrev %in% selected_traits_umap_fig3a[seq(1, length(selected_traits_umap_fig3a), by = 2)]],
    aes(label = label), size = 3, max.overlaps = Inf, nudge_x = 0.01, nudge_y = 0.01,
    box.padding = 0.6, point.padding = 0.4, segment.size = 0.4, segment.alpha = 0.8,
    show.legend = FALSE, min.segment.length = 0
  ) +
  scale_fill_manual(name = "Category", values = category_colors_all,
                    breaks = all_trait_categories, labels = all_trait_categories, drop = FALSE) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right", panel.grid = element_line(linewidth = 0.3, colour = "grey85")) +
  labs(title = "All clusters", x = "UMAP 1", y = "UMAP 2")

ggsave("figures/figureS5_umap.pdf", umap_post_clust_plot, width = 8.5, height = 6.5)
message("FigureS5 complete: figures/figureS5_umap.pdf")

write.table(umap_df, file = "results/umap-coordinates.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Figure 1 assembly ----

# Scree panel (pre-computed multiscale data)
cache2 <- readRDS(fig2_cache)
list2env(cache2, envir = environment())

# Panel a: PC1 vs PC2
pca_plot_data$label_b <- ifelse(pca_plot_data$trait %in% fixed_pca_labels$pc12, pca_plot_data$trait, NA_character_)
p_a <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, fill = phenotype_group)) +
  geom_point(shape = 21, size = 2, colour = "black", stroke = 0.4, alpha = 0.9) +
  geom_text_repel(data = subset(pca_plot_data, !is.na(label_b)), aes(label = label_b),
                  size = 3, box.padding = 0.45, point.padding = 0.2,
                  segment.size = 0.35, segment.alpha = 0.65,
                  max.overlaps = Inf, force = 2.2, force_pull = 0.08, min.segment.length = 0) +
  scale_fill_manual(values = category_colors_all, drop = FALSE) +
  theme_minimal() + labs(x = pc1_label, y = pc2_label, fill = "Category", tag = "a") +
  theme(legend.position = "none")

# Panel b: PC2 vs PC3
pca_plot_data$label_c <- ifelse(pca_plot_data$trait %in% fixed_pca_labels$pc23, pca_plot_data$trait, NA_character_)
p_b <- ggplot(pca_plot_data, aes(x = PC2, y = PC3, fill = phenotype_group)) +
  geom_point(shape = 21, size = 2, colour = "black", stroke = 0.4, alpha = 0.9) +
  geom_text_repel(data = subset(pca_plot_data, !is.na(label_c)), aes(label = label_c),
                  size = 3, box.padding = 0.45, point.padding = 0.2,
                  segment.size = 0.35, segment.alpha = 0.65,
                  max.overlaps = Inf, force = 2.2, force_pull = 0.08, min.segment.length = 0) +
  scale_fill_manual(values = category_colors_all, drop = FALSE) +
  theme_minimal() + labs(x = pc2_label, y = pc3_label, fill = "Category", tag = "b") +
  theme(legend.position = "none")

# Panel c: scree
base_theme2 <- theme_minimal(base_size = 10) +
  theme(panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey92", linewidth = 0.2),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3))

p_scree <- ggplot(scree_dt, aes(x = pc_index, y = var_explained * 100, color = scheme)) +
  geom_line(aes(group = scheme), linewidth = 0.55) +
  geom_point(size = 1.0) +
  scale_color_manual(values = scheme_colors, labels = label_map) +
  labs(x = "Principal Component", y = "Variance Explained (%)", color = NULL, tag = "c") +
  base_theme2 +
  theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", colour = "grey70", linewidth = 0.25))

# Panel d: category × PC heatmap
trait_pc_heat_dt <- fread(heat_tsv)
trait_pc_heat_dt[, best_p    := min(p_emp_high,     na.rm = TRUE), by = phenotype_group]
trait_pc_heat_dt[, max_z     := max(z_score,         na.rm = TRUE), by = phenotype_group]
trait_pc_heat_dt[, max_delta := max(delta_from_null, na.rm = TRUE), by = phenotype_group]
group_order <- unique(trait_pc_heat_dt[, .(phenotype_group, best_p, max_z, max_delta)]
                      )[order(best_p, -max_z, -max_delta), phenotype_group]
trait_pc_heat_dt[, phenotype_group := factor(phenotype_group, levels = group_order)]
trait_pc_heat_dt[, pc := factor(pc, levels = paste0("PC", 1:14))]
group_n_dt <- unique(pca_plot_data[, .(trait, phenotype_group)])
group_n_dt <- group_n_dt[!is.na(phenotype_group), .(n_traits = uniqueN(trait)), by = phenotype_group]
group_n_dt[, heatmap_label := sprintf("%s (N=%d)", phenotype_group, n_traits)]
heatmap_labels <- setNames(group_n_dt$heatmap_label, group_n_dt$phenotype_group)
rho_cap <- 0.2
trait_pc_heat_dt[, label_col  := fifelse(observed_median_abs_rho > 0.12, "white", "black")]
trait_pc_heat_dt[, pc_index   := as.integer(sub("^PC", "", as.character(pc)))]
strip_dt <- unique(trait_pc_heat_dt[, .(phenotype_group)])
strip_dt[, strip_col := unname(category_colors_all[as.character(phenotype_group)])]
ag_sunset_cols <- rev(colorRampPalette(as.character(paletteer::paletteer_d("rcartocolor::ag_Sunset")))(256))

p_d <- ggplot(trait_pc_heat_dt, aes(x = pc_index, y = phenotype_group, fill = observed_median_abs_rho)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_tile(data = strip_dt, aes(x = 0.35, y = phenotype_group), inherit.aes = FALSE,
            fill = strip_dt$strip_col, width = 0.18, height = 0.98, linewidth = 0) +
  geom_text(aes(label = sig, color = label_col), size = 3.2, na.rm = TRUE) +
  scale_fill_gradientn(colours = ag_sunset_cols, limits = c(0, rho_cap),
                       oob = scales::squish, na.value = "grey90", name = "Median |rho|") +
  scale_x_continuous(breaks = 1:14, labels = paste0("PC", 1:14),
                     expand = expansion(add = c(0.45, 0.2))) +
  scale_y_discrete(labels = heatmap_labels) +
  scale_color_identity() + coord_cartesian(clip = "off") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = rel(1.2)),
        axis.title  = element_blank(),
        plot.margin = margin(5.5, 5.5, 5.5, 12)) +
  labs(tag = "d")

# Legend panels
legend_levels   <- levels(pca_plot_data$phenotype_group)
legend_levels   <- legend_levels[!is.na(legend_levels)]
half_n          <- ceiling(length(legend_levels) / 2)
legend_groups_1 <- legend_levels[seq_len(half_n)]
legend_groups_2 <- legend_levels[(half_n + 1):length(legend_levels)]
legend_groups_2 <- legend_groups_2[!is.na(legend_groups_2)]

make_subset_legend <- function(levels_subset, title = NULL) {
  leg_dt <- data.frame(phenotype_group = factor(levels_subset, levels = levels_subset), x = 1, y = 1)
  cowplot::get_legend(
    ggplot(leg_dt, aes(x = x, y = y, fill = phenotype_group)) +
      geom_point(shape = 21, size = 2, colour = "black", stroke = 0.35) +
      scale_fill_manual(values = category_colors_all[levels_subset], drop = FALSE) +
      guides(fill = guide_legend(title = title, ncol = 1, byrow = TRUE)) +
      theme_void() +
      theme(legend.position = "right", legend.title = element_text(size = 11),
            legend.text = element_text(size = 10))
  )
}
legend_panel_1 <- cowplot::ggdraw(make_subset_legend(legend_groups_1, title = "Trait category"))
legend_panel_2 <- cowplot::ggdraw(make_subset_legend(legend_groups_2, title = NULL))
legend_panel_1 <- cowplot::plot_grid(legend_panel_1, NULL, ncol = 1, rel_heights = c(1, 0.55))
legend_panel_2 <- cowplot::plot_grid(legend_panel_2, NULL, ncol = 1, rel_heights = c(1, 0.55))

top_row    <- cowplot::plot_grid(p_a, p_b, p_scree, ncol = 3, rel_widths = c(1, 1.03, 1.02))
bottom_row <- cowplot::plot_grid(legend_panel_1, legend_panel_2, p_d, ncol = 3, rel_widths = c(0.46, 0.46, 1.58))
fig1       <- cowplot::plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(0.92, 0.84))
fig1       <- cowplot::ggdraw() + cowplot::draw_plot(fig1, x = 0.012, y = 0, width = 0.988, height = 1)

ggsave("figures/figure1.pdf", fig1, width = 11.9, height = 8.8)
message("Figure 1 complete: figures/figure1.pdf")
