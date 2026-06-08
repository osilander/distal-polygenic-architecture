#!/usr/bin/env Rscript
# Produces: figures/figureS9_pca_absolute_50k.pdf,
#           results/figureS_pca_absolute_50k_{trait_scores,scree,window_loadings_top200,
#                                             heatmap_matrix_top60}.tsv
# Pre-computed inputs: results/all-trait-50k-mean-{effects,pvals}.tsv,
#                      data/trait_abbrevs_categorised.txt
# Run: Rscript figure-scripts/run-figureS9.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(irlba)
  library(patchwork)
  library(RColorBrewer)
  library(circlize)
  library(ComplexHeatmap)
  library(cowplot)
  library(grid)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

out_fig    <- file.path("figures", "figureS9_pca_absolute_50k.pdf")
out_scores <- file.path("results", "figureS_pca_absolute_50k_trait_scores.tsv")
out_scree  <- file.path("results", "figureS_pca_absolute_50k_scree.tsv")
out_win    <- file.path("results", "figureS_pca_absolute_50k_window_loadings_top200.tsv")
out_heat   <- file.path("results", "figureS_pca_absolute_50k_heatmap_matrix_top60.tsv")
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)

# ---- Rebuild traits_t_scaled ----
message("FigureS9: rebuilding scaled trait matrix (~2-3 min)...")
source(file.path("figure-scripts", "helpers-50k-matrix.R"))
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
category_colors      <- setNames(custom_theme_colors[seq_along(all_trait_categories)], all_trait_categories)
used_categories      <- names(category_colors)

fixed_pca_labels <- list(
  pc12 = c("Anxious", "BasMetR", "BirthWt", "BMI_03", "Body10Y",
           "FFMTrunk", "FVitC_05", "HH_Inc", "HighStrg", "HipCirc", "HlthRating",
           "ImpedBody", "ImpedLeg", "Insomnia", "NervFeel",
           "Neurotic_30", "Obesity1", "SitHgt", "SmokeCur_16",
           "StandHgt", "UniverD", "VascHDxN", "Weight"),
  pc23 = c("DiaBP_A59", "EduYr", "FamSatis", "FFM_Leg", "FFMTrunk",
           "FVitC_05", "GripStrR", "Guilt", "HH_Inc", "HighStrg", "HipCirc",
           "HlthRating", "ImpedArm", "ImpedLeg",
           "MoodSwg", "NervFeel", "Neurotic_30", "PExpFl", "SitHgt", "StandHgt",
           "SysBP_A", "UniverD", "VascHDxN", "Weight", "WorryEmb"),
  pc34 = c("Anxious", "DiaBP_A59", "FFMTrunk",
           "FI3_Word", "FluIntSc", "FVitC_05",
           "Guilt", "NervFeel", "Neurotic_30",
           "SitHgt", "StandHgt", "SysBP_A", "Townsend", "ProfQual",
           "RiskTak", "SexPartn", "SmokeCur_25",
           "Speeding", "UniverD", "WorryEmb")
)

# ---- PCA on abs(traits_t_scaled) ----
X <- as.matrix(traits_t_scaled)
message("FigureS9: running PCA on abs(X) (", nrow(X), " traits × ", ncol(X), " windows)...")
X_abs  <- abs(X)
ncomp  <- min(80L, min(dim(X_abs)) - 1L)
pca_abs <- prcomp_irlba(X_abs, n = ncomp, center = FALSE, scale. = FALSE)
rownames(pca_abs$x)        <- rownames(X_abs)
rownames(pca_abs$rotation) <- colnames(X_abs)

max_pca_dim    <- 5L
export_max_pc  <- 80L
pc_cols        <- paste0("PC", seq_len(max_pca_dim))
pc_cols        <- intersect(pc_cols, colnames(pca_abs$x))
pc_export_cols <- intersect(paste0("PC", seq_len(export_max_pc)), colnames(pca_abs$x))
if (length(pc_cols) < 4L) stop("Need at least PC1-PC4; available: ", paste(colnames(pca_abs$x), collapse = ", "))

G_full   <- tcrossprod(X_abs) / max(1, nrow(X_abs) - 1L)
eig_full <- eigen(G_full, symmetric = TRUE, only.values = TRUE)$values
eig_full <- pmax(eig_full, 0); pve <- eig_full / sum(eig_full); pve_pct <- pve * 100
pc1_label <- paste0("PC1 (", sprintf("%.3f", pve_pct[1]), "%)")
pc2_label <- paste0("PC2 (", sprintf("%.3f", pve_pct[2]), "%)")
pc3_label <- paste0("PC3 (", sprintf("%.3f", pve_pct[3]), "%)")
pc4_label <- paste0("PC4 (", sprintf("%.3f", pve_pct[4]), "%)")

pca_scores <- as.data.frame(pca_abs$x); pca_scores$Trait <- rownames(pca_scores)
trait_cat_subset <- as.data.table(trait_categories)[trait_abbrev %in% pca_scores$Trait, .(trait_abbrev, phenotype_group)]
pca_plot_data <- pca_scores[, pc_cols, drop = FALSE]
pca_plot_data$trait <- pca_scores$Trait
pca_plot_data <- left_join(pca_plot_data, trait_cat_subset, by = c("trait" = "trait_abbrev"))
pca_plot_data$phenotype_group <- factor(pca_plot_data$phenotype_group, levels = used_categories)
setDT(pca_plot_data)

pc12_labels <- if (!is.null(fixed_pca_labels$pc12)) fixed_pca_labels$pc12 else character()
pc23_labels <- if (!is.null(fixed_pca_labels$pc23)) fixed_pca_labels$pc23 else character()
pc34_labels <- if (!is.null(fixed_pca_labels$pc34)) fixed_pca_labels$pc34 else character()

pca_plot_data$label_b <- ifelse(pca_plot_data$trait %in% pc12_labels, pca_plot_data$trait, NA_character_)
pca_plot_data$label_c <- ifelse(pca_plot_data$trait %in% pc23_labels, pca_plot_data$trait, NA_character_)
pca_plot_data$label_d <- ifelse(pca_plot_data$trait %in% pc34_labels, pca_plot_data$trait, NA_character_)

make_scatter <- function(xcol, ycol, label_col, xlab, ylab) {
  ggplot(pca_plot_data, aes(x = get(xcol), y = get(ycol), fill = phenotype_group)) +
    geom_point(shape = 21, size = 2, colour = "black", stroke = 0.4, alpha = 0.9) +
    geom_text_repel(data = subset(pca_plot_data, !is.na(get(label_col))), aes(label = get(label_col)),
                    size = 3, box.padding = 0.3, segment.size = 0.3, segment.alpha = 0.5,
                    max.overlaps = Inf, force = 1.5, force_pull = 0.1, min.segment.length = 0) +
    scale_fill_manual(values = category_colors, drop = FALSE) +
    theme_minimal() + labs(x = xlab, y = ylab, fill = "Category") + theme(legend.position = "none")
}

pc1_pc2_plot <- make_scatter("PC1", "PC2", "label_b", pc1_label, pc2_label)
pc2_pc3_plot <- make_scatter("PC2", "PC3", "label_c", pc2_label, pc3_label)
pc3_pc4_plot <- make_scatter("PC3", "PC4", "label_d", pc3_label, pc4_label)
message("FigureS9: scatter panels done")

# ---- Heatmap ----
color_palette <- rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100))
max_load      <- pca_plot_data[, max(abs(unlist(.SD))), .SDcols = pc_cols]
corr.lims     <- c(-max_load, max_load)
pca_plot_data[, total_abs := rowSums(abs(.SD)), .SDcols = pc_cols]
top_traits <- head(pca_plot_data[order(-pca_plot_data$total_abs), "trait"], 60)$trait
sub   <- subset(pca_plot_data, trait %in% top_traits); sub <- sub[match(top_traits, sub$trait), ]
mat   <- as.matrix(sub[, ..pc_cols]); rownames(mat) <- sub$trait
fwrite(as.data.table(mat, keep.rownames = "trait"), out_heat, sep = "\t")

col_fun       <- colorRamp2(seq(corr.lims[1], corr.lims[2], length.out = 100), color_palette)
row_hclust    <- hclust(dist(mat), method = "ward.D2")
row_annot_data <- sub[, .(trait, phenotype_group)]
row_annot_data$phenotype_group <- droplevels(row_annot_data$phenotype_group)
group_levels   <- levels(row_annot_data$phenotype_group)
group_levels   <- group_levels[group_levels %in% as.character(unique(row_annot_data$phenotype_group))]
group_colours  <- category_colors[group_levels]
row_ha <- rowAnnotation(Group = row_annot_data$phenotype_group, col = list(Group = group_colours),
                        annotation_legend_param = list(title = "Trait group"))
pc_trait_heat <- Heatmap(mat, name = "Loading", col = col_fun,
                          cluster_rows = as.dendrogram(row_hclust), cluster_columns = FALSE,
                          right_annotation = row_ha, row_dend_side = "left", border = TRUE,
                          rect_gp = gpar(col = "white", lwd = 0.5),
                          row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 10, fontface = "bold"),
                          column_title = "",
                          heatmap_legend_param = list(title = "PC loading",
                            at = c(corr.lims[1], 0, corr.lims[2]),
                            labels = c("Negative", "Neutral", "Positive"),
                            legend_height = unit(3, "cm")))
pc_trait_heat_plot <- wrap_elements(full = grid::grid.grabExpr(draw(pc_trait_heat)))
message("FigureS9: heatmap done")

first_col   <- pc1_pc2_plot / pc2_pc3_plot / pc3_pc4_plot
legend_seed <- data.frame(
  phenotype_group = factor(levels(pca_plot_data$phenotype_group), levels = levels(pca_plot_data$phenotype_group)),
  x = 1, y = 1
)
trait_legend <- cowplot::get_legend(
  ggplot(legend_seed, aes(x = x, y = y, fill = phenotype_group)) +
    geom_point(shape = 21, size = 2, colour = "black", stroke = 0.35) +
    scale_fill_manual(values = category_colors, drop = FALSE) +
    guides(fill = guide_legend(title = "Trait category", ncol = 1, byrow = TRUE)) +
    theme_void() + theme(legend.position = "right")
)
fig <- (first_col | wrap_elements(full = trait_legend) | pc_trait_heat_plot) +
  plot_layout(widths = c(1, 0.24, 1)) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 18))

ggsave(out_fig, fig, width = 14, height = 15)

# ---- Save tabular outputs ----
trait_scores <- as.data.table(pca_scores)[, c("Trait", pc_export_cols), with = FALSE]
setnames(trait_scores, "Trait", "trait_abbrev")
trait_scores <- merge(trait_scores, trait_cat_subset, by = "trait_abbrev", all.x = TRUE)
fwrite(trait_scores, out_scores, sep = "\t")
scree_df <- data.frame(PC = paste0("PC", seq_along(pve)), VarianceExplained = pve * 100, PC_num = seq_along(pve))
fwrite(as.data.table(scree_df), out_scree, sep = "\t")
fwrite(as.data.table(pca_abs$rotation)[, .(
  region_label = rownames(pca_abs$rotation), PC1, PC2, PC3, PC4
)][, abs_PC1 := abs(PC1)][order(-abs_PC1)][1:min(200L, .N)], out_win, sep = "\t")

message("FigureS9 complete: ", out_fig)
