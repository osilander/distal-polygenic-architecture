#!/usr/bin/env Rscript
# Produces: figures/figureS8_pca_abs50k_from_windows.pdf
#           results/pca_abs50k_from_windows/{trait_scores,window_loadings,singular_values,
#                                            scree_fullspectrum,window_coordinates,
#                                            heatmap_matrix_top60}.tsv
# Pre-computed inputs: data/windows_filtered/*_w50000.summary.tsv, data/trait_abbrevs_categorised.txt
# Run: Rscript figure-scripts/run-figureS8.R

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

windows_dir <- file.path("data", "windows_filtered")
trait_file  <- file.path("data", "trait_abbrevs_categorised.txt")
out_fig     <- file.path("figures", "figureS8_pca_abs50k_from_windows.pdf")
out_dir     <- file.path("results", "pca_abs50k_from_windows")
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Category colors ----
message("FigureS8: loading trait map and category colors...")
trait_map <- fread(trait_file)[, .(trait_abbrev, phenotype_group)]
allowed_traits <- unique(trait_map$trait_abbrev)
trait_categories <- fread(trait_file)
custom_theme_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
  "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
  "#E5C494", "#B3B3B3"
)
all_trait_categories <- sort(unique(na.omit(as.character(trait_categories$phenotype_group))))
all_trait_categories <- all_trait_categories[nzchar(all_trait_categories)]
category_colors      <- setNames(custom_theme_colors[seq_along(all_trait_categories)], all_trait_categories)

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

# ---- Load window files ----
files <- list.files(windows_dir,
                    pattern = "(\\.w50000\\.summary\\.tsv|_w50000\\.summary\\.tsv)$",
                    full.names = TRUE)
if (length(files) == 0L) stop("No 50kb window files found in ", windows_dir)

parse_trait <- function(f) {
  b <- basename(f)
  b <- sub("\\.w50000\\.summary\\.tsv$", "", b)
  b <- sub("_w50000\\.summary\\.tsv$", "", b)
  b
}

file_dt <- data.table(file = files, trait_abbrev = vapply(files, parse_trait, character(1)))
file_dt <- file_dt[trait_abbrev %in% allowed_traits]
if (nrow(file_dt) == 0L) stop("No 50kb files matched allowed trait list.")
file_dt[, bytes := file.size(file)]
setorder(file_dt, trait_abbrev, -bytes, file)
file_dt <- file_dt[, .SD[1], by = trait_abbrev]
message("FigureS8: using ", nrow(file_dt), " traits.")

read_one <- function(f, tr) {
  x <- fread(f, select = c("chrom", "start", "end", "median_es_abs"),
             na.strings = c("NA", "NaN", "nan", ".", ""))
  x[, `:=`(chrom = as.character(chrom), start = as.integer(start),
            end   = as.integer(end), value = as.numeric(median_es_abs), trait_abbrev = tr)]
  x <- x[is.finite(start) & is.finite(end)]
  x[, window_id := paste(chrom, start, end, sep = "_")]
  x[, .(window_id, chrom, start, end, trait_abbrev, value)]
}

long_dt <- rbindlist(lapply(seq_len(nrow(file_dt)), function(i) read_one(file_dt$file[i], file_dt$trait_abbrev[i])),
                     use.names = TRUE, fill = TRUE)
wide <- dcast(long_dt, window_id + chrom + start + end ~ trait_abbrev, value.var = "value")

wide <- wide[!(
  (chrom == "17" & start < 45000000  & end > 43500000) |
  (chrom == "6"  & start < 34000000  & end > 25000000) |
  (chrom == "12" & start < 113400000 & end > 111200000)
)]

trait_cols    <- setdiff(names(wide), c("window_id", "chrom", "start", "end"))
impute_number <- 5L
na_counts     <- rowSums(is.na(wide[, ..trait_cols]))
wide          <- wide[na_counts <= impute_number]
for (cc in trait_cols) {
  med <- median(wide[[cc]], na.rm = TRUE)
  wide[[cc]][is.na(wide[[cc]])] <- med
}

X_win_trait  <- as.matrix(wide[, ..trait_cols]); rownames(X_win_trait) <- wide$window_id
X_trait_win  <- t(X_win_trait)
X_trait_win  <- t(scale(t(X_trait_win), center = TRUE, scale = TRUE))
X_trait_win[!is.finite(X_trait_win)] <- 0
message("FigureS8: matrix dims = ", nrow(X_trait_win), " traits x ", ncol(X_trait_win), " windows")

# ---- PCA ----
pca <- prcomp_irlba(X_trait_win, n = 80, center = FALSE, scale. = FALSE)
rownames(pca$x)        <- rownames(X_trait_win)
rownames(pca$rotation) <- colnames(X_trait_win)

G_full   <- tcrossprod(X_trait_win) / max(1, nrow(X_trait_win) - 1L)
eig_full <- eigen(G_full, symmetric = TRUE, only.values = TRUE)$values
eig_full <- pmax(eig_full, 0); pve <- eig_full / sum(eig_full); pve_pct <- pve * 100

scores <- as.data.table(pca$x); scores[, trait_abbrev := rownames(pca$x)]
scores <- merge(scores, trait_map, by = "trait_abbrev", all.x = TRUE)
scores[, phenotype_group := factor(phenotype_group, levels = names(category_colors))]

loadings <- as.data.table(pca$rotation); loadings[, window_id := rownames(pca$rotation)]
svals    <- data.table(component = seq_along(pca$sdev), singular_value = pca$sdev)
scree    <- data.table(PC = seq_along(pve), var_explained = pve, var_explained_pct = pve * 100, cumulative = cumsum(pve))

fwrite(scores,   file.path(out_dir, "trait_scores.tsv"),     sep = "\t")
fwrite(loadings, file.path(out_dir, "window_loadings.tsv"),  sep = "\t")
fwrite(svals,    file.path(out_dir, "singular_values.tsv"),  sep = "\t")
fwrite(scree,    file.path(out_dir, "scree_fullspectrum.tsv"), sep = "\t")
fwrite(wide[, c("window_id", "chrom", "start", "end"), with = FALSE],
       file.path(out_dir, "window_coordinates.tsv"), sep = "\t")

pc1_label <- paste0("PC1 (", sprintf("%.3f", pve_pct[1]), "%)")
pc2_label <- paste0("PC2 (", sprintf("%.3f", pve_pct[2]), "%)")
pc3_label <- paste0("PC3 (", sprintf("%.3f", pve_pct[3]), "%)")
pc4_label <- paste0("PC4 (", sprintf("%.3f", pve_pct[4]), "%)")

pc_cols5 <- paste0("PC", 1:5)
plot_dt  <- copy(scores)
plot_dt[, label_b := ifelse(trait_abbrev %in% fixed_pca_labels$pc12, trait_abbrev, NA_character_)]
plot_dt[, label_c := ifelse(trait_abbrev %in% fixed_pca_labels$pc23, trait_abbrev, NA_character_)]
plot_dt[, label_d := ifelse(trait_abbrev %in% fixed_pca_labels$pc34, trait_abbrev, NA_character_)]

make_scatter <- function(xcol, ycol, label_col, xlab, ylab) {
  ggplot(plot_dt, aes(x = get(xcol), y = get(ycol), fill = phenotype_group)) +
    geom_point(shape = 21, size = 2, colour = "black", stroke = 0.4, alpha = 0.9) +
    geom_text_repel(data = plot_dt[!is.na(get(label_col))], aes(label = get(label_col)),
                    size = 3, box.padding = 0.3, segment.size = 0.3, segment.alpha = 0.5,
                    max.overlaps = Inf, force = 1.5, force_pull = 0.1, min.segment.length = 0) +
    scale_fill_manual(values = category_colors, drop = FALSE) +
    theme_minimal() + labs(x = xlab, y = ylab, fill = "Category") + theme(legend.position = "none")
}

pc1_pc2_plot <- make_scatter("PC1", "PC2", "label_b", pc1_label, pc2_label)
pc2_pc3_plot <- make_scatter("PC2", "PC3", "label_c", pc2_label, pc3_label)
pc3_pc4_plot <- make_scatter("PC3", "PC4", "label_d", pc3_label, pc4_label)

# ---- Heatmap ----
pc_cols_h <- intersect(pc_cols5, names(plot_dt))
plot_dt[, total_abs := rowSums(abs(.SD)), .SDcols = pc_cols_h]
top_traits <- head(plot_dt[order(-total_abs), trait_abbrev], 60)
sub  <- plot_dt[trait_abbrev %in% top_traits][match(top_traits, trait_abbrev)]
mat  <- as.matrix(sub[, ..pc_cols_h]); rownames(mat) <- sub$trait_abbrev
fwrite(as.data.table(mat, keep.rownames = "trait_abbrev"), file.path(out_dir, "heatmap_matrix_top60.tsv"), sep = "\t")

color_palette <- rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100))
max_load  <- max(abs(mat)); corr.lims <- c(-max_load, max_load)
col_fun   <- colorRamp2(seq(corr.lims[1], corr.lims[2], length.out = 100), color_palette)
row_hclust <- hclust(dist(mat), method = "ward.D2")
row_annot_data <- sub[, .(trait_abbrev, phenotype_group)]
row_annot_data$phenotype_group <- droplevels(row_annot_data$phenotype_group)
group_levels  <- levels(row_annot_data$phenotype_group)
group_levels  <- group_levels[group_levels %in% as.character(unique(row_annot_data$phenotype_group))]
group_colours <- category_colors[group_levels]
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

first_col <- pc1_pc2_plot / pc2_pc3_plot / pc3_pc4_plot
legend_seed <- data.frame(
  phenotype_group = factor(levels(plot_dt$phenotype_group), levels = levels(plot_dt$phenotype_group)),
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
message("FigureS8 complete: ", out_fig)
message("Results written to: ", out_dir)
