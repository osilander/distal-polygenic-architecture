#!/usr/bin/env Rscript
# Produces: figures/figure3.pdf, figures/figureS13_figure3_pc4_manhattans.pdf
# Pre-computed inputs required:
#   results/pca_loadings_50k.tsv, results/all-trait-50k-mean-effects.tsv,
#   results/all-trait-50k-mean-pvals.tsv, data/gencode.v19.genes.protein_coding.rds,
#   results/gwas_removed_distance_cdf_ribbon_50kb_bins.tsv,
#   results/gwas_removed_distance_cdf_ribbon_50kb_bins_pvalues.tsv
# Run: Rscript figure-scripts/run-figure3.R [flank_kb]
#   default flank_kb = 100

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(irlba)
  library(patchwork)
  library(RColorBrewer)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

pcs_main        <- c("PC1", "PC2", "PC3")
pcs_cdf         <- c("PC1", "PC2", "PC3", "PC4")
gws_log10_thresh <- -log10(5e-8)
top_perc_fig3   <- 0.001

args     <- commandArgs(trailingOnly = TRUE)
flank_kb <- if (length(args) >= 1L) as.integer(args[[1]]) else 100L
if (!is.finite(flank_kb) || flank_kb < 0L) stop("flank_kb must be a non-negative integer.")
flank_bp    <- as.integer(flank_kb * 1000L)
flank_label <- paste0(flank_kb, "kb")

out_fig <- if (flank_kb == 100L) {
  file.path("figures", "figure3.pdf")
} else {
  file.path("figures", paste0("figure3_gws_plus", flank_label, ".pdf"))
}
out_fig_pc4_supp <- if (flank_kb == 100L) {
  file.path("figures", "figureS13_figure3_pc4_manhattans.pdf")
} else {
  file.path("figures", paste0("figureS13_figure3_pc4_manhattans_plus", flank_label, ".pdf"))
}
out_removed_tsv <- if (flank_kb == 100L) {
  file.path("results", "figure3_removed_gws_plus_neighbors_windows.tsv")
} else {
  file.path("results", paste0("figure3_removed_gws_plus", flank_label, "_windows.tsv"))
}
out_top_removed_overlay_tsv <- if (flank_kb == 100L) {
  file.path("results", "figure3_top_removed_overlay_windows.tsv")
} else {
  file.path("results", paste0("figure3_top_removed_overlay_plus", flank_label, "_windows.tsv"))
}
cdf_stats_file <- file.path("results", "gwas_removed_distance_cdf_ribbon_50kb_bins.tsv")
cdf_pvals_file <- file.path("results", "gwas_removed_distance_cdf_ribbon_50kb_bins_pvalues.tsv")
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)

# ---- Load pca_loadings ----
pca_loadings_file <- file.path("results", "pca_loadings_50k.tsv")
if (!file.exists(pca_loadings_file)) stop("Missing: ", pca_loadings_file, "\nRun run-figure1.R first.")
pca_loadings <- fread(pca_loadings_file)

source(file.path("figure-scripts", "helpers-plotting.R"))
chrom_plot_palette <- colorRampPalette(brewer.pal(8, "Set1"))(22)

# ---- Rebuild traits_t_scaled ----
message("Figure 3: rebuilding scaled trait matrix (~2-3 min)...")
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
category_colors_all  <- setNames(custom_theme_colors[seq_along(all_trait_categories)], all_trait_categories)

# ---- Fixed PCA labels ----
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

# ---- Gene annotation ----
genes_rds <- file.path(project_root, "data", "gencode.v19.genes.protein_coding.rds")
if (!file.exists(genes_rds)) stop("Missing gene annotation RDS: ", genes_rds)
genes_gr <- readRDS(genes_rds)

# ---- Helpers ----
parse_region_cols <- function(dt) {
  x <- copy(as.data.table(dt))
  if (!all(c("chr", "start", "end") %in% names(x))) {
    if (!"region_label" %in% names(x)) stop("Need region_label to parse chr/start/end.")
    x[, c("chr", "start", "end") := tstrsplit(region_label, "_", fixed = TRUE)]
  }
  x[, `:=`(chr = as.integer(chr), start = as.integer(start), end = as.integer(end))]
  x[, mid := (start + end) / 2]
  x[, window_id := paste(chr, start, end, sep = "_")]
  setorder(x, chr, start, end)
  x
}

make_cum_axis <- function(dt_regions) {
  x <- copy(dt_regions); setorder(x, chr, mid)
  chr_lengths <- x[, .(chr_len = max(mid, na.rm = TRUE)), by = chr]
  setorder(chr_lengths, chr)
  chr_lengths[, chr_start := c(0, head(cumsum(chr_len), -1))]
  chr_lengths
}

collapse_adjacent_top <- function(dt_top) {
  z <- copy(dt_top)
  if (nrow(z) == 0L) return(z)
  setorder(z, chr, start, end)
  z[, prev_end := shift(end), by = chr]
  z[, is_adj_prev := !is.na(prev_end) & (start == prev_end)]
  z[, run_id := cumsum(!is_adj_prev), by = chr]
  keep <- z[, .SD[floor((.N + 1L) / 2L)], by = .(chr, run_id)]
  keep[, c("prev_end", "is_adj_prev", "run_id") := NULL]
  keep[]
}

# ---- Identify and remove GWAS-significant windows ----
message("Figure 3: reading p-value matrix...")
pval_file <- file.path("results", "all-trait-50k-mean-pvals.tsv")
if (!file.exists(pval_file)) stop("Missing p-value matrix: ", pval_file)
pv <- fread(pval_file, na.strings = c("NA", "NaN", "nan", ".", ""))
pv[, `:=`(chrom = as.integer(chrom), start = as.integer(start), end = as.integer(end))]
pv <- pv[is.finite(chrom) & is.finite(start) & is.finite(end)]
pv[, window_id := paste(chrom, start, end, sep = "_")]
setorder(pv, chrom, start, end)

trait_cols <- setdiff(names(pv), c("chrom", "start", "end", "window_id", "chr_pos"))
for (cc in trait_cols) {
  if (!is.numeric(pv[[cc]])) suppressWarnings(set(pv, j = cc, value = as.numeric(pv[[cc]])))
}

message("Figure 3: identifying GWS windows + ±", flank_kb, "kb flank...")
pmat   <- as.matrix(pv[, ..trait_cols]); storage.mode(pmat) <- "numeric"
any_gws <- apply(pmat > gws_log10_thresh, 1, function(v) any(v, na.rm = TRUE))
pv[, gws_any := any_gws]

if (flank_bp == 0L) {
  removed_dt <- pv[gws_any == TRUE, .(window_id, chrom, start, end, gws_any)]
} else {
  gws <- pv[gws_any == TRUE, .(chrom, gws_start = start, gws_end = end)]
  if (nrow(gws) == 0L) {
    removed_dt <- pv[0, .(window_id, chrom, start, end, gws_any)]
  } else {
    gws[, `:=`(flank_start = pmax(0L, gws_start - flank_bp), flank_end = gws_end + flank_bp)]
    cand <- merge(pv[, .(window_id, chrom, start, end, gws_any)],
                  gws[, .(chrom, flank_start, flank_end)], by = "chrom", allow.cartesian = TRUE)
    removed_ids <- unique(cand[start < flank_end & end > flank_start, window_id])
    removed_dt  <- pv[window_id %in% removed_ids, .(window_id, chrom, start, end, gws_any)]
  }
}
fwrite(removed_dt, out_removed_tsv, sep = "\t")
message("Figure 3: removed windows = ", nrow(removed_dt))

# ---- Top-labeled loci in original PCA (for grey overlay) ----
pl <- parse_region_cols(pca_loadings)
top_label_by_pc <- rbindlist(lapply(c(pcs_main, "PC4"), function(pc) {
  z <- copy(pl[, .(window_id, chr, start, end, loading = get(pc))])
  z[, abs_loading := abs(loading)]
  cutoff <- stats::quantile(z$abs_loading, probs = 1 - top_perc_fig3, na.rm = TRUE)
  z  <- z[abs_loading >= cutoff]
  z  <- collapse_adjacent_top(z)
  z[, pc := pc]
  z[, .(pc, window_id, chr, start, end, loading)]
}), use.names = TRUE, fill = TRUE)

overlay_removed <- merge(top_label_by_pc, removed_dt[, .(window_id)], by = "window_id", all = FALSE)
fwrite(overlay_removed, out_top_removed_overlay_tsv, sep = "\t")
message("Figure 3: Fig2-top loci removed = ", nrow(overlay_removed))

# ---- Reduced PCA ----
# Note: compute-gwas-removed-distances.R is the canonical GWAS-depleted pipeline
# and writes reduced_pca_window_loadings.tsv / reduced_pca_trait_scores.tsv.
# This script recomputes the reduced PCA independently because the Manhattan-plot
# and locus-zoom panels need the full prcomp_irlba object in memory (rotation
# matrix + scores), not just the TSV outputs.  Both implementations apply
# identical window removal (same threshold, same flank_kb) and identical scaling,
# so leading components agree.  n = 40 here vs n = 80 in compute-gwas-removed-
# distances.R; both are well above the 14 PCs displayed.
orig_trait_window <- as.matrix(traits_t_scaled)
remove_ids  <- intersect(colnames(orig_trait_window), removed_dt$window_id)
keep_cols   <- setdiff(colnames(orig_trait_window), remove_ids)
if (length(keep_cols) < 100L) stop("Reduced matrix has too few windows after screening.")
reduced_trait_window <- orig_trait_window[, keep_cols, drop = FALSE]
message("Figure 3: running PCA on reduced matrix (", nrow(reduced_trait_window), " traits × ", ncol(reduced_trait_window), " windows)...")
pca_reduced <- prcomp_irlba(reduced_trait_window, n = 40, center = FALSE, scale. = FALSE)
rownames(pca_reduced$rotation) <- colnames(reduced_trait_window)
rownames(pca_reduced$x)        <- rownames(reduced_trait_window)

# ---- Anchor-based sign orientation ----
anchor_map <- c(PC1 = "BMI_03", PC2 = "StandHgt", PC3 = "Neurotic_30", PC4 = "DiaBP_A59")
for (k in seq_along(anchor_map)) {
  trait_nm <- anchor_map[k]; pc_nm <- names(anchor_map)[k]
  if (trait_nm %in% rownames(pca_reduced$x) && pc_nm %in% colnames(pca_reduced$x)) {
    if (pca_reduced$x[trait_nm, pc_nm] < 0) {
      pca_reduced$x[, pc_nm]        <- -pca_reduced$x[, pc_nm]
      pca_reduced$rotation[, pc_nm] <- -pca_reduced$rotation[, pc_nm]
    }
  }
}

reduced_loadings <- as.data.table(pca_reduced$rotation)
reduced_loadings[, region_label := rownames(pca_reduced$rotation)]
setcolorder(reduced_loadings, c("region_label", setdiff(names(reduced_loadings), "region_label")))
reduced_loadings <- parse_region_cols(reduced_loadings)
axis_reduced     <- make_cum_axis(reduced_loadings[, .(region_label, chr, start, end, mid)])

overlay_for_pc <- function(pc_name) {
  ov <- merge(
    overlay_removed[pc == pc_name, .(window_id, chr, start, end, removed_loading = loading)],
    axis_reduced[, .(chr, chr_start)], by = "chr", all.x = TRUE
  )
  ov[, mid := (start + end) / 2]
  ov[, pos_cum := mid + chr_start]
  ov
}

# ---- Manhattan plots ----
message("Figure 3 / Panels A-C: Manhattan plots...")
pc1_res <- plot_pc_loadings(reduced_loadings, pc = "PC1", genes_dt = genes_gr)
pc2_res <- plot_pc_loadings(reduced_loadings, pc = "PC2", genes_dt = genes_gr)
pc3_res <- plot_pc_loadings(reduced_loadings, pc = "PC3", genes_dt = genes_gr)
pc4_orig_res <- plot_pc_loadings(pca_loadings, pc = "PC4", genes_dt = genes_gr)
pc4_red_res  <- plot_pc_loadings(reduced_loadings, pc = "PC4", genes_dt = genes_gr)

ov1 <- overlay_for_pc("PC1"); ov2 <- overlay_for_pc("PC2")
ov3 <- overlay_for_pc("PC3"); ov4 <- overlay_for_pc("PC4")

p_manh1 <- pc1_res$plot + geom_point(data = ov1, aes(x = pos_cum, y = removed_loading), inherit.aes = FALSE, colour = "grey55", size = 0.9, alpha = 0.9)
p_manh2 <- pc2_res$plot + geom_point(data = ov2, aes(x = pos_cum, y = removed_loading), inherit.aes = FALSE, colour = "grey55", size = 0.9, alpha = 0.9)
p_manh3 <- pc3_res$plot + geom_point(data = ov3, aes(x = pos_cum, y = removed_loading), inherit.aes = FALSE, colour = "grey55", size = 0.9, alpha = 0.9)
p_manh4_orig <- pc4_orig_res$plot
p_manh4_red  <- pc4_red_res$plot + geom_point(data = ov4, aes(x = pos_cum, y = removed_loading), inherit.aes = FALSE, colour = "grey55", size = 0.9, alpha = 0.9)

# ---- PC scatter panels ----
message("Figure 3 / Panels D-F: PC scatter plots...")
scores_red <- as.data.table(pca_reduced$x)
scores_red[, trait := rownames(pca_reduced$x)]
scores_red <- merge(scores_red, trait_categories[, .(trait_abbrev, phenotype_group)],
                    by.x = "trait", by.y = "trait_abbrev", all.x = TRUE)
scores_red[is.na(phenotype_group), phenotype_group := "Uncategorized"]
all_levels <- names(category_colors_all)
if ("Uncategorized" %in% scores_red$phenotype_group && !("Uncategorized" %in% all_levels)) {
  category_cols <- c(category_colors_all, "Uncategorized" = "#999999")
  all_levels <- c(all_levels, "Uncategorized")
} else {
  category_cols <- category_colors_all
}
scores_red[, phenotype_group := factor(phenotype_group, levels = all_levels)]

G_red_full   <- tcrossprod(reduced_trait_window) / max(1, nrow(reduced_trait_window) - 1L)
eig_red_full <- eigen(G_red_full, symmetric = TRUE, only.values = TRUE)$values
eig_red_full <- pmax(eig_red_full, 0)
pve_red_pct  <- 100 * eig_red_full / sum(eig_red_full)
pc1_lab <- sprintf("PC1 (%.2f%%)", pve_red_pct[1])
pc2_lab <- sprintf("PC2 (%.2f%%)", pve_red_pct[2])
pc3_lab <- sprintf("PC3 (%.2f%%)", pve_red_pct[3])
pc4_lab <- sprintf("PC4 (%.2f%%)", pve_red_pct[4])

scores_red[, label12 := fifelse(trait %in% fixed_pca_labels$pc12, trait, NA_character_)]
scores_red[, label23 := fifelse(trait %in% fixed_pca_labels$pc23, trait, NA_character_)]
scores_red[, label34 := fifelse(trait %in% fixed_pca_labels$pc34, trait, NA_character_)]

make_scatter <- function(xcol, ycol, labels_col, xlab, ylab) {
  ggplot(scores_red, aes(get(xcol), get(ycol), fill = phenotype_group)) +
    geom_point(shape = 21, size = 1.2, colour = "black", stroke = 0.4, alpha = 0.9) +
    geom_text_repel(data = scores_red[!is.na(get(labels_col))], aes(label = get(labels_col)),
                    size = 2.1, box.padding = 0.3, segment.size = 0.3, segment.alpha = 0.5,
                    max.overlaps = Inf, min.segment.length = 0) +
    scale_fill_manual(values = category_cols, drop = FALSE) +
    theme_minimal() + theme(legend.position = "none") +
    labs(x = xlab, y = ylab, fill = "Category")
}

p_pc12 <- make_scatter("PC1", "PC2", "label12", pc1_lab, pc2_lab)
p_pc23 <- make_scatter("PC2", "PC3", "label23", pc2_lab, pc3_lab)
p_pc34 <- make_scatter("PC3", "PC4", "label34", pc3_lab, pc4_lab)
scatter_row <- p_pc12 | p_pc23 | p_pc34

# ---- CDF row ----
if (!file.exists(cdf_stats_file) || !file.exists(cdf_pvals_file)) {
  stop("Missing GWAS-removed CDF inputs: ", cdf_stats_file, " and/or ", cdf_pvals_file)
}
cdf_dt <- fread(cdf_stats_file)[pc %in% pcs_cdf & flank == paste0("+", flank_label) & abs(top_frac - 0.01) < 1e-12]
cdf_p  <- fread(cdf_pvals_file)[pc %in% pcs_cdf & flank == paste0("+", flank_label) & abs(top_frac - 0.01) < 1e-12]

make_cdf_panel <- function(pc_name) {
  dtp <- cdf_dt[pc == pc_name]; pp <- cdf_p[pc == pc_name]
  ann_label <- if (nrow(pp) > 0L) sprintf("p = %.3g", pp$p_emp[1]) else NULL
  p <- ggplot() +
    geom_ribbon(data = dtp, aes(x = dist_bp / 1000, ymin = sim_lo, ymax = sim_hi),
                fill = "grey80", alpha = 0.5) +
    geom_line(data = dtp, aes(x = dist_bp / 1000, y = sim_med), color = "grey45", linewidth = 0.6) +
    geom_step(data = dtp, aes(x = dist_bp / 1000, y = obs_cdf), color = "#1f78b4", linewidth = 0.9) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
          panel.grid.minor = element_line(color = "grey95", linewidth = 0.2),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3)) +
    scale_x_continuous(limits = c(0, 2000), breaks = seq(0, 2000, by = 500))
  if (!is.null(ann_label)) p <- p + annotate("text", x = 1950, y = 0.05, label = ann_label, hjust = 1, vjust = 0, size = 2.9)
  p + annotate("text", x = 80, y = 0.95, label = pc_name, hjust = 0, vjust = 1, size = 3.2)
}

p_cdf1 <- make_cdf_panel("PC1") + labs(x = "Distance to nearest GWAS-peak window (kb)", y = "Cumulative fraction")
p_cdf2 <- make_cdf_panel("PC2") + theme(axis.title.y = element_blank()) + labs(x = NULL, y = NULL)
p_cdf3 <- make_cdf_panel("PC3") + theme(axis.title.y = element_blank()) + labs(x = NULL, y = NULL)
p_cdf4 <- make_cdf_panel("PC4") + theme(axis.title.y = element_blank()) + labs(x = NULL, y = NULL)
cdf_row <- p_cdf1 | p_cdf2 | p_cdf3 | p_cdf4

# ---- Assemble and save ----
fig3 <- (scatter_row / p_manh1 / p_manh2 / p_manh3 / cdf_row) +
  plot_layout(heights = c(2.20, 1, 1, 1, 1.15)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 18))

ggsave(out_fig, fig3, width = 11.9, height = 10.7)
message("Figure 3 complete: ", out_fig)

fig3_pc4_supp <- (p_manh4_orig + labs(tag = "a")) / (p_manh4_red + labs(tag = "b")) &
  theme(plot.tag = element_text(size = 18))
ggsave(out_fig_pc4_supp, fig3_pc4_supp, width = 11.5, height = 6.2)
message("FigureS13 complete: ", out_fig_pc4_supp)
