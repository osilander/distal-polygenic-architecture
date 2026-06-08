#!/usr/bin/env Rscript
# Produces: figures/figureS19_gwas_removed_pc_correlations_pc1_pc4.pdf,
#           results/gwas_removed_pc_correlations_pc1_pc4.tsv
# Pre-computed inputs: results/pca_loadings_50k.tsv,
#                      results/gwas_removed_distance_topload_50k_plus100kb/reduced_pca_window_loadings.tsv
# Run: Rscript figure-scripts/run-figureS19.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(irlba)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

pcs      <- paste0("PC", 1:4)
flank_kb <- 100L
out_tsv  <- file.path("results", "gwas_removed_pc_correlations_pc1_pc4.tsv")
out_fig  <- file.path("figures", "figureS19_gwas_removed_pc_correlations_pc1_pc4.pdf")
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)

pca_loadings_file <- file.path("results", "pca_loadings_50k.tsv")
if (!file.exists(pca_loadings_file)) stop("Missing: ", pca_loadings_file, "\nRun run-figure1.R first.")
pca_loadings <- fread(pca_loadings_file)

red_file <- file.path("results", paste0("gwas_removed_distance_topload_50k_plus", flank_kb, "kb"),
                      "reduced_pca_window_loadings.tsv")
if (!file.exists(red_file)) stop("Missing reduced loadings: ", red_file)

parse_region_cols <- function(dt) {
  x <- copy(as.data.table(dt))
  if (!all(c("chrom", "start", "end") %in% names(x))) {
    if (!"region_label" %in% names(x)) stop("Need region_label to parse chrom/start/end.")
    x[, c("chrom", "start", "end") := tstrsplit(region_label, "_", fixed = TRUE)]
  }
  x[, `:=`(chrom = as.integer(chrom), start = as.integer(start), end = as.integer(end))]
  x[, window_id := paste(chrom, start, end, sep = "_")]
  x
}

pl_full <- parse_region_cols(pca_loadings)
pl_red  <- parse_region_cols(fread(red_file))

shared_ids    <- intersect(pl_full$window_id, pl_red$window_id)
full_mat      <- as.matrix(pl_full[match(shared_ids, window_id), ..pcs])
red_mat       <- as.matrix(pl_red[match(shared_ids, window_id), ..pcs])

cor_mat       <- cor(full_mat, red_mat, use = "pairwise.complete.obs")
best_sign <- sapply(seq_len(ncol(cor_mat)), function(j) {
  best_i <- which.max(abs(cor_mat[, j]))
  sign(cor_mat[best_i, j])
})
best_sign[best_sign == 0] <- 1
red_mat_aligned   <- sweep(red_mat, 2, best_sign, `*`)
cor_mat_aligned   <- cor(full_mat, red_mat_aligned, use = "pairwise.complete.obs")

dt <- as.data.table(as.table(cor_mat_aligned))
setnames(dt, c("full_pc", "removed_pc", "correlation"))
dt[, `:=`(
  full_pc    = factor(as.character(full_pc),    levels = pcs),
  removed_pc = factor(as.character(removed_pc), levels = pcs)
)]
fwrite(dt, out_tsv, sep = "\t")

p <- ggplot(dt, aes(x = removed_pc, y = full_pc, fill = correlation)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", correlation)), size = 4) +
  scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac", midpoint = 0, limits = c(-1, 1)) +
  coord_equal() +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 11)) +
  labs(x = "GWAS-removed PCA", y = "Original PCA", fill = "r")

ggsave(out_fig, p, width = 5.6, height = 4.9)
message("Wrote: ", out_tsv)
message("FigureS19 complete: ", out_fig)
