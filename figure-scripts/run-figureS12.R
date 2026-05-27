#!/usr/bin/env Rscript
# Produces: figures/figureS12_allbyall_procrustes_heatmaps.pdf
# Pre-computed inputs required:
#   results/allbyall_cosine_matrices/allbyall_cosine_long_pc1_pc4.tsv
#   results/allbyall_cosine_matrices/allbyall_procrustes_long_k4_k8_k25_k50_k75.tsv
#   results/allbyall_cosine_matrices/entrysign_perm100_procrustes_normalized_12x12.tsv
# Run: Rscript post-process-scripts/run-figureS12.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(RColorBrewer)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

in_file          <- file.path("results", "allbyall_cosine_matrices", "allbyall_cosine_long_pc1_pc4.tsv")
in_proc_file     <- file.path("results", "allbyall_cosine_matrices", "allbyall_procrustes_long_k4_k8_k25_k50_k75.tsv")
in_norm_proc_file <- file.path("results", "allbyall_cosine_matrices", "entrysign_perm100_procrustes_normalized_12x12.tsv")
out_proc_fig     <- file.path("figures", "figureS12_allbyall_procrustes_heatmaps.pdf")
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

for (f in c(in_file, in_proc_file, in_norm_proc_file)) {
  if (!file.exists(f)) stop("Missing pre-computed input: ", f)
}

dt           <- fread(in_file)
dt_proc      <- fread(in_proc_file)
dt_norm_proc <- fread(in_norm_proc_file)

label_map <- c(
  "windows_25kb" = "25Kb", "observed_50k_signed" = "50Kb", "windows_100kb" = "100Kb",
  "windows_200kb" = "200 Kb", "windows_500kb" = "500Kb", "windows_1mb" = "1Mb",
  "permuted_within_500kb" = "500Kb rand.", "pickrell_randomized" = "LD blocks rand.",
  "entry_sign_randomized" = "Entry sign rand.", "gwas_removed_plus150kb" = "No GWAS loci +- 150Kb",
  "gwas_removed_plus100kb" = "No GWAS loci +- 100Kb",
  "abs_from_signed_matrix" = "Matrix abs. value", "abs_from_windows_filtered" = "Variant abs value"
)
conditions <- names(label_map)
ordered_conditions <- c(
  "windows_25kb", "observed_50k_signed", "windows_100kb", "windows_200kb", "windows_500kb",
  "windows_1mb", "permuted_within_500kb", "pickrell_randomized",
  "gwas_removed_plus150kb", "gwas_removed_plus100kb",
  "abs_from_signed_matrix", "abs_from_windows_filtered", "entry_sign_randomized"
)

dt      <- dt[condition_a %in% conditions & condition_b %in% conditions & !is.na(condition_a) & !is.na(condition_b)]
dt_proc <- dt_proc[condition_a %in% conditions & condition_b %in% conditions & !is.na(condition_a) & !is.na(condition_b)]

spec_cols <- colorRampPalette(brewer.pal(11, "Spectral"))(256)

to_long <- function(mat) {
  mat  <- mat[ordered_conditions, ordered_conditions, drop = FALSE]
  rlab <- label_map[rownames(mat)]; clab <- label_map[colnames(mat)]
  long <- as.data.table(as.table(mat))
  setnames(long, c("row_key", "col_key", "value"))
  long[, row := factor(label_map[as.character(row_key)], levels = rev(rlab))]
  long[, col := factor(label_map[as.character(col_key)], levels = clab)]
  long[, label_col := ifelse(value >= 0.15 & value <= 0.85, "black", "white")]
  long
}

make_proc_matrix <- function(k_val) {
  sub <- dt_proc[k == k_val, .(condition_a, condition_b, value = procrustes_similarity)]
  mat <- matrix(NA_real_, nrow = length(conditions), ncol = length(conditions),
                dimnames = list(conditions, conditions))
  for (i in seq_len(nrow(sub))) mat[sub$condition_a[i], sub$condition_b[i]] <- sub$value[i]
  mat
}

make_norm_proc_matrix <- function(k_val) {
  sub <- dt_norm_proc[k == k_val, .(condition_a, condition_b, value = proc_norm)]
  mat <- matrix(NA_real_, nrow = length(conditions), ncol = length(conditions),
                dimnames = list(conditions, conditions))
  for (i in seq_len(nrow(sub))) mat[sub$condition_a[i], sub$condition_b[i]] <- sub$value[i]
  mat
}

make_heat_from_long <- function(long, fill_name, title) {
  ggplot(long, aes(x = col, y = row, fill = value)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.2f", value), color = label_col), size = 2.05, fontface = "bold") +
    scale_color_identity() +
    scale_fill_gradientn(colours = spec_cols, limits = c(0, 1), name = fill_name) +
    theme_minimal(base_size = 9) +
    theme(text = element_text(face = "bold"), panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6.8),
          axis.text.y = element_text(size = 6.8), axis.title = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
          plot.title = element_text(size = 9, face = "bold")) +
    labs(title = title)
}

pk4_raw  <- make_heat_from_long(to_long(make_proc_matrix(4L)),        "Procrustes",            "k=4")
pk4_norm <- make_heat_from_long(to_long(make_norm_proc_matrix(4L)),   "Normalized Procrustes", "k=4 normalized")
pk8_norm <- make_heat_from_long(to_long(make_norm_proc_matrix(8L)),   "Normalized Procrustes", "k=8 normalized")
pk25_norm <- make_heat_from_long(to_long(make_norm_proc_matrix(25L)), "Normalized Procrustes", "k=25 normalized")
pk50_norm <- make_heat_from_long(to_long(make_norm_proc_matrix(50L)), "Normalized Procrustes", "k=50 normalized")
pk75_norm <- make_heat_from_long(to_long(make_norm_proc_matrix(75L)), "Normalized Procrustes", "k=75 normalized")

fig_proc <- cowplot::plot_grid(pk4_raw, pk4_norm, pk8_norm, pk25_norm, pk50_norm, pk75_norm,
                               ncol = 3, align = "hv")
ggsave(out_proc_fig, fig_proc, width = 17.7408, height = 7.92)
message("FigureS12 complete: ", out_proc_fig)
