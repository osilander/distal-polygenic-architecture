#!/usr/bin/env Rscript
# Produces: figures/figureS3_trait_pc_category_absrho_boxplot_pc1_pc8.pdf,
#           results/trait_pc_category_absrho_summary_pc1_pc8.tsv
# Pre-computed inputs: results/pca_loadings_50k.tsv, results/all-trait-50k-mean-effects.tsv,
#                      results/all-trait-50k-mean-pvals.tsv, data/trait_abbrevs_categorised.txt
# Run: Rscript post-process-scripts/run-figureS3.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

out_fig <- file.path("figures", "figureS3_trait_pc_category_absrho_boxplot_pc1_pc8.pdf")
out_tsv <- file.path("results", "trait_pc_category_absrho_summary_pc1_pc8.tsv")
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)

# ---- Load pca_loadings ----
pca_loadings_file <- file.path("results", "pca_loadings_50k.tsv")
if (!file.exists(pca_loadings_file)) stop("Missing: ", pca_loadings_file, "\nRun run-figure1.R first.")
pca_loadings <- fread(pca_loadings_file)

# ---- Rebuild traits_t_scaled ----
message("FigureS3: rebuilding scaled trait matrix (~2-3 min)...")
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
if (!("Uncategorized" %in% names(category_colors_all))) {
  category_colors_all <- c(category_colors_all, "Uncategorized" = "grey60")
}

pcs <- paste0("PC", 1:8)

loadings_dt <- as.data.table(pca_loadings)
if (!"region_label" %in% names(loadings_dt)) {
  loadings_dt[, region_label := paste(chrom, start, end, sep = "_")]
}
common_ids <- intersect(colnames(traits_t_scaled), loadings_dt$region_label)
if (length(common_ids) == 0L) stop("No shared window IDs between traits_t_scaled and pca_loadings.")

trait_mat <- as.matrix(traits_t_scaled[, common_ids, drop = FALSE])  # traits x windows
pc_mat    <- as.matrix(loadings_dt[match(common_ids, region_label), ..pcs])  # windows x 8
rownames(pc_mat) <- common_ids

row_rank_matrix <- t(apply(trait_mat, 1, rank, ties.method = "average"))
pc_rank_mat     <- apply(pc_mat, 2, rank, ties.method = "average")
if (is.null(dim(pc_rank_mat))) pc_rank_mat <- matrix(pc_rank_mat, ncol = length(pcs))
colnames(pc_rank_mat) <- pcs

row_cor_with_vec <- function(X, y) {
  n        <- ncol(X)
  row_mean <- rowMeans(X)
  y_mean   <- mean(y)
  x_ss     <- rowSums((X - row_mean)^2)
  y_ss     <- sum((y - y_mean)^2)
  cov_num  <- as.numeric(X %*% y) - n * row_mean * y_mean
  out      <- cov_num / sqrt(x_ss * y_ss)
  out[!is.finite(out)] <- NA_real_
  out
}

rho_rows   <- vector("list", length(pcs))
trait_names <- rownames(trait_mat)
for (i in seq_along(pcs)) {
  rho_rows[[i]] <- data.table(
    trait_abbrev = trait_names, pc = pcs[[i]],
    spearman_rho = row_cor_with_vec(row_rank_matrix, pc_rank_mat[, i])
  )
}

dt <- merge(
  rbindlist(rho_rows, use.names = TRUE),
  as.data.table(trait_categories)[, .(trait_abbrev, trait_full, phenotype_group)],
  by = "trait_abbrev", all.x = TRUE
)
dt <- dt[is.finite(spearman_rho)]
dt[, abs_rho := abs(spearman_rho)]
dt[is.na(phenotype_group) | phenotype_group == "", phenotype_group := "Uncategorized"]
dt[, pc := factor(pc, levels = pcs)]

ord_dt <- dt[, .(median_abs_rho = median(abs_rho, na.rm = TRUE)), by = .(pc, phenotype_group)]
setorder(ord_dt, pc, -median_abs_rho, phenotype_group)
ord_dt[, phenotype_group_within := paste(phenotype_group, pc, sep = "___")]
level_order <- ord_dt$phenotype_group_within

dt[, phenotype_group_within := paste(phenotype_group, pc, sep = "___")]
dt[, phenotype_group_within := factor(phenotype_group_within, levels = level_order)]

summary_dt <- dt[, .(
  n_traits = .N,
  median_abs_rho = median(abs_rho, na.rm = TRUE),
  q25_abs_rho = quantile(abs_rho, 0.25, na.rm = TRUE),
  q75_abs_rho = quantile(abs_rho, 0.75, na.rm = TRUE)
), by = .(pc, phenotype_group)]
fwrite(summary_dt, out_tsv, sep = "\t")

label_map <- setNames(ord_dt$phenotype_group, ord_dt$phenotype_group_within)

p <- ggplot(dt, aes(x = phenotype_group_within, y = abs_rho, fill = phenotype_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.72, color = "grey20", linewidth = 0.3) +
  geom_jitter(width = 0.14, height = 0, size = 0.5, alpha = 0.35, color = "black") +
  scale_fill_manual(values = category_colors_all, drop = FALSE) +
  scale_x_discrete(labels = function(x) unname(label_map[x])) +
  facet_wrap(~pc, nrow = 3, scales = "free_x") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1),
        legend.position = "none", strip.text = element_text(size = 11)) +
  labs(x = "Phenotype category", y = "|Trait-PC Spearman correlation|")

ggsave(out_fig, p, width = 16, height = 17.2)
message("Wrote: ", out_tsv)
message("FigureS3 complete: ", out_fig)
