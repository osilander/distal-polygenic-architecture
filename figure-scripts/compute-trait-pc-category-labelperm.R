#!/usr/bin/env Rscript
# Produces: results/trait_pc_category_absrho_labelperm_pc1_pc14.tsv  (used by run-figure1.R)
#           figures/figureS4_trait_pc_category_absrho_labelperm_pc1_pc14.pdf
# Pre-computed inputs:
#   results/pca_loadings_50k.tsv (from run-figure1.R)
#   results/all-trait-50k-mean-{effects,pvals}.tsv
#   data/trait_abbrevs_categorised.txt
# Run: Rscript post-process-scripts/compute-trait-pc-category-labelperm.R
# Runtime: ~5 min (1000 permutations)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(paletteer)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

source(file.path("post-process-scripts", "helpers-50k-matrix.R"))

set.seed(20260415)
pcs <- paste0("PC", 1:14)
n_perm <- 1000L

out_tsv <- file.path("results", "trait_pc_category_absrho_labelperm_pc1_pc14.tsv")
out_fig <- file.path("figures", "figureS4_trait_pc_category_absrho_labelperm_pc1_pc14.pdf")
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

# ---- Load inputs ----
message("Loading pca_loadings_50k.tsv...")
pca_loadings_file <- file.path("results", "pca_loadings_50k.tsv")
if (!file.exists(pca_loadings_file)) stop("Missing: ", pca_loadings_file, "\nRun run-figure1.R first.")
pca_loadings <- fread(pca_loadings_file)
if (!"region_label" %in% names(pca_loadings)) {
  pca_loadings[, region_label := paste(chrom, start, end, sep = "_")]
}

message("Rebuilding scaled trait matrix...")
all_betas <- fread(file.path("results", "all-trait-50k-mean-effects.tsv"))
all_pvals <- fread(file.path("results", "all-trait-50k-mean-pvals.tsv"))
prep <- prepare_X_50k(all_betas, all_pvals)
traits_t_scaled <- t(prep$X)  # traits x windows

trait_categories <- fread(file.path("data", "trait_abbrevs_categorised.txt"))

# ---- Compute Spearman correlations ----
loadings_dt <- as.data.table(pca_loadings)
common_ids <- intersect(colnames(traits_t_scaled), loadings_dt$region_label)
if (length(common_ids) == 0L) stop("No shared window IDs between traits_t_scaled and pca_loadings.")

trait_mat <- as.matrix(traits_t_scaled[, common_ids, drop = FALSE])  # traits x windows
pc_mat    <- as.matrix(loadings_dt[match(common_ids, region_label), ..pcs])  # windows x pcs

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
  denom    <- sqrt(x_ss * y_ss)
  out      <- cov_num / denom
  out[!is.finite(out)] <- NA_real_
  out
}

rho_rows <- vector("list", length(pcs))
trait_names <- rownames(trait_mat)
for (i in seq_along(pcs)) {
  pc <- pcs[[i]]
  rho_rows[[i]] <- data.table(
    trait_abbrev = trait_names,
    pc = pc,
    spearman_rho = row_cor_with_vec(row_rank_matrix, pc_rank_mat[, i])
  )
}

dt <- merge(
  rbindlist(rho_rows, use.names = TRUE),
  as.data.table(trait_categories)[, .(trait_abbrev, phenotype_group)],
  by = "trait_abbrev", all.x = TRUE
)
dt <- dt[is.finite(spearman_rho)]
dt[is.na(phenotype_group) | phenotype_group == "", phenotype_group := "Uncategorized"]
dt[, abs_rho := abs(spearman_rho)]
dt[, pc := factor(pc, levels = pcs)]
dt[, abs_rho_rank := frank(abs_rho, ties.method = "average"), by = pc]

obs <- dt[, .(
  n_traits = .N,
  observed_median_abs_rho = median(abs_rho, na.rm = TRUE),
  observed_mean_rank = mean(abs_rho_rank, na.rm = TRUE)
), by = .(pc, phenotype_group)]

trait_lab  <- unique(dt[, .(trait_abbrev, phenotype_group)])
all_groups <- sort(unique(trait_lab$phenotype_group))

perm_rows <- vector("list", n_perm)
for (b in seq_len(n_perm)) {
  perm_map <- copy(trait_lab)
  perm_map[, phenotype_group := sample(phenotype_group, .N, replace = FALSE)]
  z <- merge(dt[, .(trait_abbrev, pc, abs_rho_rank)], perm_map, by = "trait_abbrev", all.x = FALSE)
  perm_rows[[b]] <- z[, .(perm_mean_rank = mean(abs_rho_rank, na.rm = TRUE)),
                      by = .(pc, phenotype_group)][, perm := b]
}
perm_dt <- rbindlist(perm_rows, use.names = TRUE)

null_sum <- perm_dt[, .(
  null_mean = mean(perm_mean_rank, na.rm = TRUE),
  null_sd   = sd(perm_mean_rank, na.rm = TRUE),
  q025 = quantile(perm_mean_rank, 0.025, na.rm = TRUE),
  q975 = quantile(perm_mean_rank, 0.975, na.rm = TRUE)
), by = .(pc, phenotype_group)]

p_hi <- perm_dt[obs, on = .(pc, phenotype_group)][
  , .(p_emp_high = (1 + sum(perm_mean_rank >= observed_mean_rank, na.rm = TRUE)) / (.N + 1)),
  by = .(pc, phenotype_group)
]

res <- Reduce(function(x, y) merge(x, y, by = c("pc", "phenotype_group"), all = TRUE),
              list(obs, null_sum, p_hi))
res[, z_score    := fifelse(is.finite(null_sd) & null_sd > 0,
                            (observed_mean_rank - null_mean) / null_sd, NA_real_)]
res[, delta_from_null := observed_mean_rank - null_mean]
res[, sig := fifelse(p_emp_high < 0.001, "***",
             fifelse(p_emp_high < 0.01,  "**",
             fifelse(p_emp_high < 0.05,  "*", "")))]
res[, label_rho := sprintf("%.2f", observed_median_abs_rho)]
res[, label_sig := sig]

group_order <- res[, .(
  best_p    = min(p_emp_high, na.rm = TRUE),
  max_z     = max(z_score, na.rm = TRUE),
  max_delta = max(delta_from_null, na.rm = TRUE)
), by = phenotype_group][order(best_p, -max_z, -max_delta), phenotype_group]
res[, phenotype_group := factor(phenotype_group, levels = group_order)]
res[, pc := factor(pc, levels = pcs)]
fwrite(res[order(pc, phenotype_group)], out_tsv, sep = "\t")

rho_cap        <- 0.2
ag_sunset_cols <- rev(colorRampPalette(as.character(paletteer::paletteer_d("rcartocolor::ag_Sunset")))(256))
custom_theme_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
  "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
  "#E5C494", "#B3B3B3"
)
all_trait_categories <- sort(unique(na.omit(as.character(trait_categories$phenotype_group))))
all_trait_categories <- all_trait_categories[nzchar(all_trait_categories)]
category_colors_all  <- setNames(custom_theme_colors[seq_along(all_trait_categories)], all_trait_categories)

res[, label_col  := fifelse(observed_median_abs_rho > 0.12, "white", "black")]
res[, pc_index   := as.integer(sub("^PC", "", as.character(pc)))]

strip_dt <- unique(res[, .(phenotype_group)])
strip_dt[, strip_col := unname(category_colors_all[as.character(phenotype_group)])]

p <- ggplot(res, aes(x = pc_index, y = phenotype_group, fill = observed_median_abs_rho)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_tile(data = strip_dt, aes(x = 0.35, y = phenotype_group), inherit.aes = FALSE,
            fill = strip_dt$strip_col, width = 0.18, height = 0.98, linewidth = 0) +
  geom_text(aes(label = label_rho, color = label_col), size = 3.2, na.rm = TRUE) +
  geom_text(aes(label = label_sig, color = label_col), size = 2.6, nudge_y = 0.28, na.rm = TRUE) +
  scale_fill_gradientn(colours = ag_sunset_cols, limits = c(0, rho_cap),
                       oob = scales::squish, na.value = "grey90", name = "Median |rho|") +
  scale_color_identity() +
  scale_x_continuous(breaks = 1:14, labels = paste0("PC", 1:14),
                     expand = expansion(add = c(0.45, 0.2))) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_blank()
  )

ggsave(out_fig, p, width = 8.5, height = 6.16)
message("Wrote: ", out_tsv)
message("Wrote: ", out_fig)
