#!/usr/bin/env Rscript
# Produces: results/figure2_cache_50k.rds  (used by run-figure1.R)
# Pre-computed inputs:
#   results/svd-nulls-50k/  (from generate-svd-nulls-50k.R)
#   results/chunk-withinperm-nulls-50k/  (from prepare-chunk-withinperm-base-50k.R + analyse-chunk-withinperm-nulls-50k.R)
#   results/all-trait-50k-mean-{effects,pvals}.tsv  (from build-all-trait-matrices.R)
# Run: Rscript post-process-scripts/compute-figure2-cache.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(ggrepel)
})

setDTthreads(1L)
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(project_root, "post-process-scripts", "helpers-50k-matrix.R"))

safe_fread <- function(path, ...) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    return(fread(cmd = sprintf("gzip -dc %s", shQuote(path)), ...))
  }
  fread(path, ...)
}

obs_sing_file    <- file.path(project_root, "results", "svd-nulls-50k", "observed", "observed_50k_singular_values.tsv")
entry_flip_dir   <- file.path(project_root, "results", "svd-nulls-50k", "entry_flip")
pickrell_perm_dir <- file.path(project_root, "results", "svd-nulls-50k", "withinblock_perm")
effects_file     <- file.path(project_root, "results", "all-trait-50k-mean-effects.tsv")
pvals_file       <- file.path(project_root, "results", "all-trait-50k-mean-pvals.tsv")
trait_map_file   <- file.path(project_root, "data", "trait_abbrevs_categorised.txt")
chunk_runs_file  <- file.path(project_root, "results", "chunk-withinperm-nulls-50k", "chunk_withinperm_scree_runs_50k.tsv.gz")
out_cache        <- file.path(project_root, "results", "figure2_cache_50k.rds")

stopifnot(
  file.exists(obs_sing_file),
  file.exists(effects_file),
  file.exists(pvals_file),
  file.exists(trait_map_file),
  dir.exists(entry_flip_dir),
  dir.exists(pickrell_perm_dir),
  file.exists(chunk_runs_file)
)

custom_theme_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
  "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
  "#E5C494", "#B3B3B3"
)

anchor_map <- c(PC1 = "BMI_03", PC2 = "StandHgt", PC3 = "Neurotic_30", PC4 = "DiaBP_A59")

scheme_colors <- c(
  "observed" = custom_theme_colors[3],
  "sign_randomised" = custom_theme_colors[2],
  "withinblock_perm" = custom_theme_colors[1],
  "permuted_within_500Kb" = custom_theme_colors[4]
)

label_map <- c(
  "observed" = "observed",
  "sign_randomised" = "sign randomised",
  "withinblock_perm" = "permuted within LD blocks",
  "permuted_within_500Kb" = "permuted within 500Kb"
)

load_null_scree <- function(dir_path, prefix_regex) {
  ff <- list.files(dir_path, pattern = paste0("^", prefix_regex, "_iter[0-9]{3}_50k_singular_values\\.tsv$"), full.names = TRUE)
  if (length(ff) == 0L) stop("No matching scree files in ", dir_path, " for ", prefix_regex)
  rbindlist(lapply(ff, function(f) {
    x <- fread(f)[, .(pc, var_explained)]
    x[, iter := as.integer(sub(".*_iter([0-9]{3})_50k_singular_values\\.tsv$", "\\1", basename(f)))]
    x
  }), use.names = TRUE, fill = TRUE)
}

load_first_trait_vectors <- function(dir_path, prefix_regex) {
  ff <- sort(list.files(dir_path, pattern = paste0("^", prefix_regex, "_iter[0-9]{3}_50k_trait_vectors\\.tsv\\.gz$"), full.names = TRUE))
  if (length(ff) == 0L) stop("No trait-vector files in ", dir_path, " for ", prefix_regex)
  safe_fread(ff[1L], select = c("trait", "PC1", "PC2", "PC3"))
}

full_decomp_trait_pc123 <- function(X) {
  n <- nrow(X)
  C <- crossprod(X) / max(1, n - 1)
  eg <- eigen(C, symmetric = TRUE)
  vals <- pmax(eg$values, 0)
  V <- eg$vectors
  keep <- which(vals > 1e-12)
  if (length(keep) < 3L) stop("Fewer than 3 non-zero PCs.")
  V <- V[, keep, drop = FALSE]
  out <- as.data.table(V[, 1:3, drop = FALSE])
  setnames(out, c("PC1", "PC2", "PC3"))
  out
}

join_trait_groups <- function(dt, trait_map) {
  out <- merge(
    dt,
    trait_map[, .(trait_abbrev, phenotype_group)],
    by.x = "trait",
    by.y = "trait_abbrev",
    all.x = TRUE
  )
  out[is.na(phenotype_group), phenotype_group := "Uncategorized"]
  out
}

align_anchor <- function(dt) {
  x <- copy(dt)
  for (pc in intersect(names(anchor_map), names(x))) {
    tr <- anchor_map[[pc]]
    idx <- which(x$trait == tr)
    if (length(idx) == 0L) next
    aval <- x[[pc]][idx[1]]
    if (is.finite(aval) && aval < 0) x[[pc]] <- -x[[pc]]
  }
  x
}

extreme_traits <- function(dt, pc, n_each = 3L) {
  z <- copy(dt[is.finite(get(pc)), .(trait, val = get(pc))])
  if (nrow(z) == 0L) return(character())
  pos <- z[order(-val)][seq_len(min(n_each, .N)), trait]
  neg <- z[order(val)][seq_len(min(n_each, .N)), trait]
  unique(c(pos, neg))
}

obs <- fread(obs_sing_file)[, .(pc, var_explained)]
obs[, pc_index := as.integer(sub("^PC", "", pc))]
obs[, scheme := "observed"]
obs[, `:=`(q25 = NA_real_, q75 = NA_real_)]

sc_sign <- load_null_scree(entry_flip_dir, "entry_flip")
sc_sign_sum <- sc_sign[, .(
  var_explained = mean(var_explained, na.rm = TRUE),
  q25 = as.numeric(quantile(var_explained, 0.25, na.rm = TRUE)),
  q75 = as.numeric(quantile(var_explained, 0.75, na.rm = TRUE))
), by = pc]
sc_sign_sum[, `:=`(pc_index = as.integer(sub("^PC", "", pc)), scheme = "sign_randomised")]

sc_pick <- load_null_scree(pickrell_perm_dir, "withinblock_perm")
sc_pick_sum <- sc_pick[, .(
  var_explained = mean(var_explained, na.rm = TRUE),
  q25 = as.numeric(quantile(var_explained, 0.25, na.rm = TRUE)),
  q75 = as.numeric(quantile(var_explained, 0.75, na.rm = TRUE))
), by = pc]
sc_pick_sum[, `:=`(pc_index = as.integer(sub("^PC", "", pc)), scheme = "withinblock_perm")]

sc_chunk <- safe_fread(chunk_runs_file)
sc_chunk <- sc_chunk[scheme == "chunk500kb_withinperm", .(pc, var_explained)]
if (nrow(sc_chunk) == 0L) stop("No chunk500kb_withinperm rows in ", chunk_runs_file)
sc_chunk_sum <- sc_chunk[, .(
  var_explained = mean(var_explained, na.rm = TRUE),
  q25 = as.numeric(quantile(var_explained, 0.25, na.rm = TRUE)),
  q75 = as.numeric(quantile(var_explained, 0.75, na.rm = TRUE))
), by = pc]
sc_chunk_sum[, `:=`(pc_index = as.integer(sub("^PC", "", pc)), scheme = "permuted_within_500Kb")]

scree_dt <- rbindlist(list(
  obs[, .(pc, pc_index, scheme, var_explained, q25, q75)],
  sc_sign_sum[, .(pc, pc_index, scheme, var_explained, q25, q75)],
  sc_pick_sum[, .(pc, pc_index, scheme, var_explained, q25, q75)],
  sc_chunk_sum[, .(pc, pc_index, scheme, var_explained, q25, q75)]
), use.names = TRUE, fill = TRUE)
scree_dt <- scree_dt[pc_index <= 60]
scree_dt[, scheme := factor(scheme, levels = c("observed", "sign_randomised", "withinblock_perm", "permuted_within_500Kb"))]

trait_map <- fread(trait_map_file)
sign_pc <- load_first_trait_vectors(entry_flip_dir, "entry_flip")
sign_pc <- join_trait_groups(sign_pc, trait_map)
group_levels <- sort(unique(na.omit(as.character(trait_map$phenotype_group))))
group_levels <- group_levels[nzchar(group_levels)]
if ("Uncategorized" %in% sign_pc$phenotype_group) group_levels <- c(group_levels, "Uncategorized")
if (length(group_levels) > length(custom_theme_colors)) {
  stop("Not enough colors for phenotype categories")
}
group_cols <- setNames(custom_theme_colors[seq_along(group_levels)], group_levels)
sign_pc[, phenotype_group := factor(phenotype_group, levels = group_levels)]
sign_pc <- align_anchor(sign_pc)

pick_pc <- load_first_trait_vectors(pickrell_perm_dir, "withinblock_perm")
pick_pc <- join_trait_groups(pick_pc, trait_map)
pick_pc[, phenotype_group := factor(phenotype_group, levels = group_levels)]
pick_pc <- align_anchor(pick_pc)

labels_pc12 <- unique(c(extreme_traits(pick_pc, "PC1", 3L), extreme_traits(pick_pc, "PC2", 3L)))
labels_pc23 <- unique(c(extreme_traits(pick_pc, "PC2", 3L), extreme_traits(pick_pc, "PC3", 3L)))

eff   <- fread(effects_file)
pvals <- fread(pvals_file)
prep  <- prepare_X_50k(eff, pvals, impute_number = 5L, p_thresh = NULL)
X_500 <- permute_within_chunks(prep$X_scaled, prep$coords, chunk_size = 500000L, seed = 51001L)
chunk_pc <- full_decomp_trait_pc123(X_500)
chunk_pc[, trait := colnames(prep$X_scaled)]
chunk_pc <- join_trait_groups(chunk_pc, trait_map)
chunk_pc[, phenotype_group := factor(phenotype_group, levels = group_levels)]
chunk_pc <- align_anchor(chunk_pc)

cache <- list(
  scree_dt = scree_dt,
  sign_pc = sign_pc,
  pick_pc = pick_pc,
  chunk_pc = chunk_pc,
  labels_pc12 = labels_pc12,
  labels_pc23 = labels_pc23,
  group_levels = group_levels,
  group_cols = group_cols,
  scheme_colors = scheme_colors,
  label_map = label_map,
  custom_theme_colors = custom_theme_colors
)

saveRDS(cache, out_cache)
message("Wrote: ", out_cache)
