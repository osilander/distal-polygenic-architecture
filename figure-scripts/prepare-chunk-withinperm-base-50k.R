#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

setDTthreads(1L)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
effects_file <- file.path(project_root, "results", "all-trait-50k-mean-effects.tsv")
pvals_file <- file.path(project_root, "results", "all-trait-50k-mean-pvals.tsv")
out_dir <- file.path(project_root, "results", "chunk-withinperm-nulls-50k", "hpc_base")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(project_root, "figure-scripts", "helpers-50k-matrix.R"))

impute_number <- 5L
p_thresh <- NULL
k_vals <- c(5L, 8L, 12L)
max_k <- max(k_vals)

message("Loading and preparing canonical matrix...")
eff <- fread(effects_file)
pvals <- fread(pvals_file)
prep <- prepare_X_50k(eff, pvals, impute_number, p_thresh)
X <- prep$X
coords <- prep$coords

message("Computing observed decomposition...")
obs <- decomp_from_matrix(X)
obs_upper <- lapply(k_vals, function(k) trait_abs_cos_upper(obs$V, obs$d, k))
names(obs_upper) <- as.character(k_vals)

coords_chr <- gsub("^chr", "", as.character(coords$chrom))
coords_start <- as.integer(coords$start)
grp250 <- factor(paste(coords_chr, coords_start %/% 250000L, sep = ":"))
grp500 <- factor(paste(coords_chr, coords_start %/% 500000L, sep = ":"))
idx250 <- split(seq_len(nrow(X)), grp250)
idx500 <- split(seq_len(nrow(X)), grp500)

base <- list(
  X = X,
  coords = coords,
  trait_names = colnames(X),
  k_vals = k_vals,
  max_k = max_k,
  obs_U_top = obs$U[, seq_len(max_k), drop = FALSE],
  obs_V_top = obs$V[, seq_len(max_k), drop = FALSE],
  obs_d_top = obs$d[seq_len(max_k)],
  obs_var = obs$var,
  obs_effective_rank = effective_rank(obs$d),
  obs_upper = obs_upper,
  groups = list(
    chunk250 = list(id = as.character(grp250), idx = idx250, n = length(idx250)),
    chunk500 = list(id = as.character(grp500), idx = idx500, n = length(idx500))
  )
)

saveRDS(base, file.path(out_dir, "chunk_withinperm_base_50k.rds"), compress = "xz")

obs_dt <- rbindlist(c(
  list(data.table(metric = "pc1_var_explained", k = NA_integer_, value = obs$var[1])),
  lapply(k_vals, function(k) data.table(metric = "cumvar_k", k = k, value = sum(obs$var[seq_len(min(k, length(obs$var)))]))),
  list(data.table(metric = "effective_rank", k = NA_integer_, value = base$obs_effective_rank))
), use.names = TRUE)
fwrite(obs_dt, file.path(out_dir, "observed_metrics_50k.tsv"), sep = "\t")

meta_dt <- data.table(
  window_id = rownames(X),
  chrom = coords$chrom,
  start = coords$start,
  end = coords$end,
  chunk250_id = base$groups$chunk250$id,
  chunk500_id = base$groups$chunk500$id
)
fwrite(meta_dt, file.path(out_dir, "window_chunk_map_50k.tsv.gz"), sep = "\t")

scree_obs <- data.table(pc = paste0("PC", seq_along(obs$var)), var_explained = obs$var)
fwrite(scree_obs, file.path(out_dir, "observed_scree_50k.tsv"), sep = "\t")

message("Wrote base object and observed references to: ", out_dir)
