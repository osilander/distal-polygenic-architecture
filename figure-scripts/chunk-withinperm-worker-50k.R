#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

setDTthreads(1L)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2L) {
  stop("Usage: Rscript chunk-withinperm-worker-50k.R <scheme> <iter>")
}

scheme <- args[[1]]
iter <- as.integer(args[[2]])
if (!scheme %in% c("chunk250kb_withinperm", "chunk500kb_withinperm")) stop("Invalid scheme")
if (!is.finite(iter) || iter < 1L) stop("iter must be >= 1")

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
base_file <- file.path(project_root, "results", "chunk-withinperm-nulls-50k", "hpc_base", "chunk_withinperm_base_50k.rds")
run_dir <- file.path(project_root, "results", "chunk-withinperm-nulls-50k", "hpc_runs")
metric_dir <- file.path(run_dir, "metrics")
scree_dir <- file.path(run_dir, "scree")
dir.create(metric_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(scree_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(base_file))

base <- readRDS(base_file)
X <- base$X
k_vals <- as.integer(base$k_vals)
max_k <- as.integer(base$max_k)

decomp_from_matrix <- function(X) {
  n <- nrow(X)
  C <- crossprod(X) / max(1, n - 1)
  eg <- eigen(C, symmetric = TRUE)
  vals <- pmax(eg$values, 0)
  keep <- which(vals > 1e-12)
  vals <- vals[keep]
  V <- eg$vectors[, keep, drop = FALSE]
  sdev <- sqrt(vals)
  d <- sdev * sqrt(max(1, n - 1))
  U <- X %*% V
  U <- sweep(U, 2, d, `/`)
  list(U = U, V = V, d = d, var = vals / sum(vals))
}

mean_cos_principal_angles <- function(Uobs, Urand, k) {
  kk <- min(k, ncol(Uobs), ncol(Urand))
  A <- crossprod(Uobs[, seq_len(kk), drop = FALSE], Urand[, seq_len(kk), drop = FALSE])
  cs <- svd(A, nu = 0, nv = 0)$d
  cs <- pmin(1, pmax(0, cs))
  mean(cs)
}

trait_abs_cos_upper <- function(V, d, k) {
  kk <- min(k, ncol(V), length(d))
  T <- sweep(V[, seq_len(kk), drop = FALSE], 2, d[seq_len(kk)], `*`)
  rn <- sqrt(rowSums(T^2))
  rn[rn == 0] <- NA_real_
  T <- T / rn
  K <- abs(T %*% t(T))
  as.numeric(K[upper.tri(K, diag = FALSE)])
}

effective_rank <- function(d) {
  p <- (d^2) / sum(d^2)
  exp(-sum(p * log(p + 1e-300)))
}

permute_within_groups <- function(X, idx_by_group) {
  n <- nrow(X); p <- ncol(X)
  out <- matrix(NA_real_, nrow = n, ncol = p)
  rownames(out) <- rownames(X)
  colnames(out) <- colnames(X)
  for (j in seq_len(p)) {
    v <- X[, j]
    for (g in seq_along(idx_by_group)) {
      ix <- idx_by_group[[g]]
      if (length(ix) > 1L) v[ix] <- v[ix][sample.int(length(ix))]
    }
    out[, j] <- v
  }
  out
}

grp_obj <- if (scheme == "chunk250kb_withinperm") base$groups$chunk250 else base$groups$chunk500

seed_base <- if (scheme == "chunk250kb_withinperm") 70250L else 70500L
set.seed(seed_base + iter)
Xr <- permute_within_groups(X, grp_obj$idx)
dec <- decomp_from_matrix(Xr)
window_dir <- file.path(run_dir, "window_vectors")                                                                                                                            
dir.create(window_dir, recursive = TRUE, showWarnings = FALSE)
win_vec <- as.data.table(dec$U[, 1:4])                                                                                                                                        
setnames(win_vec, paste0("PC", 1:4))                                                                                                                                          
win_vec[, window_id := rownames(dec$U)]                                                                                                                                       
setcolorder(win_vec, c("window_id", paste0("PC", 1:4)))                                                                                                                       
fwrite(win_vec, file.path(window_dir,                                                                                                                                         
  sprintf("%s_iter%03d_window_vectors_pc1_4.tsv.gz", scheme, iter)), sep = "\t")  


rows <- list()
idx <- 1L
for (k in k_vals) {
  kk <- min(k, max_k, ncol(dec$U), ncol(dec$V), length(dec$d))
  rows[[idx]] <- data.table(
    scheme = scheme, iter = iter, metric = "window_subspace_mean_cos", k = kk,
    value = mean_cos_principal_angles(base$obs_U_top, dec$U, kk)
  ); idx <- idx + 1L

  rows[[idx]] <- data.table(
    scheme = scheme, iter = iter, metric = "trait_abs_cos_upper_corr", k = kk,
    value = suppressWarnings(cor(base$obs_upper[[as.character(kk)]], trait_abs_cos_upper(dec$V, dec$d, kk), method = "pearson"))
  ); idx <- idx + 1L

  rows[[idx]] <- data.table(
    scheme = scheme, iter = iter, metric = "cumvar_k", k = kk,
    value = sum(dec$var[seq_len(min(kk, length(dec$var)))])
  ); idx <- idx + 1L
}
rows[[idx]] <- data.table(scheme = scheme, iter = iter, metric = "pc1_var_explained", k = NA_integer_, value = dec$var[1]); idx <- idx + 1L
rows[[idx]] <- data.table(scheme = scheme, iter = iter, metric = "effective_rank", k = NA_integer_, value = effective_rank(dec$d)); idx <- idx + 1L

metrics <- rbindlist(rows, use.names = TRUE, fill = TRUE)
scree <- data.table(scheme = scheme, iter = iter, pc = paste0("PC", seq_along(dec$var)), var_explained = dec$var)

fwrite(metrics, file.path(metric_dir, sprintf("%s_iter%03d_metrics.tsv", scheme, iter)), sep = "\t")
fwrite(scree, file.path(scree_dir, sprintf("%s_iter%03d_scree.tsv.gz", scheme, iter)), sep = "\t")

message("Wrote worker outputs for ", scheme, " iter ", iter)
