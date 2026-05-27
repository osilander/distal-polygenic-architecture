#!/usr/bin/env Rscript
# compute-gwasremoved-withinperm-peaks.R
#
# Within-block permutation null for the GWAS-removed PCA loading concentration.
# Correctly permutes only the retained-window sub-matrix (not the full matrix),
# runs SVD on each permutation, then computes Gini and run-length metrics.
#
# Outputs:
#   results/gwasremoved_withinperm_peaks_iter_metrics.tsv
#   results/gwasremoved_withinperm_peaks_summary.tsv
#
# Run from project root:
#   Rscript post-process-scripts/compute-gwasremoved-withinperm-peaks.R [n_perm]

suppressPackageStartupMessages(library(data.table))

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(project_root, "post-process-scripts", "helpers-50k-matrix.R"))

args   <- commandArgs(trailingOnly = TRUE)
n_perm <- if (length(args) >= 1L) as.integer(args[[1L]]) else 100L
chunk_sizes <- c(250000L, 500000L)

N_PCS   <- 4L
TOP_PCT <- 0.01
WIN_BP  <- 50000L

# ---- Metric functions ----

gini <- function(x) {
  x <- abs(x); x <- sort(x); n <- length(x)
  if (n == 0L || sum(x) == 0) return(NA_real_)
  (2 * sum(x * seq_len(n)) - (n + 1) * sum(x)) / (n * sum(x))
}

run_lengths <- function(abs_load, window_id) {
  thresh <- quantile(abs_load, 1 - TOP_PCT, na.rm = TRUE)
  in_top <- abs_load >= thresh
  parts  <- strsplit(window_id, "_", fixed = TRUE)
  chrom  <- sapply(parts, `[[`, 1)
  start  <- as.integer(sapply(parts, `[[`, 2))
  dt <- data.table(chrom = chrom, start = start, in_top = in_top)
  setorder(dt, chrom, start)
  dt[, adj_prev := chrom == shift(chrom) & start == shift(start) + WIN_BP & shift(in_top) == TRUE]
  dt[is.na(adj_prev), adj_prev := FALSE]
  dt[, run_id := cumsum(!(in_top & adj_prev))]
  top_runs <- dt[in_top == TRUE, .N, by = run_id]
  top_runs$N * WIN_BP / 1000
}

compute_metrics <- function(load_mat, window_ids) {
  rbindlist(lapply(seq_len(ncol(load_mat)), function(k) {
    v  <- load_mat[, k]
    rl <- run_lengths(abs(v), window_ids)
    data.table(
      pc            = paste0("PC", k),
      gini          = gini(v),
      mean_run_kb   = if (length(rl) == 0L) NA_real_ else mean(rl),
      median_run_kb = if (length(rl) == 0L) NA_real_ else median(rl)
    )
  }))
}

match_pcs <- function(obs_mat, perm_mat) {
  k      <- min(ncol(obs_mat), ncol(perm_mat))
  norm   <- function(m) apply(m[, seq_len(k)], 2, function(x) x / sqrt(sum(x^2)))
  obs_n  <- norm(obs_mat); perm_n <- norm(perm_mat)
  cos_mat <- abs(t(obs_n) %*% perm_n)
  used <- integer(0); ord <- integer(k); sgn <- numeric(k)
  for (i in seq_len(k)) {
    tmp <- cos_mat[i, ]; tmp[used] <- -1
    j <- which.max(tmp); used <- c(used, j); ord[i] <- j
    sgn[i] <- if (t(obs_n[, i]) %*% perm_n[, j] >= 0) 1 else -1
  }
  sweep(perm_mat[, ord, drop = FALSE], 2, sgn, `*`)
}

# ---- Build retained-window matrix ----

message("Loading full trait matrix...")
eff   <- fread(file.path(project_root, "results", "all-trait-50k-mean-effects.tsv"))
pvals <- fread(file.path(project_root, "results", "all-trait-50k-mean-pvals.tsv"))
prep  <- prepare_X_50k(eff, pvals)
X_full   <- prep$X        # all retained windows x traits (after standard filters)
coords_full <- prep$coords

rem <- fread(file.path(project_root, "results",
  "gwas_removed_distance_topload_50k_plus100kb", "removed_windows_gwas_plus100kb.tsv"))
removed_ids <- rem$window_id

keep_idx    <- which(!(rownames(X_full) %in% removed_ids))
X_ret       <- X_full[keep_idx, , drop = FALSE]
coords_ret  <- coords_full[keep_idx]
ret_ids     <- rownames(X_ret)
message(sprintf("Retained windows for permutation: %d", nrow(X_ret)))

# ---- Observed GWAS-removed loadings ----

obs_file <- file.path(project_root, "results",
  "gwas_removed_distance_topload_50k_plus100kb", "reduced_pca_window_loadings.tsv")
obs_ld   <- fread(obs_file)
id_col   <- if ("window_id" %in% names(obs_ld)) "window_id" else "region_label"
obs_ids  <- obs_ld[[id_col]]
pc_cols  <- paste0("PC", seq_len(N_PCS))
obs_mat  <- as.matrix(obs_ld[match(ret_ids, obs_ids), ..pc_cols])

obs_metrics <- compute_metrics(obs_mat, ret_ids)
obs_metrics[, `:=`(scheme = "observed", iter = 0L)]
message("Observed GWAS-removed metrics computed.")

# ---- Permutations ----

iter_rows <- list()
for (chunk_bp in chunk_sizes) {
  scheme <- paste0("chunk", chunk_bp %/% 1000L, "kb_withinperm")
  message("Running scheme: ", scheme)
  for (i in seq_len(n_perm)) {
    Xp      <- permute_within_chunks(X_ret, coords_ret, chunk_size = chunk_bp,
                                     seed = 70000L + chunk_bp + i)
    dec     <- decomp_from_matrix(Xp)
    perm_mat <- dec$U[, seq_len(N_PCS), drop = FALSE]  # window scores (= loadings)
    aligned  <- match_pcs(obs_mat, perm_mat)
    m        <- compute_metrics(aligned, ret_ids)
    m[, `:=`(scheme = scheme, iter = i)]
    iter_rows[[length(iter_rows) + 1L]] <- m
    if (i %% 20L == 0L || i == n_perm)
      message(sprintf("  %s: %d/%d", scheme, i, n_perm))
  }
}

iter_dt <- rbindlist(iter_rows, use.names = TRUE)
all_dt  <- rbindlist(list(obs_metrics, iter_dt), use.names = TRUE, fill = TRUE)
setorder(all_dt, scheme, iter, pc)

fwrite(all_dt,
  file.path(project_root, "results", "gwasremoved_withinperm_peaks_iter_metrics.tsv"),
  sep = "\t")
message("Per-iteration metrics written.")

# ---- Summary ----

metrics_long <- c("gini", "mean_run_kb", "median_run_kb")
sum_rows <- list()
for (sc in unique(iter_dt$scheme)) {
  for (pci in paste0("PC", seq_len(N_PCS))) {
    for (mm in metrics_long) {
      null_vals <- iter_dt[scheme == sc & pc == pci, get(mm)]
      obs_val   <- obs_metrics[pc == pci, get(mm)]
      if (!length(null_vals) || !is.finite(obs_val)) next
      mu  <- mean(null_vals, na.rm = TRUE)
      sdv <- sd(null_vals,   na.rm = TRUE)
      z   <- if (is.finite(sdv) && sdv > 0) (obs_val - mu) / sdv else NA_real_
      p   <- if (mm == "gini") {
        (1 + sum(null_vals >= obs_val, na.rm = TRUE)) / (length(null_vals) + 1)
      } else {
        (1 + sum(null_vals <= obs_val, na.rm = TRUE)) / (length(null_vals) + 1)
      }
      sum_rows[[length(sum_rows) + 1L]] <- data.table(
        scheme = sc, pc = pci, metric = mm,
        observed = obs_val, null_mean = mu, null_sd = sdv,
        z_score = z, empirical_p = p, n_null = sum(is.finite(null_vals))
      )
    }
  }
}

summary_dt <- rbindlist(sum_rows, use.names = TRUE)
setorder(summary_dt, scheme, pc, metric)
fwrite(summary_dt,
  file.path(project_root, "results", "gwasremoved_withinperm_peaks_summary.tsv"),
  sep = "\t")
message("Summary written.\n")
print(summary_dt)
