#!/usr/bin/env Rscript
# compute-withinperm-peaks-50k.R
#
# Compares observed PC1-4 window loading distributions against 250kb and 500kb
# within-block permutation nulls, using three metrics:
#
#   1. Gini coefficient of |loadings|        -- inequality / concentration
#   2. Mean run length of top-1% windows     -- spatial clustering
#   3. Median run length of top-1% windows   -- spatial clustering (robust)
#
# PCs are matched between observed and each permuted run by absolute cosine
# similarity (greedy), so ordering differences between SVDs don't confound results.
#
# Outputs:
#   results/withinperm_peaks_iter_metrics.tsv      -- per-iteration metrics
#   results/withinperm_peaks_summary.tsv           -- observed vs null summary
#
# Run from project root: Rscript figure-scripts/compute-withinperm-peaks-50k.R

suppressPackageStartupMessages(library(data.table))

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
vec_dir      <- file.path(project_root, "results", "permute-window-vectors")
out_dir      <- file.path(project_root, "results")

N_PCS   <- 4L
TOP_PCT <- 0.01   # top 1% for run-length metric
WIN_BP  <- 50000L # window size

# ---- Metric functions ----

gini <- function(x) {
  x <- abs(x); x <- sort(x); n <- length(x)
  if (n == 0L || sum(x) == 0) return(NA_real_)
  (2 * sum(x * seq_len(n)) - (n + 1) * sum(x)) / (n * sum(x))
}

run_lengths <- function(abs_load, window_id) {
  # Returns vector of run lengths (in kb) for top TOP_PCT windows
  thresh  <- quantile(abs_load, 1 - TOP_PCT, na.rm = TRUE)
  in_top  <- abs_load >= thresh
  parts   <- strsplit(window_id, "_", fixed = TRUE)
  chrom   <- sapply(parts, `[[`, 1)
  start   <- as.integer(sapply(parts, `[[`, 2))
  dt <- data.table(chrom = chrom, start = start, in_top = in_top)
  setorder(dt, chrom, start)
  dt[, adj_prev := chrom == shift(chrom) & start == shift(start) + WIN_BP & shift(in_top) == TRUE]
  dt[is.na(adj_prev), adj_prev := FALSE]
  dt[, run_id := cumsum(!(in_top & adj_prev))]
  top_runs <- dt[in_top == TRUE, .N, by = run_id]
  top_runs$N * WIN_BP / 1000   # in kb
}

compute_metrics <- function(load_mat, window_ids) {
  rbindlist(lapply(seq_len(ncol(load_mat)), function(k) {
    v   <- load_mat[, k]
    rl  <- run_lengths(abs(v), window_ids)
    data.table(
      pc           = paste0("PC", k),
      gini         = gini(v),
      mean_run_kb  = if (length(rl) == 0L) NA_real_ else mean(rl),
      median_run_kb = if (length(rl) == 0L) NA_real_ else median(rl)
    )
  }))
}

# ---- PC matching: greedy absolute cosine ----

match_pcs <- function(obs_mat, perm_mat) {
  # Returns column-reordered (and sign-corrected) perm_mat aligned to obs_mat
  k      <- min(ncol(obs_mat), ncol(perm_mat))
  obs_n  <- apply(obs_mat[, seq_len(k)], 2, function(x) x / sqrt(sum(x^2)))
  perm_n <- apply(perm_mat[, seq_len(k)], 2, function(x) x / sqrt(sum(x^2)))
  cos_mat <- abs(t(obs_n) %*% perm_n)   # k x k absolute cosine

  used_perm <- integer(0); order_perm <- integer(k); signs <- numeric(k)
  for (i in seq_len(k)) {
    cos_tmp <- cos_mat[i, ]
    cos_tmp[used_perm] <- -1
    j <- which.max(cos_tmp)
    used_perm <- c(used_perm, j)
    order_perm[i] <- j
    raw_cos <- t(obs_n[, i]) %*% perm_n[, j]
    signs[i] <- if (raw_cos >= 0) 1 else -1
  }
  out <- perm_mat[, order_perm, drop = FALSE]
  sweep(out, 2, signs, `*`)
}

# ---- Load observed loadings ----

cat("Loading observed loadings ...\n")
obs_dt   <- fread(file.path(project_root, "results", "pca_loadings_50k.tsv"),
                  showProgress = FALSE)
pc_cols  <- paste0("PC", seq_len(N_PCS))
obs_mat  <- as.matrix(obs_dt[, ..pc_cols])
obs_ids  <- obs_dt$region_label   # chrom_start_end

obs_metrics <- compute_metrics(obs_mat, obs_ids)
obs_metrics[, `:=`(scheme = "observed", iter = 0L)]
cat("Observed metrics computed.\n")

# ---- Process permuted files ----

perm_files <- list.files(vec_dir, pattern = "_window_vectors_pc1_4\\.tsv\\.gz$",
                         full.names = TRUE)
cat(sprintf("Found %d permuted files\n", length(perm_files)))

iter_rows <- list()
for (f in perm_files) {
  bn     <- basename(f)
  scheme <- sub("_iter.*", "", bn)
  iter   <- as.integer(sub(".*_iter([0-9]+)_.*", "\\1", bn))

  perm_dt  <- fread(f, showProgress = FALSE)
  # align window order to observed
  perm_dt  <- perm_dt[match(obs_ids, window_id)]
  if (anyNA(perm_dt$window_id)) {
    warning("Window mismatch in ", bn, " — skipping"); next
  }
  perm_mat <- as.matrix(perm_dt[, ..pc_cols])

  aligned  <- match_pcs(obs_mat, perm_mat)
  m        <- compute_metrics(aligned, obs_ids)
  m[, `:=`(scheme = scheme, iter = iter)]
  iter_rows[[length(iter_rows) + 1L]] <- m
  if (iter %% 20L == 0L) cat(sprintf("  %s iter %d done\n", scheme, iter))
}

iter_dt <- rbindlist(iter_rows, use.names = TRUE)
all_dt  <- rbindlist(list(obs_metrics, iter_dt), use.names = TRUE, fill = TRUE)
setorder(all_dt, scheme, iter, pc)

fwrite(all_dt, file.path(out_dir, "withinperm_peaks_iter_metrics.tsv"), sep = "\t")
cat("Per-iteration metrics written.\n")

# ---- Summary: observed vs null ----

metrics_long <- c("gini", "mean_run_kb", "median_run_kb")

sum_rows <- list()
for (sc in unique(iter_dt$scheme)) {
  for (pci in paste0("PC", seq_len(N_PCS))) {
    for (mm in metrics_long) {
      null_vals <- iter_dt[scheme == sc & pc == pci, get(mm)]
      obs_val   <- obs_metrics[pc == pci, get(mm)]
      if (length(null_vals) == 0L || !is.finite(obs_val)) next
      mu  <- mean(null_vals, na.rm = TRUE)
      sdv <- sd(null_vals, na.rm = TRUE)
      z   <- if (is.finite(sdv) && sdv > 0) (obs_val - mu) / sdv else NA_real_
      # run metrics: observed < null (shorter = more focused), so use lower tail
      # gini: observed > null (more concentrated), so use upper tail
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
fwrite(summary_dt, file.path(out_dir, "withinperm_peaks_summary.tsv"), sep = "\t")
cat("Summary written.\n\n")
print(summary_dt)
