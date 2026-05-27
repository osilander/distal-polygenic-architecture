#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

setDTthreads(1L)
args <- commandArgs(trailingOnly = TRUE)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out_dir <- file.path(project_root, "results", "chunk-withinperm-nulls-50k")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
base_file <- file.path(project_root, "results", "chunk-withinperm-nulls-50k", "hpc_base", "chunk_withinperm_base_50k.rds")
source(file.path(project_root, "post-process-scripts", "helpers-50k-matrix.R"))

n_perm <- if (length(args) >= 1L) as.integer(args[[1]]) else 100L
if (!is.finite(n_perm) || n_perm < 1L) stop("n_perm must be a positive integer")
k_vals <- c(5L, 8L, 12L)
chunk_sizes <- c(250000L, 500000L)
stopifnot(file.exists(base_file))
message("Loading prepared chunk-withinperm base object...")
base <- readRDS(base_file)
X <- base$X
coords <- base$coords
obs <- list(U = base$obs_U_top, V = base$obs_V_top, d = base$obs_d_top, var = base$obs_var)
obs_upper <- lapply(k_vals, function(k) trait_abs_cos_upper(obs$V, obs$d, k))
names(obs_upper) <- as.character(k_vals)

obs_metric_rows <- list()
idx <- 1L
for (k in k_vals) {
  kk <- min(k, length(obs$var))
  obs_metric_rows[[idx]] <- data.table(scheme = "observed", iter = 0L, metric = "window_subspace_mean_cos", k = kk, value = 1); idx <- idx + 1L
  obs_metric_rows[[idx]] <- data.table(scheme = "observed", iter = 0L, metric = "trait_abs_cos_upper_corr", k = kk, value = 1); idx <- idx + 1L
  obs_metric_rows[[idx]] <- data.table(scheme = "observed", iter = 0L, metric = "cumvar_k", k = kk, value = sum(obs$var[seq_len(kk)])); idx <- idx + 1L
}
obs_metric_rows[[idx]] <- data.table(scheme = "observed", iter = 0L, metric = "pc1_var_explained", k = NA_integer_, value = obs$var[1]); idx <- idx + 1L
obs_metric_rows[[idx]] <- data.table(scheme = "observed", iter = 0L, metric = "effective_rank", k = NA_integer_, value = base$obs_effective_rank); idx <- idx + 1L
obs_metrics <- rbindlist(obs_metric_rows, use.names = TRUE, fill = TRUE)

run_rows <- list()
ridx <- 1L
scree_rows <- list()
sidx <- 1L

for (chunk_bp in chunk_sizes) {
  scheme <- paste0("chunk", chunk_bp %/% 1000L, "kb_withinperm")
  message("Running scheme: ", scheme)
  grp <- factor(paste(gsub("^chr", "", as.character(coords$chrom)), as.integer(coords$start) %/% chunk_bp, sep = ":"))
  n_chunks <- uniqueN(grp)
  fwrite(
    data.table(scheme = scheme, chunk_bp = chunk_bp, n_chunks = n_chunks, n_perm = n_perm),
    file.path(out_dir, paste0(scheme, "_runinfo_50k.tsv")),
    sep = "\t"
  )

  for (i in seq_len(n_perm)) {
    Xr <- permute_within_chunks(X, coords, chunk_size = chunk_bp, seed = 70000L + chunk_bp + i)
    dec <- decomp_from_matrix(Xr)

    # scree for mean curves
    scree_rows[[sidx]] <- data.table(
      scheme = scheme,
      iter = i,
      pc = paste0("PC", seq_along(dec$var)),
      var_explained = dec$var
    )
    sidx <- sidx + 1L

    for (k in k_vals) {
      kk <- min(k, length(dec$var), ncol(obs$U), ncol(dec$U))
      sub_cos <- mean_cos_principal_angles(obs$U, dec$U, kk)
      run_rows[[ridx]] <- data.table(scheme = scheme, iter = i, metric = "window_subspace_mean_cos", k = kk, value = sub_cos); ridx <- ridx + 1L

      tr_corr <- suppressWarnings(cor(obs_upper[[as.character(kk)]], trait_abs_cos_upper(dec$V, dec$d, kk), method = "pearson"))
      run_rows[[ridx]] <- data.table(scheme = scheme, iter = i, metric = "trait_abs_cos_upper_corr", k = kk, value = tr_corr); ridx <- ridx + 1L

      run_rows[[ridx]] <- data.table(scheme = scheme, iter = i, metric = "cumvar_k", k = kk, value = sum(dec$var[seq_len(kk)])); ridx <- ridx + 1L
    }

    run_rows[[ridx]] <- data.table(scheme = scheme, iter = i, metric = "pc1_var_explained", k = NA_integer_, value = dec$var[1]); ridx <- ridx + 1L
    run_rows[[ridx]] <- data.table(scheme = scheme, iter = i, metric = "effective_rank", k = NA_integer_, value = effective_rank(dec$d)); ridx <- ridx + 1L

    if (i %% 5L == 0L || i == n_perm) message("  ", scheme, ": ", i, "/", n_perm)
  }
}

run_metrics <- rbindlist(run_rows, use.names = TRUE, fill = TRUE)
scree_runs <- rbindlist(scree_rows, use.names = TRUE, fill = TRUE)
all_metrics <- rbindlist(list(obs_metrics, run_metrics), use.names = TRUE, fill = TRUE)

fwrite(all_metrics, file.path(out_dir, "chunk_withinperm_metrics_runs_50k.tsv"), sep = "\t")
fwrite(scree_runs, file.path(out_dir, "chunk_withinperm_scree_runs_50k.tsv.gz"), sep = "\t")

# Summaries vs observed
sum_rows <- list()
sid <- 1L
for (sc in unique(run_metrics$scheme)) {
  for (mm in unique(run_metrics$metric)) {
    sub <- run_metrics[scheme == sc & metric == mm]
    combos <- if (mm %in% c("window_subspace_mean_cos", "trait_abs_cos_upper_corr", "cumvar_k")) unique(sub[, .(k)]) else data.table(k = NA_integer_)
    for (i in seq_len(nrow(combos))) {
      kk <- combos$k[i]
      x <- if (is.na(kk)) sub$value else sub[k == kk, value]
      obs_val <- if (is.na(kk)) obs_metrics[metric == mm, value][1] else obs_metrics[metric == mm & k == kk, value][1]
      if (!length(x) || !is.finite(obs_val)) next
      mu <- mean(x, na.rm = TRUE)
      sdv <- sd(x, na.rm = TRUE)
      z <- if (is.finite(sdv) && sdv > 0) (obs_val - mu) / sdv else NA_real_
      p_emp <- if (mm == "effective_rank") {
        (1 + sum(x <= obs_val, na.rm = TRUE)) / (length(x) + 1)
      } else {
        (1 + sum(x >= obs_val, na.rm = TRUE)) / (length(x) + 1)
      }
      sum_rows[[sid]] <- data.table(
        scheme = sc, metric = mm, k = kk,
        observed_value = obs_val, null_mean = mu, null_sd = sdv, z_score = z, empirical_p = p_emp, n_null = length(x)
      )
      sid <- sid + 1L
    }
  }
}
summary_dt <- rbindlist(sum_rows, use.names = TRUE, fill = TRUE)
fwrite(summary_dt, file.path(out_dir, "chunk_withinperm_metrics_summary_50k.tsv"), sep = "\t")

# Scree means and plot
scree_mean <- scree_runs[, .(var_explained = mean(var_explained, na.rm = TRUE)), by = .(scheme, pc)]
scree_obs <- data.table(scheme = "observed", pc = paste0("PC", seq_along(obs$var)), var_explained = obs$var)
scree_cmp <- rbindlist(list(scree_obs, scree_mean), use.names = TRUE, fill = TRUE)
scree_cmp[, pc_index := as.integer(sub("^PC", "", pc))]
fwrite(scree_cmp, file.path(out_dir, "chunk_withinperm_scree_observed_vs_null_means_50k.tsv"), sep = "\t")

p <- ggplot(scree_cmp[pc_index <= 60], aes(x = pc_index, y = var_explained, color = scheme)) +
  geom_line(linewidth = 0.8) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank()) +
  labs(title = "Observed vs chunk-withinperm null scree (50k)", x = "PC", y = "Variance explained")
ggsave(file.path(out_dir, "chunk_withinperm_scree_observed_vs_null_means_50k.pdf"), p, width = 9, height = 5)

message("Done. Wrote outputs to: ", out_dir)
