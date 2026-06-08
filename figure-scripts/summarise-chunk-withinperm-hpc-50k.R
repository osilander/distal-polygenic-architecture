#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

setDTthreads(1L)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
base_dir <- file.path(project_root, "results", "chunk-withinperm-nulls-50k", "hpc_base")
run_dir <- file.path(project_root, "results", "chunk-withinperm-nulls-50k", "hpc_runs")
out_dir <- file.path(project_root, "results", "chunk-withinperm-nulls-50k")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

obs_file <- file.path(base_dir, "observed_metrics_50k.tsv")
obs_scree_file <- file.path(base_dir, "observed_scree_50k.tsv")
metric_files <- list.files(file.path(run_dir, "metrics"), pattern = "_metrics\\.tsv$", full.names = TRUE)
scree_files <- list.files(file.path(run_dir, "scree"), pattern = "_scree\\.tsv\\.gz$", full.names = TRUE)

stopifnot(file.exists(obs_file), file.exists(obs_scree_file))
if (length(metric_files) == 0L) stop("No worker metric files found in ", file.path(run_dir, "metrics"))
if (length(scree_files) == 0L) stop("No worker scree files found in ", file.path(run_dir, "scree"))

obs <- fread(obs_file)
obs_scree <- fread(obs_scree_file)[, `:=`(scheme = "observed", iter = 0L)]
run_metrics <- rbindlist(lapply(metric_files, fread), use.names = TRUE, fill = TRUE)
run_scree <- rbindlist(lapply(scree_files, fread), use.names = TRUE, fill = TRUE)

# Observed metrics for same schema.
obs_rows <- list()
i <- 1L
for (k in sort(unique(run_metrics[metric %in% c("window_subspace_mean_cos", "trait_abs_cos_upper_corr", "cumvar_k"), k]))) {
  obs_rows[[i]] <- data.table(scheme = "observed", iter = 0L, metric = "window_subspace_mean_cos", k = k, value = 1); i <- i + 1L
  obs_rows[[i]] <- data.table(scheme = "observed", iter = 0L, metric = "trait_abs_cos_upper_corr", k = k, value = 1); i <- i + 1L
  obs_rows[[i]] <- data.table(scheme = "observed", iter = 0L, metric = "cumvar_k", k = k, value = obs[metric == "cumvar_k" & k == k, value][1]); i <- i + 1L
}
obs_rows[[i]] <- data.table(scheme = "observed", iter = 0L, metric = "pc1_var_explained", k = NA_integer_, value = obs[metric == "pc1_var_explained", value][1]); i <- i + 1L
obs_rows[[i]] <- data.table(scheme = "observed", iter = 0L, metric = "effective_rank", k = NA_integer_, value = obs[metric == "effective_rank", value][1]); i <- i + 1L
obs_metrics_full <- rbindlist(obs_rows, use.names = TRUE, fill = TRUE)

all_metrics <- rbindlist(list(obs_metrics_full, run_metrics), use.names = TRUE, fill = TRUE)
all_scree <- rbindlist(list(obs_scree, run_scree), use.names = TRUE, fill = TRUE)
fwrite(all_metrics, file.path(out_dir, "chunk_withinperm_metrics_runs_50k.tsv"), sep = "\t")
fwrite(all_scree, file.path(out_dir, "chunk_withinperm_scree_runs_50k.tsv.gz"), sep = "\t")

# Summary vs observed.
sum_rows <- list()
sid <- 1L
for (sc in unique(run_metrics$scheme)) {
  for (mm in unique(run_metrics$metric)) {
    sub <- run_metrics[scheme == sc & metric == mm]
    combos <- if (mm %in% c("window_subspace_mean_cos", "trait_abs_cos_upper_corr", "cumvar_k")) unique(sub[, .(k)]) else data.table(k = NA_integer_)
    for (r in seq_len(nrow(combos))) {
      kk <- combos$k[r]
      x <- if (is.na(kk)) sub$value else sub[k == kk, value]
      obs_val <- if (is.na(kk)) obs_metrics_full[metric == mm, value][1] else obs_metrics_full[metric == mm & k == kk, value][1]
      if (length(x) == 0L || !is.finite(obs_val)) next
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

# Mean scree and plot
scree_mean <- run_scree[, .(var_explained = mean(var_explained, na.rm = TRUE)), by = .(scheme, pc)]
scree_obs <- obs_scree[, .(scheme, pc, var_explained)]
scree_cmp <- rbindlist(list(scree_obs, scree_mean), use.names = TRUE, fill = TRUE)
scree_cmp[, pc_index := as.integer(sub("^PC", "", pc))]
fwrite(scree_cmp, file.path(out_dir, "chunk_withinperm_scree_observed_vs_null_means_50k.tsv"), sep = "\t")

p_scree <- ggplot(scree_cmp[pc_index <= 60], aes(x = pc_index, y = var_explained, color = scheme)) +
  geom_line(linewidth = 0.8) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank()) +
  labs(title = "Observed vs chunk-withinperm null scree (50k)", x = "PC", y = "Variance explained")
ggsave(file.path(out_dir, "chunk_withinperm_scree_observed_vs_null_means_50k.pdf"), p_scree, width = 9, height = 5)

message("Done. Wrote summaries to: ", out_dir)
