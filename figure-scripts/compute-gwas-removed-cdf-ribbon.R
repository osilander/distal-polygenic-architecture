#!/usr/bin/env Rscript
# Produces: results/gwas_removed_distance_cdf_ribbon_50kb_bins.tsv
#           results/gwas_removed_distance_cdf_ribbon_50kb_bins_pvalues.tsv
# Pre-computed inputs:
#   results/gwas_removed_distance_topload_50k_plus100kb/ (from compute-gwas-removed-distances.R 100)
#   results/gwas_removed_distance_topload_50k_plus150kb/ (from compute-gwas-removed-distances.R 150)
# Run: Rscript figure-scripts/compute-gwas-removed-cdf-ribbon.R

suppressPackageStartupMessages({
  library(data.table)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

dirs <- list(
  `+100kb` = file.path("results", "gwas_removed_distance_topload_50k_plus100kb"),
  `+150kb` = file.path("results", "gwas_removed_distance_topload_50k_plus150kb")
)
bin_bp    <- 50000L
out_stats <- file.path("results", "gwas_removed_distance_cdf_ribbon_50kb_bins.tsv")
out_pvals <- file.path("results", "gwas_removed_distance_cdf_ribbon_50kb_bins_pvalues.tsv")
dir.create("results", recursive = TRUE, showWarnings = FALSE)

build_curves <- function(dir_path, flank_label) {
  top_file <- file.path(dir_path, "top_windows_nearest_removed_distances.tsv")
  rnd_file <- file.path(dir_path, "random_chrmatched_nearest_removed_distances.tsv.gz")
  if (!file.exists(top_file) || !file.exists(rnd_file))
    stop("Missing files in ", dir_path)

  top <- fread(top_file)[is.finite(nearest_removed_distance_bp),
                         .(pc, top_frac, nearest_removed_distance_bp)]
  rnd <- tryCatch(
    fread(rnd_file),
    error = function(e) fread(cmd = paste("gzip -dc", shQuote(rnd_file)))
  )
  rnd <- rnd[is.finite(nearest_removed_distance_bp),
             .(pc, top_frac, random_rep, nearest_removed_distance_bp)]

  max_bp <- max(c(top$nearest_removed_distance_bp, rnd$nearest_removed_distance_bp), na.rm = TRUE)
  bins   <- seq.int(0L, as.integer(ceiling(max_bp / bin_bp) * bin_bp), by = bin_bp)

  obs <- top[, {
    x <- nearest_removed_distance_bp
    data.table(dist_bp = bins,
               cdf = vapply(bins, function(b) mean(x <= b), numeric(1)))
  }, by = .(pc, top_frac)]

  rnd_curve <- rnd[, {
    x <- nearest_removed_distance_bp
    data.table(dist_bp = bins,
               cdf = vapply(bins, function(b) mean(x <= b), numeric(1)))
  }, by = .(pc, top_frac, random_rep)]

  rnd_sum <- rnd_curve[, .(
    sim_med = median(cdf, na.rm = TRUE),
    sim_lo  = quantile(cdf, 0.025, na.rm = TRUE),
    sim_hi  = quantile(cdf, 0.975, na.rm = TRUE)
  ), by = .(pc, top_frac, dist_bp)]

  obs[,     `:=`(flank = flank_label, dist_kb = dist_bp / 1000)]
  rnd_sum[, `:=`(flank = flank_label, dist_kb = dist_bp / 1000)]
  list(obs = obs, rnd = rnd_sum)
}

all_obs <- list(); all_rnd <- list(); ii <- 1L
for (nm in names(dirs)) {
  message("Building ribbon CDF for ", nm, "...")
  cur <- build_curves(dirs[[nm]], nm)
  all_obs[[ii]] <- cur$obs
  all_rnd[[ii]] <- cur$rnd
  ii <- ii + 1L
}

obs_dt <- rbindlist(all_obs, use.names = TRUE, fill = TRUE)
rnd_dt <- rbindlist(all_rnd, use.names = TRUE, fill = TRUE)
obs_dt <- obs_dt[abs(top_frac - 0.01) < 1e-12]
rnd_dt <- rnd_dt[abs(top_frac - 0.01) < 1e-12]

plot_dt <- merge(
  rnd_dt,
  obs_dt[, .(pc, top_frac, flank, dist_bp, obs_cdf = cdf)],
  by = c("pc", "top_frac", "flank", "dist_bp"),
  all.x = TRUE
)
fwrite(plot_dt, out_stats, sep = "\t")
message("Wrote: ", out_stats)

# ---- Empirical p-values ----
# Direct permutation test: observed stat = area(CDF_obs, mean_null_CDF);
# null stat for each rep r = area(CDF_rand_r, mean_null_CDF).
# p = (1 + #{null_r >= obs_stat}) / (1 + n_reps). Min p = 1/2001 ~ 0.0005.
calc_emp_p <- function(rnd_raw, obs_raw, pcv, flankv, topv) {
  ro <- rnd_raw[pc == pcv & flank == flankv & top_frac == topv]
  oo <- obs_raw[pc == pcv & flank == flankv & top_frac == topv]
  if (nrow(ro) == 0L || nrow(oo) == 0L)
    return(data.table(pc = pcv, flank = flankv, top_frac = topv,
                      obs_area_stat = NA_real_, null_area_mean = NA_real_,
                      null_area_sd = NA_real_, p_emp = NA_real_, n_null_reps = 0L))

  bins_p <- unique(c(0L, as.integer(sort(unique(
    c(ro$nearest_removed_distance_bp, oo$nearest_removed_distance_bp)
  )))))
  bins_p <- bins_p[is.finite(bins_p) & bins_p >= 0]
  if (length(bins_p) < 2L) bins_p <- c(0L, 50000L)

  rep_cdf <- ro[, .(
    dist_bp = bins_p,
    cdf     = vapply(bins_p, function(b) mean(nearest_removed_distance_bp <= b), numeric(1))
  ), by = random_rep]

  # Mean null CDF across all reps (reference)
  mean_null <- rep_cdf[, .(cdf_null = mean(cdf, na.rm = TRUE)), by = dist_bp]

  obs_cdf <- data.table(
    dist_bp = bins_p,
    cdf_obs = vapply(bins_p, function(b) mean(oo$nearest_removed_distance_bp <= b), numeric(1))
  )

  # Observed statistic: area between CDF_obs and mean null CDF
  obs_stat <- merge(obs_cdf, mean_null, by = "dist_bp")[
    , sum(abs(cdf_obs - cdf_null), na.rm = TRUE)] * bin_bp

  # Null stats: area between each rep's CDF and mean null CDF
  null_stats <- merge(rep_cdf, mean_null, by = "dist_bp", all.x = TRUE)[
    , .(stat = sum(abs(cdf - cdf_null), na.rm = TRUE) * bin_bp), by = random_rep]

  p_emp <- (1 + sum(null_stats$stat >= obs_stat, na.rm = TRUE)) / (1 + nrow(null_stats))
  data.table(pc = pcv, flank = flankv, top_frac = topv,
             obs_area_stat  = obs_stat,
             null_area_mean = mean(null_stats$stat, na.rm = TRUE),
             null_area_sd   = sd(null_stats$stat,   na.rm = TRUE),
             p_emp = p_emp, n_null_reps = nrow(null_stats))
}

rnd_raw <- rbindlist(lapply(names(dirs), function(nm) {
  f <- file.path(dirs[[nm]], "random_chrmatched_nearest_removed_distances.tsv.gz")
  x <- tryCatch(fread(f), error = function(e) fread(cmd = paste("gzip -dc", shQuote(f))))
  x[, flank := nm]; x
}), use.names = TRUE, fill = TRUE)
obs_raw <- rbindlist(lapply(names(dirs), function(nm) {
  x <- fread(file.path(dirs[[nm]], "top_windows_nearest_removed_distances.tsv"))
  x[, flank := nm]; x
}), use.names = TRUE, fill = TRUE)
obs_raw <- obs_raw[abs(top_frac - 0.01) < 1e-12]
rnd_raw <- rnd_raw[abs(top_frac - 0.01) < 1e-12]

p_rows <- rbindlist(lapply(sort(unique(obs_raw$pc)), function(pcv) {
  rbindlist(lapply(sort(unique(obs_raw$flank)), function(fv) {
    calc_emp_p(rnd_raw, obs_raw, pcv, fv, 0.01)
  }), use.names = TRUE, fill = TRUE)
}), use.names = TRUE, fill = TRUE)

fwrite(p_rows, out_pvals, sep = "\t")
message("Wrote: ", out_pvals)
message("Done.")
