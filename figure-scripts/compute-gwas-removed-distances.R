#!/usr/bin/env Rscript
# Produces: results/gwas_removed_distance_topload_50k_plus{flank_kb}kb/
#   reduced_pca_window_loadings.tsv, reduced_pca_trait_scores.tsv,
#   reduced_pca_scree.tsv, removed_windows_gwas_plus{flank_kb}kb.tsv,
#   top_windows_nearest_removed_distances.tsv,
#   random_chrmatched_nearest_removed_distances.tsv.gz,
#   summary_top_vs_random_distances.tsv
# Pre-computed inputs: results/all-trait-50k-mean-{effects,pvals}.tsv
# Run: Rscript post-process-scripts/compute-gwas-removed-distances.R [flank_kb]
#   default flank_kb = 100; run again with 150 for the +150kb variant

suppressPackageStartupMessages({
  library(data.table)
  library(irlba)
  library(GenomicRanges)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

gws_log10_thresh <- -log10(5e-8)
args     <- commandArgs(trailingOnly = TRUE)
flank_kb <- if (length(args) >= 1L) as.integer(args[[1]]) else 100L
if (!is.finite(flank_kb) || flank_kb < 0L) stop("flank_kb must be a non-negative integer.")
flank_bp   <- as.integer(flank_kb * 1000L)
top_fracs  <- c(0.01, 0.02, 0.025)
pcs        <- paste0("PC", 1:4)
n_random_reps <- 2000L
anchor_map <- c(PC1 = "BMI_03", PC2 = "StandHgt", PC3 = "Neurotic_30", PC4 = "DiaBP_A59")
set.seed(20260317)

out_dir      <- file.path("results", paste0("gwas_removed_distance_topload_50k_plus", flank_kb, "kb"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_removed  <- file.path(out_dir, paste0("removed_windows_gwas_plus", flank_kb, "kb.tsv"))
out_loadings <- file.path(out_dir, "reduced_pca_window_loadings.tsv")
out_traits   <- file.path(out_dir, "reduced_pca_trait_scores.tsv")
out_scree    <- file.path(out_dir, "reduced_pca_scree.tsv")
out_top_dist <- file.path(out_dir, "top_windows_nearest_removed_distances.tsv")
out_rand_dist <- file.path(out_dir, "random_chrmatched_nearest_removed_distances.tsv.gz")
out_summary  <- file.path(out_dir, "summary_top_vs_random_distances.tsv")

# ---- Build scaled trait matrix ----
message("GWAS distance analysis: rebuilding scaled trait matrix (~2-3 min)...")
source(file.path("post-process-scripts", "helpers-50k-matrix.R"))
all_betas <- fread(file.path("results", "all-trait-50k-mean-effects.tsv"))
all_pvals <- fread(file.path("results", "all-trait-50k-mean-pvals.tsv"))
prep <- prepare_X_50k(all_betas, all_pvals)
traits_t_scaled <- t(prep$X)  # traits x windows

# ---- Identify GWS windows ----
message("GWAS distance analysis: deriving removed loci (p<5e-8 +/-", flank_kb, "kb)...")
pv <- copy(all_pvals)
pv[, `:=`(chrom = as.integer(chrom), start = as.integer(start), end = as.integer(end))]
pv <- pv[is.finite(chrom) & is.finite(start) & is.finite(end)]
pv[, window_id := paste(chrom, start, end, sep = "_")]
trait_cols <- setdiff(names(pv), c("chrom", "start", "end", "window_id"))
for (cc in trait_cols) {
  if (!is.numeric(pv[[cc]])) suppressWarnings(set(pv, j = cc, value = as.numeric(pv[[cc]])))
}
pmat <- as.matrix(pv[, ..trait_cols]); storage.mode(pmat) <- "numeric"
pv[, gws_any := apply(pmat > gws_log10_thresh, 1, function(v) any(v, na.rm = TRUE))]

gws <- pv[gws_any == TRUE, .(chrom, gws_start = start, gws_end = end)]
if (nrow(gws) == 0L) stop("No GWAS-significant windows found; cannot continue.")
gws[, `:=`(flank_start = pmax(0L, gws_start - flank_bp), flank_end = gws_end + flank_bp)]
gws_peaks_dt <- unique(pv[gws_any == TRUE, .(window_id, chrom, start, end)])
setorder(gws_peaks_dt, chrom, start, end)
cand <- merge(
  pv[, .(window_id, chrom, start, end, gws_any)],
  gws[, .(chrom, flank_start, flank_end)],
  by = "chrom", allow.cartesian = TRUE
)
removed_ids <- unique(cand[start < flank_end & end > flank_start, window_id])
removed_dt  <- pv[window_id %in% removed_ids, .(window_id, chrom, start, end, gws_any)]
setorder(removed_dt, chrom, start, end)
fwrite(removed_dt, out_removed, sep = "\t")
message("Removed windows: ", nrow(removed_dt))

# ---- Reduced PCA ----
message("GWAS distance analysis: rerunning PCA on reduced matrix...")
X         <- as.matrix(traits_t_scaled)  # traits x windows
keep_cols <- setdiff(colnames(X), removed_ids)
if (length(keep_cols) < 100L) stop("Too few windows remain after removal.")
X_red <- X[, keep_cols, drop = FALSE]
pca_red <- prcomp_irlba(X_red, n = 80, center = FALSE, scale. = FALSE)
rownames(pca_red$rotation) <- colnames(X_red)
rownames(pca_red$x)        <- rownames(X_red)

trait_scores <- as.data.table(pca_red$x)
trait_scores[, trait_abbrev := rownames(pca_red$x)]
for (pc in names(anchor_map)) {
  if (!(pc %in% names(trait_scores))) next
  tr  <- anchor_map[[pc]]
  idx <- which(trait_scores$trait_abbrev == tr)
  if (length(idx) > 0L) {
    aval <- trait_scores[[pc]][idx[1]]
    flip <- if (is.finite(aval) && aval < 0) -1 else 1
    trait_scores[[pc]]          <- trait_scores[[pc]] * flip
    pca_red$rotation[, pc] <- pca_red$rotation[, pc] * flip
  }
}
fwrite(trait_scores, out_traits, sep = "\t")

load_dt <- as.data.table(pca_red$rotation)
load_dt[, window_id := rownames(pca_red$rotation)]
setcolorder(load_dt, c("window_id", setdiff(names(load_dt), "window_id")))
load_dt[, c("chrom", "start", "end") := tstrsplit(window_id, "_", fixed = TRUE)]
load_dt[, `:=`(chrom = as.integer(chrom), start = as.integer(start), end = as.integer(end))]
load_dt <- load_dt[is.finite(chrom) & is.finite(start) & is.finite(end)]
setorder(load_dt, chrom, start, end)
fwrite(load_dt, out_loadings, sep = "\t")

eig   <- pca_red$sdev^2
scree <- data.table(
  PC = seq_along(eig),
  var_explained_trunc     = eig / sum(eig),
  var_explained_trunc_pct = 100 * eig / sum(eig),
  cumulative_trunc        = cumsum(eig) / sum(eig)
)
fwrite(scree, out_scree, sep = "\t")

# ---- Nearest-GWAS-peak distances ----
nearest_gws_peak <- function(query_dt, gws_peak_windows_dt) {
  if (nrow(query_dt) == 0L) return(copy(query_dt))
  out  <- copy(query_dt)
  out[, `:=`(nearest_gws_peak_window = NA_character_, nearest_gws_peak_distance_bp = NA_real_)]
  gr_q <- GRanges(seqnames = as.character(out$chrom),
                  ranges = IRanges(start = out$start + 1L, end = out$end))
  gr_p <- GRanges(seqnames = as.character(gws_peak_windows_dt$chrom),
                  ranges = IRanges(start = gws_peak_windows_dt$start + 1L, end = gws_peak_windows_dt$end))
  hits <- distanceToNearest(gr_q, gr_p, ignore.strand = TRUE)
  if (length(hits) > 0L) {
    qi <- queryHits(hits); pi <- subjectHits(hits)
    out$nearest_gws_peak_window[qi]       <- gws_peak_windows_dt$window_id[pi]
    out$nearest_gws_peak_distance_bp[qi]  <- mcols(hits)$distance
  }
  out
}

sample_chr_matched <- function(pool_dt, chr_count_dt) {
  rows <- vector("list", nrow(chr_count_dt))
  for (i in seq_len(nrow(chr_count_dt))) {
    chr_i  <- chr_count_dt$chrom[i]
    n_i    <- chr_count_dt$n[i]
    cand_i <- pool_dt[chrom == chr_i]
    if (n_i > nrow(cand_i))
      stop("Not enough candidate windows for chromosome ", chr_i, ": need ", n_i, ", have ", nrow(cand_i))
    rows[[i]] <- cand_i[sample.int(nrow(cand_i), n_i, replace = FALSE)]
  }
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

message("GWAS distance analysis: precomputing nearest-GWAS-peak distances for all retained windows...")
dist_lookup <- nearest_gws_peak(
  load_dt[, .(window_id, chrom, start, end)],
  gws_peaks_dt
)[, .(window_id, nearest_gws_peak_window, nearest_gws_peak_distance_bp)]

top_rows    <- list()
rand_rows   <- list()
summary_rows <- list()
ii <- jj <- kk <- 1L

for (pc in pcs) {
  if (!(pc %in% names(load_dt))) stop("Missing ", pc, " in reduced loadings.")
  z_all <- load_dt[, .(window_id, chrom, start, end, loading = get(pc))]
  z_all <- z_all[is.finite(loading)]
  z_all <- merge(z_all, dist_lookup, by = "window_id", all.x = TRUE)
  z_all[, abs_loading := abs(loading)]
  message("Processing ", pc, "...")

  for (tf in top_fracs) {
    message("  top fraction ", tf, " with ", n_random_reps, " chr-matched random reps")
    cutoff <- as.numeric(stats::quantile(z_all$abs_loading, probs = 1 - tf, na.rm = TRUE))
    top_dt <- z_all[abs_loading >= cutoff]
    top_dt[, `:=`(pc = pc, top_frac = tf, set = "top")]
    top_rows[[ii]] <- copy(top_dt)
    ii <- ii + 1L

    chr_counts <- top_dt[, .(n = .N), by = chrom]
    pool_dt    <- z_all[!(window_id %in% top_dt$window_id)]
    if (nrow(pool_dt) == 0L) stop("Random pool empty for ", pc, " top ", tf)

    rep_list <- vector("list", n_random_reps)
    for (rr in seq_len(n_random_reps)) {
      rnd <- sample_chr_matched(pool_dt, chr_counts)
      rnd[, `:=`(pc = pc, top_frac = tf, set = "random", random_rep = rr)]
      rep_list[[rr]] <- rnd
    }
    rand_dt <- rbindlist(rep_list, use.names = TRUE, fill = TRUE)
    rand_rows[[jj]] <- rand_dt
    jj <- jj + 1L

    obs <- top_dt[is.finite(nearest_gws_peak_distance_bp), nearest_gws_peak_distance_bp]
    rand_rep_stats <- rand_dt[is.finite(nearest_gws_peak_distance_bp), .(
      median_bp = median(nearest_gws_peak_distance_bp),
      mean_bp   = mean(nearest_gws_peak_distance_bp)
    ), by = random_rep]
    obs_median <- if (length(obs) > 0L) median(obs) else NA_real_
    obs_mean   <- if (length(obs) > 0L) mean(obs)   else NA_real_

    summary_rows[[kk]] <- data.table(
      pc = pc, top_frac = tf, n_top = nrow(top_dt),
      obs_median_bp          = obs_median,
      obs_mean_bp            = obs_mean,
      rand_median_mean_bp    = mean(rand_rep_stats$median_bp, na.rm = TRUE),
      rand_median_sd_bp      = sd(rand_rep_stats$median_bp, na.rm = TRUE),
      rand_mean_mean_bp      = mean(rand_rep_stats$mean_bp, na.rm = TRUE),
      rand_mean_sd_bp        = sd(rand_rep_stats$mean_bp, na.rm = TRUE),
      emp_p_median_low = (1 + sum(rand_rep_stats$median_bp <= obs_median, na.rm = TRUE)) /
                         (1 + nrow(rand_rep_stats)),
      emp_p_mean_low   = (1 + sum(rand_rep_stats$mean_bp   <= obs_mean,   na.rm = TRUE)) /
                         (1 + nrow(rand_rep_stats))
    )
    kk <- kk + 1L
  }
}

top_all  <- rbindlist(top_rows,     use.names = TRUE, fill = TRUE)
rand_all <- rbindlist(rand_rows,    use.names = TRUE, fill = TRUE)
sum_all  <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)

# Aliases for downstream compatibility with CDF-ribbon and enrichment scripts.
# NOTE: "nearest_removed_*" is shorthand for "nearest original GWAS-significant peak window".
# The actual quantity is distance to the nearest window that had any GWS hit (gws_any == TRUE),
# computed by nearest_gws_peak() above — NOT distance to any removed window in the flanked set.
top_all[ , `:=`(nearest_removed_window      = nearest_gws_peak_window,
                nearest_removed_distance_bp  = nearest_gws_peak_distance_bp)]
rand_all[, `:=`(nearest_removed_window      = nearest_gws_peak_window,
                nearest_removed_distance_bp  = nearest_gws_peak_distance_bp)]

fwrite(top_all,  out_top_dist,  sep = "\t")
fwrite(rand_all, out_rand_dist, sep = "\t")
fwrite(sum_all,  out_summary,   sep = "\t")

message("Done. Outputs in: ", out_dir)
