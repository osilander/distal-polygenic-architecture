#!/usr/bin/env Rscript
# Produces: results/observed-nearest-patterns/
#   observed_topfracs_nearestgene_bin_enrichment.tsv
#   observed_topfracs_nearestgene_null_replicate_stats.tsv.gz
#   observed_topfracs_nearestgene_z_summary.tsv
#   (and analogous nearestgwaspeak_* files when analysis_target includes "gwas")
# Pre-computed inputs: results/pca_loadings_50k.tsv,
#                      data/gencode.v19.genes.protein_coding.rds,
#                      results/all-trait-50k-mean-pvals.tsv (for gwas target only)
# Run: Rscript post-process-scripts/compute-nearestgene-patterns.R [gene|gwas|all] [n_reps]
#   default: gene 1000

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

set.seed(20260328)

args <- commandArgs(trailingOnly = TRUE)
analysis_target <- if (length(args) >= 1L) tolower(args[[1]]) else "gene"
if (!analysis_target %in% c("all", "gene", "gwas"))
  stop("First argument must be one of: all, gene, gwas")
n_reps <- if (length(args) >= 2L) as.integer(args[[2]]) else 1000L
if (!is.finite(n_reps) || n_reps < 1L) stop("n_reps must be a positive integer")

pcs            <- c("PC1", "PC2", "PC3", "PC4")
top_fracs      <- c(0.01, 0.02, 0.05)
gwas_log10_thresh <- 7.30103

bin_breaks <- c(-Inf, 0, 50000, 100000, 250000, 500000, 750000, 1000000, Inf)
bin_labels <- c("0", "(0,50kb]", "(50,100kb]", "(100,250kb]",
                "(250,500kb]", "(500,750kb]", "(750kb,1Mb]", ">1Mb")

out_dir <- file.path("results", "observed-nearest-patterns")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_gene_z    <- file.path(out_dir, "observed_topfracs_nearestgene_z_summary.tsv")
out_gene_bins <- file.path(out_dir, "observed_topfracs_nearestgene_bin_enrichment.tsv")
out_gene_rep  <- file.path(out_dir, "observed_topfracs_nearestgene_null_replicate_stats.tsv.gz")
out_gws_z     <- file.path(out_dir, "observed_topfracs_nearestgwaspeak_z_summary.tsv")
out_gws_bins  <- file.path(out_dir, "observed_topfracs_nearestgwaspeak_bin_enrichment.tsv")
out_gws_rep   <- file.path(out_dir, "observed_topfracs_nearestgwaspeak_null_replicate_stats.tsv.gz")

get_nearest_gene_distances <- function(win_dt) {
  genes_rds <- file.path(project_root, "data", "gencode.v19.genes.protein_coding.rds")
  if (!file.exists(genes_rds)) stop("Missing genes RDS: ", genes_rds)
  genes_gr <- readRDS(genes_rds)
  if (any(grepl("^chr", seqlevels(genes_gr)))) seqlevels(genes_gr) <- sub("^chr", "", seqlevels(genes_gr))
  GenomeInfoDb::seqlevelsStyle(genes_gr) <- "NCBI"
  valid_chr <- intersect(unique(win_dt$chrom), as.character(seqlevels(genes_gr)))
  z <- copy(win_dt[chrom %in% valid_chr])
  if (nrow(z) == 0L) stop("No windows remain after chromosome harmonization.")
  window_gr <- GRanges(seqnames = z$chrom, ranges = IRanges(start = z$start + 1L, end = z$end), window_id = z$window_id)
  GenomeInfoDb::seqlevelsStyle(window_gr) <- "NCBI"
  hits <- distanceToNearest(window_gr, genes_gr, ignore.strand = TRUE)
  nearest_dist <- integer(length(window_gr))
  nearest_dist[queryHits(hits)] <- mcols(hits)$distance
  z[, dist_bp := nearest_dist]
  z[, .(window_id, chrom, start, end, dist_bp)]
}

get_nearest_gwaspeak_distances <- function(win_dt) {
  pval_file <- file.path("results", "all-trait-50k-mean-pvals.tsv")
  if (!file.exists(pval_file)) stop("Missing p-value matrix: ", pval_file)
  pv <- fread(pval_file, na.strings = c("NA", "NaN", "nan", ".", ""))
  pv[, `:=`(chrom = as.integer(chrom), start = as.integer(start), end = as.integer(end))]
  trait_cols <- setdiff(names(pv), c("chrom", "start", "end"))
  for (cc in trait_cols) if (!is.numeric(pv[[cc]])) suppressWarnings(set(pv, j = cc, value = as.numeric(pv[[cc]])))
  pmat <- as.matrix(pv[, ..trait_cols]); storage.mode(pmat) <- "numeric"
  pv[, gws_any := apply(pmat > gwas_log10_thresh, 1, function(v) any(v, na.rm = TRUE))]
  peaks <- pv[gws_any == TRUE, .(chrom, start, end)]
  if (nrow(peaks) == 0L) stop("No GWAS-significant windows found.")
  z <- copy(win_dt)
  z[, chrom_int := suppressWarnings(as.integer(chrom))]
  peaks[, chrom_int := as.integer(chrom)]
  z <- z[is.finite(chrom_int)]; peaks <- peaks[is.finite(chrom_int)]
  window_gr <- GRanges(seqnames = as.character(z$chrom_int), ranges = IRanges(start = z$start + 1L, end = z$end), window_id = z$window_id)
  peak_gr   <- GRanges(seqnames = as.character(peaks$chrom_int), ranges = IRanges(start = peaks$start + 1L, end = peaks$end))
  hits <- distanceToNearest(window_gr, peak_gr, ignore.strand = TRUE)
  nearest_dist <- integer(length(window_gr))
  nearest_dist[queryHits(hits)] <- mcols(hits)$distance
  z[, dist_bp := nearest_dist]
  z[, .(window_id, chrom, start, end, dist_bp)]
}

deadjacent_top <- function(top_dt) {
  x <- copy(top_dt)
  x[, chrom_ord := suppressWarnings(as.integer(chrom))]
  x[is.na(chrom_ord), chrom_ord := 1e6]
  setorder(x, chrom_ord, start, end)
  x[, prev_end    := shift(end), by = chrom]
  x[, is_adj_prev := !is.na(prev_end) & (start == prev_end)]
  x[, run_id      := cumsum(!is_adj_prev), by = chrom]
  keep <- x[, .SD[floor((.N + 1L) / 2L)], by = .(chrom, run_id)]
  keep[, .(window_id, chrom, start, end)]
}

sample_nonadjacent_windows_chr <- function(chr_dt, n_pick) {
  if (n_pick == 0L) return(character())
  n <- nrow(chr_dt)
  max_pick <- ceiling(n / 2)
  if (n_pick > max_pick)
    stop(sprintf("Impossible non-adjacent sample on chr%s: requested %d, max %d from %d windows.",
                 unique(chr_dt$chrom)[1], n_pick, max_pick, n))
  y <- sort(sample.int(n - n_pick + 1L, n_pick, replace = FALSE))
  chr_dt$window_id[y + seq_len(n_pick) - 1L]
}

sample_null_set <- function(universe_by_chr, chr_counts) {
  unlist(lapply(seq_along(chr_counts), function(i) {
    sample_nonadjacent_windows_chr(universe_by_chr[[names(chr_counts)[i]]], as.integer(chr_counts[[i]]))
  }), use.names = FALSE)
}

frac_in_bins <- function(d) {
  b   <- cut(d, breaks = bin_breaks, labels = bin_labels, include.lowest = TRUE, right = TRUE)
  tab <- table(b)
  out <- as.numeric(tab / sum(tab)); names(out) <- names(tab); out
}

run_pattern <- function(ld, dist_lookup) {
  z_rows <- list(); bin_rows <- list(); rep_rows <- list()
  ii <- jj <- kk <- 1L

  for (top_frac in top_fracs) {
    for (pc in pcs) {
      message("  ", pc, ": top ", sprintf("%.2f", 100 * top_frac), "% + ", n_reps, " null reps")
      x <- copy(ld[, .(window_id, chrom, start, end, loading = get(pc))])
      x <- x[is.finite(loading)]
      x[, abs_loading := abs(loading)]; setorder(x, -abs_loading)

      n_top_raw <- max(1L, ceiling(top_frac * nrow(x)))
      top_raw  <- x[seq_len(n_top_raw), .(window_id, chrom, start, end)]
      top_obs  <- deadjacent_top(top_raw)
      if (nrow(top_obs) < 2L) stop("Too few observed peaks after de-adjacency for ", pc)

      chr_counts_dt <- top_obs[, .N, by = chrom][order(as.integer(chrom))]
      chr_counts    <- setNames(chr_counts_dt$N, chr_counts_dt$chrom)

      universe <- unique(dist_lookup[, .(window_id, chrom, start, end, dist_bp)])
      universe[, chrom_ord := suppressWarnings(as.integer(chrom))]
      universe[is.na(chrom_ord), chrom_ord := 1e6]
      setorder(universe, chrom_ord, start, end); universe[, chrom_ord := NULL]
      universe_by_chr <- split(universe, by = "chrom", keep.by = FALSE)

      obs_d     <- dist_lookup[match(top_obs$window_id, window_id), dist_bp]
      obs_d     <- obs_d[is.finite(obs_d)]
      obs_zero  <- mean(obs_d == 0)
      obs_nz    <- obs_d[obs_d > 0]
      obs_med_nz <- if (length(obs_nz) > 0L) median(obs_nz) else NA_real_
      obs_bins  <- frac_in_bins(obs_d)

      null_zero  <- numeric(n_reps); null_med_nz <- numeric(n_reps)
      null_bin   <- matrix(NA_real_, nrow = n_reps, ncol = length(bin_labels), dimnames = list(NULL, bin_labels))
      for (r in seq_len(n_reps)) {
        ids <- sample_null_set(universe_by_chr, chr_counts)
        d   <- dist_lookup[match(ids, window_id), dist_bp]; d <- d[is.finite(d)]
        null_zero[r]   <- mean(d == 0)
        nz <- d[d > 0]; null_med_nz[r] <- if (length(nz) > 0L) median(nz) else NA_real_
        null_bin[r, ]  <- frac_in_bins(d)
        if (r %% 1000L == 0L) message("    rep ", r, "/", n_reps)
      }

      z_rows[[ii]] <- data.table(
        pc = pc, top_frac = top_frac, n_observed = nrow(top_obs),
        obs_zero_fraction = obs_zero,
        null_zero_mean = mean(null_zero, na.rm = TRUE),
        null_zero_sd   = sd(null_zero, na.rm = TRUE),
        z0_overlap     = (obs_zero - mean(null_zero, na.rm = TRUE)) / sd(null_zero, na.rm = TRUE),
        p0_high        = (1 + sum(null_zero >= obs_zero, na.rm = TRUE)) / (1 + sum(is.finite(null_zero))),
        obs_nonzero_median = obs_med_nz,
        null_nonzero_median_mean = mean(null_med_nz, na.rm = TRUE),
        null_nonzero_median_sd   = sd(null_med_nz, na.rm = TRUE),
        zplus_closer             = (mean(null_med_nz, na.rm = TRUE) - obs_med_nz) / sd(null_med_nz, na.rm = TRUE),
        p_nonzero_median_low  = (1 + sum(null_med_nz <= obs_med_nz, na.rm = TRUE)) / (1 + sum(is.finite(null_med_nz))),
        p_nonzero_median_high = (1 + sum(null_med_nz >= obs_med_nz, na.rm = TRUE)) / (1 + sum(is.finite(null_med_nz)))
      ); ii <- ii + 1L

      for (bb in seq_along(bin_labels)) {
        null_vec <- null_bin[, bb]; bname <- bin_labels[bb]
        null_mean_bin <- mean(null_vec, na.rm = TRUE)
        bin_rows[[jj]] <- data.table(
          pc = pc, top_frac = top_frac, bin = bname,
          obs_frac = obs_bins[[bname]], null_mean = null_mean_bin,
          null_lo  = quantile(null_vec, 0.025, na.rm = TRUE),
          null_hi  = quantile(null_vec, 0.975, na.rm = TRUE),
          enrich_diff    = obs_bins[[bname]] - null_mean_bin,
          enrich_diff_lo = obs_bins[[bname]] - quantile(null_vec, 0.975, na.rm = TRUE),
          enrich_diff_hi = obs_bins[[bname]] - quantile(null_vec, 0.025, na.rm = TRUE),
          p_emp_high     = (1 + sum(null_vec >= obs_bins[[bname]], na.rm = TRUE)) / (1 + sum(is.finite(null_vec))),
          p_emp_low      = (1 + sum(null_vec <= obs_bins[[bname]], na.rm = TRUE)) / (1 + sum(is.finite(null_vec))),
          p_emp_two_sided = {
            obs_dev <- abs(obs_bins[[bname]] - null_mean_bin)
            null_dev <- abs(null_vec - null_mean_bin)
            (1 + sum(null_dev >= obs_dev, na.rm = TRUE)) / (1 + sum(is.finite(null_dev)))
          }
        ); jj <- jj + 1L
      }

      rep_rows[[kk]] <- data.table(
        pc = pc, top_frac = top_frac, replicate = seq_len(n_reps),
        null_zero_fraction = null_zero, null_nonzero_median = null_med_nz,
        matrix(null_bin, ncol = length(bin_labels), dimnames = list(NULL, bin_labels))
      ); kk <- kk + 1L
    }
  }
  list(z    = rbindlist(z_rows,    use.names = TRUE, fill = TRUE),
       bins = rbindlist(bin_rows,  use.names = TRUE, fill = TRUE),
       reps = rbindlist(rep_rows,  use.names = TRUE, fill = TRUE))
}

# ---- Load pca_loadings ----
message("Loading pca_loadings_50k.tsv...")
pca_loadings_file <- file.path("results", "pca_loadings_50k.tsv")
if (!file.exists(pca_loadings_file)) stop("Missing: ", pca_loadings_file, "\nRun run-figure1.R first.")
pca_loadings <- fread(pca_loadings_file)

ld <- copy(as.data.table(pca_loadings))
if (!all(c("chrom", "start", "end") %in% names(ld))) {
  if (!"region_label" %in% names(ld)) stop("Need chrom/start/end or region_label in pca_loadings.")
  ld[, c("chrom", "start", "end") := tstrsplit(region_label, "_", fixed = TRUE)]
}
ld[, `:=`(chrom = as.character(chrom), start = as.integer(start), end = as.integer(end))]
ld <- ld[is.finite(start) & is.finite(end) & end > start]
ld[, chrom := gsub("^chr", "", chrom)]
ld[, window_id := paste(chrom, start, end, sep = "_")]
ld[, chrom_ord := suppressWarnings(as.integer(chrom))]
ld[is.na(chrom_ord), chrom_ord := 1e6]
setorder(ld, chrom_ord, start, end); ld[, chrom_ord := NULL]

if (analysis_target %in% c("all", "gene")) {
  message("Observed analysis (nearest gene)...")
  dist_gene <- get_nearest_gene_distances(ld[, .(window_id, chrom, start, end)])
  res_gene  <- run_pattern(ld, dist_gene)
  fwrite(res_gene$z,    out_gene_z,    sep = "\t")
  fwrite(res_gene$bins, out_gene_bins, sep = "\t")
  fwrite(res_gene$reps, out_gene_rep,  sep = "\t")
  message("Wrote: ", out_gene_bins)
}

if (analysis_target %in% c("all", "gwas")) {
  message("Observed analysis (nearest GWAS-significant window)...")
  dist_gws <- get_nearest_gwaspeak_distances(ld[, .(window_id, chrom, start, end)])
  res_gws  <- run_pattern(ld, dist_gws)
  fwrite(res_gws$z,    out_gws_z,    sep = "\t")
  fwrite(res_gws$bins, out_gws_bins, sep = "\t")
  fwrite(res_gws$reps, out_gws_rep,  sep = "\t")
  message("Wrote: ", out_gws_bins)
}

message("Done.")
