#!/usr/bin/env Rscript
# Produces: results/observed-nearest-patterns-{window_tag}-top1_2_5_band45_50/
#   observed_nearestgene_bin_enrichment_{window_tag}_top1_2_5_band45_50.tsv
#   observed_nearestgene_bandsummary_{window_tag}_top1_2_5_band45_50.tsv
#   observed_nearestgene_null_replicates_{window_tag}_top1_2_5_band45_50.tsv.gz
#   observed_nearestgene_ranked_windows_{window_tag}.tsv.gz
#   observed_nearestgwaspeak_bin_enrichment_{window_tag}_top1_2_5_band45_50.tsv
#   observed_nearestgwaspeak_bandsummary_{window_tag}_top1_2_5_band45_50.tsv
#   observed_nearestgwaspeak_null_replicates_{window_tag}_top1_2_5_band45_50.tsv.gz
#   observed_nearestgwaspeak_ranked_windows_{window_tag}.tsv.gz
# Pre-computed inputs:
#   50k: results/pca_loadings_50k.tsv (from run-figure1.R), results/all-trait-50k-mean-pvals.tsv †
#   100k: results/pca_multiscale_anchored/pca_windows_100k_anchored.tsv †,
#          results/all-trait-100k-mean-pvals.tsv †
#   data/gencode.v19.genes.protein_coding.rds (via data/ symlink)
# Run: Rscript figure-scripts/compute-nearestgene-pattern-bands.R 50k
#      Rscript figure-scripts/compute-nearestgene-pattern-bands.R 100k
# Runtime: ~60-90 min (3000 null reps × 4 PCs × 4 bands × 2 analyses)

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

args <- commandArgs(trailingOnly = TRUE)
window_tag <- if (length(args) >= 1L) args[[1]] else "50k"
if (!window_tag %in% c("50k", "100k", "gwasremoved")) stop("window_tag must be one of: 50k, 100k, gwasremoved")
bin_mode <- if (length(args) >= 2L) args[[2]] else "default"
analysis_kind <- if (length(args) >= 3L) args[[3]] else "both"
if (!bin_mode %in% c("default", "common100")) stop("bin_mode must be one of: default, common100")
if (!analysis_kind %in% c("both", "gene_only", "gwas_only")) stop("analysis_kind must be one of: both, gene_only, gwas_only")
# For gwasremoved: optional 4th arg is flank_kb (default 100); gene analysis only.
flank_kb <- if (window_tag == "gwasremoved" && length(args) >= 4L) as.integer(args[[4]]) else 100L
if (window_tag == "gwasremoved" && analysis_kind == "both") {
  message("Note: gwasremoved tag only supports gene_only analysis; skipping GWAS-peak distances.")
  analysis_kind <- "gene_only"
}
file_tag <- if (window_tag == "gwasremoved") sprintf("gwasremoved%dkb", flank_kb) else window_tag

set.seed(20260330)

pcs <- c("PC1", "PC2", "PC3", "PC4")
pc_cols <- c("PC1" = "#1f78b4", "PC2" = "#33a02c", "PC3" = "#e31a1c", "PC4" = "#984EA3")
n_reps <- 3000L

bands <- data.table(
  band_id = c("top1", "top2", "top5", "band45_50"),
  start_frac = c(0.00, 0.00, 0.00, 0.45),
  end_frac = c(0.01, 0.02, 0.05, 0.50),
  label = c("Top 1%", "Top 2%", "Top 5%", "45-50%")
)

if (bin_mode == "common100") {
  bin_breaks <- c(-Inf, 0, 100000, 200000, 300000, Inf)
  bin_labels <- c("0", "(0,100kb]", "(100,200kb]", "(200,300kb]", ">300kb")
} else if (window_tag %in% c("50k", "gwasremoved")) {
  bin_breaks <- c(-Inf, 0, 50000, 100000, 250000, Inf)
  bin_labels <- c("0", "(0,50kb]", "(50,100kb]", "(100,250kb]", ">250kb")
} else {
  bin_breaks <- c(-Inf, 0, 100000, 200000, 300000, Inf)
  bin_labels <- c("0", "(0,100kb]", "(100,200kb]", "(200,300kb]", ">300kb")
}

out_dir <- file.path("results", sprintf("observed-nearest-patterns-%s-top1_2_5_band45_50", file_tag))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

out_gene_z      <- file.path(out_dir, sprintf("observed_nearestgene_bandsummary_%s_top1_2_5_band45_50.tsv", file_tag))
out_gene_bins   <- file.path(out_dir, sprintf("observed_nearestgene_bin_enrichment_%s_top1_2_5_band45_50.tsv", file_tag))
out_gene_rep    <- file.path(out_dir, sprintf("observed_nearestgene_null_replicates_%s_top1_2_5_band45_50.tsv.gz", file_tag))
out_gene_ranked <- file.path(out_dir, sprintf("observed_nearestgene_ranked_windows_%s.tsv.gz", file_tag))
out_gws_z      <- file.path(out_dir, sprintf("observed_nearestgwaspeak_bandsummary_%s_top1_2_5_band45_50.tsv", file_tag))
out_gws_bins   <- file.path(out_dir, sprintf("observed_nearestgwaspeak_bin_enrichment_%s_top1_2_5_band45_50.tsv", file_tag))
out_gws_rep    <- file.path(out_dir, sprintf("observed_nearestgwaspeak_null_replicates_%s_top1_2_5_band45_50.tsv.gz", file_tag))
out_gws_ranked <- file.path(out_dir, sprintf("observed_nearestgwaspeak_ranked_windows_%s.tsv.gz", file_tag))
out_progress   <- file.path(out_dir, sprintf("progress_%s_top1_2_5_band45_50.log", file_tag))

log_progress <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", paste0(..., collapse = ""))
  message(msg)
  cat(msg, "\n", file = out_progress, append = TRUE)
  flush.console()
}

get_nearest_gene_distances <- function(win_dt) {
  genes_rds <- file.path(project_root, "data", "gencode.v19.genes.protein_coding.rds")
  if (!file.exists(genes_rds)) stop("Missing genes RDS: ", genes_rds)
  genes_gr <- readRDS(genes_rds)
  if (any(grepl("^chr", seqlevels(genes_gr)))) {
    seqlevels(genes_gr) <- sub("^chr", "", seqlevels(genes_gr))
  }
  GenomeInfoDb::seqlevelsStyle(genes_gr) <- "NCBI"

  valid_chr <- intersect(unique(win_dt$chrom), as.character(seqlevels(genes_gr)))
  z <- copy(win_dt[chrom %in% valid_chr])
  if (nrow(z) == 0L) stop("No windows remain after chromosome harmonization for gene distances.")

  window_gr <- GRanges(
    seqnames = z$chrom,
    ranges = IRanges(start = z$start + 1L, end = z$end),
    window_id = z$window_id
  )
  GenomeInfoDb::seqlevelsStyle(window_gr) <- "NCBI"
  nearest_all <- distanceToNearest(window_gr, genes_gr, ignore.strand = TRUE)
  nearest_dist <- integer(length(window_gr))
  nearest_dist[queryHits(nearest_all)] <- mcols(nearest_all)$distance
  z[, dist_bp := nearest_dist]
  z[, .(window_id, chrom, start, end, dist_bp)]
}

get_nearest_gwaspeak_distances <- function(win_dt) {
  pval_file <- file.path("results", sprintf("all-trait-%s-mean-pvals.tsv", window_tag))
  if (!file.exists(pval_file)) stop("Missing p-value matrix: ", pval_file)
  pv <- fread(pval_file, na.strings = c("NA", "NaN", "nan", ".", ""))
  if (!all(c("chrom", "start", "end") %in% names(pv))) {
    stop("P-value matrix missing chrom/start/end.")
  }
  pv[, `:=`(chrom = as.integer(chrom), start = as.integer(start), end = as.integer(end))]
  trait_cols <- setdiff(names(pv), c("chrom", "start", "end"))
  for (cc in trait_cols) {
    if (!is.numeric(pv[[cc]])) suppressWarnings(set(pv, j = cc, value = as.numeric(pv[[cc]])))
  }
  pmat <- as.matrix(pv[, ..trait_cols]); storage.mode(pmat) <- "numeric"
  pv[, gws_any := apply(pmat > 7.30103, 1, function(v) any(v, na.rm = TRUE))]
  peaks <- pv[gws_any == TRUE, .(chrom, start, end)]
  if (nrow(peaks) == 0L) stop("No GWAS-significant windows found in ", pval_file)

  z <- copy(win_dt)
  z[, chrom_int := suppressWarnings(as.integer(chrom))]
  peaks[, chrom_int := as.integer(chrom)]
  z <- z[is.finite(chrom_int)]
  peaks <- peaks[is.finite(chrom_int)]

  window_gr <- GRanges(
    seqnames = as.character(z$chrom_int),
    ranges = IRanges(start = z$start + 1L, end = z$end),
    window_id = z$window_id
  )
  peak_gr <- GRanges(
    seqnames = as.character(peaks$chrom_int),
    ranges = IRanges(start = peaks$start + 1L, end = peaks$end)
  )
  nearest_all <- distanceToNearest(window_gr, peak_gr, ignore.strand = TRUE)
  nearest_dist <- integer(length(window_gr))
  nearest_dist[queryHits(nearest_all)] <- mcols(nearest_all)$distance
  z[, dist_bp := nearest_dist]
  z[, .(window_id, chrom, start, end, dist_bp)]
}

deadjacent_top <- function(top_dt) {
  x <- copy(top_dt)
  x[, chrom_ord := suppressWarnings(as.integer(chrom))]
  x[is.na(chrom_ord), chrom_ord := 1e6]
  setorder(x, chrom_ord, start, end)
  x[, prev_end := shift(end), by = chrom]
  x[, is_adj_prev := !is.na(prev_end) & (start == prev_end)]
  x[, run_id := cumsum(!is_adj_prev), by = chrom]
  keep <- x[, .SD[floor((.N + 1L) / 2L)], by = .(chrom, run_id)]
  keep[, .(window_id, chrom, start, end)]
}

sample_nonadjacent_windows_chr <- function(chr_dt, n_pick) {
  n <- nrow(chr_dt)
  if (n_pick == 0L) return(character())
  max_pick <- ceiling(n / 2)
  if (n_pick > max_pick) {
    stop(sprintf("Impossible non-adjacent sample on chr%s: requested %d, max %d from %d windows.",
                 unique(chr_dt$chrom)[1], n_pick, max_pick, n))
  }
  y <- sort(sample.int(n - n_pick + 1L, n_pick, replace = FALSE))
  idx <- y + seq_len(n_pick) - 1L
  chr_dt$window_id[idx]
}

sample_null_set <- function(universe_by_chr, chr_counts) {
  sampled <- vector("list", length(chr_counts))
  chrs <- names(chr_counts)
  for (i in seq_along(chrs)) {
    chr_i <- chrs[i]
    n_i <- as.integer(chr_counts[[i]])
    chr_dt <- universe_by_chr[[chr_i]]
    if (is.null(chr_dt)) stop("Chromosome missing from universe: ", chr_i)
    sampled[[i]] <- sample_nonadjacent_windows_chr(chr_dt, n_i)
  }
  unlist(sampled, use.names = FALSE)
}

frac_in_bins <- function(d) {
  b <- cut(d, breaks = bin_breaks, labels = bin_labels, include.lowest = TRUE, right = TRUE)
  tab <- table(b)
  out <- as.numeric(tab / sum(tab))
  names(out) <- names(tab)
  out
}

build_ranked_windows <- function(ld, dist_lookup) {
  ranked_rows <- vector("list", length(pcs))
  for (i in seq_along(pcs)) {
    pc <- pcs[[i]]
    x <- copy(ld[, .(window_id, chrom, start, end, loading = get(pc))])
    x <- x[is.finite(loading)]
    x[, abs_loading := abs(loading)]
    setorder(x, -abs_loading, chrom, start, end)
    x[, rank := seq_len(.N)]
    x[, rank_frac := rank / .N]
    x[, percentile := 100 * rank_frac]
    x[, pc := pc]
    x <- dist_lookup[x, on = "window_id"]
    setcolorder(x, c("pc", "window_id", "chrom", "start", "end", "loading", "abs_loading", "rank", "rank_frac", "percentile", "dist_bp"))
    ranked_rows[[i]] <- x[]
  }
  rbindlist(ranked_rows, use.names = TRUE, fill = TRUE)
}

run_pattern <- function(ld, dist_lookup) {
  z_rows <- list()
  bin_rows <- list()
  rep_rows <- list()
  ii <- 1L
  jj <- 1L
  kk <- 1L

  universe <- unique(dist_lookup[, .(window_id, chrom, start, end, dist_bp)])
  universe[, chrom_ord := suppressWarnings(as.integer(chrom))]
  universe[is.na(chrom_ord), chrom_ord := 1e6]
  setorder(universe, chrom_ord, start, end)
  universe[, chrom_ord := NULL]
  universe_by_chr <- split(universe, by = "chrom", keep.by = FALSE)

  for (b in seq_len(nrow(bands))) {
    band <- bands[b]
    for (pc in pcs) {
      log_progress("Running ", pc, ": ", band$label, " + ", n_reps, " null reps")
      x <- copy(ld[, .(window_id, chrom, start, end, loading = get(pc))])
      x <- x[is.finite(loading)]
      x[, abs_loading := abs(loading)]
      setorder(x, -abs_loading)

      n <- nrow(x)
      start_idx <- floor(band$start_frac * n) + 1L
      end_idx <- max(start_idx, ceiling(band$end_frac * n))
      band_raw <- x[start_idx:end_idx, .(window_id, chrom, start, end)]
      band_obs <- deadjacent_top(band_raw)
      n_obs <- nrow(band_obs)
      if (n_obs < 2L) stop("Too few observed peaks after de-adjacency for ", pc, " / ", band$label)

      chr_counts_dt <- band_obs[, .N, by = chrom][order(as.integer(chrom))]
      chr_counts <- setNames(chr_counts_dt$N, chr_counts_dt$chrom)

      obs_d <- dist_lookup[match(band_obs$window_id, window_id), dist_bp]
      obs_d <- obs_d[is.finite(obs_d)]
      obs_zero <- mean(obs_d == 0)
      obs_nonzero <- obs_d[obs_d > 0]
      obs_med_nz <- if (length(obs_nonzero) > 0L) median(obs_nonzero) else NA_real_
      obs_bins <- frac_in_bins(obs_d)

      null_zero <- numeric(n_reps)
      null_med_nz <- numeric(n_reps)
      null_bin <- matrix(NA_real_, nrow = n_reps, ncol = length(bin_labels), dimnames = list(NULL, bin_labels))
      for (r in seq_len(n_reps)) {
        ids <- sample_null_set(universe_by_chr, chr_counts)
        d <- dist_lookup[match(ids, window_id), dist_bp]
        d <- d[is.finite(d)]
        null_zero[r] <- mean(d == 0)
        nz <- d[d > 0]
        null_med_nz[r] <- if (length(nz) > 0L) median(nz) else NA_real_
        null_bin[r, ] <- frac_in_bins(d)
      }

      z_rows[[ii]] <- data.table(
        pc = pc,
        band_id = band$band_id,
        band_label = band$label,
        start_frac = band$start_frac,
        end_frac = band$end_frac,
        n_observed = n_obs,
        obs_zero_fraction = obs_zero,
        null_zero_mean = mean(null_zero, na.rm = TRUE),
        null_zero_sd = sd(null_zero, na.rm = TRUE),
        z0_overlap = (obs_zero - mean(null_zero, na.rm = TRUE)) / sd(null_zero, na.rm = TRUE),
        p0_high = (1 + sum(null_zero >= obs_zero, na.rm = TRUE)) / (1 + sum(is.finite(null_zero))),
        obs_nonzero_median = obs_med_nz,
        null_nonzero_median_mean = mean(null_med_nz, na.rm = TRUE),
        null_nonzero_median_sd = sd(null_med_nz, na.rm = TRUE),
        zplus_closer = (mean(null_med_nz, na.rm = TRUE) - obs_med_nz) / sd(null_med_nz, na.rm = TRUE),
        p_nonzero_median_low = (1 + sum(null_med_nz <= obs_med_nz, na.rm = TRUE)) / (1 + sum(is.finite(null_med_nz))),
        p_nonzero_median_high = (1 + sum(null_med_nz >= obs_med_nz, na.rm = TRUE)) / (1 + sum(is.finite(null_med_nz)))
      )
      ii <- ii + 1L

      for (bb in seq_along(bin_labels)) {
        null_vec <- null_bin[, bb]
        bname <- bin_labels[bb]
        bin_rows[[jj]] <- data.table(
          pc = pc,
          band_id = band$band_id,
          band_label = band$label,
          bin = bname,
          obs_frac = obs_bins[[bname]],
          null_mean = mean(null_vec, na.rm = TRUE),
          null_lo = quantile(null_vec, 0.025, na.rm = TRUE),
          null_hi = quantile(null_vec, 0.975, na.rm = TRUE),
          enrich_diff = obs_bins[[bname]] - mean(null_vec, na.rm = TRUE),
          enrich_diff_lo = obs_bins[[bname]] - quantile(null_vec, 0.975, na.rm = TRUE),
          enrich_diff_hi = obs_bins[[bname]] - quantile(null_vec, 0.025, na.rm = TRUE),
          p_emp_high = (1 + sum(null_vec >= obs_bins[[bname]], na.rm = TRUE)) / (1 + sum(is.finite(null_vec))),
          p_emp_low = (1 + sum(null_vec <= obs_bins[[bname]], na.rm = TRUE)) / (1 + sum(is.finite(null_vec))),
          p_emp_two_sided = {
            null_mean_bin <- mean(null_vec, na.rm = TRUE)
            obs_dev <- abs(obs_bins[[bname]] - null_mean_bin)
            null_dev <- abs(null_vec - null_mean_bin)
            (1 + sum(null_dev >= obs_dev, na.rm = TRUE)) / (1 + sum(is.finite(null_dev)))
          },
          p_emp_directional = if (obs_bins[[bname]] - mean(null_vec, na.rm = TRUE) >= 0) {
            (1 + sum(null_vec >= obs_bins[[bname]], na.rm = TRUE)) / (1 + sum(is.finite(null_vec)))
          } else {
            (1 + sum(null_vec <= obs_bins[[bname]], na.rm = TRUE)) / (1 + sum(is.finite(null_vec)))
          }
        )
        jj <- jj + 1L
      }

      rep_rows[[kk]] <- data.table(
        pc = pc,
        band_id = band$band_id,
        band_label = band$label,
        replicate = seq_len(n_reps),
        null_zero_fraction = null_zero,
        null_nonzero_median = null_med_nz
      )
      kk <- kk + 1L
    }
  }

  list(
    z = rbindlist(z_rows, use.names = TRUE, fill = TRUE),
    bins = rbindlist(bin_rows, use.names = TRUE, fill = TRUE),
    reps = rbindlist(rep_rows, use.names = TRUE, fill = TRUE)
  )
}

# ---- Load window loadings ----
log_progress("Loading ", window_tag, " window loadings...")
if (window_tag == "50k") {
  pca_loadings_file <- file.path("results", "pca_loadings_50k.tsv")
  if (!file.exists(pca_loadings_file)) stop("Missing: ", pca_loadings_file, "\nRun run-figure1.R first.")
  ld <- fread(pca_loadings_file)
  if ("chrom" %in% names(ld) && !"window_id" %in% names(ld)) {
    ld[, window_id := paste(chrom, start, end, sep = "_")]
  }
  ld[, `:=`(chrom = as.character(chrom), start = as.integer(start), end = as.integer(end))]
} else if (window_tag == "gwasremoved") {
  reduced_file <- file.path("results",
    paste0("gwas_removed_distance_topload_50k_plus", flank_kb, "kb"),
    "reduced_pca_window_loadings.tsv")
  if (!file.exists(reduced_file))
    stop("Missing: ", reduced_file,
         "\nRun: Rscript figure-scripts/compute-gwas-removed-distances.R ", flank_kb)
  ld <- fread(reduced_file)
  if (!all(c("chrom", "start", "end") %in% names(ld))) {
    ld[, c("chrom", "start", "end") := tstrsplit(window_id, "_", fixed = TRUE)]
  }
  ld[, `:=`(chrom = as.character(chrom), start = as.integer(start), end = as.integer(end))]
  log_progress("Loaded ", nrow(ld), " retained (GWAS-removed) windows from ", reduced_file)
} else {
  anchored_file <- file.path("results", "pca_multiscale_anchored", "pca_windows_100k_anchored.tsv")
  if (!file.exists(anchored_file)) stop("Missing: ", anchored_file, "\nCopy from ../directional-coherence/results/pca_multiscale_anchored/")
  ld <- fread(anchored_file)
  if (!all(c("window_id", pcs) %in% names(ld))) stop("Missing expected columns in pca_windows file")
  ld[, c("chrom", "start", "end") := tstrsplit(window_id, "_", fixed = TRUE)]
  ld[, `:=`(chrom = as.character(chrom), start = as.integer(start), end = as.integer(end))]
}
ld <- ld[is.finite(start) & is.finite(end) & end > start]
ld[, chrom := gsub("^chr", "", chrom)]
if (!"window_id" %in% names(ld)) ld[, window_id := paste(chrom, start, end, sep = "_")]
ld[, chrom_ord := suppressWarnings(as.integer(chrom))]
ld[is.na(chrom_ord), chrom_ord := 1e6]
setorder(ld, chrom_ord, start, end)
ld[, chrom_ord := NULL]

if (analysis_kind %in% c("both", "gene_only")) {
  log_progress("Observed analysis (nearest gene, ", window_tag, ")...")
  dist_gene <- get_nearest_gene_distances(ld[, .(window_id, chrom, start, end)])
  ranked_gene <- build_ranked_windows(ld, dist_gene)
  fwrite(ranked_gene, out_gene_ranked, sep = "\t")
  log_progress("Wrote ", out_gene_ranked)
  res_gene <- run_pattern(ld, dist_gene)
  fwrite(res_gene$z, out_gene_z, sep = "\t")
  fwrite(res_gene$bins, out_gene_bins, sep = "\t")
  fwrite(res_gene$reps, out_gene_rep, sep = "\t")
  log_progress("Wrote ", out_gene_z)
  log_progress("Wrote ", out_gene_bins)
  log_progress("Wrote ", out_gene_rep)
}

if (analysis_kind %in% c("both", "gwas_only")) {
  log_progress("Observed analysis (nearest GWAS-significant window, ", window_tag, ")...")
  dist_gws <- get_nearest_gwaspeak_distances(ld[, .(window_id, chrom, start, end)])
  ranked_gws <- build_ranked_windows(ld, dist_gws)
  fwrite(ranked_gws, out_gws_ranked, sep = "\t")
  log_progress("Wrote ", out_gws_ranked)
  res_gws <- run_pattern(ld, dist_gws)
  fwrite(res_gws$z, out_gws_z, sep = "\t")
  fwrite(res_gws$bins, out_gws_bins, sep = "\t")
  fwrite(res_gws$reps, out_gws_rep, sep = "\t")
  log_progress("Wrote ", out_gws_z)
  log_progress("Wrote ", out_gws_bins)
  log_progress("Wrote ", out_gws_rep)
}
