#!/usr/bin/env Rscript
# Produces: figures/figureS27_gwasremoved_nearestgene_nonzero_distance_cdf_ribbon_50k_bands_plus{flank_kb}kb.pdf
#           results/gwasremoved_nearestgene_nonzero_distance_cdf_ribbon_50k_bands_plus{flank_kb}kb.tsv
#           results/gwasremoved_nearestgene_nonzero_distance_cdf_ribbon_50k_bands_plus{flank_kb}kb_pvalues.tsv
# Pre-computed inputs:
#   results/gwas_removed_distance_topload_50k_plus{flank_kb}kb/reduced_pca_window_loadings.tsv
#     (from compute-gwas-removed-distances.R)
#   data/gencode.v19.genes.protein_coding.rds
# Run: Rscript post-process-scripts/run-figureS27.R [flank_kb]
#   default flank_kb = 100
# Runtime: ~5-8 hr (10000 null replicates × 4 PCs × 4 bands)
# Scientific question: do top-loading windows remain near gene bodies after GWAS hits are removed?
# Null: windows drawn from the same GWAS-removed universe (retained windows only).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)
set.seed(20260402)

args     <- commandArgs(trailingOnly = TRUE)
flank_kb <- if (length(args) >= 1L) as.integer(args[[1]]) else 100L
if (!is.finite(flank_kb) || flank_kb < 0L) stop("flank_kb must be a non-negative integer.")
flank_label <- paste0("plus", flank_kb, "kb")

pcs     <- c("PC1", "PC2", "PC3", "PC4")
pc_cols <- c("PC1" = "#1f78b4", "PC2" = "#33a02c", "PC3" = "#e31a1c", "PC4" = "#984EA3")
bands   <- data.table(
  band_id    = c("top1", "top2", "top5", "band45_50"),
  start_frac = c(0.00, 0.00, 0.00, 0.45),
  end_frac   = c(0.01, 0.02, 0.05, 0.50),
  label      = c("Top 1%", "Top 2%", "Top 5%", "45-50%")
)
n_reps <- 10000L; bin_bp <- 50000L

out_fig   <- file.path("figures",  sprintf("figureS27_gwasremoved_nearestgene_nonzero_distance_cdf_ribbon_50k_bands_%s.pdf",   flank_label))
out_stats <- file.path("results",  sprintf("gwasremoved_nearestgene_nonzero_distance_cdf_ribbon_50k_bands_%s.tsv",            flank_label))
out_pvals <- file.path("results",  sprintf("gwasremoved_nearestgene_nonzero_distance_cdf_ribbon_50k_bands_%s_pvalues.tsv",    flank_label))
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)

# ---- Load GWAS-removed loadings ----
loadings_file <- file.path("results",
  paste0("gwas_removed_distance_topload_50k_", flank_label),
  "reduced_pca_window_loadings.tsv")
if (!file.exists(loadings_file))
  stop("Missing reduced loadings: ", loadings_file,
       "\nRun: Rscript post-process-scripts/compute-gwas-removed-distances.R ", flank_kb)
ld <- fread(loadings_file)
if (!all(c("chrom", "start", "end") %in% names(ld))) {
  ld[, c("chrom", "start", "end") := tstrsplit(window_id, "_", fixed = TRUE)]
}
ld[, `:=`(chrom = as.character(chrom), start = as.integer(start), end = as.integer(end))]
if (!"window_id" %in% names(ld)) ld[, window_id := paste(chrom, start, end, sep = "_")]
ld[, chrom := gsub("^chr", "", chrom)]
ld <- ld[is.finite(start) & is.finite(end)]

message("Loaded ", nrow(ld), " retained windows from ", loadings_file)

# ---- Load gene annotation ----
genes_rds <- file.path(project_root, "data", "gencode.v19.genes.protein_coding.rds")
if (!file.exists(genes_rds)) stop("Missing genes RDS: ", genes_rds)
genes_gr <- readRDS(genes_rds)
if (any(grepl("^chr", seqlevels(genes_gr)))) seqlevels(genes_gr) <- sub("^chr", "", seqlevels(genes_gr))
GenomeInfoDb::seqlevelsStyle(genes_gr) <- "NCBI"

# ---- Nearest-gene distance for all retained windows ----
get_nearest_gene_distances <- function(win_dt) {
  valid_chr <- intersect(unique(win_dt$chrom), as.character(seqlevels(genes_gr)))
  z <- copy(win_dt[chrom %in% valid_chr])
  window_gr <- GRanges(
    seqnames = z$chrom,
    ranges   = IRanges(start = z$start + 1L, end = z$end),
    window_id = z$window_id
  )
  GenomeInfoDb::seqlevelsStyle(window_gr) <- "NCBI"
  nearest_all  <- distanceToNearest(window_gr, genes_gr, ignore.strand = TRUE)
  nearest_dist <- integer(length(window_gr))
  nearest_dist[queryHits(nearest_all)] <- mcols(nearest_all)$distance
  z[, dist_bp := nearest_dist]
  z[, .(window_id, chrom, start, end, dist_bp)]
}

# ---- De-adjacency: keep middle window from runs of adjacent top windows ----
deadjacent_top <- function(top_dt) {
  x <- copy(top_dt)
  x[, chrom_ord := suppressWarnings(as.integer(chrom))]
  x[is.na(chrom_ord), chrom_ord := 1e6]
  setorder(x, chrom_ord, start, end)
  x[, prev_end  := shift(end),  by = chrom]
  x[, is_adj_prev := !is.na(prev_end) & (start == prev_end)]
  x[, run_id := cumsum(!is_adj_prev), by = chrom]
  keep <- x[, .SD[floor((.N + 1L) / 2L)], by = .(chrom, run_id)]
  keep[, .(window_id, chrom, start, end)]
}

# ---- Non-adjacent null sampling ----
sample_nonadjacent_windows_chr <- function(chr_dt, n_pick) {
  n <- nrow(chr_dt)
  if (n_pick == 0L) return(character())
  max_pick <- ceiling(n / 2)
  if (n_pick > max_pick) {
    warning(sprintf("Non-adjacent sample on chr%s: requested %d > max %d; clamping.",
                    unique(chr_dt$chrom)[1], n_pick, max_pick))
    n_pick <- max_pick
  }
  y <- sort(sample.int(n - n_pick + 1L, n_pick, replace = FALSE))
  chr_dt$window_id[y + seq_len(n_pick) - 1L]
}

sample_null_set <- function(universe_by_chr, chr_counts) {
  chrs <- names(chr_counts)
  unlist(lapply(seq_along(chrs), function(i) {
    sample_nonadjacent_windows_chr(universe_by_chr[[chrs[i]]], as.integer(chr_counts[[i]]))
  }), use.names = FALSE)
}

# ---- Compute distances for all retained windows ----
message("GWAS-removed nearest-gene CDF: computing distances for all retained windows...")
dist_lookup <- get_nearest_gene_distances(ld[, .(window_id, chrom, start, end)])

# Build null universe from the retained windows only
universe <- unique(dist_lookup[, .(window_id, chrom, start, end, dist_bp)])
universe[, chrom_ord := suppressWarnings(as.integer(chrom))]
universe[is.na(chrom_ord), chrom_ord := 1e6]
setorder(universe, chrom_ord, start, end)
universe[, chrom_ord := NULL]
universe_by_chr <- split(universe, by = "chrom", keep.by = FALSE)

# ---- Main loop: 4 bands × 4 PCs ----
curve_rows <- list(); pval_rows <- list(); ii <- 1L; jj <- 1L

for (band_i in seq_len(nrow(bands))) {
  band <- bands[band_i]
  for (pc in pcs) {
    message("GWAS-removed CDF: ", pc, " / ", band$label, "...")
    x <- copy(ld[, .(window_id, chrom, start, end, loading = get(pc))])
    x <- x[is.finite(loading)]
    x[, abs_loading := abs(loading)]; setorder(x, -abs_loading)
    n_all     <- nrow(x)
    start_idx <- min(max(floor(band$start_frac * n_all) + 1L, 1L), n_all)
    end_idx   <- min(max(ceiling(band$end_frac * n_all), start_idx), n_all)
    top_raw   <- x[start_idx:end_idx, .(window_id, chrom, start, end)]
    top_obs   <- deadjacent_top(top_raw)
    chr_counts_dt <- top_obs[, .N, by = chrom]
    chr_counts    <- setNames(chr_counts_dt$N, chr_counts_dt$chrom)

    obs_d  <- dist_lookup[match(top_obs$window_id, window_id), dist_bp]
    obs_d  <- obs_d[is.finite(obs_d) & obs_d > 0]
    if (length(obs_d) == 0L) {
      warning("No non-zero distances for ", pc, " / ", band$label, "; skipping.")
      next
    }
    max_bp <- max(obs_d, na.rm = TRUE)

    null_list <- vector("list", n_reps)
    for (r in seq_len(n_reps)) {
      ids <- sample_null_set(universe_by_chr, chr_counts)
      d   <- dist_lookup[match(ids, window_id), dist_bp]
      null_list[[r]] <- d[is.finite(d) & d > 0]
      if (length(null_list[[r]]) > 0L) max_bp <- max(max_bp, max(null_list[[r]]), na.rm = TRUE)
      if (r %% 1000L == 0L) message("  ", pc, " [", band$label, "] rep ", r, "/", n_reps)
    }
    bins <- seq.int(bin_bp, as.integer(ceiling(max_bp / bin_bp) * bin_bp), by = bin_bp)
    if (length(bins) == 0L) bins <- bin_bp

    obs_curve <- data.table(pc = pc, band_label = band$label, dist_bp = bins,
                             obs_cdf = vapply(bins, function(b) mean(obs_d <= b), numeric(1)))
    rnd_curve <- rbindlist(lapply(seq_len(n_reps), function(r) {
      d <- null_list[[r]]
      data.table(pc = pc, band_label = band$label, random_rep = r, dist_bp = bins,
                 cdf = vapply(bins, function(b) mean(d <= b), numeric(1)))
    }))
    rnd_sum <- rnd_curve[, .(sim_med = median(cdf, na.rm = TRUE),
                              sim_lo  = quantile(cdf, 0.025, na.rm = TRUE),
                              sim_hi  = quantile(cdf, 0.975, na.rm = TRUE)),
                          by = .(pc, band_label, dist_bp)]
    curve_rows[[ii]] <- merge(rnd_sum, obs_curve, by = c("pc", "band_label", "dist_bp"))
    ii <- ii + 1L

    # Empirical p-value: area between observed and null CDFs, vs. null-vs-null distribution
    obs_stat <- rnd_curve[, .(stat = sum(abs(cdf - obs_curve$obs_cdf), na.rm = TRUE) * bin_bp),
                           by = random_rep][, median(stat)]
    rep_ids <- seq_len(n_reps); n_pairs <- 20000L
    r1 <- sample(rep_ids, n_pairs, replace = TRUE)
    r2 <- sample(rep_ids, n_pairs, replace = TRUE)
    eq <- r1 == r2
    while (any(eq)) { r2[eq] <- sample(rep_ids, sum(eq), replace = TRUE); eq <- r1 == r2 }
    pairs <- data.table(pair_id = seq_len(n_pairs), r1 = r1, r2 = r2)
    null_stats <- merge(
      merge(pairs, rnd_curve[, .(random_rep, dist_bp, cdf1 = cdf)],
            by.x = "r1", by.y = "random_rep", all.x = TRUE, allow.cartesian = TRUE),
      rnd_curve[, .(random_rep, dist_bp, cdf2 = cdf)],
      by.x = c("r2", "dist_bp"), by.y = c("random_rep", "dist_bp"),
      all.x = TRUE, allow.cartesian = TRUE
    )[, .(stat = sum(abs(cdf1 - cdf2), na.rm = TRUE) * bin_bp), by = pair_id]
    p_emp <- (1 + sum(null_stats$stat >= obs_stat, na.rm = TRUE)) / (1 + nrow(null_stats))
    pval_rows[[jj]] <- data.table(pc = pc, band_label = band$label,
                                   obs_area_stat_median = obs_stat, p_emp = p_emp)
    jj <- jj + 1L
  }
}

plot_dt  <- rbindlist(curve_rows, use.names = TRUE)
pvals_dt <- rbindlist(pval_rows,  use.names = TRUE)
fwrite(plot_dt,  out_stats, sep = "\t")
fwrite(pvals_dt, out_pvals, sep = "\t")

# ---- Plot ----
ann_dt <- copy(pvals_dt)
ann_dt[, `:=`(x = 1950, y = 0.05, label = sprintf("p = %.3g", p_emp))]
ann_dt[, band_label := factor(band_label, levels = bands$label)]
plot_dt[, band_label := factor(band_label, levels = bands$label)]
plot_dt[, pc := factor(pc, levels = pcs)]
ann_dt[,  pc := factor(pc, levels = pcs)]

p <- ggplot() +
  geom_ribbon(data = plot_dt, aes(x = dist_bp / 1000, ymin = sim_lo, ymax = sim_hi),
              fill = "grey80", alpha = 0.5) +
  geom_line(data = plot_dt, aes(x = dist_bp / 1000, y = sim_med),
            color = "grey45", linewidth = 0.6) +
  geom_step(data = plot_dt, aes(x = dist_bp / 1000, y = obs_cdf, color = pc),
            linewidth = 0.9) +
  geom_text(data = ann_dt, aes(x = x, y = y, label = label, color = pc),
            inherit.aes = FALSE, hjust = 1, vjust = 0, size = 2.7) +
  scale_color_manual(values = pc_cols, drop = FALSE) +
  scale_x_continuous(limits = c(0, 2000), breaks = seq(0, 2000, by = 500)) +
  facet_grid(pc ~ band_label, scales = "fixed") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.2),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3)) +
  labs(title = sprintf("GWAS-removed PCA (±%d kb flanks): nearest-gene distance CDF", flank_kb),
       x = "Distance to most proximal gene body (non-zero, kb)",
       y = "Cumulative fraction")

ggsave(out_fig, p, width = 14.8, height = 10.5)
message("Wrote: ", out_fig)
message("Wrote: ", out_stats)
message("Wrote: ", out_pvals)
