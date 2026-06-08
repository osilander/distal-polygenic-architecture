#!/usr/bin/env Rscript
# Produces: results/allbyall_cosine_matrices/entrysign_perm100_procrustes_normalized.tsv
#           results/allbyall_cosine_matrices/entrysign_observed_procrustes_vs_others.tsv
#           results/allbyall_cosine_matrices/entrysign_perm100_procrustes_summary.tsv
#           results/allbyall_cosine_matrices/entrysign_perm100_procrustes_normalized_overall.tsv
# Pre-computed inputs:
#   results/all-trait-50k-mean-{effects,pvals}.tsv
#   results/pca_multiscale_anchored/pca_traits_{25k,100k,200k,500k,1m}_anchored.tsv
#   results/svd-nulls-50k/withinblock_perm/withinblock_perm_iter001_50k_trait_vectors.tsv.gz
#   results/gwas_removed_distance_topload_50k_plus{100,150}kb/reduced_pca_trait_scores.tsv
#   results/figureS_pca_absolute_50k_trait_scores.tsv
#   results/pca_abs50k_from_windows/trait_scores.tsv
# Run: Rscript figure-scripts/run-entrysign-perm-procrustes-normalized.R

suppressPackageStartupMessages({
  library(data.table)
  library(irlba)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)
source(file.path("figure-scripts", "helpers-50k-matrix.R"))

out_dir <- file.path("results", "allbyall_cosine_matrices")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_rand         <- file.path(out_dir, "entrysign_perm100_procrustes_vs_others.tsv")
out_obs          <- file.path(out_dir, "entrysign_observed_procrustes_vs_others.tsv")
out_summary      <- file.path(out_dir, "entrysign_perm100_procrustes_summary.tsv")
out_norm         <- file.path(out_dir, "entrysign_perm100_procrustes_normalized.tsv")
out_norm_overall <- file.path(out_dir, "entrysign_perm100_procrustes_normalized_overall.tsv")
n_perm  <- 100L
n_pcs   <- 80L
pcs_all <- paste0("PC", 1:75)
k_values <- c(4L, 8L, 25L, 50L, 75L)
set.seed(20260320)

anchor_map <- c(PC1 = "BMI_03", PC2 = "StandHgt", PC3 = "Neurotic_30", PC4 = "DiaBP_A59")

label_map <- c(
  "windows_25kb"             = "25Kb",
  "observed_50k_signed"      = "50Kb",
  "windows_100kb"            = "100Kb",
  "windows_200kb"            = "200 Kb",
  "windows_500kb"            = "500Kb",
  "windows_1mb"              = "1Mb",
  "permuted_within_500kb"    = "500Kb rand.",
  "pickrell_randomized"      = "LD blocks rand.",
  "entry_sign_randomized"    = "Entry sign rand.",
  "gwas_removed_plus150kb"    = "No GWAS loci +- 150Kb",
  "gwas_removed_plus100kb"   = "No GWAS loci +- 100Kb",
  "abs_from_signed_matrix"   = "Matrix abs. value",
  "abs_from_windows_filtered" = "Variant abs value"
)
condition_order <- names(label_map)
target_conds <- setdiff(condition_order, "entry_sign_randomized")  # 12 targets

load_trait_scores <- function(path, trait_col_candidates = c("trait_abbrev", "trait", "Trait")) {
  if (!file.exists(path)) stop("Missing file: ", path)
  dt <- tryCatch(
    fread(path),
    error = function(e) fread(cmd = paste("gzip -dc", shQuote(path)))
  )
  tr_col <- trait_col_candidates[trait_col_candidates %in% names(dt)][1]
  if (is.na(tr_col)) stop("No trait column in ", path)
  setnames(dt, tr_col, "trait_abbrev")
  keep <- c("trait_abbrev", intersect(pcs_all, names(dt)))
  out <- dt[, ..keep]
  for (pc in setdiff(names(out), "trait_abbrev")) set(out, j = pc, value = as.numeric(out[[pc]]))
  out
}

align_anchor <- function(dt) {
  x <- copy(dt)
  for (pc in intersect(names(anchor_map), names(x))) {
    tr <- anchor_map[[pc]]
    idx <- which(x$trait_abbrev == tr)
    if (length(idx) == 0L) next
    aval <- x[[pc]][idx[1]]
    if (is.finite(aval) && aval < 0) x[[pc]] <- -x[[pc]]
  }
  x
}

procrustes_similarity <- function(A, B) {
  if (!all(dim(A) == dim(B))) return(NA_real_)
  if (nrow(A) == 0L || ncol(A) == 0L) return(NA_real_)
  sa <- sqrt(sum(A * A)); sb <- sqrt(sum(B * B))
  if (!is.finite(sa) || !is.finite(sb) || sa == 0 || sb == 0) return(NA_real_)
  d <- svd(crossprod(A, B), nu = 0, nv = 0)$d
  ps <- sum(d) / (sa * sb)
  if (is.finite(ps)) ps <- max(0, min(1, ps))
  ps
}

pc_cols_sorted <- function(dt) {
  pcs <- grep("^PC[0-9]+$", names(dt), value = TRUE)
  pcs[order(as.integer(sub("^PC", "", pcs)))]
}

proc_pair <- function(A_dt, B_dt, k) {
  common_pcs <- intersect(pc_cols_sorted(A_dt), pc_cols_sorted(B_dt))
  if (length(common_pcs) < k) return(NA_real_)
  use <- common_pcs[seq_len(k)]
  m <- merge(
    A_dt[, c("trait_abbrev", use), with = FALSE],
    B_dt[, c("trait_abbrev", use), with = FALSE],
    by = "trait_abbrev",
    suffixes = c("_a", "_b")
  )
  if (nrow(m) == 0L) return(NA_real_)
  A <- as.matrix(m[, paste0(use, "_a"), with = FALSE])
  B <- as.matrix(m[, paste0(use, "_b"), with = FALSE])
  procrustes_similarity(A, B)
}

message("Loading 50k matrix...")
eff   <- fread(file.path("results", "all-trait-50k-mean-effects.tsv"))
pvals <- fread(file.path("results", "all-trait-50k-mean-pvals.tsv"))
prep  <- prepare_X_50k(eff, pvals)
X     <- t(prep$X)  # traits x windows

message("Computing fixed condition trait-score tables...")
p_obs <- prcomp_irlba(X, n = n_pcs, center = FALSE, scale. = FALSE)
rownames(p_obs$x) <- rownames(X)
obs <- as.data.table(p_obs$x)
obs[, trait_abbrev := rownames(p_obs$x)]
obs <- obs[, c("trait_abbrev", intersect(pcs_all, names(obs))), with = FALSE]

# 500kb within-chunk permutation (deterministic single condition).
win_ids <- colnames(X)
coords <- tstrsplit(win_ids, "_", fixed = TRUE)
coords_dt <- data.table(chrom = as.integer(coords[[1]]), start = as.integer(coords[[2]]))
coords_dt[, chunk_id := paste0(chrom, "_", floor(start / 500000L))]
idx_list <- split(seq_len(ncol(X)), coords_dt$chunk_id)
set.seed(51001)
Xp <- X
for (r in seq_len(nrow(Xp))) {
  xv <- Xp[r, ]
  for (idx in idx_list) if (length(idx) > 1L) xv[idx] <- xv[sample(idx, length(idx), replace = FALSE)]
  Xp[r, ] <- xv
}
p_500 <- prcomp_irlba(Xp, n = n_pcs, center = FALSE, scale. = FALSE)
rownames(p_500$x) <- rownames(X)
perm500 <- as.data.table(p_500$x)
perm500[, trait_abbrev := rownames(p_500$x)]
perm500 <- perm500[, c("trait_abbrev", intersect(pcs_all, names(perm500))), with = FALSE]

fixed <- list(
  windows_25kb     = load_trait_scores(file.path("results", "pca_multiscale_anchored", "pca_traits_25k_anchored.tsv")),
  observed_50k_signed = obs,
  windows_100kb    = load_trait_scores(file.path("results", "pca_multiscale_anchored", "pca_traits_100k_anchored.tsv")),
  windows_200kb    = load_trait_scores(file.path("results", "pca_multiscale_anchored", "pca_traits_200k_anchored.tsv")),
  windows_500kb    = load_trait_scores(file.path("results", "pca_multiscale_anchored", "pca_traits_500k_anchored.tsv")),
  windows_1mb      = load_trait_scores(file.path("results", "pca_multiscale_anchored", "pca_traits_1m_anchored.tsv")),
  permuted_within_500kb  = perm500,
  pickrell_randomized    = load_trait_scores(file.path("results", "svd-nulls-50k", "withinblock_perm", "withinblock_perm_iter001_50k_trait_vectors.tsv.gz")),
  gwas_removed_plus150kb  = load_trait_scores(file.path("results", "gwas_removed_distance_topload_50k_plus150kb/reduced_pca_trait_scores.tsv")),
  gwas_removed_plus100kb = load_trait_scores(file.path("results", "gwas_removed_distance_topload_50k_plus100kb", "reduced_pca_trait_scores.tsv")),
  abs_from_signed_matrix    = load_trait_scores(file.path("results", "figureS_pca_absolute_50k_trait_scores.tsv")),
  abs_from_windows_filtered = load_trait_scores(file.path("results", "pca_abs50k_from_windows", "trait_scores.tsv"))
)
fixed <- lapply(fixed, align_anchor)

message("Computing Proc_obs vs 12 targets...")
obs_rows <- list()
ii <- 1L
for (k in k_values) {
  for (cond in target_conds) {
    ps <- proc_pair(fixed[["observed_50k_signed"]], fixed[[cond]], k)
    obs_rows[[ii]] <- data.table(k = k, condition = cond, proc_obs = ps)
    ii <- ii + 1L
  }
}
obs_dt <- rbindlist(obs_rows, use.names = TRUE, fill = TRUE)
fwrite(obs_dt, out_obs, sep = "\t")

message("Running ", n_perm, " entry-wise sign-randomized permutations...")
rand_rows <- vector("list", n_perm)
for (pp in seq_len(n_perm)) {
  signs <- matrix(sample(c(-1, 1), length(X), replace = TRUE), nrow = nrow(X), ncol = ncol(X))
  Xe <- X * signs
  p_e <- prcomp_irlba(Xe, n = n_pcs, center = FALSE, scale. = FALSE)
  rownames(p_e$x) <- rownames(X)
  entry_dt <- as.data.table(p_e$x)
  entry_dt[, trait_abbrev := rownames(p_e$x)]
  entry_dt <- entry_dt[, c("trait_abbrev", intersect(pcs_all, names(entry_dt))), with = FALSE]
  entry_dt <- align_anchor(entry_dt)

  rr <- list()
  jj <- 1L
  for (k in k_values) {
    for (cond in target_conds) {
      ps <- proc_pair(entry_dt, fixed[[cond]], k)
      rr[[jj]] <- data.table(iter = pp, k = k, condition = cond, proc_rand = ps)
      jj <- jj + 1L
    }
  }
  rand_rows[[pp]] <- rbindlist(rr, use.names = TRUE, fill = TRUE)
  if (pp == 1L || pp %% 10L == 0L || pp == n_perm) {
    message("  permutation ", pp, "/", n_perm, " done")
  }
}
rand_dt <- rbindlist(rand_rows, use.names = TRUE, fill = TRUE)
fwrite(rand_dt, out_rand, sep = "\t")

summary_dt <- rand_dt[, .(
  proc_rand_mean = mean(proc_rand, na.rm = TRUE),
  proc_rand_sd   = sd(proc_rand, na.rm = TRUE),
  proc_rand_se   = sd(proc_rand, na.rm = TRUE) / sqrt(.N),
  n = .N
), by = .(k, condition)]
fwrite(summary_dt, out_summary, sep = "\t")

norm_dt <- merge(obs_dt, summary_dt, by = c("k", "condition"), all = TRUE)
norm_dt[, proc_norm_raw := (proc_obs - proc_rand_mean) / (1 - proc_rand_mean)]
norm_dt[, proc_norm := pmax(0, pmin(1, proc_norm_raw))]
norm_dt[, condition_label := label_map[condition]]
norm_dt[, condition_label := factor(condition_label, levels = label_map[target_conds])]
fwrite(norm_dt, out_norm, sep = "\t")

overall <- norm_dt[, .(
  proc_obs_mean = mean(proc_obs, na.rm = TRUE),
  proc_rand_mean_of_means = mean(proc_rand_mean, na.rm = TRUE)
), by = k]
overall[, proc_norm_overall_raw := (proc_obs_mean - proc_rand_mean_of_means) / (1 - proc_rand_mean_of_means)]
overall[, proc_norm_overall := pmax(0, pmin(1, proc_norm_overall_raw))]
fwrite(overall, out_norm_overall, sep = "\t")

message("Done.")
message("Random per-iter table:  ", out_rand)
message("Observed baseline table:", out_obs)
message("Summary table:          ", out_summary)
message("Normalized table:       ", out_norm)
message("Overall normalized:     ", out_norm_overall)
