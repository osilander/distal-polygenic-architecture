#!/usr/bin/env Rscript
# Produces: results/svd_robustness/
#   trait_subsample_pc_cosine.tsv
#   trait_subsample_procrustes.tsv
#   leave_one_out_observed.tsv      (long format: one row per category x PC)
#   leave_one_out_null.tsv          (long format: one row per category x PC x rep)
#   window_subsample_pc_cosine.tsv
#   window_subsample_procrustes.tsv
#
# Checkpointing: each analysis saves its own intermediate TSV on completion.
# If an intermediate file already exists it is loaded and that analysis is skipped.
# Delete individual files to force re-run of specific analyses.
#
# Run: Rscript post-process-scripts/compute-svd-robustness.R
# Runtime: ~60-90 min at n_perm=100, n_perm_null=30

suppressPackageStartupMessages({
  library(data.table)
  library(irlba)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)
source(file.path("post-process-scripts", "helpers-50k-matrix.R"))

# ---- Parameters ----
set.seed(20260519)
n_perm       <- 100L   # replicates for trait / window subsampling
n_perm_null  <- 30L    # size-matched null replicates per category (analysis 3)
n_irlba      <- 50L    # components in subsampled PCAs
n_irlba_null <- 20L    # components for LOO null runs (only need k<=4)
frac         <- 2/3
k_cos        <- 8L     # PCs for per-PC cosine similarity
k_cos_loo    <- 8L     # PCs for LOO (PC1-PC8)
k_proc       <- c(4L, 8L)

out_dir <- file.path("results", "svd_robustness")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Intermediate file paths (one per analysis for checkpointing)
f_cos_random  <- file.path(out_dir, "_cache_cos_random.tsv")
f_proc_random <- file.path(out_dir, "_cache_proc_random.tsv")
f_cos_strat   <- file.path(out_dir, "_cache_cos_strat.tsv")
f_proc_strat  <- file.path(out_dir, "_cache_proc_strat.tsv")
f_loo_obs     <- file.path(out_dir, "_cache_loo_obs.tsv")
f_loo_null    <- file.path(out_dir, "_cache_loo_null.tsv")
f_cos_win     <- file.path(out_dir, "_cache_cos_win_full.tsv")
f_proc_win    <- file.path(out_dir, "_cache_proc_win_full.tsv")
f_cos_gd      <- file.path(out_dir, "_cache_cos_win_gd.tsv")
f_proc_gd     <- file.path(out_dir, "_cache_proc_win_gd.tsv")

# ---- Helpers ----

procrustes_sim <- function(A, B, k) {
  k  <- min(k, ncol(A), ncol(B))
  A  <- A[, seq_len(k), drop = FALSE]
  B  <- B[, seq_len(k), drop = FALSE]
  sa <- sqrt(sum(A * A)); sb <- sqrt(sum(B * B))
  if (!is.finite(sa) || !is.finite(sb) || sa == 0 || sb == 0) return(NA_real_)
  d  <- svd(crossprod(A, B), nu = 0, nv = 0)$d
  max(0, min(1, sum(d) / (sa * sb)))
}

greedy_cosines <- function(scores_new, scores_orig, k) {
  k  <- min(k, ncol(scores_new), ncol(scores_orig))
  nc <- function(M) { n <- sqrt(colSums(M^2)); n[n == 0] <- 1; sweep(M, 2, n, "/") }
  Sn <- nc(scores_new[,  seq_len(k), drop = FALSE])
  So <- nc(scores_orig[, seq_len(k), drop = FALSE])
  cos_mat    <- crossprod(Sn, So)
  avail_new  <- seq_len(k)
  avail_orig <- seq_len(k)
  matched    <- numeric(k)   # indexed by orig PC
  for (i in seq_len(k)) {
    sub  <- abs(cos_mat[avail_new, avail_orig, drop = FALSE])
    best <- which(sub == max(sub), arr.ind = TRUE)[1, ]
    bn   <- avail_new[best[1]]; bo <- avail_orig[best[2]]
    matched[bo]  <- abs(cos_mat[bn, bo])
    avail_new    <- avail_new[avail_new   != bn]
    avail_orig   <- avail_orig[avail_orig != bo]
  }
  matched
}

run_one <- function(X_new_tw, scores_ref, k_cos, k_proc, n_comp,
                    analysis = NULL, dataset = NULL) {
  pca_r      <- prcomp_irlba(X_new_tw, n = n_comp, center = FALSE, scale. = FALSE)
  scores_new <- pca_r$x
  cosines    <- greedy_cosines(scores_new, scores_ref, k = k_cos)
  cos_dt     <- data.table(PC = paste0("PC", seq_len(k_cos)), cosine = cosines)
  if (!is.null(analysis)) cos_dt[, analysis := analysis]
  if (!is.null(dataset))  cos_dt[, dataset  := dataset]
  proc_vals  <- sapply(k_proc, function(k) procrustes_sim(scores_new, scores_ref, k))
  proc_dt    <- data.table(k = k_proc, procrustes = proc_vals)
  if (!is.null(analysis)) proc_dt[, analysis := analysis]
  if (!is.null(dataset))  proc_dt[, dataset  := dataset]
  list(cos_dt = cos_dt, proc_dt = proc_dt)
}

# ---- Load data (always needed for reference PCA) ----
message("Loading 50k matrix...")
all_betas   <- fread(file.path("results", "all-trait-50k-mean-effects.tsv"))
all_pvals   <- fread(file.path("results", "all-trait-50k-mean-pvals.tsv"))
prep        <- prepare_X_50k(all_betas, all_pvals)
X_wt        <- prep$X
trait_names <- prep$trait_names
n_traits    <- length(trait_names)
n_windows   <- nrow(X_wt)
X_tw        <- t(X_wt)
message("  ", n_traits, " traits x ", n_windows, " windows")

message("Original PCA (n=", n_irlba, ")...")
pca_orig    <- prcomp_irlba(X_tw, n = n_irlba, center = FALSE, scale. = FALSE)
scores_orig <- pca_orig$x

trait_cats <- fread(file.path("data", "trait_abbrevs_categorised.txt"))
cat_col    <- if ("phenotype_group" %in% names(trait_cats)) "phenotype_group" else names(trait_cats)[2]
trait_col  <- names(trait_cats)[1]
trait_cats <- trait_cats[, .(trait_abbrev = get(trait_col), category = get(cat_col))]
trait_cats <- trait_cats[trait_abbrev %in% trait_names & !is.na(category) & nzchar(category)]
cat_map       <- setNames(trait_cats$category, trait_cats$trait_abbrev)
cats_per_idx  <- cat_map[trait_names]
n_sample_trait <- floor(frac * n_traits)
n_sample_win   <- floor(frac * n_windows)

# ============================================================
# Analysis 1: Random trait subsampling
# ============================================================
if (file.exists(f_cos_random) && file.exists(f_proc_random)) {
  message("\nAnalysis 1: Loading cached results...")
  cosine_trait_random <- fread(f_cos_random)
  proc_trait_random   <- fread(f_proc_random)
} else {
  message("\nAnalysis 1: Random trait subsampling (", n_perm, " reps)...")
  cos_list <- vector("list", n_perm)
  pro_list <- vector("list", n_perm)
  for (r in seq_len(n_perm)) {
    if (r %% 20 == 0) message("  rep ", r, "/", n_perm)
    idx    <- sort(sample.int(n_traits, n_sample_trait, replace = FALSE))
    result <- run_one(X_tw[idx, , drop = FALSE], scores_orig[idx, , drop = FALSE],
                      k_cos, k_proc, n_irlba, analysis = "random_trait")
    cos_list[[r]] <- result$cos_dt[, rep := r]
    pro_list[[r]] <- result$proc_dt[, rep := r]
  }
  cosine_trait_random <- rbindlist(cos_list)
  proc_trait_random   <- rbindlist(pro_list)
  fwrite(cosine_trait_random, f_cos_random,  sep = "\t")
  fwrite(proc_trait_random,   f_proc_random, sep = "\t")
  message("  Saved.")
}

# ============================================================
# Analysis 2: Stratified trait subsampling
# ============================================================
if (file.exists(f_cos_strat) && file.exists(f_proc_strat)) {
  message("\nAnalysis 2: Loading cached results...")
  cosine_trait_strat <- fread(f_cos_strat)
  proc_trait_strat   <- fread(f_proc_strat)
} else {
  message("\nAnalysis 2: Stratified trait subsampling (", n_perm, " reps)...")
  cos_list <- vector("list", n_perm)
  pro_list <- vector("list", n_perm)
  for (r in seq_len(n_perm)) {
    if (r %% 20 == 0) message("  rep ", r, "/", n_perm)
    cats_factor <- ifelse(is.na(cats_per_idx), "__none__", cats_per_idx)
    idx <- sort(unlist(lapply(split(seq_len(n_traits), cats_factor), function(g)
      sample(g, max(1L, floor(frac * length(g))), replace = FALSE))))
    result <- run_one(X_tw[idx, , drop = FALSE], scores_orig[idx, , drop = FALSE],
                      k_cos, k_proc, n_irlba, analysis = "stratified_trait")
    cos_list[[r]] <- result$cos_dt[, rep := r]
    pro_list[[r]] <- result$proc_dt[, rep := r]
  }
  cosine_trait_strat <- rbindlist(cos_list)
  proc_trait_strat   <- rbindlist(pro_list)
  fwrite(cosine_trait_strat, f_cos_strat,  sep = "\t")
  fwrite(proc_trait_strat,   f_proc_strat, sep = "\t")
  message("  Saved.")
}

# ============================================================
# Analysis 3: Leave-one-category-out + size-matched null
# Long format: one row per category x PC (x rep for null)
# ============================================================
if (file.exists(f_loo_obs) && file.exists(f_loo_null)) {
  message("\nAnalysis 3: Loading cached results...")
  loo_obs_dt  <- fread(f_loo_obs)
  loo_null_dt <- fread(f_loo_null)
} else {
  message("\nAnalysis 3: Leave-one-category-out...")
  categories    <- sort(unique(trait_cats$category))
  loo_obs_list  <- vector("list", length(categories))
  loo_null_list <- vector("list", length(categories))

  for (ci in seq_along(categories)) {
    cat_i      <- categories[ci]
    remove_idx <- which(trait_names %in% trait_cats[category == cat_i, trait_abbrev])
    keep_idx   <- setdiff(seq_len(n_traits), remove_idx)
    n_removed  <- length(remove_idx)
    message("  [", ci, "/", length(categories), "] ", cat_i,
            "  (remove n=", n_removed, ")")

    result <- run_one(X_tw[keep_idx, , drop = FALSE], scores_orig[keep_idx, , drop = FALSE],
                      k_cos_loo, k_proc, n_irlba, analysis = cat_i)
    obs_cos <- result$cos_dt[, .(category = cat_i, n_removed = n_removed, PC, cosine)]
    obs_cos[, procrustes_k4 := result$proc_dt[k == 4L, procrustes]]
    obs_cos[, procrustes_k8 := result$proc_dt[k == 8L, procrustes]]
    loo_obs_list[[ci]] <- obs_cos

    null_rows <- vector("list", n_perm_null)
    for (r in seq_len(n_perm_null)) {
      null_remove <- sample.int(n_traits, n_removed, replace = FALSE)
      null_keep   <- setdiff(seq_len(n_traits), null_remove)
      res_null    <- run_one(X_tw[null_keep, , drop = FALSE],
                             scores_orig[null_keep, , drop = FALSE],
                             k_cos_loo, k_proc, n_irlba_null, analysis = cat_i)
      nr <- res_null$cos_dt[, .(category = cat_i, n_removed = n_removed, rep = r, PC, cosine)]
      nr[, procrustes_k4 := res_null$proc_dt[k == 4L, procrustes]]
      nr[, procrustes_k8 := res_null$proc_dt[k == 8L, procrustes]]
      null_rows[[r]] <- nr
    }
    loo_null_list[[ci]] <- rbindlist(null_rows)
  }

  loo_obs_dt  <- rbindlist(loo_obs_list)
  loo_null_dt <- rbindlist(loo_null_list)
  fwrite(loo_obs_dt,  f_loo_obs,  sep = "\t")
  fwrite(loo_null_dt, f_loo_null, sep = "\t")
  message("  Saved.")
}

# ============================================================
# Analysis 4: Window subsampling — full matrix
# ============================================================
if (file.exists(f_cos_win) && file.exists(f_proc_win)) {
  message("\nAnalysis 4: Loading cached results...")
  cosine_win_dt <- fread(f_cos_win)
  proc_win_dt   <- fread(f_proc_win)
} else {
  message("\nAnalysis 4: Window subsampling — full matrix (", n_perm, " reps)...")
  cos_list <- vector("list", n_perm)
  pro_list <- vector("list", n_perm)
  for (r in seq_len(n_perm)) {
    if (r %% 20 == 0) message("  rep ", r, "/", n_perm)
    win_idx <- sort(sample.int(n_windows, n_sample_win, replace = FALSE))
    result  <- run_one(X_tw[, win_idx, drop = FALSE], scores_orig,
                       k_cos, k_proc, n_irlba, dataset = "full")
    cos_list[[r]] <- result$cos_dt[, rep := r]
    pro_list[[r]] <- result$proc_dt[, rep := r]
  }
  cosine_win_dt <- rbindlist(cos_list)
  proc_win_dt   <- rbindlist(pro_list)
  fwrite(cosine_win_dt, f_cos_win,  sep = "\t")
  fwrite(proc_win_dt,   f_proc_win, sep = "\t")
  message("  Saved.")
}

# ============================================================
# Analysis 4b: Window subsampling — GWAS-depleted matrix
# ============================================================
removed_file <- file.path("results",
  "gwas_removed_distance_topload_50k_plus100kb",
  "removed_windows_gwas_plus100kb.tsv")

if (file.exists(f_cos_gd) && file.exists(f_proc_gd)) {
  message("\nAnalysis 4b: Loading cached results...")
  cosine_gd_dt <- fread(f_cos_gd)
  proc_gd_dt   <- fread(f_proc_gd)
} else if (!file.exists(removed_file)) {
  message("\nAnalysis 4b: Skipping — removed_windows file not found.")
  message("  Run: Rscript post-process-scripts/compute-gwas-removed-distances.R 100")
  cosine_gd_dt <- data.table()
  proc_gd_dt   <- data.table()
} else {
  message("\nAnalysis 4b: Window subsampling — GWAS-depleted matrix (", n_perm, " reps)...")
  removed_dt  <- fread(removed_file)
  removed_ids <- unique(removed_dt$window_id)
  harmonise   <- function(ids) {
    p <- tstrsplit(ids, "_", fixed = TRUE)
    paste(gsub("^chr", "", p[[1]]), p[[2]], p[[3]], sep = "_")
  }
  full_ids_h    <- harmonise(prep$window_id)
  removed_ids_h <- harmonise(removed_ids)
  retained_idx  <- which(!(full_ids_h %in% removed_ids_h))
  n_retained    <- length(retained_idx)
  message("  Retained windows: ", n_retained, " of ", n_windows)

  anchor_map <- c(PC1 = "BMI_03", PC2 = "StandHgt", PC3 = "Neurotic_30", PC4 = "DiaBP_A59")
  X_gd_tw    <- X_tw[, retained_idx, drop = FALSE]
  pca_gd     <- prcomp_irlba(X_gd_tw, n = n_irlba, center = FALSE, scale. = FALSE)
  scores_gd  <- pca_gd$x
  for (pc_name in names(anchor_map)) {
    tr  <- anchor_map[[pc_name]]
    idx_tr <- which(trait_names == tr)
    if (length(idx_tr) == 0L || !(pc_name %in% colnames(scores_gd))) next
    if (scores_gd[idx_tr[1], pc_name] < 0) scores_gd[, pc_name] <- -scores_gd[, pc_name]
  }

  n_sample_gd <- floor(frac * n_retained)
  message("  Subsampling ", n_sample_gd, " windows per rep...")
  cos_list <- vector("list", n_perm)
  pro_list <- vector("list", n_perm)
  for (r in seq_len(n_perm)) {
    if (r %% 20 == 0) message("  rep ", r, "/", n_perm)
    win_sub <- sort(sample.int(n_retained, n_sample_gd, replace = FALSE))
    result  <- run_one(X_gd_tw[, win_sub, drop = FALSE], scores_gd,
                       k_cos, k_proc, n_irlba, dataset = "gwas_depleted")
    cos_list[[r]] <- result$cos_dt[, rep := r]
    pro_list[[r]] <- result$proc_dt[, rep := r]
  }
  cosine_gd_dt <- rbindlist(cos_list)
  proc_gd_dt   <- rbindlist(pro_list)
  fwrite(cosine_gd_dt, f_cos_gd,  sep = "\t")
  fwrite(proc_gd_dt,   f_proc_gd, sep = "\t")
  message("  Saved.")
}

# ============================================================
# Write final combined outputs
# ============================================================
message("\nWriting final combined outputs...")
fwrite(rbindlist(list(cosine_trait_random, cosine_trait_strat)),
       file.path(out_dir, "trait_subsample_pc_cosine.tsv"), sep = "\t")
fwrite(rbindlist(list(proc_trait_random, proc_trait_strat)),
       file.path(out_dir, "trait_subsample_procrustes.tsv"), sep = "\t")
fwrite(loo_obs_dt,  file.path(out_dir, "leave_one_out_observed.tsv"), sep = "\t")
fwrite(loo_null_dt, file.path(out_dir, "leave_one_out_null.tsv"),     sep = "\t")
win_cos_all  <- rbindlist(list(cosine_win_dt,
  if (nrow(cosine_gd_dt) > 0) cosine_gd_dt else NULL), fill = TRUE)
win_proc_all <- rbindlist(list(proc_win_dt,
  if (nrow(proc_gd_dt) > 0) proc_gd_dt else NULL), fill = TRUE)
fwrite(win_cos_all,  file.path(out_dir, "window_subsample_pc_cosine.tsv"), sep = "\t")
fwrite(win_proc_all, file.path(out_dir, "window_subsample_procrustes.tsv"), sep = "\t")
message("Done.")
