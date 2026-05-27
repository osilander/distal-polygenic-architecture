#!/usr/bin/env Rscript
# Produces: results/allbyall_cosine_matrices/allbyall_cosine_long_pc1_pc4.tsv
#           results/allbyall_cosine_matrices/allbyall_procrustes_long_k4_k8_k25_k50_k75.tsv
#           results/allbyall_cosine_matrices/matrix_PC{1..4}_{cosine,procrustes}*.tsv
# Pre-computed inputs:
#   results/all-trait-50k-mean-{effects,pvals}.tsv
#   results/pca_multiscale_anchored/pca_traits_{25k,100k,200k,500k,1m}_anchored.tsv
#   results/svd-nulls-50k/withinblock_perm/withinblock_perm_iter001_50k_trait_vectors.tsv.gz
#   results/gwas_removed_distance_topload_50k_plus{100,150}kb/reduced_pca_trait_scores.tsv
#   results/figureS_pca_absolute_50k_trait_scores.tsv
#   results/pca_abs50k_from_windows/trait_scores.tsv
# Run: Rscript post-process-scripts/run-allbyall-cosine-matrices.R

suppressPackageStartupMessages({
  library(data.table)
  library(irlba)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)
source(file.path("post-process-scripts", "helpers-50k-matrix.R"))

out_dir <- file.path("results", "allbyall_cosine_matrices")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

anchor_map <- c(PC1 = "BMI_03", PC2 = "StandHgt", PC3 = "Neurotic_30", PC4 = "DiaBP_A59")
pcs <- paste0("PC", 1:4)
pcs_all <- paste0("PC", 1:75)
proc_k_values <- c(4L, 8L, 25L, 50L, 75L)

load_trait_scores <- function(path, trait_col_candidates = c("trait_abbrev", "trait", "Trait")) {
  if (!file.exists(path)) stop("Missing file: ", path)
  dt <- tryCatch(
    fread(path),
    error = function(e) fread(cmd = paste("gzip -dc", shQuote(path)))
  )
  tr_col <- trait_col_candidates[trait_col_candidates %in% names(dt)][1]
  if (is.na(tr_col)) stop("No trait column found in ", path)
  setnames(dt, tr_col, "trait_abbrev")
  miss <- setdiff(pcs, names(dt))
  if (length(miss) > 0L) stop("Missing PCs in ", path, ": ", paste(miss, collapse = ", "))
  keep <- c("trait_abbrev", intersect(pcs_all, names(dt)))
  out <- dt[, ..keep]
  for (pc in setdiff(names(out), "trait_abbrev")) set(out, j = pc, value = as.numeric(out[[pc]]))
  out
}

align_anchor <- function(dt) {
  x <- copy(dt)
  for (pc in pcs) {
    tr <- anchor_map[[pc]]
    idx <- which(x$trait_abbrev == tr)
    if (length(idx) == 0L) next
    aval <- x[[pc]][idx[1]]
    if (is.finite(aval) && aval < 0) x[[pc]] <- -x[[pc]]
  }
  x
}

cosine <- function(a, b) {
  na <- sqrt(sum(a * a)); nb <- sqrt(sum(b * b))
  if (!is.finite(na) || !is.finite(nb) || na == 0 || nb == 0) return(NA_real_)
  sum(a * b) / (na * nb)
}

procrustes_similarity <- function(A, B) {
  if (!all(dim(A) == dim(B))) return(NA_real_)
  if (nrow(A) == 0L || ncol(A) == 0L) return(NA_real_)
  sa <- sqrt(sum(A * A)); sb <- sqrt(sum(B * B))
  if (!is.finite(sa) || !is.finite(sb) || sa == 0 || sb == 0) return(NA_real_)
  d <- svd(crossprod(A, B), nu = 0, nv = 0)$d
  sum(d) / (sa * sb)
}

message("All-by-all cosine: loading 50k matrix...")
eff   <- fread(file.path("results", "all-trait-50k-mean-effects.tsv"))
pvals <- fread(file.path("results", "all-trait-50k-mean-pvals.tsv"))
prep  <- prepare_X_50k(eff, pvals)
X     <- t(prep$X)  # traits x windows
trait_ids <- rownames(X)
win_ids   <- colnames(X)

# Observed 50k signed trait scores.
p_obs <- prcomp_irlba(X, n = 80, center = FALSE, scale. = FALSE)
rownames(p_obs$x) <- rownames(X)
obs <- as.data.table(p_obs$x)
obs[, trait_abbrev := rownames(p_obs$x)]
obs_keep <- c("trait_abbrev", intersect(pcs_all, names(obs)))
obs <- obs[, ..obs_keep]

# 500kb within-chunk permutation on observed matrix (single deterministic run).
coords <- tstrsplit(win_ids, "_", fixed = TRUE)
coords_dt <- data.table(
  chrom = as.integer(coords[[1]]),
  start = as.integer(coords[[2]]),
  end   = as.integer(coords[[3]])
)
coords_dt[, chunk_id := paste0(chrom, "_", floor(start / 500000L))]
idx_list <- split(seq_len(ncol(X)), coords_dt$chunk_id)
set.seed(51001)
Xp <- X
for (r in seq_len(nrow(Xp))) {
  xv <- Xp[r, ]
  for (idx in idx_list) {
    if (length(idx) > 1L) xv[idx] <- xv[sample(idx, length(idx), replace = FALSE)]
  }
  Xp[r, ] <- xv
}
p_500perm <- prcomp_irlba(Xp, n = 80, center = FALSE, scale. = FALSE)
rownames(p_500perm$x) <- rownames(X)
perm500 <- as.data.table(p_500perm$x)
perm500[, trait_abbrev := rownames(p_500perm$x)]
perm_keep <- c("trait_abbrev", intersect(pcs_all, names(perm500)))
perm500 <- perm500[, ..perm_keep]

# Entry-wise sign randomization on observed matrix (independent flips per entry).
set.seed(51002)
sign_mat <- matrix(sample(c(-1, 1), length(X), replace = TRUE), nrow = nrow(X), ncol = ncol(X))
Xe <- X * sign_mat
p_entry <- prcomp_irlba(Xe, n = 80, center = FALSE, scale. = FALSE)
rownames(p_entry$x) <- rownames(X)
entry_rand <- as.data.table(p_entry$x)
entry_rand[, trait_abbrev := rownames(p_entry$x)]
entry_keep <- c("trait_abbrev", intersect(pcs_all, names(entry_rand)))
entry_rand <- entry_rand[, ..entry_keep]

cond <- list(
  windows_25kb = load_trait_scores(file.path("results", "pca_multiscale_anchored", "pca_traits_25k_anchored.tsv")),
  observed_50k_signed = obs,
  gwas_removed_plus100kb = load_trait_scores(file.path("results", "gwas_removed_distance_topload_50k_plus100kb", "reduced_pca_trait_scores.tsv")),
  gwas_removed_plus150kb  = load_trait_scores(file.path("results", "gwas_removed_distance_topload_50k_plus150kb/reduced_pca_trait_scores.tsv")),
  pickrell_randomized = load_trait_scores(file.path("results", "svd-nulls-50k", "withinblock_perm", "withinblock_perm_iter001_50k_trait_vectors.tsv.gz")),
  entry_sign_randomized = entry_rand,
  permuted_within_500kb = perm500,
  windows_100kb = load_trait_scores(file.path("results", "pca_multiscale_anchored", "pca_traits_100k_anchored.tsv")),
  windows_200kb = load_trait_scores(file.path("results", "pca_multiscale_anchored", "pca_traits_200k_anchored.tsv")),
  windows_500kb = load_trait_scores(file.path("results", "pca_multiscale_anchored", "pca_traits_500k_anchored.tsv")),
  windows_1mb   = load_trait_scores(file.path("results", "pca_multiscale_anchored", "pca_traits_1m_anchored.tsv")),
  abs_from_signed_matrix   = load_trait_scores(file.path("results", "figureS_pca_absolute_50k_trait_scores.tsv")),
  abs_from_windows_filtered = load_trait_scores(file.path("results", "pca_abs50k_from_windows", "trait_scores.tsv"))
)

cond <- lapply(cond, align_anchor)
conds <- names(cond)

long_rows <- list()
ii <- 1L
for (pc in pcs) {
  for (a in conds) {
    for (b in conds) {
      da <- cond[[a]][, .(trait_abbrev, v = get(pc))]
      db <- cond[[b]][, .(trait_abbrev, v = get(pc))]
      m <- merge(da, db, by = "trait_abbrev", suffixes = c("_a", "_b"))
      c_signed <- cosine(m$v_a, m$v_b)
      long_rows[[ii]] <- data.table(
        pc = pc,
        condition_a = a,
        condition_b = b,
        n_traits_overlap = nrow(m),
        cosine_signed = c_signed,
        cosine_abs = abs(c_signed),
        cosine_distance_signed = 1 - c_signed,
        cosine_distance_abs = 1 - abs(c_signed)
      )
      ii <- ii + 1L
    }
  }
}
long_dt <- rbindlist(long_rows, use.names = TRUE, fill = TRUE)
fwrite(long_dt, file.path(out_dir, "allbyall_cosine_long_pc1_pc4.tsv"), sep = "\t")

# Procrustes similarities on top-k subspaces in trait-score space.
proc_rows <- list()
jj <- 1L
for (k in proc_k_values) {
  for (a in conds) {
    for (b in conds) {
      common_pcs <- intersect(
        setdiff(names(cond[[a]]), "trait_abbrev"),
        setdiff(names(cond[[b]]), "trait_abbrev")
      )
      common_pcs <- common_pcs[grepl("^PC[0-9]+$", common_pcs)]
      common_pcs <- common_pcs[order(as.integer(sub("^PC", "", common_pcs)))]
      if (length(common_pcs) < k) {
        da <- cond[[a]][, .(trait_abbrev)]
        db <- cond[[b]][, .(trait_abbrev)]
        m <- merge(da, db, by = "trait_abbrev")
        ps <- NA_real_
        n_used <- length(common_pcs)
      } else {
        kpcs <- common_pcs[seq_len(k)]
        da <- cond[[a]][, c("trait_abbrev", kpcs), with = FALSE]
        db <- cond[[b]][, c("trait_abbrev", kpcs), with = FALSE]
        m <- merge(da, db, by = "trait_abbrev", suffixes = c("_a", "_b"))
        A <- as.matrix(m[, paste0(kpcs, "_a"), with = FALSE])
        B <- as.matrix(m[, paste0(kpcs, "_b"), with = FALSE])
        ps <- procrustes_similarity(A, B)
        if (is.finite(ps)) ps <- max(0, min(1, ps))
        n_used <- k
      }
      proc_rows[[jj]] <- data.table(
        k = k,
        condition_a = a,
        condition_b = b,
        n_traits_overlap = nrow(m),
        n_pc_used = n_used,
        procrustes_similarity = ps,
        procrustes_distance = 1 - ps
      )
      jj <- jj + 1L
    }
  }
}
proc_dt <- rbindlist(proc_rows, use.names = TRUE, fill = TRUE)
fwrite(proc_dt, file.path(out_dir, "allbyall_procrustes_long_k4_k8_k25_k50_k75.tsv"), sep = "\t")

# Write matrix tables per PC and metric.
for (pcv in pcs) {
  for (metric in c("cosine_signed", "cosine_abs", "cosine_distance_signed", "cosine_distance_abs")) {
    sub <- long_dt[pc == pcv, .(condition_a, condition_b, value = get(metric))]
    mat <- dcast(sub, condition_a ~ condition_b, value.var = "value")
    setorder(mat, condition_a)
    fwrite(mat, file.path(out_dir, paste0("matrix_", pcv, "_", metric, ".tsv")), sep = "\t")
  }
}

# Write matrix tables per k and Procrustes metric.
for (kv in sort(unique(proc_dt$k))) {
  for (metric in c("procrustes_similarity", "procrustes_distance")) {
    sub <- proc_dt[k == kv, .(condition_a, condition_b, value = get(metric))]
    mat <- dcast(sub, condition_a ~ condition_b, value.var = "value")
    setorder(mat, condition_a)
    fwrite(mat, file.path(out_dir, paste0("matrix_k", kv, "_", metric, ".tsv")), sep = "\t")
  }
}

message("All-by-all cosine analysis complete.")
message("Outputs: ", out_dir)
