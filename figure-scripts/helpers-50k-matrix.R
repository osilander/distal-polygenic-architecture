suppressPackageStartupMessages({
  library(data.table)
})

filter_effect_matrix <- function(eff, pvals, p_thresh = NULL) {
  stopifnot(all(dim(eff) == dim(pvals)))
  trait_cols <- 4:ncol(eff)
  eff_mat <- as.matrix(eff[, ..trait_cols])
  pval_mat <- as.matrix(pvals[, ..trait_cols])
  mask <- matrix(FALSE, nrow = nrow(eff_mat), ncol = ncol(eff_mat))
  if (!is.null(p_thresh)) {
    p_mask <- pval_mat > -log10(p_thresh)
    mask <- mask | p_mask
  }
  eff_mat[mask] <- NA_real_
  out <- copy(eff)
  out[, (trait_cols) := as.data.table(eff_mat)]
  out
}

prepare_X_50k <- function(eff, pvals, impute_number = 5L, p_thresh = NULL) {
  betas_filtered <- filter_effect_matrix(eff, pvals, p_thresh = p_thresh)
  trait_names <- colnames(betas_filtered)[-(1:3)]
  betas_subset <- copy(betas_filtered[, c("chrom", "start", "end", trait_names), with = FALSE])

  betas_subset <- betas_subset[
    !(
      (chrom == "17" & start < 45000000  & end > 43500000) |
      (chrom == "6"  & start < 34000000  & end > 25000000) |
      (chrom == "12" & start < 113400000 & end > 111200000)
    )
  ]

  window_id <- paste(
    betas_subset$chrom,
    formatC(as.integer(betas_subset$start), format = "f", digits = 0, big.mark = ""),
    formatC(as.integer(betas_subset$end), format = "f", digits = 0, big.mark = ""),
    sep = "_"
  )
  stopifnot(length(window_id) == nrow(betas_subset), uniqueN(window_id) == nrow(betas_subset))

  trait_dt <- betas_subset[, ..trait_names]
  keep <- rowSums(is.na(trait_dt)) <= impute_number
  trait_dt <- trait_dt[keep]
  coords <- betas_subset[keep, .(chrom, start, end)]
  window_id <- window_id[keep]

  for (cc in trait_names) {
    med <- median(trait_dt[[cc]], na.rm = TRUE)
    trait_dt[[cc]][is.na(trait_dt[[cc]])] <- med
  }
  stopifnot(!anyNA(trait_dt))

  X <- as.matrix(trait_dt) # windows x traits
  rownames(X) <- window_id
  colnames(X) <- trait_names
  X_scaled <- scale(X, center = TRUE, scale = TRUE)
  X_scaled <- as.matrix(X_scaled)
  rownames(X_scaled) <- window_id
  colnames(X_scaled) <- trait_names

  list(
    X = X_scaled,
    X_scaled = X_scaled,
    coords = coords,
    trait_names = trait_names,
    window_id = window_id,
    n_windows_before_na_filter = nrow(betas_subset),
    n_windows_after_na_filter = nrow(X_scaled)
  )
}

decomp_from_matrix <- function(X, eig_tol = 1e-12) {
  n <- nrow(X)
  C <- crossprod(X) / max(1, n - 1)
  eg <- eigen(C, symmetric = TRUE)
  vals <- pmax(eg$values, 0)
  keep <- which(vals > eig_tol)
  vals <- vals[keep]
  V <- eg$vectors[, keep, drop = FALSE]
  sdev <- sqrt(vals)
  d <- sdev * sqrt(max(1, n - 1))
  U <- X %*% V
  U <- sweep(U, 2, d, `/`)
  list(U = U, V = V, d = d, var = vals / sum(vals))
}

permute_within_chunks <- function(X, coords, chunk_size, seed) {
  set.seed(seed)
  chunk_id <- paste0(coords$chrom, "_", floor(as.integer(coords$start) / as.integer(chunk_size)))
  permute_within_groups(X, chunk_id, seed = NULL)
}

permute_within_groups <- function(X, group_id, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  idx_list <- split(seq_len(nrow(X)), group_id)
  Xp <- X
  for (j in seq_len(ncol(Xp))) {
    xv <- Xp[, j]
    for (idx in idx_list) {
      if (length(idx) > 1L) xv[idx] <- xv[sample(idx, length(idx), replace = FALSE)]
    }
    Xp[, j] <- xv
  }
  Xp
}

mean_cos_principal_angles <- function(Uobs, Urand, k) {
  kk <- min(k, ncol(Uobs), ncol(Urand))
  A <- crossprod(Uobs[, seq_len(kk), drop = FALSE], Urand[, seq_len(kk), drop = FALSE])
  cs <- svd(A, nu = 0, nv = 0)$d
  cs <- pmin(1, pmax(0, cs))
  mean(cs)
}

trait_abs_cos_upper <- function(V, d, k) {
  kk <- min(k, ncol(V), length(d))
  T <- sweep(V[, seq_len(kk), drop = FALSE], 2, d[seq_len(kk)], `*`)
  rn <- sqrt(rowSums(T^2))
  rn[rn == 0] <- NA_real_
  T <- T / rn
  K <- abs(T %*% t(T))
  as.numeric(K[upper.tri(K, diag = FALSE)])
}

effective_rank <- function(d) {
  p <- (d^2) / sum(d^2)
  exp(-sum(p * log(p + 1e-300)))
}
