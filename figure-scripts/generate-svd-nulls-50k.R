#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

setDTthreads(1L)
args <- commandArgs(trailingOnly = TRUE)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
effects_file <- file.path(project_root, "results", "all-trait-50k-mean-effects.tsv")
pvals_file <- file.path(project_root, "results", "all-trait-50k-mean-pvals.tsv")
pickrell_file <- file.path(project_root, "data", "pickrell_blocks.bed")

out_root <- file.path(project_root, "results", "svd-nulls-50k")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

n_perm <- if (length(args) >= 1L) as.integer(args[[1]]) else 100L
if (!is.finite(n_perm) || n_perm < 1L) stop("n_perm must be a positive integer.")

impute_number <- 5L
p_thresh <- NULL
k_window_save <- 30L
save_full_window <- FALSE
eig_tol <- 1e-12
source(file.path(project_root, "figure-scripts", "helpers-50k-matrix.R"))

assign_pickrell_blocks <- function(coords, pickrell_file) {
  stopifnot(file.exists(pickrell_file))
  blocks <- fread(pickrell_file, header = TRUE)
  setnames(blocks, c("chr", "start", "stop"), c("chrom", "block_start", "block_end"))
  blocks[, chrom := gsub("^chr", "", as.character(chrom))]
  blocks[, block_start := as.integer(block_start)]
  blocks[, block_end := as.integer(block_end)]
  setkey(blocks, chrom, block_start, block_end)

  wins <- copy(coords)
  wins[, chrom := gsub("^chr", "", as.character(chrom))]
  wins[, `:=`(start = as.integer(start), end = as.integer(end))]
  wins[, row_i := .I]
  wins[, mid := (start + end) %/% 2L]
  wins_mid <- wins[, .(chrom, start = mid, end = mid, row_i)]
  setkey(wins_mid, chrom, start, end)

  ov <- foverlaps(
    wins_mid, blocks,
    by.x = c("chrom", "start", "end"),
    by.y = c("chrom", "block_start", "block_end"),
    type = "within",
    nomatch = NA_integer_
  )
  setorder(ov, row_i)

  ov[, block_id := .GRP, by = .(chrom, block_start, block_end)]
  max_block <- suppressWarnings(max(ov$block_id, na.rm = TRUE))
  if (!is.finite(max_block)) max_block <- 0L
  if (anyNA(ov$block_id)) {
    ov[is.na(block_id), block_id := max_block + seq_len(sum(is.na(block_id)))]
  }
  stopifnot(nrow(ov) == nrow(coords))
  ov$block_id
}

full_decomp <- function(X, coords, trait_names, k_window_save, save_full_window, eig_tol) {
  n <- nrow(X)
  p <- ncol(X)

  C <- crossprod(X) / max(1, n - 1) # traits x traits
  eg <- eigen(C, symmetric = TRUE)
  vals <- pmax(eg$values, 0)
  V <- eg$vectors

  keep <- which(vals > eig_tol)
  if (length(keep) == 0L) stop("No positive eigenvalues in decomposition.")
  vals <- vals[keep]
  V <- V[, keep, drop = FALSE]

  sdev <- sqrt(vals)
  d <- sdev * sqrt(max(1, n - 1))
  var_exp <- vals / sum(vals)
  cum_var <- cumsum(var_exp)

  k_win <- if (isTRUE(save_full_window)) length(d) else min(k_window_save, length(d))
  V_top <- V[, seq_len(k_win), drop = FALSE]
  d_top <- d[seq_len(k_win)]
  U_top <- X %*% V_top
  U_top <- sweep(U_top, 2, d_top, `/`)

  list(
    singular_value = d,
    sdev = sdev,
    var_explained = var_exp,
    cum_var_explained = cum_var,
    V = V,
    U_top = U_top,
    k_win = k_win,
    coords = coords,
    trait_names = trait_names,
    pc_count = length(d)
  )
}

save_decomp <- function(dec, scheme_dir, run_prefix) {
  dir.create(scheme_dir, recursive = TRUE, showWarnings = FALSE)

  sing_dt <- data.table(
    pc = paste0("PC", seq_len(dec$pc_count)),
    singular_value = dec$singular_value,
    sdev = dec$sdev,
    var_explained = dec$var_explained,
    cum_var_explained = dec$cum_var_explained
  )
  fwrite(sing_dt, file.path(scheme_dir, sprintf("%s_singular_values.tsv", run_prefix)), sep = "\t")

  V_dt <- as.data.table(dec$V)
  n_pc_v <- ncol(V_dt)
  setnames(V_dt, paste0("PC", seq_len(n_pc_v)))
  V_dt[, trait := dec$trait_names]
  setcolorder(V_dt, c("trait", paste0("PC", seq_len(n_pc_v))))
  fwrite(V_dt, file.path(scheme_dir, sprintf("%s_trait_vectors.tsv.gz", run_prefix)), sep = "\t")

  U_dt <- cbind(copy(dec$coords), as.data.table(dec$U_top))
  setnames(U_dt, c("chrom", "start", "end", paste0("PC", seq_len(dec$k_win))))
  U_dt[, window_id := paste(chrom, start, end, sep = "_")]
  setcolorder(U_dt, c("window_id", "chrom", "start", "end", paste0("PC", seq_len(dec$k_win))))
  fwrite(U_dt, file.path(scheme_dir, sprintf("%s_window_vectors_top%d.tsv.gz", run_prefix, dec$k_win)), sep = "\t")
}

message("Loading 50k inputs...")
stopifnot(file.exists(effects_file), file.exists(pvals_file), file.exists(pickrell_file))
eff <- fread(effects_file)
pvals <- fread(pvals_file)

message("Preparing canonical matrix X (windows x traits)...")
prep <- prepare_X_50k(eff, pvals, impute_number = impute_number, p_thresh = p_thresh)
X_obs <- prep$X_scaled
stopifnot(identical(rownames(X_obs), prep$window_id))

message("Assigning Pickrell LD blocks to windows...")
block_id <- assign_pickrell_blocks(prep$coords, pickrell_file)
n_blocks <- uniqueN(block_id)

config <- data.table(
  dataset = "50k",
  n_perm = n_perm,
  n_windows = nrow(X_obs),
  n_traits = ncol(X_obs),
  n_blocks = n_blocks,
  impute_number = impute_number,
  p_thresh = ifelse(is.null(p_thresh), NA_real_, p_thresh),
  k_window_save = k_window_save,
  save_full_window = save_full_window,
  row_id_check_passed = identical(rownames(X_obs), prep$window_id)
)
fwrite(config, file.path(out_root, "config_50k.tsv"), sep = "\t")
fwrite(
  data.table(window_id = prep$window_id, chrom = prep$coords$chrom, start = prep$coords$start, end = prep$coords$end, block_id = block_id),
  file.path(out_root, "window_block_map_50k.tsv.gz"),
  sep = "\t"
)

message("Observed decomposition...")
obs_dir <- file.path(out_root, "observed")
obs_dec <- full_decomp(
  X = X_obs,
  coords = prep$coords,
  trait_names = prep$trait_names,
  k_window_save = k_window_save,
  save_full_window = save_full_window,
  eig_tol = eig_tol
)
save_decomp(obs_dec, obs_dir, "observed_50k")

scheme_dirs <- list(
  trait_flip = file.path(out_root, "trait_flip"),
  entry_flip = file.path(out_root, "entry_flip"),
  withinblock_perm = file.path(out_root, "withinblock_perm")
)
for (d in scheme_dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

run_scheme <- function(scheme_name) {
  message("Starting scheme: ", scheme_name)
  sdir <- scheme_dirs[[scheme_name]]
  for (i in seq_len(n_perm)) {
    if (scheme_name == "trait_flip") {
      set.seed(10000 + i)
      signs <- sample(c(-1, 1), ncol(X_obs), replace = TRUE)
      Xr <- sweep(X_obs, 2, signs, `*`)
    } else if (scheme_name == "entry_flip") {
      set.seed(20000 + i)
      S <- matrix(sample(c(-1, 1), length(X_obs), replace = TRUE), nrow = nrow(X_obs), ncol = ncol(X_obs))
      Xr <- X_obs * S
    } else if (scheme_name == "withinblock_perm") {
      Xr <- permute_within_groups(X_obs, block_id, seed = 30000L + i)
    } else {
      stop("Unknown scheme: ", scheme_name)
    }
    rownames(Xr) <- rownames(X_obs)
    colnames(Xr) <- colnames(X_obs)

    dec <- full_decomp(
      X = Xr,
      coords = prep$coords,
      trait_names = prep$trait_names,
      k_window_save = k_window_save,
      save_full_window = save_full_window,
      eig_tol = eig_tol
    )
    run_prefix <- sprintf("%s_iter%03d_50k", scheme_name, i)
    save_decomp(dec, sdir, run_prefix)

    if (i %% 5L == 0L || i == n_perm) {
      message("  ", scheme_name, ": ", i, "/", n_perm)
    }
  }
}

run_scheme("trait_flip")
run_scheme("entry_flip")
run_scheme("withinblock_perm")

message("Done. Wrote decompositions to: ", out_root)
