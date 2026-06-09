#!/usr/bin/env Rscript
# compute-mean-effect-pca.R
#
# PCA of traits from the signed mean-effect matrices (50k and 100k window sizes).
# Produces trait PC scores used by run-figureS10.R for within/between coherence analysis.
#
# Outputs:
#   results/pca_mean_effect/trait_scores_50k_mean.tsv
#   results/pca_mean_effect/trait_scores_100k_mean.tsv
#
# Run from project root: Rscript figure-scripts/compute-mean-effect-pca.R

suppressPackageStartupMessages({
  library(data.table)
  library(irlba)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)
source(file.path("figure-scripts", "helpers-50k-matrix.R"))

out_dir <- file.path("results", "pca_mean_effect")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

N_PCS <- 20L

run_mean_pca <- function(size_label) {
  eff_file  <- file.path("results", sprintf("all-trait-%s-mean-effects.tsv",  size_label))
  pval_file <- file.path("results", sprintf("all-trait-%s-mean-pvals.tsv",    size_label))
  if (!file.exists(eff_file))  stop("Missing: ", eff_file)
  if (!file.exists(pval_file)) stop("Missing: ", pval_file)

  message("Loading ", size_label, " mean-effect matrix...")
  eff  <- fread(eff_file)
  pval <- fread(pval_file)

  prep <- prepare_X_50k(eff, pval)
  X    <- prep$X   # windows x traits

  message("Running PCA (n=", N_PCS, ") on ", nrow(X), " windows x ", ncol(X), " traits...")
  set.seed(42)
  pca <- prcomp_irlba(X, n = N_PCS, center = TRUE, scale. = FALSE)
  rownames(pca$x) <- rownames(X)

  scores <- as.data.table(pca$x, keep.rownames = "trait_abbrev")
  setnames(scores, c("trait_abbrev", paste0("PC", seq_len(N_PCS))))

  out_file <- file.path(out_dir, sprintf("trait_scores_%s_mean.tsv", size_label))
  fwrite(scores, out_file, sep = "\t")
  message("Wrote: ", out_file)
}

run_mean_pca("50k")
run_mean_pca("100k")
message("compute-mean-effect-pca complete.")
