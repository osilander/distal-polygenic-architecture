#!/usr/bin/env Rscript
# Produces: results/all-trait-{50k,100k,500k,1m,...}-mean-{effects,pvals,rand-mean-effects}.tsv
#   for every window size found in data/windows_filtered/.
#   (25k is also built here; run-pca-25k-min20-anchored.R separately builds the
#    SNP-filtered 25k matrix all-trait-25k-mean-effects-snp20.tsv.)
# Pre-computed inputs: data/windows_filtered/*_w*.summary.tsv (via data/ symlink)
# Run: Rscript post-process-scripts/build-all-trait-matrices.R
# Runtime: ~10-30 min depending on number of window sizes

suppressPackageStartupMessages({
  library(data.table)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

in_dir      <- file.path("data", "windows_filtered")
results_dir <- file.path("results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(in_dir)) stop("Missing: ", in_dir)

size_to_label <- function(w) {
  if (w %% 1000000L == 0L) paste0(w %/% 1000000L, "m")
  else if (w %% 1000L == 0L) paste0(w %/% 1000L, "k")
  else as.character(w)
}

files <- list.files(in_dir, pattern = "_w[0-9]+\\.summary\\.tsv$", full.names = TRUE)
if (length(files) == 0L) stop("No window summary files found in: ", in_dir)

meta <- data.table(file = files, base = basename(files))
meta[, trait_abbrev := sub("_w[0-9]+\\.summary\\.tsv$", "", base)]
meta[, wbp := as.integer(sub("^.*_w([0-9]+)\\.summary\\.tsv$", "\\1", base))]

window_sizes <- sort(unique(meta$wbp))
message("Detected window sizes: ", paste(sapply(window_sizes, size_to_label), collapse = ", "))

required_cols <- c("chrom", "start", "end", "median_es_signed", "max_lp", "median_es_randsign")

for (w in window_sizes) {
  lab <- size_to_label(w)
  effect_out <- file.path(results_dir, sprintf("all-trait-%s-mean-effects.tsv", lab))
  pval_out   <- file.path(results_dir, sprintf("all-trait-%s-mean-pvals.tsv", lab))
  rand_out   <- file.path(results_dir, sprintf("all-trait-%s-rand-mean-effects.tsv", lab))

  sub <- meta[wbp == w]
  setorder(sub, trait_abbrev)
  n_trait <- nrow(sub)
  message("\nProcessing ", lab, " with ", n_trait, " traits")

  ref <- fread(sub$file[1], select = c("chrom", "start", "end"),
               na.strings = c("NA", ".", "NaN", "nan", ""))
  coords <- unique(ref[, .(chrom = as.character(chrom),
                            start = as.integer(start),
                            end = as.integer(end))])
  setkey(coords, chrom, start, end)
  n_win <- nrow(coords)

  E <- matrix(NA_real_, nrow = n_win, ncol = n_trait)
  P <- matrix(NA_real_, nrow = n_win, ncol = n_trait)
  R <- matrix(NA_real_, nrow = n_win, ncol = n_trait)

  for (i in seq_len(n_trait)) {
    tr <- sub$trait_abbrev[i]
    dt <- fread(sub$file[i], na.strings = c("NA", ".", "NaN", "nan", ""))
    miss <- setdiff(required_cols, names(dt))
    if (length(miss) > 0L) {
      warning("Skipping ", tr, " (missing columns: ", paste(miss, collapse = ", "), ")")
      next
    }
    dt[, `:=`(chrom = as.character(chrom), start = as.integer(start), end = as.integer(end))]
    for (cc in c("median_es_signed", "max_lp", "median_es_randsign")) {
      dt[[cc]] <- as.numeric(dt[[cc]])
    }
    merged <- merge(coords,
                    dt[, .(chrom, start, end, median_es_signed, max_lp, median_es_randsign)],
                    by = c("chrom", "start", "end"), all.x = TRUE, sort = FALSE)
    E[, i] <- merged$median_es_signed
    P[, i] <- merged$max_lp
    R[, i] <- merged$median_es_randsign
    if (i == 1L || i %% 100L == 0L || i == n_trait) {
      message("  ", lab, ": ", i, "/", n_trait, " (", tr, ")")
    }
  }

  eff_dt <- cbind(copy(coords), as.data.table(E))
  pvl_dt <- cbind(copy(coords), as.data.table(P))
  rnd_dt <- cbind(copy(coords), as.data.table(R))
  setnames(eff_dt, old = names(eff_dt)[-(1:3)], new = sub$trait_abbrev)
  setnames(pvl_dt, old = names(pvl_dt)[-(1:3)], new = sub$trait_abbrev)
  setnames(rnd_dt, old = names(rnd_dt)[-(1:3)], new = sub$trait_abbrev)

  fwrite(eff_dt, effect_out, sep = "\t", na = "NA")
  fwrite(pvl_dt, pval_out,   sep = "\t", na = "NA")
  fwrite(rnd_dt, rand_out,   sep = "\t", na = "NA")
  message("  wrote: ", basename(effect_out), ", ", basename(pval_out), ", ", basename(rand_out))

  if (w == 50000L) {
    alias <- file.path(results_dir, "all-trait-rand-mean-effects.tsv")
    fwrite(rnd_dt, alias, sep = "\t", na = "NA")
    message("  wrote alias: all-trait-rand-mean-effects.tsv")
  }
}

message("\nDone.")
