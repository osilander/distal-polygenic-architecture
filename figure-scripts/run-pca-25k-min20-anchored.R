#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(irlba)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

in_dir <- file.path("data", "windows_filtered")
trait_map_file <- file.path("data", "trait_abbrevs_categorised.txt")
out_dir <- file.path("results", "pca_multiscale_anchored")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_effect <- file.path("results", "all-trait-25k-mean-effects-snp20.tsv")
out_scores <- file.path(out_dir, "pca_traits_25k_anchored.tsv")
out_loadings <- file.path(out_dir, "pca_windows_25k_anchored.tsv")
out_scree <- file.path(out_dir, "pca_scree_25k_anchored.tsv")
out_flips <- file.path(out_dir, "pca_anchor_flips_25k_anchored.tsv")
out_meta <- file.path(out_dir, "pca_input_25k_snp20_meta.tsv")

impute_number <- 5L
min_snps <- 20L
anchor_map <- c(PC1 = "BMI_03", PC2 = "StandHgt", PC3 = "Neurotic_30", PC4 = "DiaBP_A59")
n_pcs <- 80L

message("25k PCA (min SNPs): loading trait map...")
trait_map <- fread(trait_map_file)
stopifnot(all(c("study_name", "trait_abbrev") %in% names(trait_map)))
valid_traits <- unique(trait_map$trait_abbrev)

files <- list.files(in_dir, pattern = "_w25000\\.summary\\.tsv$", full.names = TRUE)
if (length(files) == 0L) stop("No *_w25000.summary.tsv files found in ", in_dir)

meta <- data.table(file = files, base = basename(files))
meta[, trait_abbrev := sub("_w25000\\.summary\\.tsv$", "", base)]
meta <- meta[trait_abbrev %in% valid_traits]
if (nrow(meta) == 0L) stop("No mapped 25k trait files after trait-map filter.")
setorder(meta, trait_abbrev)
message("25k PCA: using ", nrow(meta), " mapped traits.")

# Use first file as canonical 25k grid coordinates.
ref <- fread(meta$file[1], select = c("chrom", "start", "end"))
ref[, `:=`(chrom = as.character(chrom), start = as.integer(start), end = as.integer(end))]
coords <- unique(ref[, .(chrom, start, end)])
setkey(coords, chrom, start, end)

E <- matrix(NA_real_, nrow = nrow(coords), ncol = nrow(meta))
colnames(E) <- meta$trait_abbrev

for (i in seq_len(nrow(meta))) {
  tr <- meta$trait_abbrev[i]
  dt <- fread(meta$file[i], na.strings = c("NA", ".", "NaN", "nan", ""))
  need <- c("chrom", "start", "end", "median_es_signed", "count_snps")
  miss <- setdiff(need, names(dt))
  if (length(miss) > 0L) stop("Missing columns in ", basename(meta$file[i]), ": ", paste(miss, collapse = ", "))

  dt[, `:=`(
    chrom = as.character(chrom),
    start = as.integer(start),
    end = as.integer(end),
    median_es_signed = as.numeric(median_es_signed),
    count_snps = as.numeric(count_snps)
  )]

  # Apply per-trait row filter requested for 25k matrix.
  dt <- dt[is.finite(median_es_signed) & is.finite(count_snps) & count_snps >= min_snps]

  merged <- merge(
    coords,
    dt[, .(chrom, start, end, median_es_signed)],
    by = c("chrom", "start", "end"),
    all.x = TRUE,
    sort = FALSE
  )
  E[, i] <- merged$median_es_signed

  if (i == 1L || i %% 100L == 0L || i == nrow(meta)) {
    message("  loaded ", i, "/", nrow(meta), " traits (", tr, ")")
  }
}

eff_dt <- cbind(copy(coords), as.data.table(E))
setnames(eff_dt, old = names(eff_dt)[-(1:3)], new = meta$trait_abbrev)
fwrite(eff_dt, out_effect, sep = "\t", na = "NA")
message("25k matrix written: ", out_effect)

# PCA preprocessing (same structure used in existing multiscale/figure scripts).
target_traits <- intersect(names(eff_dt), valid_traits)
target_traits <- setdiff(target_traits, c("chrom", "start", "end"))
betas_subset <- eff_dt[, c("chrom", "start", "end", target_traits), with = FALSE]

betas_subset <- betas_subset[
  !(
    (chrom == "17" & start < 45000000  & end > 43500000) |
      (chrom == "6"  & start < 34000000 & end > 25000000) |
      (chrom == "12" & start < 113400000 & end > 111200000)
  )
]

na_counts <- rowSums(is.na(betas_subset[, ..target_traits]))
betas_subset <- betas_subset[na_counts <= impute_number]

for (cc in target_traits) {
  med <- median(betas_subset[[cc]], na.rm = TRUE)
  betas_subset[[cc]][is.na(betas_subset[[cc]])] <- med
}

win_ids <- paste(
  betas_subset$chrom,
  formatC(as.integer(betas_subset$start), format = "f", digits = 0, big.mark = ""),
  formatC(as.integer(betas_subset$end), format = "f", digits = 0, big.mark = ""),
  sep = "_"
)

trait_mat <- as.data.frame(betas_subset[, ..target_traits])
rownames(trait_mat) <- win_ids
traits_t <- t(trait_mat) # traits x windows
traits_t_scaled <- t(scale(t(traits_t), center = TRUE, scale = TRUE))
traits_t_scaled[!is.finite(traits_t_scaled)] <- 0

message("25k PCA: running prcomp_irlba...")
pca <- prcomp_irlba(traits_t_scaled, n = n_pcs, center = FALSE, scale. = FALSE)
rownames(pca$x) <- rownames(traits_t_scaled)
rownames(pca$rotation) <- colnames(traits_t_scaled)

pc_cols <- intersect(paste0("PC", seq_len(n_pcs)), colnames(pca$x))

scores <- as.data.table(pca$x)[, ..pc_cols]
scores[, trait_abbrev := rownames(pca$x)]
setcolorder(scores, c("trait_abbrev", pc_cols))

loadings <- as.data.table(pca$rotation)[, ..pc_cols]
loadings[, window_id := rownames(pca$rotation)]
setcolorder(loadings, c("window_id", pc_cols))

flip_rows <- list()
for (pc in intersect(names(anchor_map), pc_cols)) {
  tr <- anchor_map[[pc]]
  idx <- which(scores$trait_abbrev == tr)
  if (length(idx) == 0L) {
    flip_rows[[pc]] <- data.table(pc = pc, anchor_trait = tr, flip_applied = NA_integer_, anchor_value = NA_real_)
    next
  }
  aval <- scores[[pc]][idx[1]]
  flip <- if (is.finite(aval) && aval < 0) -1 else 1
  scores[[pc]] <- scores[[pc]] * flip
  loadings[[pc]] <- loadings[[pc]] * flip
  flip_rows[[pc]] <- data.table(pc = pc, anchor_trait = tr, flip_applied = as.integer(flip == -1), anchor_value = aval)
}
flips <- rbindlist(flip_rows, use.names = TRUE, fill = TRUE)

var_expl <- (pca$sdev)^2
scree <- data.table(
  PC = seq_along(var_expl),
  var_explained = var_expl / sum(var_expl),
  var_explained_pct = (var_expl / sum(var_expl)) * 100,
  cumulative = cumsum(var_expl) / sum(var_expl)
)

meta_out <- data.table(
  min_snps = min_snps,
  impute_number = impute_number,
  n_traits = nrow(scores),
  n_windows_after_filters = nrow(loadings),
  n_windows_before_filters = nrow(eff_dt)
)

fwrite(scores, out_scores, sep = "\t")
fwrite(loadings, out_loadings, sep = "\t")
fwrite(scree, out_scree, sep = "\t")
fwrite(flips, out_flips, sep = "\t")
fwrite(meta_out, out_meta, sep = "\t")

message("Done.")
message("Trait scores: ", out_scores)
message("Window loadings: ", out_loadings)
message("Scree: ", out_scree)
message("Anchor flips: ", out_flips)
message("Input/filter meta: ", out_meta)
