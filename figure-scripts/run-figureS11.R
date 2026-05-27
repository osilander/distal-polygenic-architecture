#!/usr/bin/env Rscript
# Produces: figures/figureS11_multiscale_pca_cosine_50k_1m.pdf
#           results/pca_multiscale_anchored/{pca_traits,pca_windows,pca_scree,pca_anchor_flips}_{50k,100k,200k,500k,1m}_anchored.tsv
#           results/pca_multiscale_anchored/{cosine_heatmap_long_pc1_pc4,procrustes_heatmap_long_k4_k8_k25_k50_k75,pca_variance_audit_50k_1m}.tsv
# Pre-computed inputs:
#   results/all-trait-{50k,100k,200k,500k,1m}-mean-effects.tsv (from build-all-trait-matrices.R)
#   data/trait_abbrevs_categorised.txt
# Run: Rscript post-process-scripts/run-figureS11.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(irlba)
  library(patchwork)
  library(cowplot)
  library(grid)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

scales <- c("50k", "100k", "200k", "500k", "1m")
scale_labels <- c("50k" = "50Kb", "100k" = "100Kb", "200k" = "200Kb", "500k" = "500Kb", "1m" = "1Mb")
effect_files <- c(
  "50k"  = file.path("results", "all-trait-50k-mean-effects.tsv"),
  "100k" = file.path("results", "all-trait-100k-mean-effects.tsv"),
  "200k" = file.path("results", "all-trait-200k-mean-effects.tsv"),
  "500k" = file.path("results", "all-trait-500k-mean-effects.tsv"),
  "1m"   = file.path("results", "all-trait-1m-mean-effects.tsv")
)
anchor_map <- c(PC1 = "BMI_03", PC2 = "StandHgt", PC3 = "Neurotic_30", PC4 = "DiaBP_A59")
impute_number <- 5L
procrustes_k_values <- c(4L, 8L, 25L, 50L, 75L)

out_fig      <- file.path("figures", "figureS11_multiscale_pca_cosine_50k_1m.pdf")
out_dir      <- file.path("results", "pca_multiscale_anchored")
out_audit    <- file.path(out_dir, "pca_variance_audit_50k_1m.tsv")
out_cos_long <- file.path(out_dir, "cosine_heatmap_long_pc1_pc4.tsv")
out_proc_long <- file.path(out_dir, "procrustes_heatmap_long_k4_k8_k25_k50_k75.tsv")
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Trait categories and colors ----
trait_file <- file.path("data", "trait_abbrevs_categorised.txt")
if (!file.exists(trait_file)) stop("Missing: ", trait_file)
trait_categories <- fread(trait_file)
trait_map <- trait_categories[, .(trait_abbrev, phenotype_group)]
target_traits <- unique(trait_map$trait_abbrev)
custom_theme_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
  "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
  "#E5C494", "#B3B3B3"
)
all_trait_categories <- sort(unique(na.omit(as.character(trait_categories$phenotype_group))))
all_trait_categories <- all_trait_categories[nzchar(all_trait_categories)]
category_colors <- setNames(custom_theme_colors[seq_along(all_trait_categories)], all_trait_categories)
fixed_pca_labels <- list(
  pc12 = c("Anxious", "BasMetR", "BirthWt", "BMI_03", "Body10Y",
           "FFMTrunk", "FVitC_05", "HH_Inc", "HighStrg", "HipCirc", "HlthRating",
           "ImpedBody", "ImpedLeg", "Insomnia", "NervFeel",
           "Neurotic_30", "Obesity1", "SitHgt", "SmokeCur_16",
           "StandHgt", "UniverD", "VascHDxN", "Weight"),
  pc23 = c("DiaBP_A59", "EduYr", "FamSatis", "FFM_Leg", "FFMTrunk",
           "FVitC_05", "GripStrR", "Guilt", "HH_Inc", "HighStrg", "HipCirc",
           "HlthRating", "ImpedArm", "ImpedLeg",
           "MoodSwg", "NervFeel", "Neurotic_30", "PExpFl", "SitHgt", "StandHgt",
           "SysBP_A", "UniverD", "VascHDxN", "Weight", "WorryEmb"),
  pc34 = c("Anxious", "DiaBP_A59", "FFMTrunk",
           "FI3_Word", "FluIntSc", "FVitC_05",
           "Guilt", "NervFeel", "Neurotic_30",
           "SitHgt", "StandHgt", "SysBP_A", "Townsend", "ProfQual",
           "RiskTak", "SexPartn", "SmokeCur_25",
           "Speeding", "UniverD", "WorryEmb")
)

prep_matrix <- function(file_path, keep_traits) {
  dt <- fread(file_path)
  stopifnot(all(c("chrom", "start", "end") %in% names(dt)))
  dt[, `:=`(chrom = as.character(chrom), start = as.integer(start), end = as.integer(end))]

  trait_cols <- intersect(setdiff(names(dt), c("chrom", "start", "end")), keep_traits)
  if (length(trait_cols) < 50L) stop("Too few matched traits in ", basename(file_path), ": ", length(trait_cols))
  dt <- dt[, c("chrom", "start", "end", trait_cols), with = FALSE]

  dt <- dt[
    !(
      (chrom == "17" & start < 45000000  & end > 43500000) |
      (chrom == "6"  & start < 34000000  & end > 25000000) |
      (chrom == "12" & start < 113400000 & end > 111200000)
    )
  ]

  na_counts <- rowSums(is.na(dt[, ..trait_cols]))
  dt <- dt[na_counts <= impute_number]

  for (cc in trait_cols) {
    med <- median(dt[[cc]], na.rm = TRUE)
    dt[[cc]][is.na(dt[[cc]])] <- med
  }

  win_ids <- paste(dt$chrom, as.integer(dt$start), as.integer(dt$end), sep = "_")
  imputed_traits <- as.data.frame(dt[, ..trait_cols])
  rownames(imputed_traits) <- win_ids

  traits_t <- t(imputed_traits)
  traits_t_scaled <- t(scale(t(traits_t), center = TRUE, scale = TRUE))
  nonfinite <- !is.finite(traits_t_scaled)
  if (any(nonfinite)) {
    message("Non-finite values after scaling in ", basename(file_path), ": ", sum(nonfinite), " (set to 0)")
    traits_t_scaled[nonfinite] <- 0
  }
  list(X = traits_t_scaled, trait_names = rownames(traits_t_scaled), window_ids = colnames(traits_t_scaled))
}

align_pc_signs <- function(scores_dt, loadings_dt, anchors) {
  flips <- list()
  for (pc in names(anchors)) {
    if (!(pc %in% names(scores_dt))) next
    tr <- anchors[[pc]]
    idx <- which(scores_dt$trait_abbrev == tr)
    if (length(idx) == 0L) {
      flips[[pc]] <- data.table(pc = pc, anchor_trait = tr, flip_applied = NA_integer_, anchor_value = NA_real_)
      next
    }
    aval <- scores_dt[[pc]][idx[1]]
    flip <- if (is.finite(aval) && aval < 0) -1 else 1
    scores_dt[[pc]] <- scores_dt[[pc]] * flip
    if (pc %in% names(loadings_dt)) loadings_dt[[pc]] <- loadings_dt[[pc]] * flip
    flips[[pc]] <- data.table(pc = pc, anchor_trait = tr, flip_applied = as.integer(flip == -1), anchor_value = aval)
  }
  list(scores = scores_dt, loadings = loadings_dt, flips = rbindlist(flips, use.names = TRUE, fill = TRUE))
}

build_scatter <- function(dt, xpc, ypc, xlab, ylab, labels_keep) {
  lab_col <- paste0("lbl_", xpc, "_", ypc)
  dt[, (lab_col) := ifelse(trait_abbrev %in% labels_keep, trait_abbrev, NA_character_)]
  ggplot(dt, aes(x = get(xpc), y = get(ypc), fill = phenotype_group)) +
    geom_point(shape = 21, size = 2, colour = "black", stroke = 0.4, alpha = 0.9) +
    geom_text_repel(
      data = dt[!is.na(get(lab_col))],
      aes(label = get(lab_col)),
      size = 3,
      box.padding = 0.3,
      segment.size = 0.3,
      segment.alpha = 0.5,
      max.overlaps = Inf,
      force = 1.5,
      force_pull = 0.1,
      min.segment.length = 0
    ) +
    scale_fill_manual(values = category_colors, drop = FALSE) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
      panel.grid.minor = element_line(color = "grey94", linewidth = 0.2),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
      plot.title = element_text(size = 10, face = "plain")
    ) +
    labs(x = xlab, y = ylab)
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

make_cosine_matrix <- function(pc) {
  mat <- matrix(NA_real_, nrow = length(scales), ncol = length(scales), dimnames = list(scale_labels[scales], scale_labels[scales]))
  for (i in seq_along(scales)) {
    for (j in seq_along(scales)) {
      si <- scales[i]; sj <- scales[j]
      a <- score_vecs[[si]][, .(trait_abbrev, v = get(pc))]
      b <- score_vecs[[sj]][, .(trait_abbrev, v = get(pc))]
      m <- merge(a, b, by = "trait_abbrev", suffixes = c("_a", "_b"))
      mat[i, j] <- abs(cosine(m$v_a, m$v_b))
    }
  }
  mat
}

make_procrustes_matrix <- function(k) {
  pcs_k <- paste0("PC", seq_len(k))
  mat <- matrix(NA_real_, nrow = length(scales), ncol = length(scales), dimnames = list(scale_labels[scales], scale_labels[scales]))
  for (i in seq_along(scales)) {
    for (j in seq_along(scales)) {
      si <- scales[i]; sj <- scales[j]
      a <- score_vecs[[si]][, c("trait_abbrev", pcs_k), with = FALSE]
      b <- score_vecs[[sj]][, c("trait_abbrev", pcs_k), with = FALSE]
      m <- merge(a, b, by = "trait_abbrev", suffixes = c("_a", "_b"))
      A <- as.matrix(m[, paste0(pcs_k, "_a"), with = FALSE])
      B <- as.matrix(m[, paste0(pcs_k, "_b"), with = FALSE])
      mat[i, j] <- procrustes_similarity(A, B)
    }
  }
  mat
}

make_cosine_heat <- function(pc) {
  mat <- make_cosine_matrix(pc)
  dt <- as.data.table(as.table(mat))
  setnames(dt, c("Var1", "Var2", "cosine"))
  dt[, Var1 := factor(as.character(Var1), levels = scale_labels[scales])]
  dt[, Var2 := factor(as.character(Var2), levels = scale_labels[scales])]
  dt[, is_diag := Var1 == Var2]
  dt_diag <- dt[is_diag == TRUE]
  dt_off  <- dt[is_diag == FALSE]
  ggplot() +
    geom_tile(data = dt_diag, aes(x = Var2, y = Var1), fill = "grey80", color = "white", linewidth = 0.5) +
    geom_tile(data = dt_off, aes(x = Var2, y = Var1, fill = cosine), color = "white", linewidth = 0.5) +
    geom_text(data = dt, aes(x = Var2, y = Var1, label = sprintf("%.2f", cosine)), size = 3) +
    scale_fill_gradientn(
      colours = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
      limits = c(0, 1), name = "|cosine|"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(), axis.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 10, face = "plain")
    ) +
    labs(title = paste0(pc, " |cosine|"))
}

make_procrustes_heat <- function(k) {
  mat <- make_procrustes_matrix(k)
  dt <- as.data.table(as.table(mat))
  setnames(dt, c("Var1", "Var2", "procrustes"))
  dt[, Var1 := factor(as.character(Var1), levels = scale_labels[scales])]
  dt[, Var2 := factor(as.character(Var2), levels = scale_labels[scales])]
  dt[, is_diag := Var1 == Var2]
  dt_diag <- dt[is_diag == TRUE]
  dt_off  <- dt[is_diag == FALSE]
  ggplot() +
    geom_tile(data = dt_diag, aes(x = Var2, y = Var1), fill = "grey80", color = "white", linewidth = 0.5) +
    geom_tile(data = dt_off, aes(x = Var2, y = Var1, fill = procrustes), color = "white", linewidth = 0.5) +
    geom_text(data = dt, aes(x = Var2, y = Var1, label = sprintf("%.2f", procrustes)), size = 3) +
    scale_fill_gradientn(
      colours = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
      limits = c(0, 1), name = "Procrustes"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(), axis.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 10, face = "plain")
    ) +
    labs(title = paste0("k=", k, " Procrustes"))
}

# ---- Main loop: PCA at each scale ----
pca_by_scale  <- list()
score_vecs    <- list()
panels_left   <- list()
panels_right  <- list()
audit_rows    <- list()
aa <- 1L

for (sc in scales) {
  message("Scale ", sc, ": preprocessing...")
  prep <- prep_matrix(effect_files[[sc]], target_traits)
  message("Scale ", sc, ": PCA...")
  pca <- prcomp_irlba(prep$X, n = 80, center = FALSE, scale. = FALSE)
  rownames(pca$x) <- rownames(prep$X)
  rownames(pca$rotation) <- colnames(prep$X)

  G <- tcrossprod(prep$X) / max(1, nrow(prep$X) - 1L)
  eig_all <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
  eig_all <- pmax(eig_all, 0)
  pve_full <- eig_all / sum(eig_all)
  pve_pct  <- pve_full * 100

  score_keep <- paste0("PC", seq_len(max(procrustes_k_values)))
  scores <- as.data.table(pca$x)[, ..score_keep]
  scores[, trait_abbrev := rownames(pca$x)]
  setcolorder(scores, c("trait_abbrev", score_keep))
  scores <- merge(scores, trait_map, by = "trait_abbrev", all.x = TRUE)
  scores[, phenotype_group := factor(phenotype_group, levels = names(category_colors))]

  loadings <- as.data.table(pca$rotation)[, .(window_id = rownames(pca$rotation), PC1, PC2, PC3, PC4)]

  aligned  <- align_pc_signs(scores, loadings, anchor_map)
  scores   <- aligned$scores
  loadings <- aligned$loadings
  flips    <- aligned$flips
  flips[, scale := sc]

  scree <- data.table(
    PC = seq_along(pve_full),
    var_explained = pve_full,
    var_explained_pct = pve_full * 100,
    cumulative = cumsum(pve_full)
  )
  fwrite(scores,   file.path(out_dir, paste0("pca_traits_",       sc, "_anchored.tsv")), sep = "\t")
  fwrite(loadings, file.path(out_dir, paste0("pca_windows_",      sc, "_anchored.tsv")), sep = "\t")
  fwrite(scree,    file.path(out_dir, paste0("pca_scree_",        sc, "_anchored.tsv")), sep = "\t")
  fwrite(flips,    file.path(out_dir, paste0("pca_anchor_flips_", sc, "_anchored.tsv")), sep = "\t")

  pca_by_scale[[sc]] <- list(scores = scores, loadings = loadings, pve_full = pve_full)
  score_vecs[[sc]]   <- scores[, c("trait_abbrev", score_keep), with = FALSE]

  audit_rows[[aa]] <- data.table(
    scale = sc, pve_pc1 = pve_full[1], pve_pc2 = pve_full[2],
    pve_pc3 = pve_full[3], pve_pc4 = pve_full[4],
    n_traits = nrow(scores), n_windows = nrow(loadings)
  )
  aa <- aa + 1L

  xlab12 <- paste0("PC1 (", sprintf("%.2f", pve_pct[1]), "%)")
  ylab12 <- paste0("PC2 (", sprintf("%.2f", pve_pct[2]), "%)")
  xlab23 <- paste0("PC2 (", sprintf("%.2f", pve_pct[2]), "%)")
  ylab23 <- paste0("PC3 (", sprintf("%.2f", pve_pct[3]), "%)")

  p12 <- build_scatter(copy(scores), "PC1", "PC2", xlab12, ylab12, fixed_pca_labels$pc12) +
    labs(title = scale_labels[[sc]])
  p23 <- build_scatter(copy(scores), "PC2", "PC3", xlab23, ylab23, fixed_pca_labels$pc23)
  panels_left[[sc]]  <- p12
  panels_right[[sc]] <- p23
  message("Scale ", sc, ": scatter panels done")
}

# ---- Assemble figure ----
row_panels <- lapply(scales, function(sc) {
  cowplot::plot_grid(panels_left[[sc]], panels_right[[sc]], ncol = 2, align = "hv")
})
top_grid <- cowplot::plot_grid(plotlist = row_panels, ncol = 1, align = "v")

# ---- Write cosine/procrustes cache ----
audit_dt <- rbindlist(audit_rows, use.names = TRUE, fill = TRUE)
audit_dt[, `:=`(
  pve_pc1_pct = pve_pc1 * 100, pve_pc2_pct = pve_pc2 * 100,
  pve_pc3_pct = pve_pc3 * 100, pve_pc4_pct = pve_pc4 * 100
)]
fwrite(audit_dt, out_audit, sep = "\t")
message("Variance audit written: ", out_audit)

cos_long <- rbindlist(lapply(paste0("PC", 1:4), function(pc) {
  mat <- make_cosine_matrix(pc)
  dt  <- as.data.table(as.table(mat))
  setnames(dt, c("scale_a", "scale_b", "value"))
  dt[, `:=`(metric = "cosine_abs", component = pc)]
  dt
}), use.names = TRUE, fill = TRUE)
proc_long <- rbindlist(lapply(procrustes_k_values, function(k) {
  mat <- make_procrustes_matrix(k)
  dt  <- as.data.table(as.table(mat))
  setnames(dt, c("scale_a", "scale_b", "value"))
  dt[, `:=`(metric = "procrustes_similarity", component = paste0("k", k))]
  dt
}), use.names = TRUE, fill = TRUE)
fwrite(cos_long,  out_cos_long,  sep = "\t")
fwrite(proc_long, out_proc_long, sep = "\t")

ggsave(out_fig, top_grid, width = 14, height = 30)
message("FigureS11 complete: ", out_fig)
