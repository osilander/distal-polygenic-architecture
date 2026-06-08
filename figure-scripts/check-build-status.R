#!/usr/bin/env Rscript
# check-build-status.R
# Verifies the repo environment and reports what can and cannot be built.
# Run from the project root: Rscript figure-scripts/check-build-status.R

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

ok  <- function(msg) cat(sprintf("  [OK]      %s\n", msg))
warn <- function(msg) cat(sprintf("  [MISSING] %s\n", msg))
hdr <- function(msg) cat(sprintf("\n=== %s ===\n", msg))

errors <- 0L

# ---- 1. Data symlink ----
hdr("Data symlink")
data_dir <- file.path(project_root, "data")
if (!file.exists(data_dir)) {
  warn("data/ does not exist — symlink to source manuscript data tree is required")
  errors <- errors + 1L
} else if (!file.info(data_dir)$isdir) {
  warn("data/ exists but is not a directory or symlink")
  errors <- errors + 1L
} else {
  ok("data/ exists")
}

windows_dir <- file.path(project_root, "data", "windows_filtered")
if (!dir.exists(windows_dir)) {
  warn("data/windows_filtered/ not found — Step 0 (build-all-trait-matrices.R) will fail")
  errors <- errors + 1L
} else {
  n_win <- length(list.files(windows_dir, pattern = "\\.summary\\.tsv$"))
  if (n_win == 0L) {
    warn("data/windows_filtered/ is empty — no window summary files found")
    errors <- errors + 1L
  } else {
    ok(sprintf("data/windows_filtered/ — %d window summary files found", n_win))
  }
}

pickrell <- file.path(project_root, "data", "pickrell_blocks.bed")
if (!file.exists(pickrell)) {
  warn("data/pickrell_blocks.bed not found — generate-svd-nulls-50k.R will fail if re-run")
} else {
  ok("data/pickrell_blocks.bed")
}

genes_rds <- file.path(project_root, "data", "gencode.v19.genes.protein_coding.rds")
if (!file.exists(genes_rds)) {
  warn("data/gencode.v19.genes.protein_coding.rds not found — run-figureS-nearestgene-cdf-top1pct.R, FigureS25, and FigureS27 will fail")
  errors <- errors + 1L
} else {
  ok("data/gencode.v19.genes.protein_coding.rds")
}

# ---- 2. HPC null outputs ----
hdr("HPC null outputs (must be imported, not regenerated locally)")

hpc_items <- list(
  list(
    path = file.path(project_root, "results", "svd-nulls-50k", "observed", "observed_50k_singular_values.tsv"),
    label = "results/svd-nulls-50k/observed/  (needed by compute-figure2-cache.R and allbyall scripts)"
  ),
  list(
    path = file.path(project_root, "results", "svd-nulls-50k", "withinblock_perm", "withinblock_perm_iter001_50k_trait_vectors.tsv.gz"),
    label = "results/svd-nulls-50k/withinblock_perm/iter001  (needed by run-allbyall-cosine-matrices.R)"
  ),
  list(
    path = file.path(project_root, "results", "chunk-withinperm-nulls-50k", "chunk_withinperm_scree_runs_50k.tsv.gz"),
    label = "results/chunk-withinperm-nulls-50k/chunk_withinperm_scree_runs_50k.tsv.gz  (needed by compute-figure2-cache.R)"
  ),
  list(
    path = file.path(project_root, "results", "allbyall_cosine_matrices", "entrysign_perm100_procrustes_normalized.tsv"),
    label = "results/allbyall_cosine_matrices/  (needed by run-figureS12.R)"
  )
)

hpc_missing <- 0L
for (item in hpc_items) {
  if (!file.exists(item$path)) {
    warn(item$label)
    hpc_missing <- hpc_missing + 1L
    errors <- errors + 1L
  } else {
    ok(item$label)
  }
}

# ---- 3. Core intermediate results ----
hdr("Core intermediate results")

core_items <- list(
  list(path = file.path(project_root, "results", "all-trait-50k-mean-effects.tsv"),   label = "results/all-trait-50k-mean-effects.tsv"),
  list(path = file.path(project_root, "results", "all-trait-50k-mean-pvals.tsv"),     label = "results/all-trait-50k-mean-pvals.tsv"),
  list(path = file.path(project_root, "results", "all-trait-100k-mean-effects.tsv"),  label = "results/all-trait-100k-mean-effects.tsv"),
  list(path = file.path(project_root, "results", "all-trait-200k-mean-effects.tsv"),  label = "results/all-trait-200k-mean-effects.tsv"),
  list(path = file.path(project_root, "results", "all-trait-500k-mean-effects.tsv"),  label = "results/all-trait-500k-mean-effects.tsv"),
  list(path = file.path(project_root, "results", "all-trait-1m-mean-effects.tsv"),    label = "results/all-trait-1m-mean-effects.tsv"),
  list(path = file.path(project_root, "results", "figure2_cache_50k.rds"),            label = "results/figure2_cache_50k.rds"),
  list(path = file.path(project_root, "results", "pca_loadings_50k.tsv"),             label = "results/pca_loadings_50k.tsv"),
  list(path = file.path(project_root, "results", "trait_pc_category_absrho_labelperm_pc1_pc14.tsv"), label = "results/trait_pc_category_absrho_labelperm_pc1_pc14.tsv"),
  list(path = file.path(project_root, "results", "pca_multiscale_anchored", "pca_traits_25k_anchored.tsv"),  label = "results/pca_multiscale_anchored/pca_traits_25k_anchored.tsv"),
  list(path = file.path(project_root, "results", "pca_multiscale_anchored", "pca_traits_50k_anchored.tsv"),  label = "results/pca_multiscale_anchored/pca_traits_50k_anchored.tsv"),
  list(path = file.path(project_root, "results", "pca_multiscale_anchored", "pca_traits_100k_anchored.tsv"), label = "results/pca_multiscale_anchored/pca_traits_100k_anchored.tsv"),
  list(path = file.path(project_root, "results", "pca_multiscale_anchored", "pca_traits_200k_anchored.tsv"), label = "results/pca_multiscale_anchored/pca_traits_200k_anchored.tsv"),
  list(path = file.path(project_root, "results", "gwas_removed_distance_topload_50k_plus100kb", "reduced_pca_trait_scores.tsv"), label = "results/gwas_removed_distance_topload_50k_plus100kb/"),
  list(path = file.path(project_root, "results", "gwas_removed_distance_topload_50k_plus150kb", "reduced_pca_trait_scores.tsv"), label = "results/gwas_removed_distance_topload_50k_plus150kb/")
)

core_missing <- 0L
for (item in core_items) {
  if (!file.exists(item$path)) {
    warn(item$label)
    core_missing <- core_missing + 1L
  } else {
    ok(item$label)
  }
}

# ---- 4. Figure outputs ----
hdr("Figure outputs")

figures <- list.files(file.path(project_root, "figures"), pattern = "\\.pdf$")
main_figs <- grep("^figure[0-9]+\\.pdf$", figures, value = TRUE)
supp_figs <- grep("^figureS", figures, value = TRUE)
cat(sprintf("  Main figures: %d present\n", length(main_figs)))
cat(sprintf("  Supplementary figures: %d present\n", length(supp_figs)))
if (length(main_figs) < 4L) warn("Some main figures are missing")

expected_supp <- c(
  "figureS1_window_snp_count_histograms_50k.pdf",
  "figureS2_robustness_checks.pdf",
  "figureS3_trait_pc_category_absrho_boxplot_pc1_pc8.pdf",
  "figureS4_trait_pc_category_absrho_labelperm_pc1_pc14.pdf",
  "figureS5_umap.pdf",
  "figureS6_trait_loading_heatmap_top60_pc1_pc5.pdf",
  "figureS7_full_enrichment_dotplot_nearest_gene_top2pct.pdf",
  "figureS8_pca_abs50k_from_windows.pdf",
  "figureS9_pca_absolute_50k.pdf",
  "figureS10_within_over_between_abs_cosine_k4_k8_normalized_expanded.pdf",
  "figureS11_multiscale_pca_cosine_50k_1m.pdf",
  "figureS12_allbyall_procrustes_heatmaps.pdf",
  "figureS13_figure3_pc4_manhattans.pdf",
  "figureS14_withinperm_peaks_pc1_4.pdf",
  "figureS15_figure2_extra_zooms.pdf",
  "figureS16_observed_nearestgwaspeak_pattern_summary_50k_100k_comparison.pdf",
  "figureS17_chr9_chr17_pval_scatter.pdf",
  "figureS18_gwas_removed_trait_pc_category_absrho_labelperm_pc1_pc14.pdf",
  "figureS19_gwas_removed_pc_correlations_pc1_pc4.pdf",
  "figureS20_gwasremoved_enrichment_dotplot_nearest_gene_plus100kb_top2pct.pdf",
  "figureS21_gwasremoved_withinperm_peaks_pc1_4.pdf",
  "figureS22_figure3_extra_zooms_plus100kb.pdf",
  "figureS23_gwas_removed_subspace_scatter.pdf",
  "figureS24_observed_nearestgene_pattern_summary_50k_100k_comparison.pdf",
  "figureS25_nearestgene_nonzero_distance_cdf_ribbon_50k_bands.pdf",
  "figureS26_gwasremoved_nearestgene_pattern_summary_gwasremoved100_150kb.pdf",
  "figureS27_gwasremoved_nearestgene_nonzero_distance_cdf_ribbon_50k_bands_plus100kb.pdf"
)
missing_supp <- setdiff(expected_supp, figures)
if (length(missing_supp) > 0L) {
  for (f in missing_supp) warn(paste("Missing supplementary figure:", f))
} else {
  ok(sprintf("All %d expected supplementary figures present", length(expected_supp)))
}

# ---- 5. Summary ----
hdr("Build capability summary")

can_run_step0 <- dir.exists(windows_dir) && length(list.files(windows_dir, pattern = "\\.summary\\.tsv$")) > 0L
has_50k_matrix <- file.exists(file.path(project_root, "results", "all-trait-50k-mean-effects.tsv"))
has_pca_loadings <- file.exists(file.path(project_root, "results", "pca_loadings_50k.tsv"))
has_hpc_nulls <- hpc_missing == 0L
has_gwas_removed <- file.exists(file.path(project_root, "results", "gwas_removed_distance_topload_50k_plus100kb", "reduced_pca_trait_scores.tsv"))

cat("\n")
cat(sprintf("  %-50s %s\n", "Step 0 (build-all-trait-matrices.R):",    if (can_run_step0) "READY" else "BLOCKED — data/windows_filtered/ missing"))
cat(sprintf("  %-50s %s\n", "Steps 1-3 (permutations, 25k, GWAS rm):", if (has_50k_matrix) "READY" else "BLOCKED — run Step 0 first"))
cat(sprintf("  %-50s %s\n", "Steps 4-5 (null cache, figure2_cache):",  if (has_hpc_nulls) "READY" else "BLOCKED — import HPC null outputs first"))
cat(sprintf("  %-50s %s\n", "All main figures (figure1-4):",            if (has_pca_loadings && has_gwas_removed) "READY" else "BLOCKED — run figure1 and GWAS-removed steps first"))
cat(sprintf("  %-50s %s\n", "FigureS12 (allbyall Procrustes):",         if (has_hpc_nulls) "READY" else "BLOCKED — import HPC allbyall_cosine_matrices/ first"))

if (errors == 0L) {
  cat("\nAll checks passed. Full build is possible.\n")
} else {
  cat(sprintf("\n%d issue(s) found. See above for details.\n", errors))
}
