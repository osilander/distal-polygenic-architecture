# distal-polygenic-architecture
Post-processing pipeline for all figures in the distal polygenic architecture manuscript.
Scripts are run from the project root (`distal-polygenic-architecture/`).
> **Repo intent:** code release. `figures/`, `data/`, and `results/` are not tracked in git and will not exist on a fresh checkout — figure scripts create them on first run. Pre-computed HPC outputs and reference files are distributed via Zenodo (see below); `data/windows_filtered/` must be generated via HPC. A working copy may also contain local-only files (`*.tar.gz`, `*.pdf`, `*.txt`) that are gitignored and not part of the tracked repo surface.

---
## Quick build-status check
```r
Rscript figure-scripts/check-build-status.R
```
Verifies data symlink, HPC inputs, intermediate results, and reports what can and cannot be built.

---
## Build tiers
| Tier | What it covers | Requirement |
|------|----------------|-------------|
| **Fully local** | All figure scripts (run-figure\*.R, run-figureS\*.R) and downstream summaries | Zenodo data downloaded and symlinked into `results/` (see Setup below) |
| **Local if windows available** | Step 0 matrix build; FigureS1, FigureS8, FigureS11 | `data/windows_filtered/` (generate via `slurms/build-windows.slurm`) |
| **Blocked — HPC required** | Figure 1 scree panel | `chunk_withinperm_scree_runs_50k.tsv.gz` (generate via `slurms/chunk-withinperm-*.slurm`) |

---
## Directory structure
```
distal-polygenic-architecture/
  data/                              → not in git; populate from Zenodo + HPC (see Setup)
    svd-nulls-50k/                   → Zenodo: svd-nulls-50k.tar.gz
      entry_flip/                      100 × singular values; iter001 trait vectors
      withinblock_perm/                100 × singular values; iter001 trait vectors
      observed/                        Observed singular values and trait/window vectors
    allbyall_cosine_matrices/        → Zenodo: allbyall_cosine_matrices.tar.gz
    chunk_withinperm_base_50k.rds    → Zenodo: single file (134 MB); base object for scree null analysis
    gencode.v19.annotation.gtf.gz    → Zenodo: reference-data.tar.gz
    gencode.v19.genes.protein_coding.rds → Zenodo: reference-data.tar.gz
    pickrell_blocks.bed              → Zenodo: reference-data.tar.gz
    trait_abbrevs_categorised.txt   → Zenodo: reference-data.tar.gz
    windows_filtered/                → HPC only (not on Zenodo): slurms/build-windows.slurm
  results/                           → generated outputs (not in git)
  figures/                           → figure PDFs (generated; not in git)
  figure-scripts/                    → R scripts
  slurms/                            → HPC SLURM scripts
    build-windows.slurm                Array job: per-trait window summaries → data/windows_filtered/
    chunk-withinperm-prepare.slurm     Step 1: prepare chunk within-perm base RDS on HPC
    chunk-withinperm-run.slurm         Step 2: array job workers → results/chunk-withinperm-nulls-50k/
  README.md
  data/README.md                     → Zenodo deposit description
```

---
## Setup — populating data/ and results/

### Step A — Download from Zenodo
Download the archives from the Zenodo deposit and extract into the `data/` directory:
```bash
cd data/
tar -xzf svd-nulls-50k.tar.gz              # → data/svd-nulls-50k/{entry_flip,withinblock_perm,observed}/
tar -xzf allbyall_cosine_matrices.tar.gz   # → data/allbyall_cosine_matrices/
# chunk_withinperm_base_50k.rds is a single file; place at data/chunk_withinperm_base_50k.rds
# reference files (gencode.*, pickrell_blocks.bed, trait_abbrevs_categorised.txt) go directly in data/
# window summaries (optional; only needed for Step 0):
# tar -xzf windows_w50000.tar.gz           # → data/windows_filtered/ (repeat for each size needed)
cd ..
```

### Step B — Symlink data/ into results/
Scripts expect HPC null outputs under `results/svd-nulls-50k/`, `results/allbyall_cosine_matrices/`,
and `results/chunk-withinperm-nulls-50k/hpc_base/`:
```bash
mkdir -p results/svd-nulls-50k results/chunk-withinperm-nulls-50k/hpc_base

ln -s ../../data/svd-nulls-50k/entry_flip            results/svd-nulls-50k/entry_flip
ln -s ../../data/svd-nulls-50k/withinblock_perm      results/svd-nulls-50k/withinblock_perm
ln -s ../../data/svd-nulls-50k/observed              results/svd-nulls-50k/observed
ln -s ../data/allbyall_cosine_matrices               results/allbyall_cosine_matrices
ln -s ../../../data/chunk_withinperm_base_50k.rds    results/chunk-withinperm-nulls-50k/hpc_base/chunk_withinperm_base_50k.rds
```

### Step C — HPC-only inputs (not on Zenodo)
Two items are not in the Zenodo deposit and must be generated via HPC:

| What | How | Needed for |
|------|-----|-----------|
| `data/windows_filtered/` | `slurms/build-windows.slurm` (array job over all traits); or extract Zenodo `windows_w*.tar.gz` per window size needed | Step 0 matrix build, FigureS1/S8/S11 |
| `results/chunk-withinperm-nulls-50k/chunk_withinperm_scree_runs_50k.tsv.gz` | Re-run locally: `Rscript figure-scripts/analyse-chunk-withinperm-nulls-50k.R` (uses `data/chunk_withinperm_base_50k.rds` from Zenodo); or run on HPC via `slurms/chunk-withinperm-*.slurm` | Figure 1 scree panel |

---
## Full build order
Scripts must be run from the **project root** (`distal-polygenic-architecture/`). Listed in dependency order.
### Step 0 — Build all trait × window matrices
```r
Rscript figure-scripts/build-all-trait-matrices.R
```
Reads `data/windows_filtered/*_w*.summary.tsv` (via `data/` symlink) and writes
`results/all-trait-{25k,50k,100k,200k,500k,1m,...}-mean-{effects,pvals,rand-mean-effects}.tsv`
for every window size present.

---
### Step 1 — Pre-compute trait-PC category permutation table (needed by Figure 1)
```r
Rscript figure-scripts/compute-trait-pc-category-labelperm.R
```
**Outputs:** `results/trait_pc_category_absrho_labelperm_pc1_pc14.tsv`,
`figures/figureS4_trait_pc_category_absrho_labelperm_pc1_pc14.pdf`
**Requires:** Step 0 complete, `results/pca_loadings_50k.tsv`.
**Appearance:** tiles filled by median |ρ| (ag_Sunset palette, cap 0.2), median |ρ| value shown per tile, significance stars above; category colour strip on left. Matches Figure 1 heatmap style.

---
### Step 2 — 25k PCA (needed by FigureS10 and S12)
```r
Rscript figure-scripts/run-pca-25k-min20-anchored.R
```
**Outputs:** `results/pca_multiscale_anchored/pca_traits_25k_anchored.tsv` + window loadings, scree, flips.
Also writes `results/all-trait-25k-mean-effects-snp20.tsv` (SNP≥20 filter).
**Requires:** Step 0 complete.

---
### Step 3 — GWAS-removed distance analyses (needed by Figures 3, S13, S18–S23)
```r
Rscript figure-scripts/compute-gwas-removed-distances.R 100
Rscript figure-scripts/compute-gwas-removed-distances.R 150
Rscript figure-scripts/compute-gwas-removed-cdf-ribbon.R
```

---
### Step 4 — SVD null permutations (needed by Figure 1 scree, FigureS12)
The SVD null outputs (`results/svd-nulls-50k/`) are provided via Zenodo (`hpc-nulls.tar.gz`) and
symlinked into `results/` in Setup Step B. If re-running from scratch:
```bash
Rscript figure-scripts/generate-svd-nulls-50k.R 100
```
The chunk within-permutation scree runs are generated separately. Two options:
```bash
# Option A — HPC array job (faster):
sbatch slurms/chunk-withinperm-prepare.slurm          # produces hpc_base/ RDS
sbatch --array=1-200 slurms/chunk-withinperm-run.slurm  # 100 iters × 2 chunk sizes
Rscript figure-scripts/summarise-chunk-withinperm-hpc-50k.R  # collect into .tsv.gz

# Option B — local re-run from Zenodo base RDS (slower, no HPC needed):
Rscript figure-scripts/analyse-chunk-withinperm-nulls-50k.R
```
**Outputs:** `results/svd-nulls-50k/`, `results/chunk-withinperm-nulls-50k/chunk_withinperm_scree_runs_50k.tsv.gz`

---
### Step 5 — Figure 1 scree cache (needed by Figure 1)
```r
Rscript figure-scripts/compute-figure2-cache.R
```
**Outputs:** `results/figure2_cache_50k.rds`
**Requires:** Step 0 complete, Step 4 outputs present.

---
### Step 6 — All-by-all cosine / Procrustes matrices (needed by FigureS12)
These are HPC-generated. Copy from the source repo rather than re-running locally.
If re-running is needed:
```bash
Rscript figure-scripts/run-allbyall-cosine-matrices.R
Rscript figure-scripts/run-entrysign-perm-procrustes-normalized.R
```
**Outputs:** `results/allbyall_cosine_matrices/`
**Requires:** Steps 0, 2, 3, 4 complete; FigureS8 and FigureS9 complete.

---
## Figure generation

---
### FigureS1 — Window SNP count distributions
```r
Rscript figure-scripts/run-figureS1.R
```
**Outputs:** `figures/figureS1_window_snp_count_histograms_50k.pdf`,
`results/window_snp_count_histograms_50k.tsv`
**Requires:** `data/windows_filtered/*_w50000.summary.tsv` (via symlink).

---
### FigureS2 — SVD robustness checks
```r
Rscript figure-scripts/run-figureS2.R
```
**Outputs:** `figures/figureS2_robustness_checks.pdf`
**Requires:** `results/svd_robustness/*.tsv` (from `compute-svd-robustness.R`).
**Appearance:** Five-panel figure. (a) Per-PC cosine similarity to original decomposition across 100 random and stratified 2/3 trait subsampling replicates. (b) Procrustes similarity for trait subsampling (k = 4, 8). (c) Per-PC cosine similarity for 100 window subsampling replicates on the full and GWAS-depleted matrices. (d) Procrustes similarity for window subsampling. (e) Leave-one-category-out cosine similarity for PC1–PC8, sorted by PC1 disruption; observed values colour-coded using the ag_Sunset palette (warm = low cosine = structurally important category); grey bars show size-matched null (30 replicates, median + 95% interval).
To recompute the robustness TSVs from scratch:
```r
Rscript figure-scripts/compute-svd-robustness.R
```

---
### FigureS3 — Trait-category × PC Spearman ρ boxplot
```r
Rscript figure-scripts/run-figureS3.R
```
**Outputs:** `figures/figureS3_trait_pc_category_absrho_boxplot_pc1_pc8.pdf`,
`results/trait_pc_category_absrho_summary_pc1_pc8.tsv`
**Requires:** Step 0 complete, `results/pca_loadings_50k.tsv`.

---
### FigureS4 — Trait-category × PC label-permutation heatmap (standalone)
```r
Rscript figure-scripts/compute-trait-pc-category-labelperm.R
```
**Outputs:** `figures/figureS4_trait_pc_category_absrho_labelperm_pc1_pc14.pdf`
(also produced as part of Step 1; see above)

---
### Figure 1 + FigureS5 + FigureS6 — PCA scatters, UMAP, trait heatmap
```r
Rscript figure-scripts/run-figure1.R
```
**Outputs:**
- `figures/figure1.pdf`
- `figures/figureS5_umap.pdf`
- `figures/figureS6_trait_loading_heatmap_top60_pc1_pc5.pdf`
- `results/pca_loadings_50k.tsv`, `results/umap-coordinates.tsv`
**Requires:** Steps 0–1 complete; `results/figure2_cache_50k.rds` **†**.

---
### FigureS7 — Full enrichment dotplot (nearest gene, unfiltered PCA)
```r
Rscript figure-scripts/run-figureS7.R
```
**Outputs:** `figures/figureS7_full_enrichment_dotplot_nearest_gene_top2pct.pdf`
**Requires:** `results/full_enrichment/enrichment_results.tsv` (from `compute-full-gene-enrichment.R`); `GOSemSim`, `org.Hs.eg.db`.
**Note:** GO terms are redundancy-filtered using Wang semantic similarity (cutoff 0.7) via GOSemSim before plotting. Top 2% band, GO:BP/MF and GO:CC panels, intersection size ≥ 15.

---
### FigureS8 — Absolute-value PCA from raw window files
```r
Rscript figure-scripts/run-figureS8.R
```
**Outputs:** `figures/figureS8_pca_abs50k_from_windows.pdf`,
`results/pca_abs50k_from_windows/{trait_scores,window_loadings,singular_values,scree_fullspectrum,window_coordinates,heatmap_matrix_top60}.tsv`
**Requires:** `data/windows_filtered/*_w50000.summary.tsv` (via symlink).

---
### FigureS9 — Absolute-value PCA from mean effects
```r
Rscript figure-scripts/run-figureS9.R
```
**Outputs:** `figures/figureS9_pca_absolute_50k.pdf`,
`results/figureS_pca_absolute_50k_{trait_scores,scree,window_loadings_top200,heatmap_matrix_top60}.tsv`
**Requires:** Step 0 complete.

---
### FigureS10 — Within/over/between cosine coherence
```r
Rscript figure-scripts/run-figureS10.R
```
**Outputs:** `figures/figureS10_within_over_between_abs_cosine_k4_k8_normalized_expanded.pdf`
**Requires:** FigureS9 complete, FigureS11 complete, Step 2 complete.

---
### FigureS11 — Multiscale PCA scatter + cosine/Procrustes heatmaps
```r
Rscript figure-scripts/run-figureS11.R
```
**Outputs:** `figures/figureS11_multiscale_pca_cosine_50k_1m.pdf`,
`results/pca_multiscale_anchored/{pca_traits,pca_windows,pca_scree,pca_anchor_flips}_{50k,100k,200k,500k,1m}_anchored.tsv`
**Requires:** Step 0 complete.

---
### FigureS12 — All-by-all Procrustes heatmaps
```r
Rscript figure-scripts/run-figureS12.R
```
**Outputs:** `figures/figureS12_allbyall_procrustes_heatmaps.pdf`
**Requires:** `results/allbyall_cosine_matrices/` **†**.

---
### Figure 2 + FigureS15 — Manhattan plots and locus zooms
```r
Rscript figure-scripts/run-figure2.R
```
**Outputs:** `figures/figure2.pdf`, `figures/figureS15_figure2_extra_zooms.pdf`,
`results/Supplementary_PC_Top_Loadings_figure2.xlsx`
**Requires:** `results/pca_loadings_50k.tsv`.

---
### Figure 3 + FigureS13 — GWAS-removed PCA
```r
Rscript figure-scripts/run-figure3.R
```
**Outputs:** `figures/figure3.pdf`, `figures/figureS13_figure3_pc4_manhattans.pdf`
**Requires:** Step 3 complete, `results/pca_loadings_50k.tsv`.

---
### FigureS14 — Within-permutation loading concentration (withinperm peaks)
Generate the input tables first, then plot:
```r
Rscript figure-scripts/compute-withinperm-peaks-50k.R
Rscript figure-scripts/run-figureS14.R
```
**Outputs:** `figures/figureS14_withinperm_peaks_pc1_4.pdf`
**Intermediates:** `results/withinperm_peaks_iter_metrics.tsv`, `results/withinperm_peaks_summary.tsv`
**Requires:** `results/pca_loadings_50k.tsv`; chunk within-perm window vectors in `results/chunk-withinperm-nulls-50k/hpc_runs/window_vectors/` (produced by `slurms/chunk-withinperm-run.slurm`).

---
### FigureS16 + FigureS24 — Nearest-peak and nearest-gene pattern comparison (50k vs 100k)
```r
Rscript figure-scripts/compute-nearestgene-pattern-bands.R 50k
Rscript figure-scripts/compute-nearestgene-pattern-bands.R 100k
Rscript figure-scripts/run-figureS16-and-S24.R
```
**Outputs:** `figures/figureS16_observed_nearestgwaspeak_pattern_summary_50k_100k_comparison.pdf`,
`figures/figureS24_observed_nearestgene_pattern_summary_50k_100k_comparison.pdf`
**Note:** uses patchwork so each window-size row has its own axes.
**Requires:** `results/pca_loadings_50k.tsv`, `results/pca_multiscale_anchored/pca_windows_100k_anchored.tsv`
(generated by `run-figureS11.R`), `results/all-trait-100k-mean-pvals.tsv` (from Step 0).

---
### FigureS17 — Chr9 and Chr17 per-trait p-value scatter
```r
Rscript figure-scripts/run-figureS17.R
```
**Outputs:** `figures/figureS17_chr9_chr17_pval_scatter.pdf`
**Requires:** `results/all-trait-50k-mean-pvals.tsv`.
**Shows:** Per-trait −log10(p) across a 1 Mb window centred on the chr9 (~17.15 Mb) and chr17 (~21.4 Mb) top PC1 loci. Background dots = all traits; purple dots/labels = trait with maximum −log10(p) per window.

---
### FigureS18 — GWAS-removed trait-category × PC heatmap
```r
Rscript figure-scripts/run-figureS18.R
```
**Outputs:** `figures/figureS18_gwas_removed_trait_pc_category_absrho_labelperm_pc1_pc14.pdf`
**Requires:** Step 3 complete.
**Appearance:** tiles filled by median |ρ| (ag_Sunset palette), median |ρ| value shown per tile, significance stars above; category colour strip on left. Two panels (+/-100 kb, +/-150 kb).

---
### FigureS19 — GWAS-removed PC correlation heatmap
```r
Rscript figure-scripts/run-figureS19.R
```
**Outputs:** `figures/figureS19_gwas_removed_pc_correlations_pc1_pc4.pdf`
**Requires:** Step 3 complete, `results/pca_loadings_50k.tsv`.

---
### FigureS20 — GWAS-removed GO enrichment dotplot
```r
Rscript figure-scripts/run-figureS20.R [flank_kb]
```
**Outputs:** `figures/figureS20_gwasremoved_enrichment_dotplot_nearest_gene_plus{flank_kb}kb_top2pct.pdf`
**Requires:** `results/gwas_removed_enrichment/enrichment_results_plus{flank_kb}kb.tsv` (from `compute-gwas-removed-gene-enrichment.R`); `GOSemSim`, `org.Hs.eg.db`.
**Note:** GO terms are redundancy-filtered using Wang semantic similarity (cutoff 0.7) via GOSemSim before plotting. Top 2% band, GO:BP/MF and GO:CC panels, intersection size ≥ 15.

---
### FigureS21 — GWAS-removed within-permutation loading concentration
```r
Rscript figure-scripts/compute-gwasremoved-withinperm-peaks.R [n_perm]
Rscript figure-scripts/run-figureS21.R
```
**Outputs:** `figures/figureS21_gwasremoved_withinperm_peaks_pc1_4.pdf`
**Inputs:** `results/gwasremoved_withinperm_peaks_iter_metrics.tsv`, `results/gwasremoved_withinperm_peaks_summary.tsv`
**Requires:** Step 3 complete; `results/all-trait-50k-mean-{effects,pvals}.tsv`.

---
### FigureS22 — Figure 3 extra locus zooms (±100 kb)
```r
Rscript figure-scripts/run-figureS22.R
```
**Outputs:** `figures/figureS22_figure3_extra_zooms_plus100kb.pdf`
**Requires:** Step 3 complete, `results/pca_loadings_50k.tsv`.

---
### FigureS23 — GWAS-removed subspace scatter
```r
Rscript figure-scripts/run-figureS23.R
```
**Outputs:** `figures/figureS23_gwas_removed_subspace_scatter.pdf`,
`results/figureS23_gwas_removed_subspace_scatter.tsv`
**Requires:** Step 3 complete, `results/pca_loadings_50k.tsv`.

---
### Figure 4 — Nearest-gene enrichment panels
Run upstream scripts first, then assemble the figure:
```r
Rscript figure-scripts/compute-nearestgene-patterns.R gene 10000
Rscript figure-scripts/compute-genicclass-enrichment.R
Rscript figure-scripts/compute-nearestgene-cdf-top1pct.R
Rscript figure-scripts/run-figure4.R
```
**Outputs:** `figures/figure4.pdf`
**Also writes** (upstream scripts):
`results/observed-nearest-patterns/observed_topfracs_nearestgene_bin_enrichment.tsv`,
`results/observed-nearest-patterns/observed_topfracs_nearestgene_null_replicate_stats.tsv.gz`,
`results/genicclass_rank_enrichment_50k/cumulative_genicclass_enrichment_curves_pc1_pc4.tsv`,
`results/nearestgene_nonzero_distance_cdf_ribbon_50k_top1pct.tsv`,
`results/nearestgene_nonzero_distance_cdf_ribbon_50k_top1pct_pvalues.tsv`
**Requires:** `results/pca_loadings_50k.tsv`, `data/gencode.v19.genes.protein_coding.rds`,
`data/gencode.v19.annotation.gtf.gz`.

---
### FigureS25 — Nearest-gene distance CDF ribbons (bands)
```r
Rscript figure-scripts/run-figureS25.R
```
**Outputs:** `figures/figureS25_nearestgene_nonzero_distance_cdf_ribbon_50k_bands.pdf`,
`results/nearestgene_nonzero_distance_cdf_ribbon_50k_bands.tsv`,
`results/nearestgene_nonzero_distance_cdf_ribbon_50k_bands_pvalues.tsv`
**Requires:** `results/pca_loadings_50k.tsv`, `data/gencode.v19.genes.protein_coding.rds`.

---
### FigureS26 — GWAS-removed nearest-gene pattern (100 kb and 150 kb flanks stacked)
```r
Rscript figure-scripts/run-figureS26.R
```
**Outputs:** `figures/figureS26_gwasremoved_nearestgene_pattern_summary_gwasremoved100_150kb.pdf`
**Note:** stacks 100 kb and 150 kb removal flanks as separate rows.
**Requires:** `results/observed-nearest-patterns-gwasremoved{100,150}kb-top1_2_5_band45_50/`
(from `compute-nearestgene-pattern-bands.R gwasremoved default gene_only {100|150}`).

---
### FigureS27 — GWAS-removed nearest-gene distance CDF ribbons (bands)
```r
Rscript figure-scripts/run-figureS27.R [flank_kb]
```
**Outputs:** `figures/figureS27_gwasremoved_nearestgene_nonzero_distance_cdf_ribbon_50k_bands_plus{flank_kb}kb.pdf`,
`results/gwasremoved_nearestgene_nonzero_distance_cdf_ribbon_50k_bands_plus{flank_kb}kb.tsv`,
`results/gwasremoved_nearestgene_nonzero_distance_cdf_ribbon_50k_bands_plus{flank_kb}kb_pvalues.tsv`
Default flank = 100 kb.
**Requires:** Step 3 complete, `data/gencode.v19.genes.protein_coding.rds`.

---
## R package dependencies
```r
install.packages(c("data.table", "dplyr", "ggplot2", "ggrepel", "irlba",
                   "patchwork", "RColorBrewer", "scales", "cowplot",
                   "uwot", "dbscan", "paletteer"))
BiocManager::install(c("ComplexHeatmap", "circlize", "GenomicRanges",
                       "GenomeInfoDb", "IRanges"))
```
