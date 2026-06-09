#!/usr/bin/env Rscript
# Produces: results/full_enrichment/gene_sets.tsv
#           results/full_enrichment/enrichment_results.tsv
# Pre-computed inputs:
#   results/pca_loadings_50k.tsv
#   data/gencode.v19.genes.protein_coding.rds
# Run: Rscript figure-scripts/compute-full-gene-enrichment.R
# Requires: clusterProfiler, org.Hs.eg.db, msigdbr  (all local, no internet needed)
# Databases: GO (BP, MF, CC), MSigDB Hallmarks, MSigDB C2 (curated pathways)
# Parallel to compute-gwas-removed-gene-enrichment.R but uses the full (unfiltered) PCA loadings.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

pcs         <- paste0("PC", 1:14)
bands       <- c(0.02)   # only Top 2% is used by run-figureS7.R
band_labels <- c("0.02" = "Top 2%")
dist_modes  <- c("nearest_gene")
p_cutoff    <- 0.05

out_dir <- file.path("results", "full_enrichment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_genesets <- file.path(out_dir, "gene_sets.tsv")
out_enrich   <- file.path(out_dir, "enrichment_results.tsv")

loadings_file <- file.path("results", "pca_loadings_50k.tsv")
if (!file.exists(loadings_file))
  stop("Missing: ", loadings_file, "\nRun: Rscript figure-scripts/run-figure1.R first.")
genes_rds <- file.path(project_root, "data", "gencode.v19.genes.protein_coding.rds")
if (!file.exists(genes_rds)) stop("Missing genes RDS: ", genes_rds)

# ---- Load data ----
message("Loading loadings and gene annotation...")
ld <- fread(loadings_file)
# pca_loadings_50k.tsv has PC columns first, then region_label, chrom, start, end
if (!all(c("chrom", "start", "end") %in% names(ld))) {
  if ("region_label" %in% names(ld)) {
    ld[, c("chrom", "start", "end") := tstrsplit(region_label, "_", fixed = TRUE)]
  } else {
    stop("Cannot find coordinate columns in ", loadings_file)
  }
}
ld[, `:=`(chrom = as.character(chrom), start = as.integer(start), end = as.integer(end))]
ld[, chrom := gsub("^chr", "", chrom)]
if (!"window_id" %in% names(ld)) ld[, window_id := paste(chrom, start, end, sep = "_")]
ld <- ld[is.finite(start) & is.finite(end)]
message("Retained windows: ", nrow(ld))

genes_gr <- readRDS(genes_rds)
if (any(grepl("^chr", seqlevels(genes_gr)))) seqlevels(genes_gr) <- sub("^chr", "", seqlevels(genes_gr))
suppressWarnings(GenomeInfoDb::seqlevelsStyle(genes_gr) <- "NCBI")

# ---- L2 norm (cross-PC importance) ----
pc_mat <- as.matrix(ld[, ..pcs]); pc_mat[!is.finite(pc_mat)] <- 0
ld[, l2_norm := sqrt(rowSums(pc_mat^2))]

# ---- De-adjacency ----
deadjacent <- function(win_dt) {
  x <- copy(win_dt)
  x[, chrom_ord := suppressWarnings(as.integer(chrom))]
  x[is.na(chrom_ord), chrom_ord := 1e6]
  setorder(x, chrom_ord, start, end)
  x[, prev_end    := shift(end), by = chrom]
  x[, is_adj_prev := !is.na(prev_end) & (start == prev_end)]
  x[, run_id      := cumsum(!is_adj_prev), by = chrom]
  keep <- x[, .SD[floor((.N + 1L) / 2L)], by = .(chrom, run_id)]
  keep[, c("chrom_ord", "prev_end", "is_adj_prev", "run_id") := NULL]
  keep[]
}

get_top_windows <- function(score_col, frac) {
  x <- copy(ld[, .(window_id, chrom, start, end, score = get(score_col))])
  x <- x[is.finite(score)]
  x[, abs_score := abs(score)]
  setorder(x, -abs_score)
  deadjacent(x[seq_len(max(1L, ceiling(frac * nrow(x))))])
}

# ---- Gene annotation for all retained windows ----
message("Building gene annotation...")
build_gene_annotation <- function(win_dt) {
  valid_chr <- intersect(unique(win_dt$chrom), as.character(seqlevels(genes_gr)))
  z <- copy(win_dt[chrom %in% valid_chr])
  wgr <- GRanges(seqnames = z$chrom,
                 ranges = IRanges(start = z$start + 1L, end = z$end),
                 window_id = z$window_id)
  suppressWarnings(GenomeInfoDb::seqlevelsStyle(wgr) <- "NCBI")

  hits_ov <- findOverlaps(wgr, genes_gr, ignore.strand = TRUE)
  ov_dt <- data.table(
    window_id   = mcols(wgr)$window_id[queryHits(hits_ov)],
    gene_symbol = mcols(genes_gr)$symbol[subjectHits(hits_ov)],
    dist_bp     = 0L
  )
  ov_dt <- ov_dt[!is.na(gene_symbol) & nzchar(gene_symbol)]

  nn <- distanceToNearest(wgr, genes_gr, ignore.strand = TRUE)
  nn_dt <- data.table(
    window_id   = mcols(wgr)$window_id[queryHits(nn)],
    gene_symbol = mcols(genes_gr)$symbol[subjectHits(nn)],
    dist_bp     = mcols(nn)$distance
  )
  nn_dt <- nn_dt[!is.na(gene_symbol) & nzchar(gene_symbol)]
  no_ov_ids <- setdiff(nn_dt$window_id, ov_dt$window_id)
  rbindlist(list(ov_dt, nn_dt[window_id %in% no_ov_ids]), use.names = TRUE)
}

full_ann <- build_gene_annotation(ld[, .(window_id, chrom, start, end)])
bg_nearest_gene <- unique(full_ann[, gene_symbol])
message("Background: ", length(bg_nearest_gene), " genes (nearest gene, all retained windows)")

# ---- Symbol → Entrez conversion ----
sym2entrez <- function(symbols) {
  m <- mapIds(org.Hs.eg.db, keys = symbols, keytype = "SYMBOL",
              column = "ENTREZID", multiVals = "first")
  m[!is.na(m)]
}

# ---- MSigDB gene sets (local) ----
message("Loading MSigDB gene sets...")
msigdbr_safe <- function(...) {
  args <- list(...)
  tryCatch(
    do.call(msigdbr, c(list(species = "Homo sapiens"), args)),
    error = function(e) {
      names(args)[names(args) == "collection"]    <- "category"
      names(args)[names(args) == "subcollection"] <- "subcategory"
      do.call(msigdbr, c(list(species = "Homo sapiens"), args))
    }
  )
}
msig_h  <- suppressWarnings(msigdbr_safe(collection = "H"))
msig_c2 <- suppressWarnings(msigdbr_safe(collection = "C2", subcollection = "CP"))
msig_term2gene_h  <- data.table(gs_name     = msig_h$gs_name,
                                entrez_gene = as.character(msig_h$entrez_gene))
msig_term2gene_c2 <- data.table(gs_name     = msig_c2$gs_name,
                                entrez_gene = as.character(msig_c2$entrez_gene))

# ---- Build gene sets table ----
message("Building gene sets...")
geneset_rows <- list(); ii <- 1L
score_cols   <- c(setNames(pcs, pcs), c("cross_pc" = "l2_norm"))

for (sc_name in names(score_cols)) {
  for (frac in bands) {
    band_lbl <- band_labels[[as.character(frac)]]
    top_win  <- get_top_windows(score_cols[[sc_name]], frac)
    ann      <- merge(top_win, full_ann, by = "window_id", all.x = TRUE)
    for (dm in dist_modes) {
      genes_dm <- unique(ann[!is.na(gene_symbol) & nzchar(gene_symbol), gene_symbol])
      genes_dm <- genes_dm[!is.na(genes_dm) & nzchar(genes_dm)]
      geneset_rows[[ii]] <- data.table(
        pc = sc_name, band = band_lbl, band_frac = frac, dist_mode = dm,
        n_windows = nrow(top_win), n_genes = length(genes_dm), gene_symbol = genes_dm
      )
      message(sprintf("  %s | %s | %s: %d windows → %d genes",
                      sc_name, band_lbl, dm, nrow(top_win), length(genes_dm)))
      ii <- ii + 1L
    }
  }
}
genesets_dt <- rbindlist(geneset_rows, use.names = TRUE, fill = TRUE)
fwrite(genesets_dt, out_genesets, sep = "\t")
message("Wrote: ", out_genesets)

# ---- Enrichment runner ----
run_enrichment <- function(gene_symbols, bg_symbols, label) {
  rows <- list()
  query_ez <- sym2entrez(gene_symbols)
  bg_ez    <- sym2entrez(bg_symbols)
  if (length(query_ez) < 3L) return(NULL)

  for (ont in c("BP", "MF", "CC")) {
    res <- tryCatch(
      suppressMessages(enrichGO(
        gene          = query_ez,
        universe      = bg_ez,
        OrgDb         = org.Hs.eg.db,
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = p_cutoff,
        qvalueCutoff  = 0.2,
        readable      = TRUE
      )),
      error = function(e) NULL
    )
    if (!is.null(res) && nrow(res@result) > 0L) {
      r <- as.data.table(res@result)[p.adjust < p_cutoff]
      if (nrow(r) > 0L) { r[, source := paste0("GO:", ont)]; rows <- c(rows, list(r)) }
    }
  }

  res_h <- tryCatch(
    suppressMessages(enricher(
      gene          = query_ez,
      universe      = bg_ez,
      TERM2GENE     = as.data.frame(msig_term2gene_h),
      pAdjustMethod = "BH",
      pvalueCutoff  = p_cutoff,
      qvalueCutoff  = 0.2
    )),
    error = function(e) NULL
  )
  if (!is.null(res_h) && nrow(res_h@result) > 0L) {
    r <- as.data.table(res_h@result)[p.adjust < p_cutoff]
    if (nrow(r) > 0L) { r[, source := "MSigDB_Hallmarks"]; rows <- c(rows, list(r)) }
  }

  res_c2 <- tryCatch(
    suppressMessages(enricher(
      gene          = query_ez,
      universe      = bg_ez,
      TERM2GENE     = as.data.frame(msig_term2gene_c2),
      pAdjustMethod = "BH",
      pvalueCutoff  = p_cutoff,
      qvalueCutoff  = 0.2
    )),
    error = function(e) NULL
  )
  if (!is.null(res_c2) && nrow(res_c2@result) > 0L) {
    r <- as.data.table(res_c2@result)[p.adjust < p_cutoff]
    if (nrow(r) > 0L) { r[, source := "MSigDB_C2"]; rows <- c(rows, list(r)) }
  }

  if (length(rows) == 0L) return(NULL)
  out <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  if ("Description" %in% names(out) && !"term_name" %in% names(out))
    setnames(out, "Description", "term_name")
  if ("Count" %in% names(out) && !"intersection_size" %in% names(out))
    setnames(out, "Count", "intersection_size")
  if ("pvalue" %in% names(out) && !"p_value" %in% names(out))
    setnames(out, "pvalue", "p_value")
  out
}

# ---- Run all enrichments ----
message("Running enrichments...")
enrich_rows <- list(); jj <- 1L

for (sc_name in names(score_cols)) {
  for (frac in bands) {
    band_lbl <- band_labels[[as.character(frac)]]
    for (dm in dist_modes) {
      genes_dm <- genesets_dt[pc == sc_name & band == band_lbl & dist_mode == dm, gene_symbol]
      bg_dm    <- bg_nearest_gene
      message(sprintf("  Enriching %s | %s | %s (%d genes)...",
                      sc_name, band_lbl, dm, length(genes_dm)))
      res <- run_enrichment(genes_dm, bg_dm, label = paste(sc_name, band_lbl, dm))
      if (!is.null(res)) {
        res[, `:=`(pc = sc_name, band = band_lbl, band_frac = frac,
                   dist_mode = dm, n_query_genes = length(genes_dm))]
        enrich_rows[[jj]] <- res; jj <- jj + 1L
        message(sprintf("    → %d significant terms", nrow(res)))
      } else {
        message("    → no significant terms")
      }
    }
  }
}

if (length(enrich_rows) > 0L) {
  enrich_dt <- rbindlist(enrich_rows, use.names = TRUE, fill = TRUE)
  fwrite(enrich_dt, out_enrich, sep = "\t")
  message("Wrote: ", out_enrich, " (", nrow(enrich_dt), " significant terms total)")
} else {
  message("No significant enrichments found.")
  fwrite(data.table(), out_enrich, sep = "\t")
}

message("compute-full-gene-enrichment complete.")
