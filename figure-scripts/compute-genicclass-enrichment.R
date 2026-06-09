#!/usr/bin/env Rscript
# Produces: results/genicclass_rank_enrichment_50k/
#   cumulative_genicclass_enrichment_curves_pc1_pc4.tsv  (input for run-figure4.R)
# Pre-computed inputs: results/pca_loadings_50k.tsv,
#                      data/gencode.v19.annotation.gtf.gz
# Run: Rscript figure-scripts/compute-genicclass-enrichment.R
# Runtime: ~5-10 min (100 null reps × 4 PCs × 4 annotation classes × many fracs)

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

pcs         <- c("PC1", "PC2", "PC3", "PC4")
n_reps      <- 100L
promoter_bp <- 2000L
set.seed(20260329)

gtf_file <- file.path("data", "gencode.v19.annotation.gtf.gz")
if (!file.exists(gtf_file)) stop("Missing GTF: ", gtf_file)

out_dir     <- file.path("results", "genicclass_rank_enrichment_50k")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_cum     <- file.path(out_dir, "cumulative_genicclass_enrichment_curves_pc1_pc4.tsv")

sample_chr_matched <- function(z_dt, chr_counts) {
  out_idx <- integer()
  for (i in seq_len(nrow(chr_counts))) {
    chr_i <- chr_counts$chr[i]; n_i <- chr_counts$N[i]
    cand  <- which(z_dt$chr == chr_i)
    if (length(cand) < n_i) stop("Insufficient candidates for chr ", chr_i)
    out_idx <- c(out_idx, sample(cand, n_i, replace = FALSE))
  }
  out_idx
}

# ---- Load pca_loadings ----
message("Loading pca_loadings_50k.tsv...")
pca_loadings_file <- file.path("results", "pca_loadings_50k.tsv")
if (!file.exists(pca_loadings_file)) stop("Missing: ", pca_loadings_file, "\nRun run-figure1.R first.")
pca_loadings <- fread(pca_loadings_file)

ld <- as.data.table(pca_loadings)
if ("chrom" %in% names(ld) && !"chr" %in% names(ld)) setnames(ld, "chrom", "chr")
if (!all(c("chr", "start", "end") %in% names(ld))) {
  if (!"region_label" %in% names(ld)) stop("Need chr/start/end or region_label in pca_loadings.")
  ld[, c("chr", "start", "end") := tstrsplit(region_label, "_", fixed = TRUE)]
}
ld[, `:=`(chr = gsub("^chr", "", as.character(chr)), start = as.integer(start), end = as.integer(end))]
ld <- ld[is.finite(start) & is.finite(end) & end > start]
ld[, window_id := paste(chr, start, end, sep = "_")]
ld[, chr_ord := suppressWarnings(as.integer(chr))]
ld[is.na(chr_ord), chr_ord := 1e6]
setorder(ld, chr_ord, start, end); ld[, chr_ord := NULL]

# ---- Parse GTF ----
message("Reading GENCODE v19 GTF...")
extract_attr <- function(x, key) {
  pat <- paste0(key, ' "([^"]+)"')
  m   <- regexec(pat, x)
  res <- regmatches(x, m)
  vapply(res, function(z) if (length(z) >= 2L) z[2] else NA_character_, character(1))
}
gtf <- fread(cmd = paste("gzip -dc", shQuote(gtf_file)), sep = "\t", header = FALSE, quote = "",
             data.table = TRUE, showProgress = FALSE)
setnames(gtf, c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))
gtf <- gtf[!startsWith(seqname, "#")]
gtf[, seqname       := gsub("^chr", "", seqname)]
gtf[, gene_type      := extract_attr(attribute, "gene_type")]
gtf[, transcript_type := extract_attr(attribute, "transcript_type")]

genes <- gtf[feature == "gene" & gene_type == "protein_coding",
             .(chr = seqname, start = as.integer(start), end = as.integer(end), strand = strand)]
exons <- gtf[feature == "exon" & (gene_type == "protein_coding" | transcript_type == "protein_coding"),
             .(chr = seqname, start = as.integer(start), end = as.integer(end))]
if (nrow(genes) == 0L || nrow(exons) == 0L) stop("Failed to derive genes/exons from GTF.")

# ---- Build annotation classes ----
message("Building promoter/exon/intron/intergenic classes...")
genes_gr    <- GRanges(seqnames = genes$chr, ranges = IRanges(start = genes$start, end = genes$end), strand = genes$strand)
tss_pos     <- ifelse(as.character(strand(genes_gr)) == "-", end(genes_gr), start(genes_gr))
promoter_gr <- reduce(GRanges(seqnames = seqnames(genes_gr),
                               ranges = IRanges(start = pmax(1L, tss_pos - promoter_bp), end = tss_pos + promoter_bp)))
gene_body_gr <- reduce(GRanges(seqnames = genes$chr, ranges = IRanges(start = genes$start, end = genes$end)))
exon_gr      <- reduce(GRanges(seqnames = exons$chr, ranges = IRanges(start = exons$start, end = exons$end)))
intron_gr    <- setdiff(gene_body_gr, exon_gr, ignore.strand = TRUE)

win_gr       <- GRanges(seqnames = ld$chr, ranges = IRanges(start = ld$start + 1L, end = ld$end))
ov_promoter  <- overlapsAny(win_gr, promoter_gr,  ignore.strand = TRUE)
ov_exon      <- overlapsAny(win_gr, exon_gr,      ignore.strand = TRUE)
ov_intron    <- overlapsAny(win_gr, intron_gr,    ignore.strand = TRUE)

class_levels <- c("promoter", "exonic", "intronic", "intergenic")
ld[, class := "intergenic"]
ld[ov_intron,   class := "intronic"]
ld[ov_exon,     class := "exonic"]
ld[ov_promoter, class := "promoter"]
for (nm in class_levels) ld[, (nm) := class == nm]
message("Class counts:")
print(ld[, .N, by = class][match(class_levels, class)])

# ---- Rank curves ----
message("Computing rank curves with chromosome-matched null ribbons...")
bin_edges  <- seq(0, 0.10, by = 0.01)
bin_starts <- head(bin_edges, -1); bin_ends <- tail(bin_edges, -1)
cum_fracs  <- c(seq(0.01, 0.10, by = 0.01), 0.15, 0.20, 0.30, 0.40, 0.50, 0.75, 1.00)

curve_rows <- list(); curve_rows_cum <- list(); rank_rows <- list()
ii <- jj <- kk <- 1L

for (pc in pcs) {
  message("  Ranking and bins for ", pc, "...")
  z <- ld[, c("window_id", "chr", "start", "end", class_levels), with = FALSE]
  z[, loading := ld[[pc]]]; z <- z[is.finite(loading)]
  z[, abs_loading := abs(loading)]; setorder(z, -abs_loading)
  z[, rank := seq_len(.N)]; z[, rank_frac := rank / .N]

  for (nm in class_levels) {
    z_ann <- z[, .(window_id, chr, start, end, loading, abs_loading, rank, rank_frac, overlap = get(nm))]
    z_ann[, `:=`(pc = pc, annotation = nm)]
    rank_rows[[ii]] <- z_ann; ii <- ii + 1L

    base_overlap <- mean(z_ann$overlap)
    n <- nrow(z_ann)

    for (tf in cum_fracs) {
      n_top  <- max(1L, ceiling(tf * n))
      top_dt <- z_ann[seq_len(n_top)]
      chr_counts <- top_dt[, .N, by = chr]
      obs    <- mean(top_dt$overlap)
      enrich <- if (isTRUE(base_overlap > 0)) obs / base_overlap else NA_real_
      null_overlap <- numeric(n_reps)
      for (r in seq_len(n_reps)) null_overlap[r] <- mean(z_ann$overlap[sample_chr_matched(z_ann, chr_counts)])
      curve_rows_cum[[kk]] <- data.table(
        pc = pc, annotation = nm, top_frac = tf, top_pct = as.integer(round(100 * tf)), n_top = n_top,
        observed_overlap = obs, genomewide_overlap = base_overlap, enrichment_ratio = enrich,
        null_overlap_mean = mean(null_overlap, na.rm = TRUE),
        null_overlap_lo   = quantile(null_overlap, 0.025, na.rm = TRUE),
        null_overlap_hi   = quantile(null_overlap, 0.975, na.rm = TRUE),
        null_enrichment_mean = if (isTRUE(base_overlap > 0)) mean(null_overlap, na.rm = TRUE) / base_overlap else NA_real_,
        null_enrichment_lo   = if (isTRUE(base_overlap > 0)) quantile(null_overlap, 0.025, na.rm = TRUE) / base_overlap else NA_real_,
        null_enrichment_hi   = if (isTRUE(base_overlap > 0)) quantile(null_overlap, 0.975, na.rm = TRUE) / base_overlap else NA_real_
      ); kk <- kk + 1L
    }

    for (bb in seq_along(bin_starts)) {
      start_idx <- floor(bin_starts[bb] * n) + 1L; end_idx <- ceiling(bin_ends[bb] * n)
      if (end_idx < start_idx) next
      bin_dt     <- z_ann[start_idx:end_idx]
      chr_counts <- bin_dt[, .N, by = chr]
      obs    <- mean(bin_dt$overlap)
      enrich <- if (isTRUE(base_overlap > 0)) obs / base_overlap else NA_real_
      null_overlap <- numeric(n_reps)
      for (r in seq_len(n_reps)) null_overlap[r] <- mean(z_ann$overlap[sample_chr_matched(z_ann, chr_counts)])
      curve_rows[[jj]] <- data.table(
        pc = pc, annotation = nm,
        bin_start_frac = bin_starts[bb], bin_end_frac = bin_ends[bb],
        bin_start_pct  = 100 * bin_starts[bb], bin_end_pct = 100 * bin_ends[bb],
        bin_mid_pct    = 100 * (bin_starts[bb] + bin_ends[bb]) / 2,
        n_bin = length(start_idx:end_idx),
        observed_overlap = obs, genomewide_overlap = base_overlap, enrichment_ratio = enrich,
        null_overlap_mean = mean(null_overlap, na.rm = TRUE),
        null_overlap_lo   = quantile(null_overlap, 0.025, na.rm = TRUE),
        null_overlap_hi   = quantile(null_overlap, 0.975, na.rm = TRUE),
        null_enrichment_mean = if (isTRUE(base_overlap > 0)) mean(null_overlap, na.rm = TRUE) / base_overlap else NA_real_,
        null_enrichment_lo   = if (isTRUE(base_overlap > 0)) quantile(null_overlap, 0.025, na.rm = TRUE) / base_overlap else NA_real_,
        null_enrichment_hi   = if (isTRUE(base_overlap > 0)) quantile(null_overlap, 0.975, na.rm = TRUE) / base_overlap else NA_real_
      ); jj <- jj + 1L
    }
  }
}

fwrite(rbindlist(curve_rows_cum, use.names = TRUE, fill = TRUE), out_cum, sep = "\t")
message("Wrote: ", out_cum)
message("Done.")
