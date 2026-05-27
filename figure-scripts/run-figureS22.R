#!/usr/bin/env Rscript
# Produces: figures/figureS22_figure3_extra_zooms_plus100kb.pdf
# Pre-computed inputs:
#   results/pca_loadings_50k.tsv,
#   data/gencode.v19.genes.protein_coding.rds,
#   results/gwas_removed_distance_topload_50k_plus100kb/reduced_pca_window_loadings.tsv,
#   results/gwas_removed_distance_topload_50k_plus100kb/removed_windows_gwas_plus100kb.tsv
# Run: Rscript post-process-scripts/run-figureS22.R [flank_kb]
#   default flank_kb = 100

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

args     <- commandArgs(trailingOnly = TRUE)
flank_kb <- if (length(args) >= 1L) as.integer(args[[1]]) else 100L
if (!is.finite(flank_kb) || flank_kb < 0L) stop("flank_kb must be a non-negative integer.")
flank_label <- paste0(flank_kb, "kb")

out_pdf <- if (flank_kb == 100L) {
  file.path("figures", "figureS22_figure3_extra_zooms_plus100kb.pdf")
} else {
  file.path("figures", paste0("figureS22_figure3_extra_zooms_plus", flank_label, ".pdf"))
}
in_file      <- file.path("results", paste0("gwas_removed_distance_topload_50k_plus", flank_label), "reduced_pca_window_loadings.tsv")
removed_file <- file.path("results", paste0("gwas_removed_distance_topload_50k_plus", flank_label), paste0("removed_windows_gwas_plus", flank_label, ".tsv"))
genes_rds    <- file.path(project_root, "data", "gencode.v19.genes.protein_coding.rds")
for (f in c(in_file, removed_file, genes_rds)) {
  if (!file.exists(f)) stop("Missing input: ", f)
}
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

# ---- Load pca_loadings ----
pca_loadings_file <- file.path("results", "pca_loadings_50k.tsv")
if (!file.exists(pca_loadings_file)) stop("Missing: ", pca_loadings_file, "\nRun run-figure1.R first.")
pca_loadings <- fread(pca_loadings_file)

pc_cols <- c("PC1" = "#1f78b4", "PC2" = "#33a02c", "PC3" = "#e31a1c", "PC4" = "#984ea3")

make_gene_track_simple <- function(genes_gr, chr, start_bp, end_bp) {
  g <- genes_gr[as.character(seqnames(genes_gr)) == as.character(chr) &
                  start(genes_gr) <= end_bp & end(genes_gr) >= start_bp]
  if (length(g) == 0L) return(NULL)
  dt <- as.data.table(as.data.frame(g))
  dt[, `:=`(gene_start_mb = start / 1e6, gene_end_mb = end / 1e6, gene_mid_mb = (start + end) / 2e6,
            symbol = as.character(symbol))]
  setorder(dt, start, end)
  lane_ends <- numeric()
  dt[, lane := {
    s <- start; assigned <- NA_integer_
    if (length(lane_ends) == 0L) { lane_ends <<- end; assigned <- 1L
    } else {
      for (ii in seq_along(lane_ends)) {
        if (s > lane_ends[ii]) { lane_ends[ii] <<- end; assigned <- ii; break }
      }
      if (is.na(assigned)) { lane_ends <<- c(lane_ends, end); assigned <- length(lane_ends) }
    }
    assigned
  }, by = seq_len(nrow(dt))]
  dt[]
}

make_layered_locus_plot <- function(loadings_dt, full_loadings_dt, removed_dt, genes_gr,
                                    chr, start_bp, end_bp, tag, pad_bp = 975000L) {
  zstart <- max(0L, as.integer(start_bp - pad_bp))
  zend   <- as.integer(end_bp + pad_bp)
  dt <- as.data.table(loadings_dt)[
    chrom == as.integer(chr) & start >= zstart & end <= zend,
    .(chrom, start, end, PC1, PC2, PC3, PC4)
  ]
  if (nrow(dt) == 0L) stop("No windows found for reduced locus plot: chr", chr, ":", zstart, "-", zend)
  dt[, pos_mb := (start + end) / 2e6]
  long_dt <- melt(dt, id.vars = c("chrom", "start", "end", "pos_mb"),
                  measure.vars = c("PC1", "PC2", "PC3", "PC4"),
                  variable.name = "pc", value.name = "loading")

  full_dt <- as.data.table(full_loadings_dt)[
    chrom == as.integer(chr) & start >= zstart & end <= zend,
    .(window_id = paste(chrom, start, end, sep = "_"), chrom, start, end, PC1, PC2, PC3, PC4)
  ]
  present_ids    <- unique(paste(dt$chrom, dt$start, dt$end, sep = "_"))
  missing_full_dt <- full_dt[!window_id %in% present_ids]
  missing_long_dt <- if (nrow(missing_full_dt) > 0L) {
    missing_full_dt <- merge(missing_full_dt, removed_dt[, .(window_id, gws_any)], by = "window_id", all.x = TRUE)
    missing_full_dt[is.na(gws_any), gws_any := FALSE]
    missing_full_dt[, pos_mb := (start + end) / 2e6]
    melt(missing_full_dt, id.vars = c("window_id", "chrom", "start", "end", "pos_mb", "gws_any"),
         measure.vars = c("PC1", "PC2", "PC3", "PC4"), variable.name = "pc", value.name = "loading")
  } else {
    data.table(window_id = character(), chrom = integer(), start = integer(), end = integer(),
               pos_mb = numeric(), gws_any = logical(), pc = character(), loading = numeric())
  }

  gene_track <- make_gene_track_simple(genes_gr, chr = chr, start_bp = zstart, end_bp = zend)
  xlim_mb    <- c(zstart / 1e6, zend / 1e6)
  x_breaks   <- pretty(xlim_mb, n = 5)

  p <- ggplot(long_dt, aes(x = pos_mb, y = loading, colour = pc)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) +
    geom_point(data = missing_long_dt, inherit.aes = FALSE,
               aes(x = pos_mb, y = loading, colour = pc, fill = gws_any),
               shape = 21, size = 1.1, stroke = 0.6, alpha = 0.95) +
    geom_point(size = 0.95, alpha = 0.85) +
    scale_color_manual(values = pc_cols, drop = FALSE) +
    scale_fill_manual(values = c(`TRUE` = "grey65", `FALSE` = "white"), guide = "none") +
    scale_x_continuous(name = paste0("Chr", chr, " position (Mb)"), breaks = x_breaks, limits = xlim_mb) +
    ylab("Loading") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text(margin = margin(t = 10)),
          plot.margin = margin(5.5, 5.5, 16, 5.5)) +
    labs(tag = tag, colour = NULL)

  y_abs_max <- max(abs(c(long_dt$loading, missing_long_dt$loading)), na.rm = TRUE)
  if (!is.finite(y_abs_max) || y_abs_max <= 0) y_abs_max <- 0.01

  if (!is.null(gene_track) && nrow(gene_track) > 0L) {
    lane_height <- (2 * y_abs_max) * 0.06
    ymin_base   <- -y_abs_max - lane_height * (max(gene_track$lane, na.rm = TRUE) + 1.1)
    ylim_abs    <- max(y_abs_max * 1.05, abs(ymin_base - 0.2 * lane_height))
    p <- p +
      geom_rect(data = gene_track, inherit.aes = FALSE,
                aes(xmin = gene_start_mb, xmax = gene_end_mb,
                    ymin = ymin_base + (lane - 1) * lane_height, ymax = ymin_base + lane * lane_height),
                fill = "grey80", alpha = 0.9, colour = "grey55", linewidth = 0.2) +
      geom_text_repel(data = gene_track, inherit.aes = FALSE,
                      aes(x = gene_mid_mb, y = ymin_base + (lane - 0.5) * lane_height, label = symbol),
                      size = 2.2, seed = 1, min.segment.length = 0, max.overlaps = Inf,
                      box.padding = 0.25, point.padding = 0.1, segment.alpha = 0.5, segment.size = 0.2, colour = "grey25") +
      coord_cartesian(xlim = xlim_mb, ylim = c(-ylim_abs, ylim_abs), clip = "off")
  } else {
    p <- p + coord_cartesian(xlim = xlim_mb, ylim = c(-1.05 * y_abs_max, 1.05 * y_abs_max), clip = "off")
  }
  p
}

genes_gr          <- readRDS(genes_rds); seqlevelsStyle(genes_gr) <- "NCBI"
reduced_loadings  <- fread(in_file)
full_loadings     <- as.data.table(pca_loadings)
full_loadings[, `:=`(chrom = as.integer(chrom), start = as.integer(start), end = as.integer(end))]
removed_windows   <- fread(removed_file)

loci <- list(
  list(tag = "a", gene = "KCTD3",   chr = 1L,  start = 215650000L, end = 215700000L),
  list(tag = "b", gene = "HS6ST1",  chr = 2L,  start = 129300000L, end = 129350000L),
  list(tag = "c", gene = "KMT2C",   chr = 7L,  start = 151800000L, end = 151850000L),
  list(tag = "d", gene = "AGTPBP1", chr = 9L,  start = 88300000L,  end = 88350000L),
  list(tag = "e", gene = "DIP2C",   chr = 10L, start = 400000L,    end = 450000L),
  list(tag = "f", gene = "C17orf51",chr = 17L, start = 21450000L,  end = 21500000L)
)

plots <- lapply(loci, function(loc) {
  make_layered_locus_plot(reduced_loadings, full_loadings, removed_windows, genes_gr,
                          chr = loc$chr, start_bp = loc$start, end_bp = loc$end,
                          tag = loc$tag, pad_bp = 975000L) + ggtitle(loc$gene)
})
for (ii in seq_along(plots)) {
  if (ii < length(plots)) plots[[ii]] <- plots[[ii]] + theme(legend.position = "none")
}

fig <- wrap_plots(plots, ncol = 3) &
  theme(plot.tag = element_text(size = 18), plot.title = element_text(size = 11, face = "plain"))

ggsave(out_pdf, fig, width = 14, height = 6.8)
message("FigureS22 complete: ", out_pdf)
