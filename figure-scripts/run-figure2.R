#!/usr/bin/env Rscript
# Produces: figures/figure2.pdf (Manhattan + locus zooms),
#           figures/figureS15_figure2_extra_zooms.pdf,
#           results/Supplementary_PC_Top_Loadings_figure2.xlsx (.tsv fallback)
# Pre-computed inputs: results/pca_loadings_50k.tsv, data/gencode.v19.genes.protein_coding.rds
# Run: Rscript post-process-scripts/run-figure2.R [--supp-only]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(RColorBrewer)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(grid)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

args <- commandArgs(trailingOnly = TRUE)
supp_only <- "--supp-only" %in% args

fig2_pdf          <- file.path("figures", "figure2.pdf")
fig2_supp_pdf     <- file.path("figures", "figureS15_figure2_extra_zooms.pdf")
fig2_table_xlsx   <- file.path("results", "Supplementary_PC_Top_Loadings_figure2.xlsx")
fig2_table_tsv    <- file.path("results", "Supplementary_PC_Top_Loadings_figure2.tsv")
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

include_zoom_regions <- TRUE
pc_cols <- c("PC1" = "#1f78b4", "PC2" = "#33a02c", "PC3" = "#e31a1c", "PC4" = "#984ea3")

# ---- Load pca_loadings ----
message("Figure 2: loading PCA loadings...")
pca_loadings_file <- file.path("results", "pca_loadings_50k.tsv")
if (!file.exists(pca_loadings_file)) stop("Missing: ", pca_loadings_file, "\nRun run-figure1.R first.")
pca_loadings <- fread(pca_loadings_file)

source(file.path("post-process-scripts", "helpers-plotting.R"))

chrom_plot_palette <- colorRampPalette(brewer.pal(8, "Set1"))(22)

# ---- Load gene annotation ----
genes_rds <- file.path(project_root, "data", "gencode.v19.genes.protein_coding.rds")
if (!file.exists(genes_rds)) stop("Missing gene annotation RDS: ", genes_rds)
message("Figure 2: loading gene annotation...")
genes_gr <- readRDS(genes_rds)

# ---- Layered locus zoom helper ----
make_layered_locus_plot <- function(
    loadings_dt, genes_gr, chr, start_bp, end_bp, tag,
    pad_bp = 500000L, suppress_symbols = character(), xlim_mb = NULL, y_limits = NULL
) {
  normalize_symbol <- function(x) {
    x <- as.character(x)
    gsub("[\u2010\u2011\u2012\u2013\u2014\u2015\u2212]", "-", x, perl = TRUE)
  }
  zstart <- max(0, as.integer(start_bp - pad_bp))
  zend   <- as.integer(end_bp + pad_bp)
  dt <- as.data.table(loadings_dt)[
    chrom == as.character(chr) & start >= zstart & end <= zend,
    .(chrom, start, end, PC1, PC2, PC3, PC4)
  ]
  if (nrow(dt) == 0L) stop("No windows found for locus plot: chr", chr, ":", zstart, "-", zend)
  dt[, pos_mb := (start + end) / 2e6]
  long_dt <- melt(dt, id.vars = c("chrom", "start", "end", "pos_mb"),
                  measure.vars = c("PC1", "PC2", "PC3", "PC4"),
                  variable.name = "pc", value.name = "loading")

  gene_track <- make_gene_track(genes_gr, chr = chr, start_bp = zstart, end_bp = zend, pos_start_chr = 0)
  label_gene_track <- gene_track
  if (!is.null(label_gene_track) && length(suppress_symbols) > 0L) {
    label_gene_track[, symbol_norm := normalize_symbol(symbol)]
    label_gene_track <- label_gene_track[!symbol_norm %in% normalize_symbol(suppress_symbols)]
    label_gene_track[, symbol_norm := NULL]
  }

  x_breaks <- pretty(if (is.null(xlim_mb)) range(long_dt$pos_mb, na.rm = TRUE) else xlim_mb, n = 5)

  p <- ggplot(long_dt, aes(x = pos_mb, y = loading, colour = pc)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) +
    geom_point(size = 1.0, alpha = 0.85) +
    scale_color_manual(values = pc_cols, drop = FALSE) +
    scale_x_continuous(name = paste0("Chr", chr, " position (Mb)"), breaks = x_breaks, limits = xlim_mb) +
    ylab("Loading") +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "bottom", legend.title = element_blank(),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(margin = margin(t = 10)),
      plot.margin = margin(5.5, 5.5, 16, 5.5)
    ) +
    labs(tag = tag, colour = NULL)

  if (!is.null(gene_track) && nrow(gene_track) > 0L) {
    y_abs_max <- max(abs(long_dt$loading), na.rm = TRUE)
    if (!is.finite(y_abs_max) || y_abs_max <= 0) y_abs_max <- 0.01
    lane_height <- (2 * y_abs_max) * 0.06
    if (!is.finite(lane_height) || lane_height <= 0) lane_height <- 0.01
    ymin_base <- -y_abs_max - lane_height * (max(gene_track$lane, na.rm = TRUE) + 1.1)
    ylim_abs  <- max(y_abs_max * 1.05, abs(ymin_base - 0.2 * lane_height))
    p <- p +
      geom_rect(data = gene_track, inherit.aes = FALSE,
                aes(xmin = gene_start_mb, xmax = gene_end_mb,
                    ymin = ymin_base + (lane - 1) * lane_height,
                    ymax = ymin_base + lane * lane_height),
                fill = "grey80", alpha = 0.9, colour = "grey55", linewidth = 0.2) +
      coord_cartesian(xlim = xlim_mb,
                      ylim = if (is.null(y_limits)) c(-ylim_abs, ylim_abs) else y_limits,
                      clip = "off")
    if (!is.null(label_gene_track) && nrow(label_gene_track) > 0L) {
      p <- p + geom_text_repel(
        data = label_gene_track, inherit.aes = FALSE,
        aes(x = gene_mid_mb, y = ymin_base + (lane - 0.5) * lane_height, label = symbol),
        size = 2.2, seed = 1, min.segment.length = 0, max.overlaps = Inf,
        box.padding = 0.25, point.padding = 0.1, segment.alpha = 0.5, segment.size = 0.2, colour = "grey25"
      )
    }
  } else {
    y_abs_max <- max(abs(long_dt$loading), na.rm = TRUE)
    if (!is.finite(y_abs_max) || y_abs_max <= 0) y_abs_max <- 0.01
    p <- p + coord_cartesian(xlim = xlim_mb,
                              ylim = if (is.null(y_limits)) c(-1.05 * y_abs_max, 1.05 * y_abs_max) else y_limits,
                              clip = "off")
  }
  p
}

# ---- Manhattan panels ----
message("Figure 2 / Panel A: PC1 Manhattan...")
res_pc1 <- plot_pc_loadings(pca_loadings, pc = "PC1", genes_dt = genes_gr)
p_pc1   <- res_pc1$plot; table_pc1 <- res_pc1$table

message("Figure 2 / Panel B: PC2 Manhattan...")
res_pc2 <- plot_pc_loadings(pca_loadings, pc = "PC2", genes_dt = genes_gr)
p_pc2   <- res_pc2$plot; table_pc2 <- res_pc2$table

message("Figure 2 / Panel C: PC3 Manhattan...")
res_pc3 <- plot_pc_loadings(pca_loadings, pc = "PC3", genes_dt = genes_gr)
p_pc3   <- res_pc3$plot; table_pc3 <- res_pc3$table

message("Figure 2 / Table: PC4 top loadings...")
res_pc4 <- plot_pc_loadings(pca_loadings, pc = "PC4", genes_dt = genes_gr)
table_pc4 <- res_pc4$table

# ---- Zoom panels ----
if (isTRUE(include_zoom_regions)) {
  message("Figure 2: building zoomed loci panels...")

  p_d <- make_layered_locus_plot(pca_loadings, genes_gr, chr = 16, start_bp = 53800000L, end_bp = 53850000L, tag = "d", pad_bp = 975000L)
  p_e <- make_layered_locus_plot(pca_loadings, genes_gr, chr = 11, start_bp = 112850000L, end_bp = 112900000L, tag = "e", pad_bp = 975000L)
  p_f <- make_layered_locus_plot(pca_loadings, genes_gr, chr = 9,  start_bp = 17150000L, end_bp = 17200000L, tag = "f", pad_bp = 975000L)

  figure2 <- (p_pc1 + labs(tag = "a")) /
    (p_pc2 + labs(tag = "b")) /
    (p_pc3 + labs(tag = "c")) /
    ((p_d + theme(legend.position = "none")) |
     (p_e + theme(legend.position = "none")) |
     (p_f + theme(legend.position = "bottom"))) +
    plot_layout(heights = c(1, 1, 1, 1)) &
    theme(plot.tag = element_text(size = 18))

  p_sa <- make_layered_locus_plot(pca_loadings, genes_gr, chr = 2,  start_bp = 0L,         end_bp = 1000000L,  tag = "a", suppress_symbols = "FAM110C", y_limits = c(-0.045, 0.045))
  p_sb <- make_layered_locus_plot(pca_loadings, genes_gr, chr = 8,  start_bp = 56500000L,  end_bp = 58000000L, tag = "b", pad_bp = 250000L)
  p_sc <- make_layered_locus_plot(pca_loadings, genes_gr, chr = 13, start_bp = 46950000L,  end_bp = 47000000L, tag = "c", pad_bp = 975000L)
  p_sd <- make_layered_locus_plot(pca_loadings, genes_gr, chr = 18, start_bp = 57500000L,  end_bp = 59000000L, tag = "d", pad_bp = 250000L)
  p_se <- make_layered_locus_plot(pca_loadings, genes_gr, chr = 17, start_bp = 28900000L,  end_bp = 28950000L, tag = "e", pad_bp = 975000L,
    suppress_symbols = c("RAB11FIP4", "CTD-2370N5.3", "EVI2B", "D29D5", "TMIGD1", "TBC1D29",
                         "AC003101.1", "ATAD5", "ANKRD13B", "SSH2", "RNF135"))
  p_sf <- make_layered_locus_plot(pca_loadings, genes_gr, chr = 4,  start_bp = 17950000L,  end_bp = 18000000L, tag = "f", pad_bp = 975000L)

  figure2_supp <- wrap_plots(
    list(p_sa + theme(legend.position = "none"),
         p_sb + theme(legend.position = "none"),
         p_sc + theme(legend.position = "none"),
         p_sd + theme(legend.position = "none"),
         p_se + theme(legend.position = "none"),
         p_sf + theme(legend.position = "bottom")),
    ncol = 2
  ) & theme(plot.tag = element_text(size = 18))
} else {
  figure2 <- (p_pc1 + labs(tag = "a")) /
    (p_pc2 + labs(tag = "b")) /
    (p_pc3 + labs(tag = "c")) &
    theme(plot.tag = element_text(size = 18))
}

# ---- Save figures ----
if (!supp_only) {
  ggsave(fig2_pdf, figure2, height = 12, width = 12)
  message("Figure 2 complete: ", fig2_pdf)
}
if (isTRUE(include_zoom_regions)) {
  ggsave(fig2_supp_pdf, figure2_supp, height = 14.4, width = 9.5)
  message("FigureS15 complete: ", fig2_supp_pdf)
}

# ---- Supplementary top-loadings tables ----
if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "PC1"); openxlsx::writeData(wb, "PC1", table_pc1)
  openxlsx::addWorksheet(wb, "PC2"); openxlsx::writeData(wb, "PC2", table_pc2)
  openxlsx::addWorksheet(wb, "PC3"); openxlsx::writeData(wb, "PC3", table_pc3)
  openxlsx::addWorksheet(wb, "PC4"); openxlsx::writeData(wb, "PC4", table_pc4)
  openxlsx::saveWorkbook(wb, fig2_table_xlsx, overwrite = TRUE)
  message("Supplementary tables: ", fig2_table_xlsx)
} else {
  supp_dt <- rbindlist(list(
    data.table(pc = "PC1", table_pc1),
    data.table(pc = "PC2", table_pc2),
    data.table(pc = "PC3", table_pc3),
    data.table(pc = "PC4", table_pc4)
  ), use.names = TRUE, fill = TRUE)
  fwrite(supp_dt, fig2_table_tsv, sep = "\t")
  message("openxlsx not installed; wrote TSV fallback: ", fig2_table_tsv)
}
