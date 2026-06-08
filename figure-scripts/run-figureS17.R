#!/usr/bin/env Rscript
# Produces: figures/figureS17_chr9_chr17_pval_scatter.pdf
# Pre-computed inputs: results/all-trait-50k-mean-pvals.tsv
# Run: Rscript figure-scripts/run-figureS17.R
# Shows per-trait -log10(p) across a 1 Mb window centered on:
#   Chr9  ~16.65-17.65 Mb (top-loaded PC1 locus)
#   Chr17 ~20.9-21.9 Mb  (top-loaded PC1 locus)
# Background points: all traits; purple dots: trait with maximum -log10(p) per window.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
  library(patchwork)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

in_pvals <- file.path("results", "all-trait-50k-mean-pvals.tsv")
out_fig  <- file.path("figures", "figureS17_chr9_chr17_pval_scatter.pdf")
if (!file.exists(in_pvals)) stop("Missing input: ", in_pvals, "\nRun build-all-trait-matrices.R first.")
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

pv <- fread(in_pvals)
pv[, chrom := gsub('"', '', chrom)]

chrom_plot_palette <- colorRampPalette(brewer.pal(8, "Set1"))(22)
col_bg9  <- chrom_plot_palette[9]
col_bg17 <- chrom_plot_palette[17]
col_max  <- "#984EA3"

make_scatter <- function(chrom_id, center_start, chrom_col, title_str) {
  win_starts <- seq(center_start - 10L * 50000L, center_start + 9L * 50000L, by = 50000L)
  sub        <- pv[chrom == chrom_id & start %in% win_starts]
  trait_cols <- setdiff(names(sub), c("chrom", "start", "end"))
  long <- melt(sub, id.vars = c("chrom", "start", "end"),
               measure.vars = trait_cols, variable.name = "trait", value.name = "log10p")
  long[, pos_mb := (start + end) / 2e6]
  long[, is_max := log10p == max(log10p, na.rm = TRUE), by = start]

  ggplot(long, aes(x = pos_mb, y = log10p)) +
    geom_point(data = long[is_max == FALSE], size = 0.35, alpha = 0.20, colour = chrom_col) +
    geom_point(data = long[is_max == TRUE],  size = 1.0,  alpha = 0.90, colour = col_max) +
    geom_text_repel(data = long[is_max == TRUE],
                    aes(label = trait), size = 1.8, colour = col_max,
                    min.segment.length = 0, segment.size = 0.2, segment.alpha = 0.5,
                    box.padding = 0.2, max.overlaps = Inf, seed = 42) +
    theme_minimal(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey90", linewidth = 0.3)) +
    labs(x = paste0("Chr", chrom_id, " position (Mb)"),
         y = expression(-log[10](p)), title = title_str)
}

p1 <- make_scatter("9",  17150000L, col_bg9,  "Chr9 ~16.65-17.65 Mb")
p2 <- make_scatter("17", 21400000L, col_bg17, "Chr17 ~20.9-21.9 Mb")

fig <- p1 / p2 & theme(plot.title = element_text(size = 9))

ggsave(out_fig, fig, width = 5.6, height = 5.6)
message("FigureS17 complete: ", out_fig)
