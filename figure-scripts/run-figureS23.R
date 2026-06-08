#!/usr/bin/env Rscript
# Produces: figures/figureS23_gwas_removed_subspace_scatter.pdf,
#           results/figureS23_gwas_removed_subspace_scatter.tsv
# Pre-computed inputs:
#   results/pca_loadings_50k.tsv,
#   data/gencode.v19.genes.protein_coding.rds,
#   results/gwas_removed_distance_topload_50k_plus100kb/reduced_pca_window_loadings.tsv
# Run: Rscript figure-scripts/run-figureS23.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

pc_cols  <- c("PC1-4" = "#1f78b4", "PC1-8" = "#984EA3")
flank_kb <- 100L
out_tsv  <- file.path("results", "figureS23_gwas_removed_subspace_scatter.tsv")
out_fig  <- file.path("figures", "figureS23_gwas_removed_subspace_scatter.pdf")
genes_rds <- file.path(project_root, "data", "gencode.v19.genes.protein_coding.rds")
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)

# ---- Load pca_loadings ----
pca_loadings_file <- file.path("results", "pca_loadings_50k.tsv")
if (!file.exists(pca_loadings_file)) stop("Missing: ", pca_loadings_file, "\nRun run-figure1.R first.")
pca_loadings <- fread(pca_loadings_file)

# ---- Load gene annotation ----
if (!file.exists(genes_rds)) stop("Missing genes RDS: ", genes_rds)
genes_gr <- readRDS(genes_rds)
if (any(grepl("^chr", seqlevels(genes_gr)))) seqlevels(genes_gr) <- sub("^chr", "", seqlevels(genes_gr))
GenomeInfoDb::seqlevelsStyle(genes_gr) <- "NCBI"

red_file <- file.path("results", paste0("gwas_removed_distance_topload_50k_plus", flank_kb, "kb"),
                      "reduced_pca_window_loadings.tsv")
if (!file.exists(red_file)) stop("Missing reduced loadings: ", red_file)

parse_region_cols <- function(dt) {
  x <- copy(as.data.table(dt))
  if (!all(c("chrom", "start", "end") %in% names(x))) {
    if (!"region_label" %in% names(x)) stop("Need region_label to parse chrom/start/end.")
    x[, c("chrom", "start", "end") := tstrsplit(region_label, "_", fixed = TRUE)]
  }
  x[, `:=`(chrom = as.integer(chrom), start = as.integer(start), end = as.integer(end))]
  x[, window_id := paste(chrom, start, end, sep = "_")]
  x
}

pl_full <- parse_region_cols(pca_loadings)
pl_red  <- parse_region_cols(fread(red_file))

shared_ids <- intersect(pl_full$window_id, pl_red$window_id)
full_sub   <- pl_full[match(shared_ids, window_id)]
red_sub    <- pl_red[match(shared_ids, window_id)]

make_norm_dt <- function(k) {
  pcs_k    <- paste0("PC", seq_len(k))
  full_norm <- sqrt(rowSums(as.matrix(full_sub[, ..pcs_k])^2))
  red_norm  <- sqrt(rowSums(as.matrix(red_sub[, ..pcs_k])^2))
  dt <- data.table(panel = paste0("PC1-", k), window_id = shared_ids,
                   full_norm = full_norm, removed_norm = red_norm)
  dt[, `:=`(full_norm_z    = as.numeric(scale(full_norm)),
             removed_norm_z = as.numeric(scale(removed_norm)))]
  dt[, full_rank    := frank(-full_norm_z,    ties.method = "average")]
  dt[, removed_rank := frank(-removed_norm_z, ties.method = "average")]
  dt[, rank_change  := removed_rank - full_rank]
  dt[, rank_change_abs := abs(rank_change)]
  n_top <- max(1L, ceiling(nrow(dt) * 0.01))
  dt[, full_top1pct    := full_rank    <= n_top]
  dt[, removed_top1pct := removed_rank <= n_top]
  dt[, delta_abs := abs(removed_norm_z - full_norm_z)]
  dt
}

plot_dt <- rbindlist(list(make_norm_dt(4), make_norm_dt(8)), use.names = TRUE)

window_meta <- unique(pl_full[match(shared_ids, window_id), .(window_id, chrom, start, end)])
window_gr   <- GRanges(seqnames = as.character(window_meta$chrom),
                       ranges = IRanges(start = window_meta$start + 1L, end = window_meta$end),
                       window_id = window_meta$window_id)
GenomeInfoDb::seqlevelsStyle(window_gr) <- "NCBI"
hits <- findOverlaps(window_gr, genes_gr, ignore.strand = TRUE)
gene_lab_dt <- data.table(
  window_id  = mcols(window_gr)$window_id[queryHits(hits)],
  gene_label = mcols(genes_gr)$symbol[subjectHits(hits)]
)
gene_lab_dt <- gene_lab_dt[!is.na(gene_label) & gene_label != "",
                            .(gene_label = paste(unique(gene_label), collapse = ", ")), by = window_id]
window_meta <- merge(window_meta, gene_lab_dt, by = "window_id", all.x = TRUE)
window_meta[is.na(gene_label) | gene_label == "",
            gene_label := sprintf("chr%s:%0.3f-%0.3fMb", chrom, start / 1e6, end / 1e6)]

cor_summary <- plot_dt[, .(correlation = cor(full_norm_z, removed_norm_z, use = "pairwise.complete.obs")), by = panel]
rank_cor_summary <- plot_dt[, .(correlation = cor(full_rank, removed_rank, use = "pairwise.complete.obs")), by = panel]
plot_dt <- merge(plot_dt, cor_summary, by = "panel", all.x = TRUE)
setnames(rank_cor_summary, "correlation", "rank_correlation")
plot_dt <- merge(plot_dt, rank_cor_summary, by = "panel", all.x = TRUE)
plot_dt <- merge(plot_dt, window_meta[, .(window_id, gene_label)], by = "window_id", all.x = TRUE)
fwrite(plot_dt, out_tsv, sep = "\t")

label_dev_dt   <- plot_dt[order(-delta_abs), head(.SD, 10), by = panel]
rank_plot_dt   <- copy(plot_dt)
label_rank_dt  <- plot_dt[removed_top1pct == TRUE][order(-rank_change_abs, removed_rank), head(.SD, 10), by = panel]
lims           <- range(c(plot_dt$full_norm_z, plot_dt$removed_norm_z), finite = TRUE)
rank_cap       <- 600L
rank_plot_dt[, `:=`(full_rank_cap = pmin(full_rank, rank_cap), removed_rank_cap = pmin(removed_rank, rank_cap))]
label_rank_dt[, `:=`(full_rank_cap = pmin(full_rank, rank_cap), removed_rank_cap = pmin(removed_rank, rank_cap))]
rank_lims      <- c(rank_cap, 1)

p_dev <- ggplot(plot_dt, aes(x = full_norm_z, y = removed_norm_z)) +
  geom_point(aes(color = panel), size = 0.8, alpha = 0.38, stroke = 0) +
  geom_abline(intercept = 0, slope = 1, colour = "grey65", linewidth = 0.4, linetype = "dashed") +
  geom_text(data = unique(plot_dt[, .(panel, correlation)]),
            aes(x = lims[1] + 0.06 * diff(lims), y = lims[2] - 0.08 * diff(lims),
                label = sprintf("r = %.3f", correlation), color = panel),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.1) +
  geom_text_repel(data = label_dev_dt, aes(label = gene_label), size = 1.25, max.overlaps = Inf,
                  box.padding = 0.25, point.padding = 0.1, segment.size = 0.25, segment.alpha = 0.5,
                  min.segment.length = 0, color = "black", show.legend = FALSE) +
  scale_color_manual(values = pc_cols, drop = FALSE) +
  coord_cartesian(xlim = lims, ylim = lims) +
  facet_wrap(~panel, ncol = 2) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(), legend.position = "none",
        strip.text = element_text(size = 11)) +
  labs(x = "Original subspace loading norm (z-score)", y = "GWAS-removed subspace loading norm (z-score)")

make_rank_panel <- function(panel_name, x_lab = NULL, y_lab = NULL) {
  dt_main  <- rank_plot_dt[panel == panel_name]
  dt_lab   <- label_rank_dt[panel == panel_name]
  dt_inset <- plot_dt[panel == panel_name]
  r_lab    <- unique(dt_main$rank_correlation)
  if (length(r_lab) == 0L || !is.finite(r_lab[1])) r_lab <- NA_real_

  p_main <- ggplot(dt_main, aes(x = full_rank_cap, y = removed_rank_cap)) +
    geom_point(color = pc_cols[[panel_name]], size = 0.8, alpha = 0.38, stroke = 0) +
    geom_abline(intercept = 0, slope = 1, colour = "grey65", linewidth = 0.4, linetype = "dashed") +
    geom_text_repel(data = dt_lab, aes(x = full_rank_cap, y = removed_rank_cap, label = gene_label),
                    size = 1.25, max.overlaps = Inf, box.padding = 0.25, point.padding = 0.1,
                    segment.size = 0.25, segment.alpha = 0.5, min.segment.length = 0, color = "black") +
    scale_x_reverse(limits = rank_lims) + scale_y_reverse(limits = rank_lims) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
          panel.grid.minor = element_blank(), plot.title = element_text(size = 11, hjust = 0.5)) +
    labs(title = panel_name, x = x_lab, y = y_lab)

  p_inset <- ggplot(dt_inset, aes(x = full_rank, y = removed_rank)) +
    geom_point(color = pc_cols[[panel_name]], size = 0.35, alpha = 0.35, stroke = 0) +
    annotate("rect", xmin = 600, xmax = 1, ymin = 600, ymax = 1,
             fill = NA, color = "black", linewidth = 0.25) +
    geom_abline(intercept = 0, slope = 1, colour = "grey70", linewidth = 0.25, linetype = "dashed") +
    annotate("text", x = 50, y = 11000, label = sprintf("r = %.3f", r_lab[1]),
             hjust = 1, vjust = 0, size = 2.2, color = "black") +
    scale_x_reverse(breaks = c(12000, 1), labels = c("12000", "1")) +
    scale_y_reverse(breaks = c(12000, 1), labels = c("12000", "1")) +
    theme_minimal(base_size = 6) +
    theme(panel.grid = element_blank(), axis.title = element_blank(),
          axis.text = element_text(size = 5), axis.ticks = element_line(linewidth = 0.2),
          plot.margin = margin(1, 1, 1, 1),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3))

  p_main + inset_element(p_inset, left = 0.70, bottom = 0.09, right = 0.95, top = 0.35)
}

p_rank4 <- make_rank_panel("PC1-4", x_lab = "Original importance rank", y_lab = "GWAS-removed importance rank")
p_rank8 <- make_rank_panel("PC1-8", x_lab = "Original importance rank", y_lab = NULL)
p <- p_dev / (p_rank4 | p_rank8)

ggsave(out_fig, p, width = 8.4, height = 8.8)
message("Wrote: ", out_tsv)
message("FigureS23 complete: ", out_fig)
