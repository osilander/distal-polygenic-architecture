#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})
setDTthreads(1L)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
in_dir <- file.path(project_root, "data", "windows_filtered")
out_fig <- file.path(project_root, "figures", "figureS1_window_snp_count_histograms_50k.pdf")
out_tsv <- file.path(project_root, "results", "window_snp_count_histograms_50k.tsv")

files <- Sys.glob(file.path(in_dir, "*_w50000.summary.tsv"))
if (length(files) == 0L) {
  stop("No 50 kb filtered window summary files found in ", in_dir)
}

message("Reading ", length(files), " filtered 50 kb trait files...")
all_pairs <- rbindlist(lapply(files, function(f) {
  dt <- fread(f, select = c("chrom", "start", "end", "count_snps"))
  dt[, trait := sub("_w50000.summary.tsv$", "", basename(f))]
  dt[, chrom := as.character(chrom)]
  dt[, count_snps := as.numeric(count_snps)]
  dt[is.finite(count_snps), .(chrom, start, end, trait, count_snps)]
}), use.names = TRUE)

all_pairs <- all_pairs[!(
  (chrom == "17" & start < 45000000 & end > 43500000) |
    (chrom == "6" & start < 34000000 & end > 25000000) |
    (chrom == "12" & start < 113400000 & end > 111200000)
)]

window_summ <- all_pairs[, .(
  median_count_snps = median(count_snps, na.rm = TRUE),
  min_count_snps = min(count_snps, na.rm = TRUE)
), by = .(chrom, start, end)]

plot_dt <- rbindlist(list(
  all_pairs[, .(panel = "All retained trait-window pairs", count_snps)],
  window_summ[, .(panel = "Median across traits per genomic window", count_snps = median_count_snps)],
  window_summ[, .(panel = "Minimum across traits per genomic window", count_snps = min_count_snps)]
), use.names = TRUE)

panel_order <- c(
  "All retained trait-window pairs",
  "Median across traits per genomic window",
  "Minimum across traits per genomic window"
)
plot_dt[, panel := factor(panel, levels = panel_order)]

panel_stats <- plot_dt[, .(
  n = .N,
  median_count = median(count_snps, na.rm = TRUE)
), by = panel]
panel_stats[, label := sprintf("N = %s\nmedian = %s", comma(n), comma(round(median_count)))]
panel_stats[, x := Inf]
panel_stats[, y := Inf]

fwrite(plot_dt, out_tsv, sep = "\t")

p <- ggplot(plot_dt, aes(x = count_snps)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 80,
    fill = "#8da0cb",
    color = "white",
    linewidth = 0.2
  ) +
  geom_vline(
    data = panel_stats,
    aes(xintercept = median_count),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey25"
  ) +
  geom_text(
    data = panel_stats,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1.03,
    vjust = 1.2,
    size = 3.2
  ) +
  facet_wrap(~ panel, ncol = 1, scales = "free_y") +
  scale_x_continuous(labels = comma, limits = c(0, 800), oob = squish) +
  labs(
    x = "SNPs per 50 kb window",
    y = "Density"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.margin = margin(8, 20, 8, 8)
  )

ggsave(out_fig, p, width = 8.4, height = 8.8)
message("Wrote: ", out_fig)
message("Wrote: ", out_tsv)
