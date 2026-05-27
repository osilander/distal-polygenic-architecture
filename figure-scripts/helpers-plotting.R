regions_to_central_gene_local <- function(regions_dt, genes_gr) {
  if (any(grepl("^chr", seqlevels(genes_gr)))) {
    seqlevels(genes_gr) <- sub("^chr", "", seqlevels(genes_gr))
  }
  if ("chrom" %in% names(regions_dt)) {
    regions_dt[, chrom := sub("^chr", "", chrom)]
  }
  GenomeInfoDb::seqlevelsStyle(genes_gr) <- "NCBI"

  regions_gr <- GRanges(
    seqnames = regions_dt$chrom,
    ranges = IRanges(start = regions_dt$start, end = regions_dt$end)
  )
  seqlevelsStyle(regions_gr) <- "NCBI"

  overlaps <- findOverlaps(regions_gr, genes_gr, ignore.strand = TRUE)
  dt <- data.table(
    Region_ID        = queryHits(overlaps),
    region_chr       = as.character(seqnames(regions_gr)[queryHits(overlaps)]),
    region_start     = start(regions_gr)[queryHits(overlaps)],
    region_end       = end(regions_gr)[queryHits(overlaps)],
    gene_chr         = as.character(seqnames(genes_gr)[subjectHits(overlaps)]),
    gene_start       = start(genes_gr)[subjectHits(overlaps)],
    gene_end         = end(genes_gr)[subjectHits(overlaps)],
    strand           = as.character(strand(genes_gr)[subjectHits(overlaps)]),
    ensembl_gene_id  = mcols(genes_gr)$ensembl_gene_id[subjectHits(overlaps)],
    symbol           = mcols(genes_gr)$symbol[subjectHits(overlaps)],
    gene_biotype     = mcols(genes_gr)$gene_biotype[subjectHits(overlaps)]
  )

  dt[, biotype_category := ifelse(gene_biotype == "protein_coding", "coding", "noncoding")]
  dt[, region_center := (region_start + region_end) / 2]
  dt[, gene_center   := (gene_start   + gene_end)   / 2]
  dt[, dist_to_center := abs(region_center - gene_center)]
  dt <- dt[order(Region_ID, dist_to_center)]
  central_genes <- dt[, .SD[1], by = Region_ID]

  missing_ids <- setdiff(seq_len(nrow(regions_dt)), unique(central_genes$Region_ID))
  if (length(missing_ids) > 0) {
    intergenic_rows <- regions_dt[missing_ids, .(
      Region_ID       = .I,
      region_chr      = chrom,
      region_start    = start,
      region_end      = end,
      gene_chr        = chrom,
      gene_start      = NA_integer_,
      gene_end        = NA_integer_,
      strand          = NA_character_,
      ensembl_gene_id = NA_character_,
      symbol          = paste0(chrom, "_", start, "_", end),
      gene_biotype    = "intergenic",
      biotype_category = "noncoding",
      region_center   = (start + end) / 2,
      gene_center     = NA_real_,
      dist_to_center  = NA_real_
    )]
    central_genes <- rbindlist(list(central_genes, intergenic_rows), use.names = TRUE, fill = TRUE)
  }

  central_genes[, location    := symbol]
  central_genes[, coding_region := ifelse(biotype_category == "coding", "gene", "noncoding")]
  central_genes[, .(Region_ID, region_chr, region_start, region_end,
                    gene_chr, gene_start, gene_end, strand,
                    ensembl_gene_id, symbol, gene_biotype,
                    biotype_category, dist_to_center, location, coding_region)]
}

make_gene_track <- function(genes_gr, chr, start_bp, end_bp, pos_start_chr = 0) {
  chr_raw  <- as.character(chr)
  chr_ucsc <- if (!grepl("^chr", chr_raw)) paste0("chr", chr_raw) else chr_raw
  seqlevelsStyle(genes_gr) <- "UCSC"
  g <- as.data.table(genes_gr)
  g[, seqname_chr := as.character(seqnames(genes_gr))]
  genes_chr <- g[seqname_chr %in% c(chr_raw, chr_ucsc)]
  g_overlap <- genes_chr[start <= end_bp & end >= start_bp]
  if (nrow(g_overlap) == 0) return(NULL)
  g_overlap[, gene_start_mb := start / 1e6]
  g_overlap[, gene_end_mb   := end   / 1e6]
  g_overlap[, gene_mid_mb   := (start + end) / 2e6]
  setorder(g_overlap, gene_start_mb)
  g_overlap[, lane := NA_integer_]
  for (i in seq_len(nrow(g_overlap))) {
    this_start <- g_overlap$gene_start_mb[i]
    this_end   <- g_overlap$gene_end_mb[i]
    if (is.na(this_start) || is.na(this_end)) {
      g_overlap$lane[i] <- ifelse(i == 1, 1L, max(g_overlap$lane, na.rm = TRUE) + 1L)
      next
    }
    placed <- FALSE
    for (l in seq_len(max(c(1, g_overlap$lane), na.rm = TRUE))) {
      lane_genes       <- g_overlap[lane == l & !is.na(gene_end_mb)]
      lane_genes_clean <- lane_genes[!is.na(gene_end_mb)]
      last_end <- if (nrow(lane_genes_clean) > 0) {
        setorder(lane_genes_clean, gene_end_mb)
        lane_genes_clean$gene_end_mb[nrow(lane_genes_clean)]
      } else {
        -Inf
      }
      if (!is.na(last_end) && this_start > last_end) {
        g_overlap$lane[i] <- l
        placed <- TRUE
        break
      }
    }
    if (!placed) g_overlap$lane[i] <- max(c(1, g_overlap$lane), na.rm = TRUE) + 1L
  }
  g_overlap[, .(symbol, start, end, gene_start_mb, gene_end_mb, gene_mid_mb, lane)]
}

plot_pc_loadings <- function(pca_loadings, pc = "PC1", pdf_file = NULL,
                             zoom_chr = NULL, zoom_start = NULL, zoom_end = NULL,
                             genes_dt = NULL, top_perc = 0.001, use_abs = FALSE) {
  dt <- copy(pca_loadings)
  if (!all(c("start", "end", "chr") %in% names(dt))) {
    dt[, c("chr", "start", "end") := tstrsplit(region_label, "_")]
    dt[, `:=`(chr = as.integer(chr), start = as.numeric(start), end = as.numeric(end))]
  }
  dt[, mid := (start + end) / 2]
  setorder(dt, chr, mid)
  chr_lengths <- dt[, .(chr_len = max(mid)), by = chr]
  chr_lengths[, chr_start := c(0, head(cumsum(chr_len), -1))]
  dt <- merge(dt, chr_lengths, by = "chr", all.x = TRUE)
  dt[, pos_cum := mid + chr_start]
  axis_df <- chr_lengths[, .(chr, center = chr_start + chr_len / 2)]

  if (!pc %in% names(dt)) stop("Column ", pc, " not found in loadings data table.")
  dt[, loading     := get(pc)]
  dt[, abs_loading := abs(loading)]
  dt[, top         := abs_loading >= quantile(abs_loading, 1 - top_perc)]

  if (!is.null(genes_dt)) {
    dt[, chr := as.character(chr)]
    regions_dt <- dt[top == TRUE, .(chrom = chr, start, end, region_label)]
    annots <- regions_to_central_gene_local(regions_dt, genes_dt)
    dt <- merge(dt, annots[, .(region_chr, region_start, region_end, symbol)],
                by.x = c("chr", "start", "end"),
                by.y = c("region_chr", "region_start", "region_end"), all.x = TRUE)
  }

  is_zoom <- !is.null(zoom_chr)
  if (is_zoom) {
    dt <- dt[chr == zoom_chr & start >= zoom_start & end <= zoom_end]
  }

  if (!"symbol" %in% names(dt)) dt[, symbol := NA_character_]
  dt[, region_label := as.character(region_label)]
  dt[, symbol       := as.character(symbol)]
  parts     <- tstrsplit(dt$region_label, "_")
  chr_raw   <- parts[[1]]
  start_raw <- suppressWarnings(as.numeric(parts[[2]]))
  end_raw   <- suppressWarnings(as.numeric(parts[[3]]))
  dt[, pretty_label := sprintf("Chr%s:%.3f-%.3fMb", chr_raw, start_raw/1e6, end_raw/1e6)]
  dt[is.na(symbol) | symbol == "" | symbol == region_label |
       grepl("^[0-9._eE\\+\\-]+$", symbol), symbol := pretty_label]
  dt[, plot_label := symbol]

  collapse_adjacent_labels <- function(x_top) {
    z <- copy(x_top)
    if (nrow(z) == 0L) return(z)
    setorder(z, chr, start, end)
    z[, prev_end    := shift(end), by = chr]
    z[, is_adj_prev := !is.na(prev_end) & (start == prev_end)]
    z[, run_id      := cumsum(!is_adj_prev), by = chr]
    keep <- z[, .SD[floor((.N + 1L) / 2L)], by = .(chr, run_id)]
    keep[, c("prev_end", "is_adj_prev", "run_id") := NULL]
    keep[]
  }

  gene_track <- NULL
  if (is_zoom && !is.null(genes_dt)) {
    gene_track <- make_gene_track(genes_dt, chr = zoom_chr, start_bp = zoom_start,
                                  end_bp = zoom_end,
                                  pos_start_chr = axis_df$chr_start[axis_df$chr == zoom_chr])
  }

  if (is_zoom) {
    dt[, pos_mb  := (start + end) / 2e6]
    mb_range  <- range(dt$pos_mb, na.rm = TRUE)
    pretty_n  <- 6L
    mb_breaks <- pretty(mb_range, n = pretty_n)
    while (length(mb_breaks) > 5L && pretty_n > 2L) {
      pretty_n  <- pretty_n - 1L
      mb_breaks <- pretty(mb_range, n = pretty_n)
    }
    x_var   <- "pos_mb"
    x_scale <- scale_x_continuous(name = paste0("Chr", zoom_chr, " position (Mb)"),
                                   breaks = mb_breaks, labels = sprintf("%.2f", mb_breaks))
    zoom_lines <- geom_vline(xintercept = mb_breaks, colour = "grey85", linewidth = 0.3)
  } else {
    x_var   <- "pos_cum"
    x_scale <- scale_x_continuous(breaks = axis_df$center, labels = axis_df$chr, name = "Chromosome")
    zoom_lines <- NULL
  }

  dt_labels    <- if (is_zoom) dt[0] else collapse_adjacent_labels(dt[top == TRUE])
  colour_scale <- if (is_zoom) NULL else scale_colour_manual(values = chrom_plot_palette)
  plot_col     <- if (use_abs) "abs_loading" else "loading"

  if (is_zoom) {
    p <- ggplot(dt, aes(x = .data[[x_var]], y = get(plot_col))) +
      zoom_lines +
      geom_point(size = 0.70, alpha = 0.70, colour = "#984EA3") +
      geom_hline(yintercept = 0, colour = "grey60") +
      geom_point(data = dt[top == TRUE], size = 1.1, alpha = 0.95, colour = "#984EA3") +
      x_scale + theme_minimal(base_size = 10) + theme(legend.position = "none") +
      ylab(paste0(pc, if (use_abs) " |loading|" else " loading"))
  } else {
    p <- ggplot(dt, aes(x = .data[[x_var]], y = get(plot_col), colour = as.factor(chr))) +
      zoom_lines +
      geom_point(size = 0.5, alpha = 0.7) +
      geom_hline(yintercept = 0, colour = "grey60") +
      geom_point(data = dt[top == TRUE], size = 1) +
      colour_scale + x_scale + theme_minimal(base_size = 10) +
      theme(legend.position = "none") +
      ylab(paste0(pc, if (use_abs) " |loading|" else " loading"))
  }

  if (nrow(dt_labels) > 0L) {
    p <- p + geom_text_repel(data = dt_labels, aes(label = plot_label),
                              size = 2, min.segment.length = 0, max.overlaps = 50, colour = "black")
  }

  if (!is.null(gene_track)) {
    lane_height <- (max(dt[[plot_col]]) - min(dt[[plot_col]])) * 0.03
    ymin_base   <- min(dt[[plot_col]], na.rm = TRUE) * 0.95
    p <- p +
      geom_rect(data = gene_track, inherit.aes = FALSE,
                aes(xmin = gene_start_mb, xmax = gene_end_mb,
                    ymin = ymin_base + (lane - 1) * lane_height,
                    ymax = ymin_base + lane * lane_height),
                fill = "grey80", alpha = 0.9, colour = "grey55", linewidth = 0.2) +
      geom_text_repel(data = gene_track, inherit.aes = FALSE,
                      aes(x = gene_mid_mb, y = ymin_base + (lane - 1) * lane_height, label = symbol),
                      seed = 1, min.segment.length = 0, max.overlaps = Inf,
                      box.padding = 0.35, point.padding = 0.20,
                      force = 3.0, force_pull = 0.10, segment.alpha = 0.6,
                      segment.size = 0.25, direction = "both", size = 2, colour = "grey25")
  }

  if (!is.null(pdf_file)) { pdf(pdf_file, width = 12, height = 6); print(p); dev.off() }

  table_dt <- dt[top == TRUE, .(chr, start, end, loading, abs_loading, symbol)]
  list(plot = p, table = table_dt[order(-abs_loading)])
}

filter_effect_matrix <- function(eff, pvals, p_thresh = NULL) {
  helper_path <- file.path("post-process-scripts", "helpers-50k-matrix.R")
  helper_env  <- new.env(parent = parent.frame())
  sys.source(helper_path, envir = helper_env)
  helper_env$filter_effect_matrix(eff = eff, pvals = pvals, p_thresh = p_thresh)
}
