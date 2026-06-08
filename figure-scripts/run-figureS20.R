#!/usr/bin/env Rscript
# Produces: figures/figureS20_gwasremoved_enrichment_dotplot_nearest_gene_plus{flank_kb}kb_top2pct.pdf
# Pre-computed inputs:
#   results/gwas_removed_enrichment/enrichment_results_plus{flank_kb}kb.tsv
# Run: Rscript figure-scripts/run-figureS20.R [flank_kb]
# Shows top 2% band; GO only; intersection_size >= 15; BP/MF and CC as separate panels.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(GOSemSim)
  library(org.Hs.eg.db)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)

args     <- commandArgs(trailingOnly = TRUE)
flank_kb <- if (length(args) >= 1L) as.integer(args[[1]]) else 100L
flank_label <- paste0("plus", flank_kb, "kb")

in_file <- file.path("results", "gwas_removed_enrichment",
                     sprintf("enrichment_results_%s.tsv", flank_label))
if (!file.exists(in_file))
  stop("Missing: ", in_file,
       "\nRun: Rscript figure-scripts/compute-gwas-removed-gene-enrichment.R ", flank_kb)

dir.create("figures", recursive = TRUE, showWarnings = FALSE)

# ---- Load & filter ----
dt_all <- fread(in_file)
dt_all <- dt_all[grepl("^GO:", source) & intersection_size >= 15L & pc != "cross_pc"]
if (nrow(dt_all) == 0L) { message("No qualifying results to plot."); quit(save = "no") }

# Save unfiltered copy for the left (unsimplified) panel
dt_orig <- copy(dt_all)

# ---- Semantic-similarity-based redundancy reduction (GOSemSim) ----
# For each (source, pc, band, dist_mode) group within GO:BP/MF/CC,
# greedily keep the most significant term and discard any term with
# Wang similarity >= 0.7 to the kept term.
simplify_go_dt <- function(dt) {
  source_to_ont <- c("GO:BP" = "BP", "GO:MF" = "MF", "GO:CC" = "CC")
  go_sources    <- intersect(unique(dt$source), names(source_to_ont))
  if (length(go_sources) == 0L) return(dt)

  go_data_cache <- list()

  keep_rows <- list()
  for (src in go_sources) {
    ont <- source_to_ont[[src]]
    # Build/reuse godata object
    if (is.null(go_data_cache[[ont]])) {
      go_data_cache[[ont]] <- tryCatch(
        # GOSemSim >= 2.28 uses 'annoDb'; fall back to 'OrgDb' for older versions
        GOSemSim::godata(annoDb = "org.Hs.eg.db", ont = ont, computeIC = FALSE),
        error = function(e) tryCatch(
          GOSemSim::godata("org.Hs.eg.db", ont = ont, computeIC = FALSE),
          error = function(e2) NULL
        )
      )
    }
    gd <- go_data_cache[[ont]]
    if (is.null(gd)) {
      # Fallback: keep all rows for this source
      keep_rows <- c(keep_rows, list(dt[source == src]))
      next
    }

    sub <- dt[source == src]
    groups <- unique(sub[, .(pc, band, dist_mode)])
    for (g in seq_len(nrow(groups))) {
      grp <- groups[g]
      grp_sub <- sub[pc == grp$pc & band == grp$band & dist_mode == grp$dist_mode]
      if (nrow(grp_sub) <= 1L) {
        keep_rows <- c(keep_rows, list(grp_sub))
        next
      }
      setorder(grp_sub, p.adjust)
      ids     <- grp_sub$ID
      kept    <- logical(length(ids))
      removed <- logical(length(ids))
      for (i in seq_along(ids)) {
        if (removed[i]) next
        kept[i] <- TRUE
        if (i < length(ids)) {
          sims <- tryCatch(
            GOSemSim::mgoSim(ids[i], ids[(i + 1):length(ids)],
                             semData = gd, measure = "Wang", combine = NULL),
            error = function(e) rep(NA_real_, length(ids) - i)
          )
          if (!is.null(sims)) {
            redundant <- which(!is.na(sims) & sims >= 0.7)
            if (length(redundant) > 0L)
              removed[(i + 1):length(ids)][redundant] <- TRUE
          }
        }
      }
      keep_rows <- c(keep_rows, list(grp_sub[kept]))
    }
  }

  # Re-attach non-GO rows (should be none here, but be safe)
  non_go <- dt[!source %in% go_sources]
  rbindlist(c(keep_rows, list(non_go)), use.names = TRUE, fill = TRUE)
}

dt_simp <- tryCatch({
  message("Applying GOSemSim redundancy reduction (cutoff = 0.7)...")
  before_n <- nrow(dt_orig)
  result   <- simplify_go_dt(dt_orig)
  message(sprintf("  Reduced from %d to %d terms.", before_n, nrow(result)))
  result
}, error = function(e) {
  warning("GOSemSim simplification failed, skipping: ", conditionMessage(e))
  copy(dt_orig)
})

# ---- Term abbreviation ----
abbrev_go <- function(x) {
  x <- gsub("plasma membrane bounded cell projection", "PM-bounded cell proj.", x)
  x <- gsub("plasma membrane bounded cell projection", "PM-bounded cell proj.", x)
  x <- gsub("cell-cell adhesion via plasma-membrane adhesion molecules",
             "cell-cell adhesion (PM)", x)
  x <- gsub("homophilic cell adhesion via plasma membrane adhesion molecules",
             "homophilic adhesion (PM)", x)
  x <- gsub("plasma membrane region", "plasma memb. region", x)
  x <- gsub("plasma membrane", "PM", x)
  x <- gsub("photoreceptor cell cilium", "photoreceptor cilium", x)
  x <- gsub("postsynaptic specialization membrane", "postsynaptic spec. memb.", x)
  x <- gsub("postsynaptic specialization", "postsynaptic spec.", x)
  x <- gsub("postsynaptic density membrane", "postsynaptic dens. memb.", x)
  x <- gsub("positive regulation of", "+reg.", x)
  x <- gsub("regulation of synapse structure or activity", "reg. synapse struct./act.", x)
  x <- gsub("regulation of", "reg.", x)
  x <- gsub("neuron projection", "neuron proj.", x)
  x <- gsub("cell projection", "cell proj.", x)
  x <- gsub("morphogenesis involved in neuron diff\\.", "morphogen./neuron diff.", x)
  x <- gsub("morphogenesis", "morphogen.", x)
  x <- gsub("differentiation", "diff.", x)
  x <- gsub("development", "dev.", x)
  x <- gsub("organization", "org.", x)
  x <- gsub("signaling pathway", "signaling", x)
  x <- gsub("nervous system", "NS", x)
  x <- gsub("generation of neurons", "neuron generation", x)
  x <- gsub(" +", " ", trimws(x))
}

# ---- Shared preparation: GO:MF → BP panel, factor PC, -log10(p), abbreviate ----
pc_levels <- paste0("PC", 1:14)
pc_labels <- setNames(pc_levels, pc_levels)
prep_dt <- function(dt) {
  dt[source == "GO:MF", source := "GO:BP"]
  dt[, panel := source]
  dt[, pc    := factor(pc, levels = pc_levels, labels = pc_labels[pc_levels])]
  dt[, log10p    := -log10(p.adjust)]
  dt[, term_short := abbrev_go(term_name)]
  dt
}
dt_orig <- prep_dt(dt_orig)
dt_simp <- prep_dt(dt_simp)

# ---- Shared colour scale: use full (unsimplified) range so both panels are comparable ----
log10p_max <- ceiling(max(dt_orig$log10p, na.rm = TRUE))
log10p_min <- 1
ylgnbu     <- brewer.pal(9, "YlGnBu")

# ---- Build one dotplot panel ----
make_panel <- function(sub, top_n = 20L, panel_label = "", col_title = NULL) {
  if (is.null(sub) || nrow(sub) == 0L) return(NULL)

  # Top top_n terms per PC by p-value, take union
  top_terms <- sub[, .SD[order(p_value)][seq_len(min(.N, top_n))],
                   by = pc][, unique(term_short)]
  sub <- sub[term_short %in% top_terms]

  # Row order: by min p-value across PCs
  term_order <- sub[, .(min_p = min(p_value)), by = term_short]
  setorder(term_order, min_p)
  sub[, term_short := factor(term_short, levels = rev(term_order$term_short))]

  ggplot(sub, aes(x = pc, y = term_short, size = intersection_size, fill = log10p)) +
    geom_point(shape = 21, colour = "black", stroke = 0.4, alpha = 0.9) +
    scale_size_continuous(range = c(2, 8), name = "Genes",
                          breaks = c(15, 25, 35, 45)) +
    scale_fill_gradientn(
      colours = ylgnbu,
      limits  = c(log10p_min, log10p_max),
      name    = expression(-log[10](p[adj]))
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x      = element_text(size = 9),
      axis.text.y      = element_text(size = 8),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.4),
      legend.position  = "right",
      plot.tag         = element_text(size = 10, face = "bold"),
      plot.title       = element_text(size = 9, face = "bold", colour = "grey30")
    ) +
    labs(x = NULL, y = NULL, tag = panel_label, title = col_title)
}

# Helper: number of displayed terms in a panel (0 if panel is NULL)
panel_n <- function(p) {
  if (is.null(p) || is.null(p$data) || nrow(p$data) == 0L) return(0L)
  length(levels(p$data$term_short))
}

# ---- Dot plots: left = unsimplified, right = simplified ----
bands_show <- c("Top 2%")
band_slug  <- c("Top 1%" = "top1pct", "Top 2%" = "top2pct")

for (band_show in bands_show) {
  dt_s <- dt_simp[band == band_show]
  if (nrow(dt_s) == 0L) next

  for (dm in unique(dt_s$dist_mode)) {

    p_bp <- make_panel(dt_s[dist_mode == dm & panel == "GO:BP"], panel_label = "A")
    p_cc <- make_panel(dt_s[dist_mode == dm & panel == "GO:CC"], panel_label = "B")

    if (is.null(p_bp) && is.null(p_cc)) next

    n_bp <- panel_n(p_bp)
    n_cc <- panel_n(p_cc)

    if (!is.null(p_bp) && !is.null(p_cc)) {
      fig <- p_bp / p_cc +
        plot_layout(heights = c(n_bp, n_cc), guides = "collect") &
        theme(legend.position = "right")
    } else {
      fig <- if (!is.null(p_bp)) p_bp else p_cc
    }

    n_pcs  <- length(unique(as.character(dt_s[dist_mode == dm, pc])))
    height <- max(5,  0.28 * (n_bp + n_cc) + 2)
    width  <- max(4,  0.55 * n_pcs + 2.5)

    fname <- file.path("figures",
      sprintf("figureS20_gwasremoved_enrichment_dotplot_%s_%s_%s.pdf",
              dm, flank_label, band_slug[[band_show]]))
    ggsave(fname, fig, width = width, height = height)
    message("Wrote: ", fname)
  }
}

message("run-figureS-gwasremoved-enrichment complete.")
