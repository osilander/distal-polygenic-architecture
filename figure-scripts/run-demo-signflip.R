#!/usr/bin/env Rscript
# Demonstrates that trait-sign flipping only reflects flipped traits in PC space.
# Produces: figures/demo_signflip_2x2.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(irlba)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(project_root, "figure-scripts", "helpers-50k-matrix.R"))
source(file.path(project_root, "figure-scripts", "helpers-trait-labels.R"))

# ── Category colour palette (identical to run-figure1.R) ─────────────────────
cat_dt <- fread("results/umap-coordinates.tsv",
                select = c("trait_abbrev", "phenotype_group"))
setnames(cat_dt, "trait_abbrev", "trait")

custom_theme_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
  "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
  "#E5C494", "#B3B3B3"
)
all_cats <- sort(unique(na.omit(cat_dt$phenotype_group)))
all_cats <- all_cats[nzchar(all_cats)]
cat_cols <- setNames(custom_theme_colors[seq_along(all_cats)], all_cats)

# ── Load and prepare matrix ───────────────────────────────────────────────────
message("Loading matrix...")
eff   <- fread("results/all-trait-50k-mean-effects.tsv")
pvals <- fread("results/all-trait-50k-mean-pvals.tsv")
prep  <- prepare_X_50k(eff, pvals)
X_wt  <- prep$X
trait_names <- prep$trait_names
X_tw <- t(X_wt)   # traits × windows

# ── Full SVD to identify extreme traits ──────────────────────────────────────
message("Full SVD (403 traits)...")
pca_full <- prcomp_irlba(X_tw, n = 4, center = FALSE, scale. = FALSE)
sc_full  <- pca_full$x[, 1:2]

orient <- function(sc, names_v, anchors) {
  for (anc in anchors) {
    nm <- anc[[1]]; j <- anc[[2]]
    idx <- match(nm, names_v)
    if (!is.na(idx) && sc[idx, j] < 0) sc[, j] <- -sc[, j]
  }
  sc
}
anchors <- list(list("BMI_03", 1L), list("StandHgt", 2L))
sc_full  <- orient(sc_full, trait_names, anchors)

# ── Select 40 traits spread in PC1/PC2 ────────────────────────────────────────
dist2   <- sqrt(sc_full[, 1]^2 + sc_full[, 2]^2)
top40   <- order(dist2, decreasing = TRUE)[1:40]
X_sub   <- X_tw[top40, ]
nm_sub  <- trait_names[top40]

# Display labels: use name_map where available, else raw abbrev; cap at 16 chars
lab_sub <- ifelse(nm_sub %in% names(name_map), name_map[nm_sub], nm_sub)
lab_sub <- substr(lab_sub, 1, 16)

# Category for each of the 40 traits
cat_sub <- cat_dt$phenotype_group[match(nm_sub, cat_dt$trait)]
cat_sub[is.na(cat_sub)] <- "Other"

# ── Base 40-trait SVD ─────────────────────────────────────────────────────────
message("Base 40-trait SVD...")
pca_base <- prcomp_irlba(X_sub, n = 4, center = FALSE, scale. = FALSE)
sc_base  <- pca_base$x[, 1:2]
sc_base  <- orient(sc_base, nm_sub, anchors)

# ── Choose traits to flip: 10 most extreme non-anchor traits ─────────────────
can_flip   <- which(!(nm_sub %in% c("BMI_03", "StandHgt")))
flip_order <- can_flip[order(sqrt(sc_base[can_flip, 1]^2 +
                                  sc_base[can_flip, 2]^2), decreasing = TRUE)]
flip10 <- flip_order[1:10]
flip5  <- flip10[1:5]
flip1  <- flip10[1]

# ── Panel builder ─────────────────────────────────────────────────────────────
make_panel <- function(X, flip_rows = integer(0), sc_ref,
                       nm_v, lab_v, cat_v, title) {
  X_mod <- X
  if (length(flip_rows)) X_mod[flip_rows, ] <- -X_mod[flip_rows, ]

  pca <- prcomp_irlba(X_mod, n = 4, center = FALSE, scale. = FALSE)
  sc  <- pca$x[, 1:2]

  # Align sign to base using unflipped traits
  unfl <- setdiff(seq_len(nrow(X)), flip_rows)
  for (j in 1:2)
    if (cor(sc[unfl, j], sc_ref[unfl, j]) < 0) sc[, j] <- -sc[, j]

  df <- data.frame(
    PC1     = sc[, 1],
    PC2     = sc[, 2],
    label   = lab_v,
    cat     = cat_v,
    flipped = seq_len(nrow(X)) %in% flip_rows,
    stringsAsFactors = FALSE
  )

  # Build in two layers: unflipped first, then flipped on top
  p <- ggplot(df, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey80") +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey80") +
    # unflipped: filled circle, category colour
    geom_point(data = subset(df, !flipped),
               aes(fill = cat), shape = 21, colour = "white",
               size = 2.2, stroke = 0.3) +
    # flipped: category fill + bold black outline ring
    geom_point(data = subset(df, flipped),
               aes(fill = cat), shape = 21, colour = "black",
               size = 3.0, stroke = 1.2) +
    # labels via repel — unflipped subtle, flipped bold
    geom_text_repel(data = subset(df, !flipped),
                    aes(label = label, colour = cat),
                    size = 1.9, max.overlaps = 30,
                    segment.size = 0.2, segment.colour = "grey70",
                    box.padding = 0.2, point.padding = 0.1,
                    min.segment.length = 0.3, seed = 42,
                    show.legend = FALSE) +
    geom_text_repel(data = subset(df, flipped),
                    aes(label = label, colour = cat),
                    size = 2.1, fontface = "bold", max.overlaps = 50,
                    segment.size = 0.3, segment.colour = "grey40",
                    box.padding = 0.3, point.padding = 0.15,
                    min.segment.length = 0.2, seed = 42,
                    show.legend = FALSE) +
    scale_fill_manual(values = cat_cols, name = "Category",
                      drop = TRUE, na.value = "grey60") +
    scale_colour_manual(values = cat_cols, guide = "none",
                        na.value = "grey40") +
    scale_x_continuous(expand = expansion(mult = 0.22)) +
    scale_y_continuous(expand = expansion(mult = 0.22)) +
    labs(title = title, x = "PC1 trait score", y = "PC2 trait score") +
    theme_bw(base_size = 9) +
    theme(panel.grid       = element_blank(),
          plot.title       = element_text(face = "bold", size = 9),
          axis.text        = element_text(size = 7),
          legend.position  = "none")

  p
}

# ── Four panels ───────────────────────────────────────────────────────────────
message("Building panels...")

pA <- make_panel(X_sub, integer(0), sc_base, nm_sub, lab_sub, cat_sub,
                 "a   Base: 40 traits, no sign flips")
pB <- make_panel(X_sub, flip1,  sc_base, nm_sub, lab_sub, cat_sub,
                 paste0("b   Flip 1 trait (", lab_sub[flip1], ")"))
pC <- make_panel(X_sub, flip5,  sc_base, nm_sub, lab_sub, cat_sub,
                 "c   Flip 5 traits")
pD <- make_panel(X_sub, flip10, sc_base, nm_sub, lab_sub, cat_sub,
                 "d   Flip 10 traits")

# ── Shared legend ─────────────────────────────────────────────────────────────
# Extract from panel A with legend turned on
legend_p <- make_panel(X_sub, integer(0), sc_base, nm_sub, lab_sub, cat_sub, "") +
  theme(legend.position  = "right",
        legend.key.size  = unit(0.35, "cm"),
        legend.text      = element_text(size = 7),
        legend.title     = element_text(size = 8, face = "bold"))
legend_grob <- cowplot::get_legend(legend_p)

# Combine
suppressPackageStartupMessages(library(cowplot))
panels <- (pA | pB) / (pC | pD)

fig <- cowplot::plot_grid(
  panels, legend_grob,
  ncol = 2, rel_widths = c(1, 0.18)
)
fig <- cowplot::ggdraw(fig) +
  cowplot::draw_label(
    "Trait-sign flipping reflects only the flipped traits; unflipped traits are stationary\n40 traits selected by extremity of projection in the full 403-trait SVD  |  outlined point = flipped trait",
    x = 0.5, y = 0.99, vjust = 1, hjust = 0.5,
    size = 8.5, fontface = "bold", colour = "grey20"
  )

out <- file.path("figures", "demo_signflip_2x2.pdf")
ggsave(out, fig, width = 13, height = 10.5)
message("Saved: ", out)
