# =============================================================================
# Forest plot: Effect on abundance response to EH (standardized coefficient)
# With SVG taxon icons on y-axis labels AND custom top legend
# =============================================================================
#
# Dependencies:
#   ggplot2, dplyr, ggtext, glue, magick, grid, gtable, cowplot, ggimage
# =============================================================================

# ── Load libraries ────────────────────────────────────────────────────────────
library(ggplot2)
library(dplyr)
library(ggtext)
library(glue)
library(magick)      # SVG → raster conversion and colorization
library(grid)        # low-level grob construction for custom legend
library(gtable)      # arranging grobs in a table layout
library(cowplot)     # stacking legend + plot via plot_grid()
library(ggimage)     # geom_image() for y-axis icons
library(rsvg)

# ── 0. Global parameters ─────────────────────────────────────────────────────

## -- Data directory --
data_dir <- "data_clean"

## -- Input CSV file names --
tree_csv   <- file.path(data_dir, "tree_data_10km.csv")
bird_csv   <- file.path(data_dir, "bird_data_10km.csv")
mammal_csv <- file.path(data_dir, "mammal_data_5km.csv")

## -- SVG icon file paths --
## Potrace-generated: black (#000000) silhouette on white background, no alpha.
icon_paths <- c(
  "Tree"   = file.path(data_dir, "tree.svg"),
  "Bird"   = file.path(data_dir, "bird.svg"),
  "Mammal" = file.path(data_dir, "mammal.svg")
)

## -- Icon rendering size (px height) for top legend --
icon_legend_px <- 120L

## -- Icon rendering size (px height) for y-axis icons --
## Higher resolution for crisp rendering at larger display sizes
icon_yaxis_px <- 200L

## -- Per-taxon display sizes for y-axis icons --
## geom_image() scales each icon within a box of this fraction of the panel.
## Different SVG aspect ratios require different sizes to look balanced:
##   Tree   (1536×1523, ~square)   → appears small at default; needs boost
##   Bird   (1536×1030, landscape) → renders wider; medium size works
##   Mammal (773×1536,  portrait)  → renders taller; medium size works
icon_yaxis_sizes <- c(
  "Tree"   = 0.13,
  "Bird"   = 0.12,
  "Mammal" = 0.06
)

## -- Aspect ratio for geom_image() --
## Controls width/height ratio of the icon bounding box.
## 1.6 works well for mixing landscape/portrait icons on the same axis.
icon_asp <- 1.6

## -- X position for y-axis icons (data coordinates) --
## Placed left of the data range; clip = "off" allows margin rendering.
icon_x_position <- -0.48

## -- Custom legend layout (cm) --
icon_w_cm   <- 0.55    # width of icon cell in top legend
icon_h_cm   <- 0.55    # height of legend row
text_size   <- 12      # font size (pt) for taxon labels in legend
gap_cm      <- 0.8     # gap between taxon groups in legend

## -- Plot output dimensions (inches) --
plot_w_in <- 10
plot_h_in <- 12

## -- Significance threshold --
p_threshold <- 0.05

## -- X-axis breaks --
x_breaks <- seq(-0.4, 0.7, by = 0.1)

# ── 1. Define mappings ────────────────────────────────────────────────────────

taxon_colors <- c("Tree"   = "#62b187",
                  "Bird"   = "#65afc1",
                  "Mammal" = "#e8964c")

category_map <- c(
  # Trees
  "SLA"                  = "Acquisitive traits",
  "leafnitro"            = "Acquisitive traits",
  "maxHt"                = "Conservative traits",
  "Seed.dry.mass"        = "Conservative traits",
  "PFT_binary"           = "Conservative traits",
  # Birds
  "Kipps.Distance"       = "Dispersal",
  "Hand.Wing.Index"      = "Dispersal",
  "range_km2"            = "Range / niche breadth",
  "vol_total"            = "Range / niche breadth",
  "vol_core"             = "Range / niche breadth",
  "vol_periphery"        = "Range / niche breadth",
  "Habitat.Density"      = "Habitat type",
  "Trophic.Niche"        = "Diet",
  "Fldg2"                = "Life history",
  "PCo1"                 = "Habitat type",
  "PCo1m"                = "Life history",
  "PCo3"                 = "Life history",
  "Primary.Lifestyle"    = "Lifestyle",
  # Mammals
  "ForStrat.Value"       = "Habitat type",
  "Terrestriality"       = "Habitat type",
  "HabitatBreadth"       = "Range / niche breadth",
  "InterbirthInterval_d" = "Life history",
  "Precip_Mean_mm"       = "Climatic niche",
  "Temp_Mean_01degC"     = "Climatic niche",
  "AET_Mean_mm"          = "Climatic niche",
  "PET_Mean_mm"          = "Climatic niche"
)

label_map <- c(
  # Trees
  "SLA"                                  = "Specific leaf area",
  "leafnitro"                            = "Leaf nitrogen",
  "maxHt"                                = "Maximum height",
  "Seed.dry.mass"                        = "Seed dry mass",
  "PFT_binaryEvergreen"                  = "Evergreen (vs deciduous)",
  
  # Birds
  "Kipps.Distance"                       = "Kipp\u2019s distance",
  "Hand.Wing.Index"                      = "Hand-wing index",
  "range_km2"                            = "Range size (km\u00b2)",
  "vol_total"                            = "Niche vol. total",
  "vol_core"                             = "Niche vol. core",
  "vol_periphery"                        = "Niche vol. periphery",
  "Habitat.DensityOpen"                  = "Open habitat (vs dense)",
  "Habitat.DensitySemi-open"             = "Semi-open habitat (vs dense)",
  "Trophic.NicheGranivore"               = "Granivore (vs predator)",
  "Fldg2"                                = "Fledging period",
  "PCo1"                                 = "PC1: small body / dense habitat",
  "PCo3"                                 = "PC3: fast life history",
  "Primary.LifestyleInsessorial"         = "Insessorial (vs aerial)",
  
  # Mammals
  "PCo1m"                                = "PC1: small body / fast life history",
  "ForStrat.ValueScansorial"             = "Scansorial (vs ground)",
  "TerrestrialityAbove-ground dwelling"  = "Above-ground (vs ground)",
  "HabitatBreadth2"                      = "Habitat breadth",
  "InterbirthInterval_d"                 = "Interbirth interval",
  "Precip_Mean_mm"                       = "Mean precipitation",
  "Temp_Mean_01degC"                     = "Mean temperature",
  "AET_Mean_mm"                          = "Actual evapotranspiration",
  "PET_Mean_mm"                          = "Potential evapotranspiration"
)

category_levels <- c(
  "Acquisitive traits",
  "Conservative traits",
  "Life history",
  "Range / niche breadth",
  "Dispersal",
  "Habitat type",
  "Diet",
  "Lifestyle",
  "Climatic niche"
)

cat_bg <- c(
  "Acquisitive traits"     = "white",
  "Conservative traits"    = "white",
  "Life history"           = "white",
  "Range / niche breadth"  = "white",
  "Dispersal"              = "white",
  "Habitat type"           = "white",
  "Diet"                   = "white",
  "Lifestyle"              = "white",
  "Climatic niche"         = "white"
)

# ── 2. Read data ──────────────────────────────────────────────────────────────
tree_data   <- read.csv(tree_csv,   stringsAsFactors = FALSE)
bird_data   <- read.csv(bird_csv,   stringsAsFactors = FALSE)
mammal_data <- read.csv(mammal_csv, stringsAsFactors = FALSE)

# ── 3. Colorize SVG helper ───────────────────────────────────────────────────
# Aim:   tint the black silhouette to `color` with transparent background.
#
# Mechanism (luminance-based masking for potrace SVGs):
#   Potrace SVGs = black paths on opaque white (no alpha channel).
#   1. Rasterize → black silhouette on white
#   2. Convert to grayscale: silhouette=0 (black), bg=255 (white)
#   3. Negate (invert): silhouette=255 (white), bg=0 (black)
#   4. Solid-color canvas with taxon color
#   5. CopyOpacity: mask luminance → alpha
#      white(255)=opaque → colored silhouette; black(0)=transparent → bg

colorize_svg <- function(svg_path, color, size_px) {
  img  <- magick::image_read_svg(svg_path, height = size_px)
  info <- magick::image_info(img)
  mask <- img %>%
    magick::image_convert(colorspace = "Gray") %>%
    magick::image_negate()
  canvas <- magick::image_blank(
    width  = info$width,
    height = info$height,
    color  = color
  )
  magick::image_composite(canvas, mask, operator = "CopyOpacity")
}

# ── 4. Save colorized icons as temp PNG files ─────────────────────────────────
# Aim:   ggimage::geom_image() needs file paths (not grobs or base64).
#        We save each colorized icon to a temporary PNG file.
#        These persist for the R session.

icon_tmpfiles <- vapply(
  names(taxon_colors),
  function(tx) {
    img <- colorize_svg(icon_paths[tx], taxon_colors[tx], icon_yaxis_px)
    tmp <- tempfile(fileext = paste0("_", tolower(tx), ".png"))
    magick::image_write(img, path = tmp, format = "png")
    tmp
  },
  character(1)
)

# ── 5. Process taxon function ─────────────────────────────────────────────────
process_taxon <- function(df, taxon_name, p_thr = p_threshold) {
  df %>%
    dplyr::filter(p < p_thr) %>%
    dplyr::mutate(
      taxon    = taxon_name,
      category = dplyr::case_when(
        level %in% names(category_map) ~ category_map[level],
        trait %in% names(category_map) ~ category_map[trait],
        TRUE ~ "Other"
      ),
      display_label = dplyr::case_when(
        level %in% names(label_map) ~ label_map[level],
        trait %in% names(label_map) ~ label_map[trait],
        TRUE ~ level
      )
    ) %>%
    dplyr::select(taxon, category, display_label, estimate, lwr, upr)
}

# ── 6. Combine all taxa ──────────────────────────────────────────────────────
trait_df <- dplyr::bind_rows(
  process_taxon(tree_data,   "Tree"),
  process_taxon(bird_data,   "Bird"),
  process_taxon(mammal_data, "Mammal")
) %>%
  dplyr::mutate(
    taxon    = factor(taxon,    levels = names(taxon_colors)),
    category = factor(category, levels = category_levels),
    icon_path = icon_tmpfiles[as.character(taxon)]
  ) %>%
  # Count rows per panel to compensate icon size for small panels
  dplyr::group_by(category) %>%
  dplyr::mutate(n_rows = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    # Base size from per-taxon lookup
    icon_size_base = icon_yaxis_sizes[as.character(taxon)],
    # Scale factor: panels with fewer rows get proportionally larger icons.
    # Reference panel size = 5 rows (typical). A 2-row panel gets 5/2 = 2.5× boost.
    # Cap the multiplier to avoid oversized icons in single-row panels.
    panel_scale = pmin(5 / n_rows, 4),
    icon_size   = icon_size_base * panel_scale
  )

# ── 7. Build y-axis labels ───────────────────────────────────────────────────
trait_df <- trait_df %>%
  dplyr::mutate(
    y_label_rich = display_label
  )

# Disambiguate duplicates with invisible zero-width spaces
dup_idx <- duplicated(trait_df$y_label_rich) | 
  duplicated(trait_df$y_label_rich, fromLast = TRUE)
if (any(dup_idx)) {
  # For each duplicate group, append increasing numbers of U+200B
  trait_df <- trait_df %>%
    dplyr::group_by(y_label_rich) %>%
    dplyr::mutate(
      y_label_rich = paste0(y_label_rich, strrep("\u200B", dplyr::row_number() - 1))
    ) %>%
    dplyr::ungroup()
}

# ── 8. Pre-compute y-axis order ──────────────────────────────────────────────
level_order <- trait_df %>%
  dplyr::arrange(category, taxon, estimate) %>%
  dplyr::pull(y_label_rich) %>%
  rev()

trait_df <- trait_df %>%
  dplyr::mutate(
    y_label_rich = factor(y_label_rich, levels = level_order)
  )

trait_df <- trait_df %>%
  dplyr::mutate(
    y_label_rich = factor(y_label_rich, levels = level_order)
  )

# ── 9. Category label data ───────────────────────────────────────────────────
cat_label_df <- trait_df %>%
  dplyr::group_by(category) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::select(category, y_label_rich)

# ── 10. Prepare colorized icon grobs for top legend ──────────────────────────
make_icon_grob <- function(svg_path, color, size_px = icon_legend_px) {
  img <- colorize_svg(svg_path, color, size_px)
  grid::rasterGrob(as.raster(img), interpolate = TRUE)
}

icon_grobs <- mapply(
  make_icon_grob,
  svg_path = icon_paths[names(taxon_colors)],
  color    = taxon_colors,
  SIMPLIFY = FALSE
)

# ── 11. Build custom icon legend grob ─────────────────────────────────────────
build_icon_legend <- function(icon_grobs, taxon_colors,
                              icon_w, icon_h, text_size, gap_cm) {
  taxa <- names(icon_grobs)
  n    <- length(taxa)
  grob_list   <- list()
  width_specs <- list()
  for (i in seq_along(taxa)) {
    tx <- taxa[i]
    grob_list   <- c(grob_list, list(icon_grobs[[tx]]))
    width_specs <- c(width_specs, list(unit(icon_w, "cm")))
    tg <- grid::textGrob(
      tx,
      x    = unit(0.15, "cm"),
      just = "left",
      gp   = grid::gpar(fontsize = text_size, col = taxon_colors[tx])
    )
    grob_list   <- c(grob_list, list(tg))
    tw <- unit(nchar(tx) * 0.25 + 0.4, "cm")
    width_specs <- c(width_specs, list(tw))
    if (i < n) {
      grob_list   <- c(grob_list, list(grid::nullGrob()))
      width_specs <- c(width_specs, list(unit(gap_cm, "cm")))
    }
  }
  widths <- do.call(unit.c, width_specs)
  gt <- gtable::gtable(widths = widths, heights = unit(icon_h, "cm"))
  for (j in seq_along(grob_list)) {
    gt <- gtable::gtable_add_grob(gt, grob_list[[j]], t = 1, l = j)
  }
  gt
}

legend_grob <- build_icon_legend(
  icon_grobs   = icon_grobs,
  taxon_colors = taxon_colors,
  icon_w       = icon_w_cm,
  icon_h       = icon_h_cm,
  text_size    = text_size,
  gap_cm       = gap_cm
)

# ── 12. Build the main forest plot ────────────────────────────────────────────
# coord_cartesian(clip = "off") allows geom_image icons placed at
# x = icon_x_position (left of axis) to render into the margin.

plot_main <- ggplot2::ggplot(
  trait_df,
  ggplot2::aes(x = estimate, y = y_label_rich, color = taxon)
) +
  # Panel background
  ggplot2::geom_rect(
    ggplot2::aes(fill = category),
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = Inf,
    alpha = 0.08,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  ggplot2::scale_fill_manual(values = cat_bg, guide = "none") +
  # White vertical grid lines
  ggplot2::geom_vline(
    xintercept = x_breaks,
    color = "white", linewidth = 0.5
  ) +
  # Zero reference line — dashed, grey10
  ggplot2::geom_vline(
    xintercept = 0,
    color = "grey10", linewidth = 0.5, linetype = "dashed"
  ) +
  # Confidence interval bars (no end-caps)
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = lwr, xmax = upr),
    linewidth = 0.7,
    height    = 0
  ) +
  # Point estimates — white fill, taxon-colored border
  ggplot2::geom_point(
    shape  = 21,
    size   = 3.5,
    stroke = 1,
    fill   = "white"
  ) +
  # Y-axis taxon icons via ggimage
  # size mapped per-taxon via aes(size = icon_size) + scale_size_identity()
  ggimage::geom_image(
    ggplot2::aes(
      x     = icon_x_position,
      y     = y_label_rich,
      image = icon_path,
      size  = icon_size
    ),
    asp         = icon_asp,
    inherit.aes = FALSE
  ) +
  ggplot2::scale_size_identity() +
  # Category label top-right
  ggplot2::geom_text(
    data        = cat_label_df,
    ggplot2::aes(x = Inf, y = Inf, label = category),
    inherit.aes = FALSE,
    hjust       = 1.05,
    vjust       = 1.3,
    size        = 4,
    color       = "black"
  ) +
  ggplot2::facet_grid(
    category ~ .,
    scales = "free_y",
    space  = "free_y"
  ) +
  ggplot2::scale_color_manual(values = taxon_colors, guide = "none") +
  ggplot2::scale_x_continuous(
    breaks = x_breaks,
    labels = function(x) sprintf("%.1f", x),
    expand = ggplot2::expansion(mult = 0.02)
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(
    x = "Standardized effects on responses of species abundance to EH",
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    strip.text         = ggplot2::element_blank(),
    strip.background   = ggplot2::element_blank(),
    panel.spacing.y    = unit(0, "cm"),
    panel.border       = ggplot2::element_rect(
      color = "black", fill = NA, linewidth = 0.4
    ),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.background   = ggplot2::element_blank(),
    axis.text.y  = ggtext::element_markdown(size = 10, hjust = 1),
    axis.ticks.y = ggplot2::element_blank(),
    axis.line.y  = ggplot2::element_blank(),
    axis.line.x  = ggplot2::element_line(color = "black"),
    axis.text.x  = ggplot2::element_text(size = 10),
    axis.title.x = ggplot2::element_text(
      size = 12, margin = ggplot2::margin(t = 8)
    ),
    legend.position = "none",
    plot.margin = ggplot2::margin(8, 16, 8, 8)
  )

# ── 13. Combine legend + plot ─────────────────────────────────────────────────
final_plot <- cowplot::plot_grid(
  legend_grob,
  plot_main,
  ncol        = 1,
  rel_heights = c(0.04, 1)
)

# ── 14. Display / save ───────────────────────────────────────────────────────
suppressWarnings(print(final_plot))


ggsave("Figure3.png", plot = final_plot, width = 8, height = 6, dpi = 300)

