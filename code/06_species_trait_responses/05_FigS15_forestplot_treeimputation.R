# =============================================================================
# 05_FigS15_forestplot_treeimputation.R
# Fig. S15 (trees): standardized trait effects on species-level abundance
# responses to EH, raw vs phylogenetically imputed traits, across five scales
# (one panel per scale; filled points P < 0.05, hollow non-significant).
#
# Source   : EmbeddingBio/analysis/sfigure15_forestplot_treeimputation.R
#            (paths only changed)
# Inputs   : data/clean/tree_data_{1,5,10,20,50}km.csv          raw traits
#            data/clean/tree_data_{1,5,10,20,50}km_imputed.csv  imputed traits
#              (output of tree_beta_trait() in functions/fun_beta_trait.R)
#            data/clean/tree.svg
# Outputs  : output/figures/FigureS15_forestplot_treeimputation.png
#
# =============================================================================

source(here::here("code", "00_config.R"))
# ── Libraries ────────────────────────────────────────────────────────────────
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(magick)
library(grid)

# ── 0. Parameters ────────────────────────────────────────────────────────────
data_dir <- DIR_CLEAN   # tree_data_<scale>(_imputed).csv and tree.svg
out_png  <- p_fig("FigureS15_forestplot_treeimputation.png")

p_threshold  <- 0.05
scale_levels <- c("1km", "5km", "10km", "20km", "50km")

# Raw vs imputed file names, keyed by scale
raw_files <- c(
  "1km"  = "tree_data_1km.csv",
  "5km"  = "tree_data_5km.csv",
  "10km" = "tree_data_10km.csv",
  "20km" = "tree_data_20km.csv",
  "50km" = "tree_data_50km.csv"
)
imp_files <- c(
  "1km"  = "tree_data_1km_imputed.csv",
  "5km"  = "tree_data_5km_imputed.csv",
  "10km" = "tree_data_10km_imputed.csv",
  "20km" = "tree_data_20km_imputed.csv",
  "50km" = "tree_data_50km_imputed.csv"
)

# Two soft, contrasting colours for the comparison
type_levels <- c("Raw", "Imputed")
type_colors <- c("Raw" = "#6FA8C7", "Imputed" = "#E0884D")  # soft blue / soft orange

# Tree icon (colorized as in the manuscript figure)
tree_icon_path  <- file.path(data_dir, "tree.svg")
tree_icon_color <- "#62b187"
icon_render_px  <- 220L

# ── 1. Trait labels ──────────────────────────────────────────────────────────
label_map <- c(
  "swd"                 = "Wood density",
  "maxHt"               = "Maximum height",
  "leafnitro"           = "Leaf nitrogen",
  "SLA"                 = "Specific leaf area",
  "Seed.dry.mass"       = "Seed dry mass",
  "PFT_binaryEvergreen" = "Evergreen (vs deciduous)",
  "range_km2"           = "Range size (km²)",
  "vol_total"           = "Niche vol. total",
  "vol_core"            = "Niche vol. core",
  "vol_periphery"       = "Niche vol. periphery"
)

make_label <- function(trait, level) {
  if (level %in% names(label_map)) return(unname(label_map[[level]]))
  if (trait %in% names(label_map)) return(unname(label_map[[trait]]))
  str_replace_all(level, "[._]", " ")
}

# ── 2. Read + combine raw and imputed across scales ──────────────────────────
read_set <- function(files, type_lab) {
  map_dfr(names(files), function(sc) {
    fp <- file.path(data_dir, files[[sc]])
    if (!file.exists(fp)) { warning("Missing file: ", fp); return(NULL) }
    df <- read.csv(fp, stringsAsFactors = FALSE)
    df$scale <- sc
    df$type  <- type_lab
    df
  })
}

dat <- bind_rows(read_set(raw_files, "Raw"),
                 read_set(imp_files, "Imputed")) %>%
  mutate(
    scale = factor(scale, levels = scale_levels),
    type  = factor(type,  levels = type_levels),
    display_label = pmap_chr(list(trait, level), make_label),
    sig     = !is.na(p) & p < p_threshold,
    fillcol = ifelse(sig, unname(type_colors[as.character(type)]), "white")
  )

message(sprintf("Rows: %d (Raw + Imputed, %d scales, %d traits)",
                nrow(dat), dplyr::n_distinct(dat$scale),
                dplyr::n_distinct(dat$display_label)))

# ── 3. Order traits on the y-axis (by raw median estimate) ───────────────────
order_tbl <- dat %>%
  filter(type == "Raw") %>%
  group_by(display_label) %>%
  summarise(med_est = median(estimate), .groups = "drop") %>%
  arrange(med_est)

dat <- dat %>%
  mutate(display_label = factor(display_label, levels = order_tbl$display_label))

# ── 4. Colorized tree icon grob ──────────────────────────────────────────────
colorize_svg <- function(svg_path, color, size_px) {
  img  <- image_read_svg(svg_path, height = size_px)
  info <- image_info(img)
  mask <- img %>% image_convert(colorspace = "Gray") %>% image_negate()
  canvas <- image_blank(width = info$width, height = info$height, color = color)
  image_composite(canvas, mask, operator = "CopyOpacity")
}
tree_img  <- colorize_svg(tree_icon_path, tree_icon_color, icon_render_px)
tree_info <- image_info(tree_img)
tree_asp  <- tree_info$width / tree_info$height
tree_grob <- rasterGrob(as.raster(tree_img), interpolate = TRUE,
                        width  = unit(1.3, "cm"),
                        height = unit(1.3 / tree_asp, "cm"))

# ── 5. Forest plot ───────────────────────────────────────────────────────────
dodge <- position_dodge(width = 0.6)
x_min <- min(dat$lwr, na.rm = TRUE)
x_max <- max(dat$upr, na.rm = TRUE)
x_pad <- 0.06 * (x_max - x_min)

p <- ggplot(dat,
            aes(x = estimate, y = display_label, color = type, group = type)) +
  geom_vline(xintercept = 0, color = "grey10",
             linewidth = 0.5, linetype = "dashed") +
  geom_errorbar(aes(xmin = lwr, xmax = upr),
                orientation = "y", width = 0, linewidth = 0.6,
                position = dodge) +
  geom_point(aes(fill = fillcol), shape = 21, size = 2.6, stroke = 0.9,
             position = dodge) +
  facet_grid(. ~ scale) +
  scale_color_manual(values = type_colors, name = "Trait data") +
  scale_fill_identity() +
  scale_x_continuous(limits = c(x_min - x_pad, x_max + x_pad),
                     labels = function(x) sprintf("%.2f", x)) +
  labs(x = "Standardized effect on species' abundance response to EH",
       y = NULL,
       caption = "Filled points: p < 0.05  •  Hollow points: non-significant") +
  guides(color = guide_legend(override.aes = list(size = 3.4), nrow = 1)) +
  theme_classic(base_size = 12) +
  theme(
    strip.text       = element_text(face = "bold", size = 11),
    strip.background  = element_rect(fill = "grey94", color = NA),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.4),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    panel.spacing.x   = unit(0.25, "cm"),
    axis.text.y  = element_text(size = 9),
    axis.text.x  = element_text(size = 8),
    axis.title.x = element_text(size = 12, margin = margin(t = 8)),
    legend.position  = "top",
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 10),
    plot.caption     = element_text(size = 9, color = "grey30", hjust = 0),
    plot.margin = margin(10, 16, 10, 8)
  )

# ── 6. Draw, with tree icon in the top-left corner ───────────────────────────
plot_w_in <- 15
plot_h_in <- max(6, dplyr::n_distinct(dat$display_label) * 0.42 + 2.5)

png(out_png, width = plot_w_in, height = plot_h_in, units = "in", res = 300,
    bg = "white")
grid::grid.newpage()
grid::grid.draw(ggplotGrob(p))
# tree icon overlaid at the top-left corner
grid::pushViewport(grid::viewport(
  x = unit(1.5, "cm"), y = unit(1, "npc") - unit(0.9, "cm"),
  width = unit(1.4, "cm"), height = unit(1.4, "cm")
))
grid::grid.draw(tree_grob)
grid::popViewport()
dev.off()

message("Saved: ", out_png, sprintf("  (%.0f x %.1f in)", plot_w_in, plot_h_in))
