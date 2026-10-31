# =============================================================================
# 04_FigS14_forestplot_multiscale.R
# Fig. S14: significant (P < 0.05) standardized trait effects on species-level
# abundance responses to EH across ALL spatial scales, for trees, birds and
# mammals (estimates coloured by scale, one facet per taxon, taxon icons in the
# facet strips; style follows 02_Figure3_forestplot.R).
#
# Source   : EmbeddingBio/analysis/sfigure14_forestplot_multiplescale.R
#            (paths only changed)
# Inputs   : data/clean/tree_data_{1,5,10,20,50}km.csv
#            data/clean/bird_data_{400m,1,5,10,20}km.csv
#            data/clean/mammal_data_{1,5,10,20,50}km.csv
#              (columns: trait, level, estimate, lwr, upr, p - output of
#               tree/bird/mammal_beta_trait() in functions/fun_beta_trait.R)
#            data/clean/{tree,bird,mammal}.svg   taxon silhouettes
# Outputs  : output/figures/FigureS14_forestplot_multiscale.png
#

# =============================================================================

source(here::here("code", "00_config.R"))
# ── Libraries ────────────────────────────────────────────────────────────────
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(magick)   # SVG -> raster + colorization
library(grid)     # rasterGrob / gtable manipulation
library(gtable)

# ── 0. Parameters ────────────────────────────────────────────────────────────
data_dir <- DIR_CLEAN   # *_data_<scale>.csv and taxon SVG icons
out_png  <- p_fig("FigureS14_forestplot_multiscale.png")

p_threshold <- 0.05

# Input files per taxon, keyed by spatial scale
file_spec <- list(
  Tree = c(
    "1km"  = "tree_data_1km.csv",
    "5km"  = "tree_data_5km.csv",
    "10km" = "tree_data_10km.csv",
    "20km" = "tree_data_20km.csv",
    "50km" = "tree_data_50km.csv"
  ),
  Bird = c(
    "400m" = "bird_data_400m.csv",
    "1km"  = "bird_data_1km.csv",
    "5km"  = "bird_data_5km.csv",
    "10km" = "bird_data_10km.csv",
    "20km" = "bird_data_20km.csv"
  ),
  Mammal = c(
    "1km"  = "mammal_data_1km.csv",
    "5km"  = "mammal_data_5km.csv",
    "10km" = "mammal_data_10km.csv",
    "20km" = "mammal_data_20km.csv",
    "50km" = "mammal_data_50km.csv"
  )
)

# Global ordering of spatial scales (small -> large), SOFT cool->warm palette
scale_levels <- c("400m", "1km", "5km", "10km", "20km", "50km")
scale_colors <- c(
  "400m" = "#6FA8C7",  # soft blue
  "1km"  = "#7DC0A7",  # soft teal
  "5km"  = "#AFD27E",  # soft green
  "10km" = "#F2D06B",  # soft yellow
  "20km" = "#EFA96B",  # soft orange
  "50km" = "#E08A95"   # soft rose
)

taxon_levels <- c("Tree", "Bird", "Mammal")

# Taxon icon colours (as in the manuscript figure)
taxon_colors <- c("Tree"   = "#62b187",
                  "Bird"   = "#65afc1",
                  "Mammal" = "#e8964c")

icon_paths <- c(
  "Tree"   = file.path(data_dir, "tree.svg"),
  "Bird"   = file.path(data_dir, "bird.svg"),
  "Mammal" = file.path(data_dir, "mammal.svg")
)
icon_render_px <- 220L   # render resolution for crisp strip icons

# ── 1. Label maps ─────────────────────────────────────────────────────────────
# Continuous traits (level == trait): pretty names
cont_label_map <- c(
  # Trees
  "swd"            = "Wood density",
  "maxHt"          = "Maximum height",
  "leafnitro"      = "Leaf nitrogen",
  "SLA"            = "Specific leaf area",
  "Seed.dry.mass"  = "Seed dry mass",
  # Birds
  "Beak.Length_Culmen" = "Beak length (culmen)",
  "Beak.Length_Nares"  = "Beak length (nares)",
  "Beak.Width"     = "Beak width",
  "Beak.Depth"     = "Beak depth",
  "Tarsus.Length"  = "Tarsus length",
  "Wing.Length"    = "Wing length",
  "Kipps.Distance" = "Kipp’s distance",
  "Hand.Wing.Index" = "Hand-wing index",
  "Tail.Length"    = "Tail length",
  "Mass"           = "Body mass",
  "Clutch_Max"     = "Clutch size (max)",
  "Incu2"          = "Incubation period",
  # Mammals
  "AET_Mean_mm"          = "Actual evapotranspiration",
  "PET_Mean_mm"          = "Potential evapotranspiration",
  "Precip_Mean_mm"       = "Mean precipitation",
  "Temp_Mean_01degC"     = "Mean temperature",
  "HuPopDen_Mean_n.km2"  = "Human population density",
  "InterbirthInterval_d" = "Interbirth interval",
  "SexualMaturityAge_d"  = "Sexual maturity age",
  "HabitatBreadth"       = "Habitat breadth",
  # Shared niche / range
  "range_km2"      = "Range size (km²)",
  "vol_total"      = "Niche vol. total",
  "vol_core"       = "Niche vol. core",
  "vol_periphery"  = "Niche vol. periphery"
)

# Explicit factor-level labels: "<level> (vs <reference>)"
level_label_map <- c(
  # Trees
  "PFT_binaryEvergreen"                  = "Evergreen (vs deciduous)",
  # Birds - Habitat density (ref: dense)
  "Habitat.DensityOpen"                  = "Open habitat (vs dense)",
  "Habitat.DensitySemi-open"             = "Semi-open habitat (vs dense)",
  # Birds - Migration (ref: sedentary)
  "MigrationMigratory"                   = "Migratory (vs sedentary)",
  "MigrationPartially migratory"         = "Partially migratory (vs sedentary)",
  # Birds - Primary lifestyle (ref: aerial)
  "Primary.LifestyleAquatic"             = "Aquatic (vs aerial)",
  "Primary.LifestyleGeneralist"          = "Generalist (vs aerial)",
  "Primary.LifestyleInsessorial"         = "Insessorial (vs aerial)",
  "Primary.LifestyleTerrestrial"         = "Terrestrial (vs aerial)",
  # Birds - Trophic level (ref: herbivore)
  "Trophic.LevelCarnivore"               = "Carnivore (vs herbivore)",
  "Trophic.LevelOmnivore"                = "Omnivore (vs herbivore)",
  "Trophic.LevelScavenger"               = "Scavenger (vs herbivore)",
  # Birds - Trophic niche (ref: aquatic predator)
  "Trophic.NicheFrugivore"               = "Frugivore (vs aquatic predator)",
  "Trophic.NicheGranivore"               = "Granivore (vs aquatic predator)",
  "Trophic.NicheHerbivore aquatic"       = "Aquatic herbivore (vs aquatic predator)",
  "Trophic.NicheHerbivore terrestrial"   = "Terrestrial herbivore (vs aquatic predator)",
  "Trophic.NicheInvertivore"             = "Invertivore (vs aquatic predator)",
  "Trophic.NicheNectarivore"             = "Nectarivore (vs aquatic predator)",
  "Trophic.NicheOmnivore"                = "Omnivore (vs aquatic predator)",
  "Trophic.NicheScavenger"               = "Scavenger (vs aquatic predator)",
  "Trophic.NicheVertivore"               = "Vertivore (vs aquatic predator)",
  # Mammals - Forage stratum (ref: ground)
  "ForStrat.ValueArboreal"               = "Arboreal (vs ground)",
  "ForStrat.ValueScansorial"             = "Scansorial (vs ground)",
  # Mammals - Terrestriality (ref: ground dwelling)
  "TerrestrialityAbove-ground dwelling"  = "Above-ground dwelling (vs ground dwelling)",
  # Mammals - Activity cycle (ref: nocturnal)
  "ActivityCycleMixed"                   = "Mixed activity (vs nocturnal)",
  "ActivityCycleDiurnal only"            = "Diurnal only (vs nocturnal)",
  # Mammals - Habitat breadth (drop the trailing "2")
  "HabitatBreadth2"                      = "Habitat breadth"
)

# Per-taxon PCo labels (same level string, different meaning by taxon)
pco_label_map <- list(
  Bird = c(
    "PCo1" = "PCo1: small body / dense habitat",
    "PCo2" = "PCo2: dispersal capacity",
    "PCo3" = "PCo3: fast life history"
  ),
  Mammal = c(
    "PCo1" = "PCo1: small body / fast life history"
  )
)

make_label <- function(taxon, trait, level) {
  # 1. PCo: taxon-specific
  if (trait %in% c("PCo1", "PCo2", "PCo3")) {
    ov <- pco_label_map[[taxon]]
    if (!is.null(ov) && level %in% names(ov)) return(unname(ov[[level]]))
    return(level)
  }
  # 2. explicit factor-level label
  if (level %in% names(level_label_map)) return(unname(level_label_map[[level]]))
  # 3. continuous trait
  if (identical(level, trait) && trait %in% names(cont_label_map))
    return(unname(cont_label_map[[trait]]))
  # 4. unmapped factor level: strip trait prefix, keep level only
  if (!identical(level, trait) && startsWith(level, trait)) {
    return(str_trim(sub(paste0("^", trait), "", level)))
  }
  # 5. fallback cleanup
  str_replace_all(level, "[._]", " ")
}

# ── 2. Read + combine all taxa / scales, keep significant only ───────────────
read_one <- function(path, taxon, scale_lab) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  df$taxon <- taxon
  df$scale <- scale_lab
  df
}

all_rows <- map_dfr(names(file_spec), function(tx) {
  spec <- file_spec[[tx]]
  map_dfr(names(spec), function(sc) {
    fp <- file.path(data_dir, spec[[sc]])
    if (!file.exists(fp)) { warning("Missing file: ", fp); return(NULL) }
    read_one(fp, tx, sc)
  })
})

sig_df <- all_rows %>%
  filter(!is.na(p), p < p_threshold) %>%
  mutate(
    taxon = factor(taxon, levels = taxon_levels),
    scale = factor(scale, levels = scale_levels),
    display_label = pmap_chr(list(as.character(taxon), trait, level), make_label)
  )

message(sprintf("Significant rows (p < %.2f): %d", p_threshold, nrow(sig_df)))

# ── 3. y-axis keys + ordering (within each taxon facet) ──────────────────────
sig_df <- sig_df %>% mutate(ykey = paste(taxon, display_label, sep = "___"))

order_tbl <- sig_df %>%
  group_by(taxon, ykey) %>%
  summarise(med_est = median(estimate), .groups = "drop") %>%
  arrange(taxon, med_est)

sig_df <- sig_df %>% mutate(ykey = factor(ykey, levels = order_tbl$ykey))
y_label_fun <- function(x) sub("^.*___", "", x)

# ── 4. Colorize SVG icons (luminance-mask method, as in figure3) ─────────────
colorize_svg <- function(svg_path, color, size_px) {
  img  <- image_read_svg(svg_path, height = size_px)
  info <- image_info(img)
  mask <- img %>% image_convert(colorspace = "Gray") %>% image_negate()
  canvas <- image_blank(width = info$width, height = info$height, color = color)
  image_composite(canvas, mask, operator = "CopyOpacity")
}

icon_grobs <- lapply(taxon_levels, function(tx) {
  img <- colorize_svg(icon_paths[[tx]], taxon_colors[[tx]], icon_render_px)
  info <- image_info(img)
  asp  <- info$width / info$height          # width / height
  w0   <- 0.95                              # icon width in cm
  grobTree(
    rectGrob(gp = gpar(col = NA, fill = "grey94")),
    rasterGrob(as.raster(img), interpolate = TRUE,
               width  = unit(w0, "cm"),
               height = unit(w0 / asp, "cm"))
  )
})
names(icon_grobs) <- taxon_levels

# ── 5. Forest plot ───────────────────────────────────────────────────────────
dodge <- position_dodge(width = 0.6)
x_min <- min(sig_df$lwr, na.rm = TRUE)
x_max <- max(sig_df$upr, na.rm = TRUE)
x_pad <- 0.05 * (x_max - x_min)

p <- ggplot(sig_df,
            aes(x = estimate, y = ykey, color = scale, group = scale)) +
  geom_vline(xintercept = 0, color = "grey10",
             linewidth = 0.5, linetype = "dashed") +
  geom_errorbar(aes(xmin = lwr, xmax = upr),
                orientation = "y", width = 0, linewidth = 0.6,
                position = dodge) +
  geom_point(shape = 21, size = 2.6, stroke = 0.9,
             fill = "white", position = dodge) +
  facet_grid(taxon ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = y_label_fun) +
  scale_color_manual(values = scale_colors, drop = FALSE, name = "Spatial scale") +
  scale_x_continuous(limits = c(x_min - x_pad, x_max + x_pad),
                     labels = function(x) sprintf("%.2f", x)) +
  labs(x = "Standardized effect on species' abundance response to EH", y = NULL) +
  guides(color = guide_legend(override.aes = list(size = 3.2), nrow = 1)) +
  theme_classic(base_size = 12) +
  theme(
    strip.text         = element_blank(),   # replaced by icon grobs below
    strip.background   = element_blank(),
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.4),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    panel.spacing.y    = unit(0.25, "cm"),
    axis.text.y  = element_text(size = 9),
    axis.text.x  = element_text(size = 10),
    axis.title.x = element_text(size = 12, margin = margin(t = 8)),
    legend.position  = "top",
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 10),
    plot.margin = margin(10, 16, 10, 8)
  )

# ── 6. Replace right facet strips with colorized taxon icons ─────────────────
g <- ggplotGrob(p)

# Right-side strips are named strip-r-1 (top) .. strip-r-3 (bottom),
# matching taxon_levels order (Tree, Bird, Mammal).
strip_names <- g$layout$name[grepl("^strip-r", g$layout$name)]
strip_names <- strip_names[order(g$layout$t[match(strip_names, g$layout$name)])]

for (i in seq_along(strip_names)) {
  nm  <- strip_names[i]
  idx <- which(g$layout$name == nm)
  g$grobs[[idx]] <- icon_grobs[[ taxon_levels[i] ]]
  # widen the strip column so the icon is clearly visible
  g$widths[g$layout$l[idx]] <- unit(1.1, "cm")
}

# ── 7. Save ──────────────────────────────────────────────────────────────────
n_labels  <- dplyr::n_distinct(sig_df$ykey)
plot_h_in <- max(8, n_labels * 0.30 + 2)
plot_w_in <- 11

png(out_png, width = plot_w_in, height = plot_h_in, units = "in", res = 300,
    bg = "white")
grid::grid.newpage()
grid::grid.draw(g)
dev.off()

message("Saved: ", out_png, sprintf("  (%.0f x %.1f in)", plot_w_in, plot_h_in))
