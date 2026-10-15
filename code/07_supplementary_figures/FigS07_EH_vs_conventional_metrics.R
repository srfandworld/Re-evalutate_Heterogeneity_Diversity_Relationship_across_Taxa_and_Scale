# =============================================================================
# FigS07_EH_vs_conventional_metrics.R
# Fig. S7: maps of the seven conventional heterogeneity components (SD of
# elevation, soil pH, temperature, precipitation; NLCD Shannon diversity; SD of
# tree canopy cover; SD of canopy height) next to AEF-based EH
#
# Source   : EmbeddingBio/analysis/sfigure_EH_covariate_map.R (paths only changed)
#            TODO: confirm this script produced Fig. S7 (file name says
#            "EH covariate map"; it may also cover parts of Fig. S6).
# Inputs   : data/raw/gee_exports/MeanCosineDist_5km_2020.shp
#            HDR_EXT_DIR/environmental_data/{soilgrids, Elevation_US_5070.tif,
#              tcc_US_5070.tif, NLCD_US_5070.tif, CHM_SD_CONUS_5070.tif}
# TODO(reviewer C3): no source code for the EH validation analysis; the
#   temperature/precipitation SD layers are not read here - locate them.
# =============================================================================

source(here::here("code", "00_config.R"))

#sfigure_EH characterization
library(scales)
library(sf)
library(terra)
library(patchwork)
library(exactextractr)

###map extent and graticule###
{
  # ── 1. Map extent ─────────────────────────────────────────────────────────────
  xlim_alb <- c(-2400000, 2300000)
  ylim_alb <- c(200000,   3200000)
  
  extent_poly <- st_bbox(
    c(xmin = xlim_alb[1], xmax = xlim_alb[2],
      ymin = ylim_alb[1], ymax = ylim_alb[2]),
    crs = st_crs(a)
  ) |> st_as_sfc()
  
  # ── 2. Graticule clipped inside box ──────────────────────────────────────────
  graticule <- st_graticule(
    lon = seq(-130, -60, by = 10),
    lat = seq(25, 55, by = 5),
    crs = st_crs(4326)
  ) |>
    st_transform(st_crs(a)) |>
    st_intersection(extent_poly)
  
  # ── 3. Build box border lines for intersection ────────────────────────────────
  # Bottom edge and left edge as sf linestrings in Albers
  bottom_edge <- st_linestring(matrix(
    c(xlim_alb[1], ylim_alb[1],
      xlim_alb[2], ylim_alb[1]), ncol = 2, byrow = TRUE
  )) |> st_sfc(crs = st_crs(a))
  
  left_edge <- st_linestring(matrix(
    c(xlim_alb[1], ylim_alb[1],
      xlim_alb[1], ylim_alb[2]), ncol = 2, byrow = TRUE
  )) |> st_sfc(crs = st_crs(a))
  
  # ── 4. Find exact intersection points of graticule with box edges ─────────────
  # Full graticule (unclipped) for intersection — need lines that reach the border
  graticule_full <- st_graticule(
    lon = seq(-130, -60, by = 10),
    lat = seq(25, 55, by = 5),
    crs = st_crs(4326)
  ) |> st_transform(st_crs(a))
  
  # Meridians (lon lines) intersect the bottom edge → lon labels
  meridians <- graticule_full |>
    filter(type == "E", degree %in% seq(-120, -70, by = 10))
  
  lon_pts <- st_intersection(meridians, bottom_edge) |>
    mutate(label = paste0(abs(degree), "°W"))
  
  # Parallels (lat lines) intersect the left edge → lat labels  
  parallels <- graticule_full |>
    filter(type == "N", degree %in% seq(30, 50, by = 5))
  
  lat_pts <- st_intersection(parallels, left_edge) |>
    mutate(label = paste0(degree, "°N"))
}

# ── Data ───────────────────────────────────────────────────────────────────────
states    <- sf::st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))
EH_US_5km <- st_read(p_gee("MeanCosineDist_5km_2020.shp")) |>
  st_transform(4326)
soilph       <- rast(p_env("soilgrids/phh2o_0-5cm_mean_1000.tif"))
elev         <- rast(p_env("Elevation_US_5070.tif"))
MAT_proj     <- rast(p_clean("MAT_usa_5070.tif"))
MAP_proj     <- rast(p_clean("MAP_usa_5070.tif"))
tcc_proj <- rast(p_env("tcc_US_5070.tif"))
NLCD_proj    <- rast(p_env("NLCD_US_5070.tif"))
chm_std <- rast(p_env("CHM_SD_CONUS_5070.tif"))

# ── Extract std per grid cell ─────────────────────────────────────────────
EH_US_5km$soilph_std <- exact_extract(soilph, EH_US_5km, fun = "stdev")
EH_US_5km$elev_std <- exact_extract(elev, EH_US_5km, fun = "stdev")
EH_US_5km$MAT_std <- exact_extract(MAT_proj, EH_US_5km, fun = "stdev")
EH_US_5km$MAP_std <- exact_extract(MAP_proj, EH_US_5km, fun = "stdev")
EH_US_5km$tcc_std <- exact_extract(tcc_proj, EH_US_5km, fun = "stdev")
EH_US_5km$lc_shannon <- exact_extract(NLCD_proj, EH_US_5km, function(values, coverage_fractions) {
  # Drop NA pixels
  valid <- !is.na(values)
  vals <- values[valid]
  wts  <- coverage_fractions[valid]
  
  if (length(vals) == 0) return(NA_real_)
  
  # Sum coverage-weighted area per land cover class
  class_area <- tapply(wts, vals, sum)
  
  # Convert to proportions
  p <- class_area / sum(class_area)
  
  # Shannon diversity index: H = -sum(p * ln(p))
  -sum(p * log(p))
})
EH_US_5km$chm_std <- exact_extract(chm_std, EH_US_5km, fun = "mean")


# ── Reproject to Albers & mask to CONUS ───────────────────────────────────────
EH_vect   <- vect(EH_US_5km)  |> project(st_crs(a)$wkt)
states_vect <- vect(states)   |> project(st_crs(a)$wkt)
EH_vect   <- EH_vect[states_vect, ]   # spatial subset instead of mask

# ── Rasterize chm_std ─────────────────────────────────────────────────────────
template  <- rast(EH_vect, resolution = 5000)   # 5 km in metres (Albers)
soilph_std_rast <- rasterize(EH_vect, template, field = "soilph_std", fun = "mean") |>
  clamp(lower = 0, values = TRUE)
elev_std_rast <- rasterize(EH_vect, template, field = "elev_std", fun = "mean") |>
  clamp(lower = 0, values = TRUE)
MAT_std_rast <- rasterize(EH_vect, template, field = "MAT_std", fun = "mean") |>
  clamp(lower = 0, values = TRUE)
MAP_std_rast <- rasterize(EH_vect, template, field = "MAP_std", fun = "mean") |>
  clamp(lower = 0, values = TRUE)
tcc_std_rast <- rasterize(EH_vect, template, field = "tcc_std", fun = "mean") |>
  clamp(lower = 0, values = TRUE)
lc_shannon_rast <- rasterize(EH_vect, template, field = "lc_shannon", fun = "mean") |>
  clamp(lower = 0, values = TRUE)
chm_std_rast <- rasterize(EH_vect, template, field = "chm_std", fun = "mean") |>
  clamp(lower = 0, values = TRUE)


# ── Plot df ───────────────────────────────────────────────────────────────────
df <- as.data.frame(chm_std_rast, xy = TRUE) |>
  setNames(c("x", "y", "std")) |>
  dplyr::filter(!is.na(std))

vmin <- quantile(df$std, 0.01)
vmax <- quantile(df$std, 0.99)


ggplot() +
  geom_raster(data = df, aes(x = x, y = y, fill = std)) +  
  scale_fill_gradientn(
    colors   = viridis::viridis(100),
    limits   = c(vmin, vmax),
    oob      = squish,
    na.value = NA,
    name     = "Standard deviation of canopy height",  #Standard deviation of tree canopy cover
    guide    = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = unit(10, "cm"),
      barheight      = unit(0.4, "cm"),
      ticks.colour   = "black",
      frame.colour   = "black"
    )
  ) +
  
  geom_sf(data = states, fill = NA, color = "white", linewidth = 0.3) +
  geom_sf(data = graticule, color = "grey10", linewidth = 0.25,
          linetype = "dashed") +
  geom_sf_text(data = lon_pts, aes(label = label),
               size = 4, color = "black",
               vjust = 1.8, fontface = "plain") +
  geom_sf_text(data = lat_pts, aes(label = label),
               size = 4, color = "black",
               hjust = 1.2, fontface = "plain") +
  
  coord_sf(
    crs    = st_crs(a),
    xlim   = xlim_alb,
    ylim   = ylim_alb,
    expand = FALSE,
    clip   = "off"
  ) +
  
  theme_minimal() +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid       = element_blank(),
    panel.background = element_rect(fill = "white"),
    axis.title       = element_blank(),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.margin    = margin(t = 10, b = 2),
    legend.title     = element_text(size = 15),
    legend.text      = element_text(size = 12),
    plot.margin      = margin(t = 5, r = 15, b = 20, l = 50)
  )

