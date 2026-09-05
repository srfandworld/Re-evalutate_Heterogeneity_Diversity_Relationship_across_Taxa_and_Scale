# =============================================================================
# Spatial models across taxa and scales
# Re-evaluating the Environmental Heterogeneity-Diversity Relationship

# =============================================================================

# ---- Libraries --------------------------------------------------------------
library(tidyverse)    # dplyr, tidyr, stringr, purrr, readr
library(sf)
library(lwgeom)       # st_startpoint()  (bird route start points)
library(spdep)        # knn2nb(), knearneigh(), nb2listw()
library(spatialreg)   # spautolm()
library(glue)

# =============================================================================
# 1. Load the per-scale tables into named objects
#    (produced by the FIA / BBS / mammal prep scripts)
# =============================================================================
scales_tree   <- c("1km", "5km", "10km", "20km", "50km")
scales_bird   <- c("400m", "1km", "5km", "10km", "20km")
scales_mammal <- c("1km", "5km", "10km", "20km", "50km")

for (s in scales_tree)   assign(glue("FIA_table_{s}"),    read.csv(glue("data_clean/FIA_table_{s}.csv")))
for (s in scales_bird)   assign(glue("BBS_table_{s}"),    read.csv(glue("data_clean/BBS_table_{s}.csv")))
for (s in scales_mammal) assign(glue("mammal_table_{s}"), read.csv(glue("data_clean/mammal_table_{s}.csv")))

# Geometry sources for the spatial weights (read once, reused)
routes_sf_bird     <- st_read("data/bbsroutes_shp/bbsrtsl020_filterd.shp", quiet = TRUE)
centroids_sf_mammal <- read_sf("data_clean/camera_array_centroids.shp") %>%
  rename(Camera_Trap_Array = Cmr_T_A)

# =============================================================================
# 2. Spatial-weights builders (k-nearest-neighbour, row-standardised)
# =============================================================================

# Birds: filter route lines to `dat`, union per route, reorder to `dat`, then
# use each route's START POINT as its coordinate.
build_lw_bird <- function(dat, routes_sf, k = 8) {
  dat_key <- str_c(dat$RTENO, dat$RTENAME, sep = "_")
  
  geom <- routes_sf %>%
    mutate(RTE = str_c(RTENO, RTENAME, sep = "_")) %>%
    filter(RTE %in% dat_key) %>%
    group_by(RTE) %>%
    summarise(geometry = st_union(geometry), .groups = "drop") %>%
    arrange(match(RTE, dat_key))
  
  stopifnot(nrow(geom) == nrow(dat))   # listw must align with data row-for-row
  coords <- st_coordinates(lwgeom::st_startpoint(geom))
  nb <- knn2nb(knearneigh(coords, k = k))
  nb2listw(nb, style = "W")
}

# Mammals: `dat` is one row per array; join centroids in `dat` order and use
# the centroid coordinates directly.
build_lw_mammal <- function(dat, centroids_sf, k = 8) {
  stopifnot(all(dat$Camera_Trap_Array %in% centroids_sf$Camera_Trap_Array))
  
  geom <- dat %>%
    dplyr::select(Camera_Trap_Array) %>%
    left_join(centroids_sf, by = "Camera_Trap_Array") %>%
    st_as_sf()
  
  stopifnot(nrow(geom) == nrow(dat))
  coords <- st_coordinates(geom)
  nb <- knn2nb(knearneigh(coords, k = k))
  nb2listw(nb, style = "W")
}

# =============================================================================
# 3. Modelling functions
# =============================================================================

# Trees: ordinary least squares (no spatial term).
sar_modeling_tree <- function(scale, y, table_env = parent.frame()) {
  dat <- get(glue("FIA_table_{scale}"), envir = table_env) %>% drop_na()
  
  form <- if (y == "richness") {
    as.formula(glue(
      "log(richness_{scale}/effort_{scale}) ~ EH_{scale}_t + I(EH_{scale}_t^2) + ",
      "MAT_{scale}_t + MAP_{scale}_t + soilph_{scale}_t + soilcec_{scale}_t + elev_mean_{scale}_t"))
  } else {
    as.formula(glue(
      "log(MSA_new_{scale}) ~ EH_{scale}_t + I(EH_{scale}_t^2) + ",
      "MAT_{scale}_t + MAP_{scale}_t + soilph_{scale}_t + soilcec_{scale}_t + elev_mean_{scale}_t"))
  }
  lm(form, data = dat)   # explicit return
}

# Birds: simultaneous autoregressive model.
sar_modeling_bird <- function(scale, y, routes_sf = routes_sf_bird, k = 8,
                              table_env = parent.frame()) {
  dat <- get(glue("BBS_table_{scale}"), envir = table_env) %>% drop_na()
  lw  <- build_lw_bird(dat, routes_sf, k = k)
  
  form <- as.formula(glue(
    "{y} ~ EH_{scale}_t + I(EH_{scale}_t^2) + MAT_{scale}_t + MAP_{scale}_t + ",
    "Tsea_{scale}_t + Psea_{scale}_t + elev_{scale}_t + PopDen_{scale}_t"))
  spautolm(formula = form, data = dat, listw = lw, family = "SAR")
}

# Mammals: aggregate to one row per array, then SAR.
sar_modeling_mammal <- function(scale, y, centroids_sf = centroids_sf_mammal, k = 8,
                                table_env = parent.frame()) {
  dat <- get(glue("mammal_table_{scale}"), envir = table_env) %>%
    drop_na() %>%
    group_by(Camera_Trap_Array) %>%
    summarise(across(c(richness_Estimator, abundance, MSA, contains("_t")),
                     ~ mean(.x, na.rm = TRUE)),
              .groups = "drop")
  lw <- build_lw_mammal(dat, centroids_sf, k = k)
  
  form <- as.formula(glue(
    "{y} ~ EH_{scale}_t + I(EH_{scale}_t^2) + elev_{scale}_t + PopDen_{scale}_t + ",
    "NPP_{scale}_t + Prec_{scale}_t + Temp_{scale}_t"))
  spautolm(formula = form, data = dat, listw = lw, family = "SAR")
}

# =============================================================================
# 4. Run all models (named by scale) and save
# =============================================================================

# ---- Trees ----
sar_tree_richness <- setNames(lapply(scales_tree, sar_modeling_tree, y = "richness"),  scales_tree)
sar_tree_MSA      <- setNames(lapply(scales_tree, sar_modeling_tree, y = "log(MSA)"),  scales_tree)
write_rds(sar_tree_richness, "D:/BiodiversityEmbedding/data/sar_tree_richness.rds")
write_rds(sar_tree_MSA,      "D:/BiodiversityEmbedding/data/sar_tree_MSA.rds")

# ---- Birds ----
sar_bird_richness <- setNames(lapply(scales_bird, sar_modeling_bird, y = "richness"),  scales_bird)
sar_bird_MSA      <- setNames(lapply(scales_bird, sar_modeling_bird, y = "log(MSA)"),  scales_bird)
write_rds(sar_bird_richness, "D:/BiodiversityEmbedding/data/sar_bird_richness.rds")
write_rds(sar_bird_MSA,      "D:/BiodiversityEmbedding/data/sar_bird_MSA.rds")

# ---- Mammals ----
sar_mammal_richness <- setNames(lapply(scales_mammal, sar_modeling_mammal, y = "richness_Estimator"), scales_mammal)
sar_mammal_MSA      <- setNames(lapply(scales_mammal, sar_modeling_mammal, y = "log(MSA)"),            scales_mammal)
write_rds(sar_mammal_richness, "D:/BiodiversityEmbedding/data/sar_mammal_richness.rds")
write_rds(sar_mammal_MSA,      "D:/BiodiversityEmbedding/data/sar_mammal_MSA.rds")