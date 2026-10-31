# =============================================================================
# fun_hdr_sar_models.R - model data, spatial weights and SAR model functions for
# the HDR models (richness / MSA ~ EH + EH^2 + covariates)
#
# Source   : Code/Figure2/01_3taxa_modeling.R, sections 1-3 (2026-10-30
#            version; paths only changed). Sourcing this file:
#   - builds the per-scale MODEL data FIA_table_<s>, BBS_table_<s>,
#     mammal_table_<s> (join EH_delta + A57 EH, keep EH_delta > 0, standardize
#     EH and covariates within each scale's sample, drop NAs in modelled columns
#     only) and writes them to data/clean/model_data/*_model_<s>.csv
#   - defines build_lw_{tree,bird,mammal}() and sar_modeling_{tree,bird,mammal}()
# The same data and models are used by Fig. 2, Tables S2-S4, Figs. S9-S11 and
# Figs. S17-S19.
# Used by  : 05_HDR_models/01_fit_HDR_models.R
# =============================================================================

source(here::here("code", "00_config.R"))

library(dplyr)
library(tidyr)        # drop_na()
library(stringr)      # str_c()
library(sf)
library(lwgeom)       # st_startpoint()  (bird route start points)
library(spdep)        # knn2nb(), knearneigh(), nb2listw()
library(spatialreg)   # spautolm()
library(glue)

# =============================================================================
# 1. Build the per-scale model data into named objects
#    Same data as the alternative-EH figures (sfigure17-19_physicallyVar_*.R),
#    Tables S2-S3 (tableS2_S3_HDR.R) and Figs. S9-S11 (figureS9_S11.R):
#    - join the mean-independent elevation range (EH_delta) and A57 EH, keep
#      EH_delta > 0
#    - standardize EH and all covariates WITHIN each scale's sample
#    - drop rows with missing values in the modelled columns only
#    Run with the working directory = EmbeddingBio.
# =============================================================================
scales_tree   <- c("1km", "5km", "10km", "20km", "50km")
scales_bird   <- c("400m", "1km", "5km", "10km", "20km")
scales_mammal <- c("1km", "5km", "10km", "20km", "50km")

# Trees (FIA plots; richness = coverage-based qD from iNEXT)
build_model_data_tree <- function(scale) {
  eledelta <- read.csv(glue(p_gee("FIA_embeddings/FIA_EH_delta/FIA_Elevation_EHdelta_{scale}.csv"))) %>%
    dplyr::select(pltID, EH_delta)
  a57 <- st_read(glue(p_gee("FIA_A57_EH/FIA_A57_EH_{scale}.shp")), quiet = TRUE) %>%
    st_drop_geometry() %>% dplyr::select(pltID, A57_sd)

  read.csv(glue(p_clean("FIA_table_{scale}_iNEXT.csv"))) %>%
    dplyr::select(focal_pltID, qD, contains(glue("_{scale}"))) %>%
    left_join(eledelta, by = join_by(focal_pltID == pltID)) %>%
    left_join(a57,      by = join_by(focal_pltID == pltID)) %>%
    filter(EH_delta > 0) %>%
    mutate("EH_{scale}_t"        := scale(.data[[glue("EH_{scale}")]])[, 1],
           EH_elevrange_t        = scale(log(.data[[glue("elev_range_{scale}")]]))[, 1],
           EH_delta_t            = scale(log(EH_delta))[, 1],
           A57_sd_t              = scale(A57_sd)[, 1],
           "soilph_{scale}_t"    := scale(log(.data[[glue("soilph_{scale}")]]))[, 1],
           "soilcec_{scale}_t"   := scale(log(.data[[glue("soilcec_{scale}")]]))[, 1],
           "elev_mean_{scale}_t" := scale(sqrt(.data[[glue("elev_mean_{scale}")]]))[, 1],
           "MAP_{scale}_t"       := scale(sqrt(.data[[glue("MAP_{scale}")]]))[, 1],
           "MAT_{scale}_t"       := scale(.data[[glue("MAT_{scale}")]])[, 1]) %>%
    drop_na(qD, all_of(glue("MSA_new_{scale}")), ends_with("_t"))
}

# Birds (BBS routes). NOTE: EH_delta (per BBSRTSL020) and A57 (per RTENO/RTENAME)
# are joined by row order (cbind), as in sfigure18_physicallyVar_BBS.R
BBS_table <- read.csv(p_clean("BBS_table.csv"))
build_model_data_bird <- function(scale) {
  eledelta <- read.csv(glue(p_gee("BBS_EH_delta/BBS_Elevation_EHdelta_{scale}.csv"))) %>%
    dplyr::select(BBSRTSL020, elev_range, EH_delta)
  a57 <- st_read(glue(p_gee("BBS_A57_EH/BBS_A57_EH_{scale}.shp")), quiet = TRUE) %>%
    st_drop_geometry() %>% dplyr::select(RTENAME, RTENO, A57_sd)
  merged <- cbind(eledelta, a57) %>%
    dplyr::select(RTENAME, RTENO, elev_range, EH_delta, A57_sd) %>%
    group_by(RTENAME, RTENO) %>%
    summarise("elev_range_{scale}" := mean(elev_range),
              EH_delta = mean(EH_delta),
              A57_sd   = mean(A57_sd), .groups = "drop")

  BBS_table %>%
    dplyr::select(route, RTENAME, RTENO, richness, MSA, contains(glue("_{scale}"))) %>%
    left_join(merged, by = c("RTENAME", "RTENO")) %>%
    filter(EH_delta > 0) %>%
    mutate("EH_{scale}_t"     := scale(.data[[glue("EH_{scale}")]])[, 1],
           EH_elevrange_t     = scale(log(.data[[glue("elev_range_{scale}")]]))[, 1],
           EH_delta_t         = scale(log(EH_delta))[, 1],
           A57_sd_t           = scale(A57_sd)[, 1],
           "MAT_{scale}_t"    := scale(.data[[glue("MAT_{scale}")]])[, 1],
           "MAP_{scale}_t"    := scale(.data[[glue("MAP_{scale}")]])[, 1],
           "Tsea_{scale}_t"   := scale(.data[[glue("Tsea_{scale}")]])[, 1],
           "Psea_{scale}_t"   := scale(.data[[glue("Psea_{scale}")]])[, 1],
           "elev_{scale}_t"   := scale(sqrt(.data[[glue("elev_{scale}")]]))[, 1],
           "PopDen_{scale}_t" := scale(log1p(.data[[glue("PopDen_{scale}")]]))[, 1]) %>%
    drop_na(richness, MSA, ends_with("_t"))
}

# Mammals (camera-trap arrays): standardize at the array-year level, then
# average to one row per array
mammal_table <- read.csv(p_clean("mammal_table.csv"))
build_model_data_mammal <- function(scale) {
  eledelta <- read.csv(glue(p_gee("mammal_EH_delta/mammal_Elevation_EHdelta_{scale}.csv"))) %>%
    dplyr::select(Camera_Trap_Array, elev_range, EH_delta)
  a57 <- st_read(glue(p_gee("mammal_A57_EH/mammal_A57_EH_{scale}.shp")), quiet = TRUE) %>%
    st_drop_geometry() %>% dplyr::select(Camera_Trap_Array = Camera_Tra, A57_sd)
  merged <- left_join(eledelta, a57, by = "Camera_Trap_Array")
  raw_cols <- c(glue("EH_{scale}"), "elev_range", "EH_delta", "A57_sd")

  mammal_table %>%
    dplyr::select(year_CTA, Year, Camera_Trap_Array,
                  richness_Estimator, MSA, contains(glue("_{scale}"))) %>%
    left_join(merged, by = "Camera_Trap_Array") %>%
    filter(EH_delta > 0) %>%
    mutate("EH_{scale}_t"     := scale(.data[[glue("EH_{scale}")]])[, 1],
           EH_elevrange_t     = scale(log(elev_range))[, 1],
           EH_delta_t         = scale(log(EH_delta))[, 1],
           A57_sd_t           = scale(A57_sd)[, 1],
           "elev_{scale}_t"   := scale(sqrt(.data[[glue("elev_{scale}")]]))[, 1],
           "PopDen_{scale}_t" := scale(log1p(.data[[glue("PopDen_{scale}")]]))[, 1],
           "NPP_{scale}_t"    := scale(.data[[glue("NPP_{scale}")]])[, 1],
           "Prec_{scale}_t"   := scale(sqrt(.data[[glue("Prec_{scale}")]]))[, 1],
           "Temp_{scale}_t"   := scale(.data[[glue("Temp_{scale}")]])[, 1]) %>%
    group_by(Camera_Trap_Array) %>%
    summarise(across(c(richness_Estimator, MSA, ends_with("_t"), all_of(raw_cols)),
                     ~ mean(.x, na.rm = TRUE)),
              .groups = "drop") %>%
    drop_na(richness_Estimator, MSA, ends_with("_t"))
}

for (s in scales_tree)   assign(glue("FIA_table_{s}"),    build_model_data_tree(s))
for (s in scales_bird)   assign(glue("BBS_table_{s}"),    build_model_data_bird(s))
for (s in scales_mammal) assign(glue("mammal_table_{s}"), build_model_data_mammal(s))

# save the model data (Figure2.R reads them for the turning-point back-transform)
model_data_dir <- p_clean("model_data")
dir.create(model_data_dir, showWarnings = FALSE)
for (s in scales_tree)   write.csv(get(glue("FIA_table_{s}")),    file.path(model_data_dir, glue("FIA_model_{s}.csv")),    row.names = FALSE)
for (s in scales_bird)   write.csv(get(glue("BBS_table_{s}")),    file.path(model_data_dir, glue("BBS_model_{s}.csv")),    row.names = FALSE)
for (s in scales_mammal) write.csv(get(glue("mammal_table_{s}")), file.path(model_data_dir, glue("mammal_model_{s}.csv")), row.names = FALSE)

# Geometry sources for the spatial weights (read once, reused)
plots_sf_tree <- st_read(p_clean("FIA_plots.shp")) %>%
  st_transform(crs = 5070)
routes_sf_bird     <- st_read(p_clean("bbsroutes_shp/bbsrtsl020_filterd.shp"), quiet = TRUE)
centroids_sf_mammal <- read_sf(p_clean("camera_array_centroids.shp")) %>%
  rename(Camera_Trap_Array = Cmr_T_A)

# =============================================================================
# 2. Spatial-weights builders (k-nearest-neighbour, row-standardised)
# =============================================================================
# fia trees
build_lw_tree <- function(dat, plots_sf, k = 8) {

  
  geom <- dat %>%
    dplyr::select(focal_pltID) %>%
    left_join(plots_sf, by = join_by(focal_pltID==pltID)) %>%
    st_as_sf()
  
  stopifnot(nrow(geom) == nrow(dat))   # listw must align with data row-for-row
  coords <- st_coordinates(geom)
  nb <- knn2nb(knearneigh(coords, k = k))
  nb2listw(nb, style = "W")
}

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

# The model data (section 1) are already complete for the modelled columns, so
# they are used as they are (no further drop_na / aggregation): the rows match
# the alternative-EH figures and Tables S2-S3 exactly.

# Trees: simultaneous autoregressive model (FIA plot locations, EPSG:5070).
sar_modeling_tree <- function(scale, y, plots_sf = plots_sf_tree, k = 8,
                              table_env = parent.frame()) {
  dat <- get(glue("FIA_table_{scale}"), envir = table_env)
  lw  <- build_lw_tree(dat, plots_sf, k = k)
  
  form <- if (y == "richness") {
    as.formula(glue(
      "log(qD) ~ EH_{scale}_t + I(EH_{scale}_t^2) + ",
      "MAT_{scale}_t + MAP_{scale}_t + soilph_{scale}_t + soilcec_{scale}_t + elev_mean_{scale}_t"))
  } else {
    as.formula(glue(
      "log(MSA_new_{scale}) ~ EH_{scale}_t + I(EH_{scale}_t^2) + ",
      "MAT_{scale}_t + MAP_{scale}_t + soilph_{scale}_t + soilcec_{scale}_t + elev_mean_{scale}_t"))
  }
  spautolm(formula = form,
                        data = dat,
                        listw = lw, 
                        method = "LU",
                        family = "SAR")
}

# Birds: simultaneous autoregressive model.
sar_modeling_bird <- function(scale, y, routes_sf = routes_sf_bird, k = 8,
                              table_env = parent.frame()) {
  dat <- get(glue("BBS_table_{scale}"), envir = table_env)
  lw  <- build_lw_bird(dat, routes_sf, k = k)
  
  form <- as.formula(glue(
    "{y} ~ EH_{scale}_t + I(EH_{scale}_t^2) + MAT_{scale}_t + MAP_{scale}_t + ",
    "Tsea_{scale}_t + Psea_{scale}_t + elev_{scale}_t + PopDen_{scale}_t"))
  spautolm(formula = form, data = dat, listw = lw, family = "SAR")
}

# Mammals: one row per array (aggregated in section 1), then SAR.
sar_modeling_mammal <- function(scale, y, centroids_sf = centroids_sf_mammal, k = 8,
                                table_env = parent.frame()) {
  dat <- get(glue("mammal_table_{scale}"), envir = table_env)
  lw <- build_lw_mammal(dat, centroids_sf, k = k)
  
  form <- as.formula(glue(
    "{y} ~ EH_{scale}_t + I(EH_{scale}_t^2) + elev_{scale}_t + PopDen_{scale}_t + ",
    "NPP_{scale}_t + Prec_{scale}_t + Temp_{scale}_t"))
  spautolm(formula = form, data = dat, listw = lw, family = "SAR")
}

