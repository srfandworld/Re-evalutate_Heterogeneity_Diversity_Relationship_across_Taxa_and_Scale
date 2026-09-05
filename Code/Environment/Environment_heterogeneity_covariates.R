# =============================================================================
# Environmental variables and heterogeneity (EH) covariates
# Re-evaluating the Environmental Heterogeneity-Diversity Relationship
#
# Builds per-unit environmental tables for three taxa:
#   FIA_env  -> FIA tree plots            -> data_clean/FIA_env.csv
#   BBS_env  -> BBS bird routes           -> data_clean/BBS_env.csv
#   CTA_env  -> mammal camera-trap arrays -> data_clean/Camera_Trap_Array_env.csv
#
# Each table = environmental heterogeneity (AEF embedding dissimilarity) +
# climate, soil, elevation, productivity, land cover, and population covariates
# extracted within multiple buffer radii.
# =============================================================================

# ---- Libraries --------------------------------------------------------------
library(terra)
library(sf)
library(dplyr)
library(tidyr)       # pivot_longer(), drop_na()
library(stringr)     # str_replace(), str_c()
library(purrr)       # map2(), reduce()
library(readr)       # read_rds()          <- still missing in the last version
library(spatialEco)

# ---- Configuration ----------------------------------------------------------
# Point `base_dir` at the raw external rasters/shapefiles; everything else is
# relative to the project root. (Replaces the repeated "D:/..." paths.)
base_dir  <- "D:/BiodiversityEmbedding"   # raw external data root  (EDIT ME)
data_dir  <- "data"                        # inputs shipped with the project
clean_dir <- "data_clean"                  # processed inputs + outputs

target_crs <- 5070   # EPSG:5070, CONUS Albers Equal Area
env_path <- function(...) file.path(base_dir, "environmental_data", ...)

# Crop a raster to the (current-CRS) boundary, then reproject to EPSG:5070
reproject_to_crs <- function(r, boundary) {
  r %>% crop(boundary) %>% project(paste0("EPSG:", target_crs), method = "bilinear")
}

# =============================================================================
# Helper functions
# =============================================================================

# Mean of a continuous raster within circular buffers around points
calculate_buffer_means <- function(points, raster, buffer_distances, var_name, id_col) {
  results <- data.frame(id = points[[id_col]])
  names(results)[1] <- id_col
  for (buf in buffer_distances) {
    cat("Processing buffer:", buf, "m for", var_name, "\n")
    buffers <- st_buffer(points, dist = buf)
    means   <- extract(raster, buffers, fun = "mean", na.rm = TRUE, touches = TRUE)
    results[[paste0(var_name, "_", buf, "m")]] <- means[, 2]
  }
  results
}

# Modal (dominant) class of a categorical raster within buffers
calculate_buffer_mode <- function(points, raster, buffer_distances, var_name, id_col) {
  results <- data.frame(id = points[[id_col]])
  names(results)[1] <- id_col
  for (buf in buffer_distances) {
    cat("Processing buffer:", buf, "m for", var_name, "\n")
    buffers <- st_buffer(points, dist = buf)
    modes   <- extract(raster, buffers, fun = modal, na.rm = TRUE, touches = FALSE)
    results[[paste0(var_name, "_", buf, "m")]] <- modes[, 2]
  }
  results
}

# =============================================================================
# Boundary
# =============================================================================
us_boundary <- st_read(file.path(base_dir, "tl_2024_us_state", "tl_2024_us_state.shp"))

# =============================================================================
# SHARED RASTER LAYERS (loaded once, reused across all three taxa)
# =============================================================================

# ---- Climate: WorldClim 2 bioclim (1970-2000, 30 arc-sec) -------------------
MAT      <- rast(env_path("wc2.1_30s_bio", "wc2.1_30s_bio_1.tif"))   # BIO1  mean annual temp
MAP      <- rast(env_path("wc2.1_30s_bio", "wc2.1_30s_bio_12.tif"))  # BIO12 mean annual precip
TempSeas <- rast(env_path("wc2.1_30s_bio", "wc2.1_30s_bio_4.tif"))   # BIO4  temp seasonality
TempMax  <- rast(env_path("wc2.1_30s_bio", "wc2.1_30s_bio_5.tif"))   # BIO5  max temp warmest month
TempMin  <- rast(env_path("wc2.1_30s_bio", "wc2.1_30s_bio_6.tif"))   # BIO6  min temp coldest month
PrecWet  <- rast(env_path("wc2.1_30s_bio", "wc2.1_30s_bio_13.tif"))  # BIO13 precip wettest month
PrecDry  <- rast(env_path("wc2.1_30s_bio", "wc2.1_30s_bio_14.tif"))  # BIO14 precip driest month
PrecSeas <- rast(env_path("wc2.1_30s_bio", "wc2.1_30s_bio_15.tif"))  # BIO15 precip seasonality

us_boundary <- st_transform(us_boundary, crs(MAT))
MAT_proj      <- reproject_to_crs(MAT,      us_boundary)
MAP_proj      <- reproject_to_crs(MAP,      us_boundary)
TempSeas_proj <- reproject_to_crs(TempSeas, us_boundary)
TempMax_proj  <- reproject_to_crs(TempMax,  us_boundary)
TempMin_proj  <- reproject_to_crs(TempMin,  us_boundary)
PrecWet_proj  <- reproject_to_crs(PrecWet,  us_boundary)
PrecDry_proj  <- reproject_to_crs(PrecDry,  us_boundary)
PrecSeas_proj <- reproject_to_crs(PrecSeas, us_boundary)

# ---- Soil: SoilGrids 2.0 (0-5 cm) + SMAP soil moisture ----------------------
SoilPH    <- rast(env_path("soilgrids", "phh2o_0-5cm_mean_1000.tif"))
SOC       <- rast(env_path("soilgrids", "soc_0-5cm_mean_1000.tif"))
SoilNitro <- rast(env_path("soilgrids", "nitrogen_0-5cm_mean_1000.tif"))
Soilcec   <- rast(env_path("soilgrids", "cec_0-5cm_mean_1000.tif"))
SoilMoist <- rast(env_path("SMAP_soilmoisture_9km", "SMAP_mean_2017_2023_world.tif"))

us_boundary <- st_transform(us_boundary, crs(SoilMoist))
SoilPH_proj    <- reproject_to_crs(SoilPH,    us_boundary)
SOC_proj       <- reproject_to_crs(SOC,       us_boundary)
SoilNitro_proj <- reproject_to_crs(SoilNitro, us_boundary)
Soilcec_proj   <- reproject_to_crs(Soilcec,   us_boundary)
SoilMoist_proj <- reproject_to_crs(SoilMoist, us_boundary)

# ---- Elevation & slope: SRTM 1 km -------------------------------------------
Elev <- do.call(mosaic, lapply(
  list.files(env_path("STRM_ele"), pattern = "Elevation_1km", full.names = TRUE), rast))
Elev_proj <- reproject_to_crs(Elev, us_boundary)

# ---- Net primary production: MODIS MOD17A3 (2020) ---------------------------
NPP <- rast(env_path("NPP", "npp_2020.tif"))
us_vect  <- vect(st_transform(us_boundary, crs(NPP)))
NPP_proj <- NPP %>% crop(us_vect) %>% mask(us_vect)

# ---- Land cover: NLCD (30 m, already EPSG:5070), reclassified ---------------
NLCD_proj <- rast(env_path("NLCD_US_5070.tif"))   # https://www.usgs.gov/node/279743
rcl <- matrix(c(
  11, 11,  1,   # Open water
  21, 21,  2,   # Developed (all intensities)
  22, 22,  2,
  23, 23,  2,
  24, 24,  2,
  31, 31,  3,   # Barren land
  41, 41,  4,   # Forest (deciduous / evergreen / mixed)
  42, 42,  4,
  43, 43,  4,
  52, 52,  5,   # Shrub / scrub
  71, 71,  7,   # Grassland / herbaceous
  81, 81, 81,   # Pasture / hay
  82, 82, 82,   # Cultivated crops
  90, 90,  9,   # Wetlands (woody / emergent herbaceous)
  95, 95,  9
), ncol = 3, byrow = TRUE)
NLCD_grouped <- terra::classify(NLCD_proj, rcl = rcl, others = NA)

# ---- Human population: GPW v4 (2020, 30 arc-sec) ----------------------------
gpw <- terra::rast(env_path(
  "gpw-v4-population-density-rev11_2020_30_sec_tif",
  "gpw_v4_population_density_rev11_2020_30_sec.tif"))
gpw_proj <- reproject_to_crs(gpw, us_boundary)

# =============================================================================
# TREES - FIA plot environmental variables
# =============================================================================
FIA_plots <- st_read(file.path(clean_dir, "FIA_plots.shp")) %>% st_transform(crs = 5070)

# Ecoregion (CEC Level 2 + Level 1)
EcoReg_l2 <- st_read(file.path(base_dir, "na_cec_eco_l2", "NA_CEC_Eco_Level2.shp")) %>%
  st_transform(crs = 5070)
FIA_plots_ecoreg <- st_join(
  FIA_plots, EcoReg_l2[, c("NA_L2CODE", "NA_L2NAME", "NA_L1CODE", "NA_L1NAME")],
  join = st_within, left = TRUE
) %>% st_drop_geometry()

# Environmental heterogeneity (AEF embedding dissimilarity), per buffer
FIA_EH_1km  <- st_read(file.path(data_dir, "FIA_embeddings", "FIA_Dissimilarity_1km_2017_2023_mean_TEST.shp"))  %>% st_drop_geometry() %>% dplyr::select(pltID, EH_1km  = mean_cosin)
FIA_EH_5km  <- st_read(file.path(data_dir, "FIA_embeddings", "FIA_Dissimilarity_5km_2017_2023_mean_TEST.shp"))  %>% st_drop_geometry() %>% dplyr::select(pltID, EH_5km  = mean_cosin)
FIA_EH_10km <- st_read(file.path(data_dir, "FIA_embeddings", "FIA_Dissimilarity_10km_2017_2023_mean_TEST.shp")) %>% st_drop_geometry() %>% dplyr::select(pltID, EH_10km = mean_cosin)
FIA_EH_20km <- st_read(file.path(data_dir, "FIA_embeddings", "FIA_Dissimilarity_20km_2017_2023_mean_TEST.shp")) %>% st_drop_geometry() %>% dplyr::select(pltID, EH_20km = mean_cosin)
FIA_EH_50km <- read.csv(file.path(data_dir, "FIA_embeddings", "FIA_Dissimilarity_50km_2017_2023_mean_TEST.csv")) %>% dplyr::select(pltID, EH_50km = mean_cosin)

FIA_EH <- FIA_plots_ecoreg %>%
  left_join(FIA_EH_1km,  by = "pltID") %>%
  left_join(FIA_EH_5km,  by = "pltID") %>%
  left_join(FIA_EH_10km, by = "pltID") %>%
  left_join(FIA_EH_20km, by = "pltID") %>%
  left_join(FIA_EH_50km, by = "pltID")

fia_buffer <- c(1000, 5000, 10000, 20000, 50000, 110000)

# Climate
fia_MAT <- calculate_buffer_means(FIA_plots, MAT_proj, fia_buffer, "MAT", "pltID")
fia_MAP <- calculate_buffer_means(FIA_plots, MAP_proj, fia_buffer, "MAP", "pltID")
climate_FIA_buffer <- cbind(fia_MAT, fia_MAP[, -1])

# Soil
fia_soilph    <- calculate_buffer_means(FIA_plots, SoilPH_proj,    fia_buffer, "soilph",    "pltID")
fia_soc       <- calculate_buffer_means(FIA_plots, SOC_proj,       fia_buffer, "soc",       "pltID")
fia_soilnitro <- calculate_buffer_means(FIA_plots, SoilNitro_proj, fia_buffer, "soilnitro", "pltID")
fia_soilcec   <- calculate_buffer_means(FIA_plots, Soilcec_proj,   fia_buffer, "soilcec",   "pltID")
fia_soilmoist <- calculate_buffer_means(FIA_plots, SoilMoist_proj, fia_buffer, "soilmoist", "pltID")
soil_FIA_buffer <- cbind(fia_soilph, fia_soc[, -1], fia_soilnitro[, -1],
                         fia_soilcec[, -1], fia_soilmoist[, -1])

# Elevation (pre-computed per-buffer stats shapefiles)
fia_buffer_names <- c("1km", "5km", "10km", "20km", "50km", "110km")
ele_FIA_buffer <- map2(
  lapply(list(1, 5, 10, 20, 50, 110), function(x)
    st_read(env_path("FIA_Elevation_Stats", paste0("FIA_Elevation_", x, "km.shp"))) %>%
      st_drop_geometry()),
  fia_buffer_names,
  function(df, suffix) {
    df %>%
      dplyr::select(-buffer_dis, -starts_with("elev_max"), -starts_with("elev_min")) %>%
      rename_with(~ paste0(., "_", suffix), -pltID)
  }
) %>% reduce(left_join, by = "pltID") %>% dplyr::select(pltID, everything())

# GPP (pre-computed buffered shapefile)
GPP_FIA_buffer <- st_read(env_path("GPP", "FIA_plots_bufferedGPP.shp")) %>%
  st_drop_geometry() %>%
  dplyr::select(pltID, starts_with("GPP_"), -GPP_0km)

# Dominant land cover
fia_lc <- calculate_buffer_mode(FIA_plots, NLCD_grouped,
                                c(1000, 5000, 10000, 20000, 50000), "dominant_lc", "pltID")

# Assemble
FIA_env <- FIA_EH %>%
  left_join(soil_FIA_buffer,    by = "pltID") %>%
  left_join(ele_FIA_buffer,     by = "pltID") %>%
  left_join(GPP_FIA_buffer,     by = "pltID") %>%
  left_join(climate_FIA_buffer, by = "pltID") %>%
  left_join(fia_lc,             by = "pltID") %>%
  rename_with(~ str_replace(., "1000m",   "1km")) %>%
  rename_with(~ str_replace(., "5000m",   "5km")) %>%
  rename_with(~ str_replace(., "10000m",  "10km")) %>%
  rename_with(~ str_replace(., "20000m",  "20km")) %>%
  rename_with(~ str_replace(., "50000m",  "50km")) %>%
  rename_with(~ str_replace(., "110000m", "110km"))

write.csv(FIA_env, file.path(clean_dir, "FIA_env.csv"), row.names = FALSE)

# =============================================================================
# BIRDS - BBS route environmental variables
# =============================================================================
# a = route-ID table row-aligned to the EH shapefiles (attached via cbind)
# b = route lookup with columns RTENO, RTENAME, route (first 3 columns)
a <- read_rds(file.path(data_dir, "a.rds"))
b <- read_rds(file.path(data_dir, "b.rds"))

# Read one BBS EH layer and collapse to one value per route.
# NOTE: cbind(a) aligns IDs by ROW ORDER - it is only correct if `a` is in the
#       exact same order as each EH shapefile. A key-based join is safer if the
#       shapefiles carry a route ID.
read_bbs_eh <- function(path, eh_name) {
  st_read(path) %>%
    st_transform(crs = 5070) %>%
    dplyr::select(RTENAME, !!eh_name := mean_cosin) %>%
    cbind(a) %>% dplyr::select(-RTENAME.1) %>%
    st_drop_geometry() %>%
    left_join(b[, 1:3], by = c("RTENO", "RTENAME")) %>%
    drop_na() %>%
    group_by(RTENAME, RTENO, route) %>%
    summarise(!!eh_name := mean(.data[[eh_name]]), .groups = "drop")
}

bbs_eh_dir <- file.path(data_dir, "BBS_embeddings")
BBS_EH_400m <- read_bbs_eh(file.path(bbs_eh_dir, "BBS_Dissimilarity_400m_2017_2023_mean.shp"), "EH_400m")
BBS_EH_1km  <- read_bbs_eh(file.path(bbs_eh_dir, "BBS_Dissimilarity_1km_2017_2023_mean.shp"),  "EH_1km")
BBS_EH_5km  <- read_bbs_eh(file.path(bbs_eh_dir, "BBS_Dissimilarity_5km_2017_2023_mean.shp"),  "EH_5km")
BBS_EH_10km <- read_bbs_eh(file.path(bbs_eh_dir, "BBS_Dissimilarity_10km_2017_2023_mean.shp"), "EH_10km")
BBS_EH_20km <- read_bbs_eh(file.path(bbs_eh_dir, "BBS_Dissimilarity_20km_2017_2023_mean.shp"), "EH_20km")

BBS_EH <- BBS_EH_400m %>%
  left_join(BBS_EH_1km,  by = c("RTENO", "RTENAME", "route")) %>%
  left_join(BBS_EH_5km,  by = c("RTENO", "RTENAME", "route")) %>%
  left_join(BBS_EH_10km, by = c("RTENO", "RTENAME", "route")) %>%
  left_join(BBS_EH_20km, by = c("RTENO", "RTENAME", "route"))

# Route geometries
bbsrtsl020 <- st_read(file.path(data_dir, "bbsroutes_shp", "bbsrtsl020.shp")) %>%
  st_transform(target_crs)

bbs_buffers <- c(400, 1000, 5000, 10000, 20000, 50000)

# Climate + population
bbs_MAT  <- calculate_buffer_means(bbsrtsl020, MAT_proj,      bbs_buffers, "MAT",    "BBSRTSL020")
bbs_MAP  <- calculate_buffer_means(bbsrtsl020, MAP_proj,      bbs_buffers, "MAP",    "BBSRTSL020")
bbs_Tsea <- calculate_buffer_means(bbsrtsl020, TempSeas_proj, bbs_buffers, "Tsea",   "BBSRTSL020")
bbs_Psea <- calculate_buffer_means(bbsrtsl020, PrecSeas_proj, bbs_buffers, "Psea",   "BBSRTSL020")
bbs_gpw  <- calculate_buffer_means(bbsrtsl020, gpw_proj,      bbs_buffers, "PopDen", "BBSRTSL020")

# Elevation (buffer means, negatives clamped to 0)
bbs_elev <- calculate_buffer_means(bbsrtsl020, Elev_proj, bbs_buffers, "elev", "BBSRTSL020") %>%
  mutate(across(-BBSRTSL020, ~ ifelse(.x < 0, 0, .x)))

# Combine covariates, average per route
bbs_env <- bbs_MAT %>%
  left_join(bbs_MAP,  by = "BBSRTSL020") %>%
  left_join(bbs_Tsea, by = "BBSRTSL020") %>%
  left_join(bbs_Psea, by = "BBSRTSL020") %>%
  left_join(bbs_gpw,  by = "BBSRTSL020") %>%
  left_join(bbs_elev, by = "BBSRTSL020") %>%
  rename_with(~ gsub("1000m",  "1km",  .x)) %>%
  rename_with(~ gsub("5000m",  "5km",  .x)) %>%
  rename_with(~ gsub("10000m", "10km", .x)) %>%
  rename_with(~ gsub("20000m", "20km", .x)) %>%
  rename_with(~ gsub("50000m", "50km", .x)) %>%
  group_by(RTENO, RTENAME) %>%
  summarise(across(-BBSRTSL020, ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  left_join(b[, 1:3], by = c("RTENO", "RTENAME")) %>%
  drop_na(route) %>%
  relocate(route, .after = RTENAME)

# Dominant land cover
BBS_lc <- read.csv(env_path("BBS_dominant_lc.csv")) %>%
  left_join(bbsrtsl020 %>% st_drop_geometry() %>% dplyr::select(BBSRTSL020, RTENAME, RTENO),
            by = "BBSRTSL020") %>%
  group_by(RTENAME, RTENO) %>%
  summarise(dominant_lc_400m = first(dominant_lc_400m), .groups = "drop")
BBS_lc_unique <- BBS_lc %>% distinct(RTENAME, RTENO, .keep_all = TRUE)

# Assemble
BBS_env <- BBS_EH %>%
  left_join(bbs_env,       by = c("RTENO", "RTENAME", "route")) %>%
  left_join(BBS_lc_unique, by = c("RTENAME", "RTENO")) %>%
  rename_with(~ str_replace(., "1000m",  "1km")) %>%
  rename_with(~ str_replace(., "5000m",  "5km")) %>%
  rename_with(~ str_replace(., "10000m", "10km")) %>%
  rename_with(~ str_replace(., "20000m", "20km"))

write.csv(BBS_env, file.path(clean_dir, "BBS_env.csv"), row.names = FALSE)

# =============================================================================
# MAMMALS - camera-trap array (CTA) environmental variables
# =============================================================================
CTA_EH <- read.csv(file.path(clean_dir, "Camera_Trap_Array_EH_2019-2023.csv"))

camera_array_centroids <- read_sf(file.path(clean_dir, "camera_array_centroids.shp")) %>%
  rename(Camera_Trap_Array = Cmr_T_A)
CTA_plots <- st_transform(camera_array_centroids, target_crs)

cta_buffers <- c(1000, 5000, 10000, 20000, 50000)

# Population density
CTA_gpw <- calculate_buffer_means(CTA_plots, gpw_proj, cta_buffers, "PopDen", "Camera_Trap_Array") %>%
  rename(PopDen_1km = PopDen_1000m, PopDen_5km = PopDen_5000m, PopDen_10km = PopDen_10000m,
         PopDen_20km = PopDen_20000m, PopDen_50km = PopDen_50000m)

# Elevation
CTA_elev <- calculate_buffer_means(CTA_plots, Elev_proj, cta_buffers, "elev", "Camera_Trap_Array") %>%
  rename(elev_1km = elev_1000m, elev_5km = elev_5000m, elev_10km = elev_10000m,
         elev_20km = elev_20000m, elev_50km = elev_50000m)

# Climate + NPP by year and buffer (pre-computed)
CTA_MAP_MAT_NPP <- read.csv(file.path(clean_dir, "CTA_MAP_MAT_NPP.csv"))

# Ecoregion (CEC Level 1), with manual fixes for coastal/island arrays
EcoReg_l1 <- st_read(file.path(base_dir, "na_cec_eco_l1", "NA_CEC_Eco_Level1.shp")) %>%
  st_transform(crs = 5070) %>% st_make_valid()
CTA_plots_ecoreg <- st_join(
  CTA_plots, EcoReg_l1[, c("NA_L1CODE", "NA_L1NAME")], join = st_within, left = TRUE
) %>% st_drop_geometry() %>%
  mutate(
    NA_L1CODE = case_when(
      Camera_Trap_Array == "Dare"          ~ "11",
      Camera_Trap_Array == "KeyDeerRefuge" ~ "15",
      Camera_Trap_Array == "PointReyes"    ~ "11",
      TRUE ~ NA_L1CODE),
    NA_L1NAME = case_when(
      Camera_Trap_Array == "Dare"          ~ "MEDITERRANEAN CALIFORNIA",
      Camera_Trap_Array == "KeyDeerRefuge" ~ "TROPICAL WET FORESTS",
      Camera_Trap_Array == "PointReyes"    ~ "MEDITERRANEAN CALIFORNIA",
      TRUE ~ NA_L1NAME)
  )

# Dominant land cover
CTA_lc <- calculate_buffer_mode(CTA_plots, NLCD_grouped, cta_buffers,
                                "dominant_lc", "Camera_Trap_Array")

# Assemble (CTA table is per array-year)
CTA_env <- CTA_EH %>%
  left_join(CTA_elev,         by = "Camera_Trap_Array") %>%
  left_join(CTA_gpw,          by = "Camera_Trap_Array") %>%
  left_join(CTA_MAP_MAT_NPP,  by = c("Year", "Camera_Trap_Array")) %>%
  left_join(CTA_plots_ecoreg, by = "Camera_Trap_Array") %>%
  left_join(CTA_lc,           by = "Camera_Trap_Array") %>%
  mutate(year_CTA = str_c(Year, Camera_Trap_Array, sep = "_")) %>%
  relocate(year_CTA, .before = Year) %>%
  rename_with(~ str_replace(., "1000m",   "1km")) %>%
  rename_with(~ str_replace(., "5000m",   "5km")) %>%
  rename_with(~ str_replace(., "10000m",  "10km")) %>%
  rename_with(~ str_replace(., "20000m",  "20km")) %>%
  rename_with(~ str_replace(., "50000m",  "50km")) %>%
  rename_with(~ str_replace(., "110000m", "110km"))

write.csv(CTA_env, file.path(clean_dir, "Camera_Trap_Array_env.csv"), row.names = FALSE)