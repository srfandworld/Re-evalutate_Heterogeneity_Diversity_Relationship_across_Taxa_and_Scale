# =============================================================================
# 05_niche_space_PCA.R
# 14-variable environmental stack -> CONUS PCA (niche space) -> KDE bandwidth
#
# Source   : stack       - EmbeddingBio/analysis/FIA tree traits.R, lines 111-124
#            PCA + H     - EmbeddingBio/analysis/niche_breadth_usingRangeMap.R,
#                          lines 1-65 and 76-110
# Inputs   : rasters under HDR_EXT_DIR/environmental_data (WorldClim 2 BIO1,4,5,6,
#            12,13,14,15; SoilGrids pH/SOC/N; SMAP; SRTM elevation + slope)
# Outputs  : data/clean/envstack_for_FIAniche.tif   (~10 GB)
#            data/clean/pca_us_climate.rds, env_us_scaled.rds,
#            scale_center.rds, scale_sd.rds, H_global.rds
#
# TODO(reviewer M5): the source built the stack from objects defined in
#   Env_Vars.R, where SoilPH_proj/SOC_proj/SoilNitro_proj are EPSG:5070 while
#   MAT is the global EPSG:4326 WorldClim layer. They are transcribed literally
#   below; confirm which soil rasters were in the saved
#   envstack_for_FIAniche.tif (Fig. S5 shows global extent) before re-running.
# TODO(reviewer M5): Fig. S5 (global maps of the 14 layers) has no plotting
#   script; see 07_supplementary_figures/FigS05_niche_env_layers.R.
# =============================================================================

source(here::here("code", "00_config.R"))

library(terra)
library(sf)
library(dplyr)
library(ks)

# ---- 1. Build the 14-layer stack (skip if already on disk) ------------------
stack_file <- p_clean("envstack_for_FIAniche.tif")
if (!file.exists(stack_file)) {
  wc <- function(i) rast(p_env("wc2.1_30s_bio", paste0("wc2.1_30s_bio_", i, ".tif")))
  MAT <- wc(1); MAP <- wc(12); TempSeas <- wc(4); TempMax <- wc(5); TempMin <- wc(6)
  PrecWet <- wc(13); PrecDry <- wc(14); PrecSeas <- wc(15)

  us_boundary_soil <- st_read(p_ext("tl_2024_us_state", "tl_2024_us_state.shp"))
  SoilMoist <- rast(p_env("SMAP_soilmoisture_9km", "SMAP_mean_2017_2023_world.tif"))
  us_boundary_soil <- st_transform(us_boundary_soil, crs(SoilMoist))
  to_proj <- function(r) r %>% crop(us_boundary_soil) %>% project(paste0("EPSG:", TARGET_CRS), method = "bilinear")
  SoilPH_proj    <- to_proj(rast(p_env("soilgrids", "phh2o_0-5cm_mean_1000.tif")))
  SOC_proj       <- to_proj(rast(p_env("soilgrids", "soc_0-5cm_mean_1000.tif")))
  SoilNitro_proj <- to_proj(rast(p_env("soilgrids", "nitrogen_0-5cm_mean_1000.tif")))

  Elev  <- do.call(mosaic, lapply(list.files(p_env("STRM_ele"), pattern = "Elevation_1km", full.names = TRUE), rast))
  Slope <- do.call(mosaic, lapply(list.files(p_env("STRM_ele"), pattern = "Slope_1km",     full.names = TRUE), rast))

  env_stack <- c(MAT, MAP, TempSeas, TempMax, TempMin, PrecWet, PrecDry, PrecSeas,
                 resample(crop(SoilPH_proj,    ext(MAT)), MAT, method = "bilinear"),
                 resample(crop(SOC_proj,       ext(MAT)), MAT, method = "bilinear"),
                 resample(crop(SoilNitro_proj, ext(MAT)), MAT, method = "bilinear"),
                 resample(crop(SoilMoist,      ext(MAT)), MAT, method = "bilinear"),
                 resample(crop(Elev,           ext(MAT)), MAT, method = "bilinear"),
                 resample(crop(Slope,          ext(MAT)), MAT, method = "bilinear"))
  names(env_stack) <- c("MAT","MAP","Tseas","Tmax","Tmin","Pmax","Pmin","Pseas",
                        "SoilPH","SOC","SoilNitro","SoilMoist","Elev","Slope")
  writeRaster(env_stack, stack_file)
}

vars <- c("MAT","MAP","Tseas","Tmax","Tmin","Pmax","Pmin","Pseas",
          "SoilPH","SOC","SoilNitro","SoilMoist","Elev","Slope")

# ---- 2. Coarsen to ~10 km and clip to CONUS ---------------------------------
env_stack        <- terra::rast(stack_file)
env_stack_coarse <- terra::aggregate(env_stack, fact = 12, fun = "mean")

us_boundary <- st_read(p_ext("tl_2024_us_state", "tl_2024_us_state.shp")) %>%
  st_union() %>% st_transform(4326)
us_vect <- terra::vect(us_boundary)
env_us  <- terra::mask(terra::crop(env_stack_coarse, us_vect), us_vect)

# ---- 3. Scale and run PCA ---------------------------------------------------
env_df <- as.data.frame(env_us, na.rm = TRUE)
colnames(env_df) <- vars
env_df_clean <- env_df %>% filter(complete.cases(.))

env_scaled <- scale(env_df_clean)

n_cells <- nrow(env_scaled)
if (n_cells > 100000) {
  set.seed(42)
  pca_input <- env_scaled[sample(n_cells, 100000), ]
} else {
  pca_input <- env_scaled
}
pca_us <- prcomp(pca_input, center = FALSE, scale. = FALSE)   # already scaled

saveRDS(pca_us,     p_clean("pca_us_climate.rds"))
saveRDS(env_scaled, p_clean("env_us_scaled.rds"))
saveRDS(attr(scale(env_df_clean), "scaled:center"), p_clean("scale_center.rds"))
saveRDS(attr(scale(env_df_clean), "scaled:scale"),  p_clean("scale_sd.rds"))

# ---- 4. Global KDE bandwidth in PC1-3 space (shared by all species) ---------
us_pc_df <- as.data.frame(predict(pca_us, newdata = env_scaled)[, 1:3])
colnames(us_pc_df) <- c("PC1", "PC2", "PC3")
set.seed(42)
sample_idx <- sample(nrow(us_pc_df), min(10000, nrow(us_pc_df)))
H_global   <- Hpi(as.matrix(us_pc_df[sample_idx, ]))
saveRDS(H_global, p_clean("H_global.rds"))
