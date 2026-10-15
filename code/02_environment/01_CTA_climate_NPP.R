# =============================================================================
# 01_CTA_climate_NPP.R
# Year-specific climate and NPP around camera-trap arrays (mammals)
#
# Source   : EmbeddingBio/analysis/Env_Vars.R (lines ~365-385). The original
#            repeated the same block once per buffer; it is written here as a
#            loop over buffers with identical logic.
# Inputs   : data/raw/gee_exports/camera_array_centroids_NPP_Climate_1_5_10_20_50km/
#              camera_array_centroids_NPP_Climate_{1,5,10,20,50}km.shp
#              Earth Engine export: per-array, per-year buffer means of
#              temperature, precipitation and NPP (2019-2023).
#              TODO(reviewer M5): name the EE climate/NPP collections and link
#              the EE script in the /gee repository.
# Outputs  : data/clean/CTA_MAP_MAT_NPP.csv   one row per Year x Camera_Trap_Array
# =============================================================================

source(here::here("code", "00_config.R"))

library(sf)
library(dplyr)
library(tidyr)
library(purrr)

cta_clim_dir <- p_gee("camera_array_centroids_NPP_Climate_1_5_10_20_50km")

read_cta_climate <- function(scale) {
  st_read(file.path(cta_clim_dir, paste0("camera_array_centroids_NPP_Climate_", scale, ".shp")),
          quiet = TRUE) %>%
    st_drop_geometry() %>%
    rename(Camera_Trap_Array = Cmr_T_A) %>%
    pivot_longer(
      cols      = -Camera_Trap_Array,
      names_to  = c(".value", "Year"),
      names_sep = "_"
    ) %>%
    mutate(Year = as.integer(Year)) %>%
    relocate(Year, .before = Camera_Trap_Array) %>%
    rename_with(~ paste0(.x, "_", scale), .cols = -c(Year, Camera_Trap_Array))
}

CTA_MAP_MAT_NPP <- SCALES_MAMMAL %>%
  map(read_cta_climate) %>%
  reduce(left_join, by = c("Year", "Camera_Trap_Array"))

write.csv(CTA_MAP_MAT_NPP, p_clean("CTA_MAP_MAT_NPP.csv"), row.names = FALSE)
