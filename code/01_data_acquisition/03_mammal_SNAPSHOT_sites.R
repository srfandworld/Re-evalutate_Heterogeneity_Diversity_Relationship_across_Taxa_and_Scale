# =============================================================================
# 03_mammal_SNAPSHOT_sites.R
# SNAPSHOT USA camera sites -> camera-trap-array (CTA) centroids
#
# Source   : Code/Mammal/01_Mammal_data_download_cleaning.R, section 1
# Inputs   : data/raw/Mammals/ssusa_finalsequences.csv
#            data/raw/Mammals/ssusa_finaldeployments.csv
#              SNAPSHOT USA 2019-2023 (https://doi.org/10.1111/geb.13941);
#              TODO(reviewer M5): add the download URL / access date
# Outputs  : data/clean/camera_array_centroids.shp   -> 02_environment
#            data/clean/CTA_habitat.csv               dominant habitat per array
# Next     : 02_environment/*, then 03_community_metrics/04_mammal_community.R
# =============================================================================

source(here::here("code", "00_config.R"))

library(tidyverse)
library(sf)

ssusa_dep <- read.csv(p_raw("Mammals", "ssusa_finaldeployments.csv"))

# Camera site points
mammal_sites <- ssusa_dep %>%
  distinct(Camera_Trap_Array, Site_Name, Deployment_ID, Latitude, Longitude) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)

# Dominant habitat per camera-trap array
CTA_habitat <- ssusa_dep %>%
  count(Camera_Trap_Array, Habitat) %>%
  group_by(Camera_Trap_Array) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(Camera_Trap_Array, Habitat)
write.csv(CTA_habitat, p_clean("CTA_habitat.csv"), row.names = FALSE)

# Array centroids (union of member sites), with habitat
# TODO(reviewer M1): 262 unique arrays here (3 outside CONUS); manuscript says 263.
camera_array_centroids <- mammal_sites %>%
  group_by(Camera_Trap_Array) %>%
  summarise(geometry = st_centroid(st_union(geometry))) %>%
  left_join(CTA_habitat, by = "Camera_Trap_Array")
write_sf(camera_array_centroids, p_clean("camera_array_centroids.shp"))
