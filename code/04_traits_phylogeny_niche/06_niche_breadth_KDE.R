# =============================================================================
# 06_niche_breadth_KDE.R
# Niche breadth (3-D KDE volumes) and range size for trees, birds, mammals
#
# Source   : birds/mammals - EmbeddingBio/analysis/niche_breadth_usingRangeMap.R,
#                            lines 262-341 (run once per taxon by swapping range_sf)
#            BBS join      - EmbeddingBio/analysis/BBS_traits_phylogeny.R, lines 464-474
#            trees         - EmbeddingBio/analysis/niche_breadth_usingOccurence.R,
#                            lines 11-44, 196-235, 307-331
#            QC plots dropped; logic unchanged.
# Inputs   : data/clean/envstack_for_FIAniche.tif, pca_us_climate.rds,
#            scale_center.rds, scale_sd.rds, H_global.rds   (05_niche_space_PCA.R)
#            data/clean/{bird,mammal}_species_range_us.shp   (04_species_ranges.R)
#            data/clean/tpa_all5.rds, BBS_species2.csv
# Outputs  : data/clean/niche_results/Mammal_niche_breadth_3d.csv
#            data/clean/niche_results/Bird_niche_breadth_3d.csv
#            data/clean/niche_results/BBSBird_niche_breadth_3d.csv  (+ aou)
#            data/clean/niche_results/FIA_niche_breadth_3d.csv
# Feeds    : Fig. 3, Table S5 (vol_total, vol_core, vol_periphery, range_km2)
# NOTE     : range size is range-polygon area for birds/mammals but the convex
#            hull of FIA occurrences for trees; say so in the Methods.
# =============================================================================

source(here::here("code", "00_config.R"))
source(here::here("code", "functions", "fun_niche_breadth.R"))

library(tidyr)

env_stack_coarse <- terra::aggregate(terra::rast(p_clean("envstack_for_FIAniche.tif")),
                                     fact = 12, fun = "mean")
pca_us       <- readRDS(p_clean("pca_us_climate.rds"))
scale_center <- readRDS(p_clean("scale_center.rds"))
scale_sd     <- readRDS(p_clean("scale_sd.rds"))
H_global     <- readRDS(p_clean("H_global.rds"))

safe_val <- function(x, default = NA) {
  if (is.null(x) || length(x) == 0) return(default)
  return(x)
}

# =============================================================================
# 1. Birds and mammals: range-map niche breadth + range-polygon area
# =============================================================================
run_range_niche <- function(range_sf) {
  species_list <- unique(range_sf$species)
  niche_results <- lapply(species_list, function(sp) {
    calculate_niche_3d(
      species_name      = sp,
      range_sf          = range_sf,
      env_stack_coarse  = env_stack_coarse,
      pca_model         = pca_us,
      sc_center         = scale_center,
      sc_sd             = scale_sd,
      H                 = H_global
    )
  })

  niche_df <- bind_rows(lapply(niche_results, function(x) {
    if (is.null(x)) return(NULL)
    data.frame(
      species       = safe_val(x$species,       NA_character_),
      vol_total     = safe_val(x$vol_total,     NA_real_),
      vol_core      = safe_val(x$vol_core,      NA_real_),
      vol_periphery = safe_val(x$vol_periphery, NA_real_),
      n_cells       = safe_val(x$n_cells,       NA_integer_),
      n_obs         = safe_val(x$n_obs,         NA_integer_),
      status        = safe_val(x$status,        "error"),
      error_msg     = safe_val(x$error,         NA_character_)
    )
  })) %>% arrange(desc(vol_total))

  range_df <- st_transform(range_sf, TARGET_CRS) %>%
    mutate(range_km2 = as.numeric(set_units(st_area(.), km^2))) %>%
    st_drop_geometry() %>%
    dplyr::select(species, range_km2)

  niche_df %>%
    left_join(range_df, by = "species") %>%
    filter(status == "success") %>%
    drop_na(vol_total, range_km2)
}

mammal_niche <- run_range_niche(st_read(p_clean("mammal_species_range_us.shp")))
write.csv(mammal_niche, p_clean("niche_results", "Mammal_niche_breadth_3d.csv"), row.names = FALSE)

bird_niche <- run_range_niche(st_read(p_clean("bird_species_range_us.shp")))
write.csv(bird_niche, p_clean("niche_results", "Bird_niche_breadth_3d.csv"), row.names = FALSE)

# Attach BBS AOU codes
BBS_species2 <- read.csv(p_clean("BBS_species2.csv"))
bbs_bird_niche <- bird_niche %>%
  left_join(BBS_species2 %>% select(aou, species),
            by = c("species" = "species"), suffix = c("", "_sp"))
write.csv(bbs_bird_niche, p_clean("niche_results", "BBSBird_niche_breadth_3d.csv"))

# =============================================================================
# 2. Trees: occurrence-based niche breadth + convex-hull range size
# =============================================================================
tpa_all5 <- readRDS(p_clean("tpa_all5.rds"))
dt0 <- tpa_all5 %>%
  dplyr::select(pltID, SCIENTIFIC_NAME, BAA, geometry) %>%
  filter(!is.na(pltID), !is.na(SCIENTIFIC_NAME), !is.na(BAA))

dt0_coords <- dt0 %>%
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2]) %>%
  st_drop_geometry() %>%
  dplyr::select(species = SCIENTIFIC_NAME, lon, lat) %>%
  distinct()

species_count <- dt0_coords %>%
  group_by(species) %>%
  summarise(n_occ = n(), .groups = "drop") %>%
  filter(n_occ >= 50)                         # minimum occurrences threshold
dt0_filtered <- dt0_coords %>% filter(species %in% species_count$species)

species_list <- unique(dt0_filtered$species)
niche_results_occ <- lapply(species_list, function(sp) {
  calculate_niche_3d_occ(
    species_name     = sp,
    occ_data         = dt0_filtered,
    env_stack_coarse = env_stack_coarse,
    pca_model        = pca_us,
    sc_center        = scale_center,
    sc_sd            = scale_sd,
    H                = H_global
  )
})

niche_df_occ <- bind_rows(lapply(niche_results_occ, function(x) {
  if (is.null(x)) return(NULL)
  data.frame(
    species       = safe_val(x$species,       NA_character_),
    vol_total     = safe_val(x$vol_total,     NA_real_),
    vol_core      = safe_val(x$vol_core,      NA_real_),
    vol_periphery = safe_val(x$vol_periphery, NA_real_),
    n_occ         = safe_val(x$n_occ,         NA_integer_),
    n_cells       = safe_val(x$n_cells,       NA_integer_),
    status        = safe_val(x$status,        "error"),
    error_msg     = safe_val(x$error,         NA_character_)
  )
})) %>% drop_na(vol_total) %>% arrange(desc(vol_total))

range_convex <- bind_rows(lapply(species_list, calculate_range_convexhull,
                                 occ_data = dt0_filtered))

niche_results <- niche_df_occ %>%
  dplyr::select(species, vol_total, vol_core, vol_periphery, n_occ) %>%
  left_join(range_convex %>% dplyr::select(species, range_km2, n_points), by = "species")
write.csv(niche_results, p_clean("niche_results", "FIA_niche_breadth_3d.csv"), row.names = FALSE)
