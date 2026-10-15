# =============================================================================
# 03_bird_community.R
# BBS route-level richness / abundance / MSA (2017-2023) and per-scale tables
#
# Source   : Code/Bird/01_BBS_data_download_cleanning.R, sections 6-8
# Inputs   : data/clean/BBS_organized_counts.rds   (01_data_acquisition/02)
#            data/clean/bbsroutes_shp/bbs_startpoints_filterd.shp
#            data/clean/BBS_env.csv                (02_environment/02)
# Outputs  : data/clean/BBS_data.csv     route x species x year counts
#            data/clean/BBS_table.csv
#            data/clean/BBS_table_{400m,1km,5km,10km,20km}.csv  (model input)
#
# TODO(reviewer M1): 584 species reported; harmonised AOU set elsewhere is 520.
#   Reconcile the AOU codes counted here (forms/hybrids/unidentified?).
# TODO(reviewer, Table S1 / "Transformations"): MAP is NOT sqrt-transformed for
#   birds below, and PopDen uses log1p(). Align the Methods or the code.
# =============================================================================

source(here::here("code", "00_config.R"))

library(tidyverse)
library(sf)

organized_data_bcr     <- readRDS(p_clean("BBS_organized_counts.rds"))
routes_sf_bcr_filtered <- st_read(p_clean("bbsroutes_shp", "bbs_startpoints_filterd.shp"))
# =============================================================================
# 6. BBS diversity (2017-2023)
# =============================================================================
# Average counts per route x species across surveyed years to damp annual
# variation (observer identity, day-to-day bird behaviour, etc.)
BBS_data <- organized_data_bcr %>%
  filter(year > 2016) %>%
  filter(route %in% unique(routes_sf_bcr_filtered$route)) %>%
  group_by(route, aou) %>%
  # TODO(reviewer M6): nyear counts years the SPECIES was detected, not years
  # the ROUTE was surveyed. Confirm which denominator is intended for MSA.
  mutate(nyear              = n_distinct(year),
         annual_mean_count  = sum(count) / nyear) %>%
  ungroup()
write.csv(BBS_data, p_clean("BBS_data.csv"), row.names = FALSE)
BBS_data <- read.csv(p_clean("BBS_data.csv"))

BBS_diversity <- BBS_data %>%
  dplyr::select(route, bcr, aou, annual_mean_count) %>%
  distinct() %>%
  group_by(aou) %>%
  mutate(max_abundance_species = max(annual_mean_count)) %>%
  ungroup() %>%
  group_by(route, bcr) %>%
  summarise(
    richness  = n_distinct(aou),
    abundance = sum(annual_mean_count / max_abundance_species) / n_distinct(aou),
    MSA       = sum(annual_mean_count) / n_distinct(aou),
    .groups = "drop"
  )   # FIX: removed dplyr::collect() (no lazy/remote source here)

# =============================================================================
# 7. Join environmental covariates -> BBS_table
# =============================================================================
BBS_env <- read.csv(p_clean("BBS_env.csv"))
BBS_table <- BBS_diversity %>%
  left_join(BBS_env, by = "route")
write.csv(BBS_table, p_clean("BBS_table.csv"), row.names = FALSE)
BBS_table <- read.csv(p_clean("BBS_table.csv"))

# Dominant land-cover codes -> labels
BBS_table <- BBS_table %>%
  mutate(across(contains("dominant_lc"), ~ factor(
    .x,
    levels = LC_LEVELS,
    labels = LC_LABELS
  )))

# =============================================================================
# 8. Per-scale standardised tables (400 m, 1, 5, 10, 20 km)
# =============================================================================
BBS_table_20km <- BBS_table %>%
  dplyr::select(route, bcr, richness, abundance, MSA, RTENO, RTENAME, contains("_20km")) %>%
  mutate(EH_20km_t     = scale(EH_20km)[, 1],
         MAT_20km_t    = scale(MAT_20km)[, 1],
         MAP_20km_t    = scale(MAP_20km)[, 1],
         Tsea_20km_t   = scale(Tsea_20km)[, 1],
         Psea_20km_t   = scale(Psea_20km)[, 1],
         elev_20km_t   = scale(sqrt(elev_20km))[, 1],
         PopDen_20km_t = scale(log1p(PopDen_20km))[, 1])

BBS_table_10km <- BBS_table %>%
  dplyr::select(route, bcr, richness, abundance, MSA, RTENO, RTENAME, contains("_10km")) %>%
  mutate(EH_10km_t     = scale(EH_10km)[, 1],
         MAT_10km_t    = scale(MAT_10km)[, 1],
         MAP_10km_t    = scale(MAP_10km)[, 1],
         Tsea_10km_t   = scale(Tsea_10km)[, 1],
         Psea_10km_t   = scale(Psea_10km)[, 1],
         elev_10km_t   = scale(sqrt(elev_10km))[, 1],
         PopDen_10km_t = scale(log1p(PopDen_10km))[, 1])

BBS_table_5km <- BBS_table %>%
  dplyr::select(route, bcr, richness, abundance, MSA, RTENO, RTENAME, contains("_5km")) %>%
  mutate(EH_5km_t     = scale(EH_5km)[, 1],
         MAT_5km_t    = scale(MAT_5km)[, 1],
         MAP_5km_t    = scale(MAP_5km)[, 1],
         Tsea_5km_t   = scale(Tsea_5km)[, 1],
         Psea_5km_t   = scale(Psea_5km)[, 1],
         elev_5km_t   = scale(sqrt(elev_5km))[, 1],
         PopDen_5km_t = scale(log1p(PopDen_5km))[, 1])

BBS_table_1km <- BBS_table %>%
  dplyr::select(route, bcr, richness, abundance, MSA, RTENO, RTENAME, contains("_1km")) %>%
  mutate(EH_1km_t     = scale(EH_1km)[, 1],
         MAT_1km_t    = scale(MAT_1km)[, 1],
         MAP_1km_t    = scale(MAP_1km)[, 1],
         Tsea_1km_t   = scale(Tsea_1km)[, 1],
         Psea_1km_t   = scale(Psea_1km)[, 1],
         elev_1km_t   = scale(sqrt(elev_1km))[, 1],
         PopDen_1km_t = scale(log1p(PopDen_1km))[, 1])

BBS_table_400m <- BBS_table %>%
  dplyr::select(route, bcr, richness, abundance, MSA, RTENO, RTENAME, contains("_400m")) %>%
  mutate(EH_400m_t     = scale(EH_400m)[, 1],
         MAT_400m_t    = scale(MAT_400m)[, 1],
         MAP_400m_t    = scale(MAP_400m)[, 1],
         Tsea_400m_t   = scale(Tsea_400m)[, 1],
         Psea_400m_t   = scale(Psea_400m)[, 1],
         elev_400m_t   = scale(sqrt(elev_400m))[, 1],
         PopDen_400m_t = scale(log1p(PopDen_400m))[, 1])

# ---- Write out the five per-scale tables ------------------------------------
BBS_tables_by_scale <- list(
  "400m" = BBS_table_400m,
  "1km"  = BBS_table_1km,
  "5km"  = BBS_table_5km,
  "10km" = BBS_table_10km,
  "20km" = BBS_table_20km
)

iwalk(BBS_tables_by_scale, function(tbl, scale) {
  write.csv(tbl, p_clean(paste0("BBS_table_", scale, ".csv")),
            row.names = FALSE)
})
