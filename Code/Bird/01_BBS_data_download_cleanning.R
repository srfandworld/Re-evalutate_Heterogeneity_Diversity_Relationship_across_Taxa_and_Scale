# =============================================================================
# BBS birds: species composition, route matching, diversity, and per-scale tables
# Re-evaluating the Environmental Heterogeneity-Diversity Relationship
#
# Sections:
#   1. Load BBS data (bbsBayes2) and build the species lookup
#   2. Observer handling + organise BCR route-year counts
#   3. Route start-points (sf)
#   4. Match BBS line routes (bbsrtsl020) to BCR routes -> lookups `a` and `b`
#   5. Filter the route shapefiles to matched routes
#   6. BBS diversity (richness / abundance / MSA) per route
#   7. Join environmental covariates -> BBS_table
#   8. Per-scale standardised tables (400 m, 1, 5, 10, 20 km) + export
#
# NOTE: this script WRITES data/a.rds and data/b.rds, which the environmental-
#       covariate script reads. Run this first.
# =============================================================================

# ---- Libraries --------------------------------------------------------------
library(tidyverse)   # dplyr, tidyr, stringr, purrr, readr, ggplot2, tibble
library(bbsBayes2)   # load_bbs_data()
library(sf)
library(lwgeom)      # st_startpoint(), st_endpoint()  (was loaded mid-script)


dir.create("data_clean", showWarnings = FALSE)
dir.create("data/bbsroutes_shp", showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. Load BBS data and species lookup
# =============================================================================
dat <- load_bbs_data(level = "state", release = 2024, sample = FALSE, quiet = FALSE)

species_hash <- dat$species %>%
  mutate(ScientificName = paste0(genus, " ", species))
write.csv(species_hash, "data_clean/BBS_species_hash.csv", row.names = FALSE)

# =============================================================================
# 2. Observers + organised BCR route-year counts
# =============================================================================
# Observer IDs / rpid for single annual surveys
observers_df <- dat$routes %>%
  mutate(route       = str_c(country_num, state_num, route, sep = "-"),
         observer_id = as.factor(obs_n)) %>%
  filter(rpid == 101) %>%
  select(year, route, route_name, observer_id)

# Counts with observer ID and a first-year-effect flag
organized_data_bcr <- dat$birds %>%
  filter(rpid == 101) %>%
  mutate(route = str_c(country_num, state_num, route, sep = "-")) %>%
  dplyr::select(year, route, bcr, aou, species_total) %>%
  left_join(observers_df, by = c("year", "route")) %>%
  na.omit() %>%
  dplyr::select(year, route, bcr, observer_id, aou, count = species_total) %>%
  group_by(observer_id) %>%
  mutate(first_year_effect = if_else(year == min(year), 1, 0)) %>%
  ungroup()

# =============================================================================
# 3. Route start-points (sf)
# =============================================================================
routes_sf_bcr <- dat$routes %>%
  mutate(route = str_c(country_num, state_num, route, sep = "-")) %>%
  filter(route %in% unique(organized_data_bcr$route)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  dplyr::select(route, bcr) %>%
  distinct()
write_sf(routes_sf_bcr, "data/bbsroutes_shp/bbs_startpoints.shp")  # write_sf overwrites

# =============================================================================
# 4. Match BBS line routes to BCR routes -> lookups `a` and `b`
# =============================================================================
bbsrtsl020 <- st_read("data/bbsroutes_shp/bbsrtsl020.shp") %>%
  dplyr::select(RTENO, RTENAME) %>%
  st_transform(crs = 5070)
# QC: do RTENO + RTENAME uniquely identify rows?
n_distinct(st_drop_geometry(bbsrtsl020)[, 1:2]) == nrow(bbsrtsl020)

routes_sf_bcr <- st_read("data/bbsroutes_shp/bbs_startpoints.shp") %>%
  st_transform(crs = 5070)

# For each line route, match whichever endpoint (start or end) is nearer to a
# BCR start-point, and record that route + distance.
start_pts <- st_startpoint(bbsrtsl020)       # extract endpoints
end_pts   <- st_endpoint(bbsrtsl020)

nearest_start <- st_nearest_feature(start_pts, routes_sf_bcr)   # find nearest BCR point
nearest_end   <- st_nearest_feature(end_pts,   routes_sf_bcr)

pt_start <- routes_sf_bcr[nearest_start, ]   # matched points
pt_end   <- routes_sf_bcr[nearest_end, ]

dist_start <- st_distance(start_pts, pt_start, by_element = TRUE)   # distances
dist_end   <- st_distance(end_pts,   pt_end,   by_element = TRUE)

use_start <- dist_start <= dist_end          # choose the closer end
nearest_route <- ifelse(use_start, routes_sf_bcr$route[nearest_start],
                        routes_sf_bcr$route[nearest_end])
nearest_dist  <- ifelse(use_start, dist_start, dist_end)   # note: ifelse drops units (metres)

result <- bbsrtsl020 %>%
  mutate(route = nearest_route, dist = nearest_dist)   # attach back to the route lines

# For each RTENO+RTENAME keep the closest match (< 1 km), then one route per line
b <- result %>%
  st_drop_geometry() %>%
  group_by(RTENO, RTENAME) %>% arrange(dist) %>% slice(1) %>% ungroup() %>%
  filter(dist < 1000) %>%
  group_by(route) %>% arrange(dist) %>% slice(1) %>% ungroup()

# `a`: BBSRTSL020 <-> RTENO/RTENAME crosswalk (from a GEE export)
export_table <- read.csv("data/export_table.csv")
n_distinct(export_table$BBSRTSL020)
# QC (interactive): row-order check against BBS_EH_50km, which is built in the
# covariate script - left commented since it isn't defined here.
# export_table$RTENAME == BBS_EH_50km$RTENAME
a <- export_table %>% dplyr::select(RTENO, RTENAME, BBSRTSL020)

write_rds(a, "data/a.rds")
write_rds(b, "data/b.rds")
a <- read_rds("data/a.rds")
b <- read_rds("data/b.rds")

# =============================================================================
# 5. Filter the route shapefiles to matched routes
# =============================================================================
bbsrtsl020_filtered <- bbsrtsl020 %>%
  filter(paste(RTENO, RTENAME) %in% paste(b$RTENO, b$RTENAME))
routes_sf_bcr_filtered <- routes_sf_bcr %>%
  filter(route %in% b$route)

write_sf(bbsrtsl020_filtered,    "data/bbsroutes_shp/bbsrtsl020_filterd.shp")
write_sf(routes_sf_bcr_filtered, "data/bbsroutes_shp/bbs_startpoints_filterd.shp")

routes_sf_bcr_filtered <- st_read("data/bbsroutes_shp/bbs_startpoints_filterd.shp")
bbsrtsl020_filtered    <- st_read("data/bbsroutes_shp/bbsrtsl020_filterd.shp")

# =============================================================================
# 6. BBS diversity (2017-2023)
# =============================================================================
# Average counts per route x species across surveyed years to damp annual
# variation (observer identity, day-to-day bird behaviour, etc.)
BBS_data <- organized_data_bcr %>%
  filter(year > 2016) %>%
  filter(route %in% unique(routes_sf_bcr_filtered$route)) %>%
  group_by(route, aou) %>%
  mutate(nyear              = n_distinct(year),
         annual_mean_count  = sum(count) / nyear) %>%
  ungroup()
write.csv(BBS_data, "data_clean/BBS_data.csv", row.names = FALSE)
BBS_data <- read.csv("data_clean/BBS_data.csv")

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
BBS_env <- read.csv("data_clean/BBS_env.csv")
BBS_table <- BBS_diversity %>%
  left_join(BBS_env, by = "route")
write.csv(BBS_table, "data_clean/BBS_table.csv", row.names = FALSE)
BBS_table <- read.csv("data_clean/BBS_table.csv")

# Dominant land-cover codes -> labels
BBS_table <- BBS_table %>%
  mutate(across(contains("dominant_lc"), ~ factor(
    .x,
    levels = c(1, 2, 3, 4, 5, 7, 9, 81, 82),
    labels = c("Open water", "Developed", "Barren land", "Forest",
               "Shrub/Scrub", "Grassland/Herbaceous", "Wetlands",
               "Pasture/Hay", "Cultivated crops")
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
  write.csv(tbl, file.path("data_clean", paste0("BBS_table_", scale, ".csv")),
            row.names = FALSE)
})