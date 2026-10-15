# =============================================================================
# 02_bird_BBS_routes.R
# Load BBS counts, match route lines to BBS route start points
#
# Source   : Code/Bird/01_BBS_data_download_cleanning.R, sections 1-5
# Inputs   : bbsBayes2::load_bbs_data(release = 2024)  (downloaded on first use)
#            data/raw/bbsroutes_shp/bbsrtsl020.shp      BBS route lines
#              TODO(reviewer M5): document the source/URL of bbsrtsl020.shp
#            data/raw/gee_exports/export_table.csv     BBSRTSL020 <-> RTENO/RTENAME
#              crosswalk exported from Earth Engine with the BBS EH layers
# Outputs  : data/clean/BBS_species_hash.csv
#            data/clean/BBS_organized_counts.rds   route-year-species counts
#            data/clean/a.rds, data/clean/b.rds    route lookups -> 02_environment
#            data/clean/bbsroutes_shp/bbs_startpoints(_filterd).shp,
#            data/clean/bbsroutes_shp/bbsrtsl020_filterd.shp
# Next     : 02_environment/02_environmental_covariates.R, then
#            03_community_metrics/03_bird_community.R
# =============================================================================

source(here::here("code", "00_config.R"))

library(tidyverse)
library(bbsBayes2)   # load_bbs_data()
library(sf)
library(lwgeom)      # st_startpoint(), st_endpoint()

route_dir <- p_raw("bbsroutes_shp")     # raw route lines (bbsrtsl020.shp)
route_out <- p_clean("bbsroutes_shp")   # derived start points / filtered lines
dir.create(route_out, showWarnings = FALSE)

# ---- 1. Load BBS data and species lookup ------------------------------------
dat <- load_bbs_data(level = "state", release = 2024, sample = FALSE, quiet = FALSE)

species_hash <- dat$species %>%
  mutate(ScientificName = paste0(genus, " ", species))
write.csv(species_hash, p_clean("BBS_species_hash.csv"), row.names = FALSE)

# ---- 2. Observers + organised BCR route-year counts -------------------------
observers_df <- dat$routes %>%
  mutate(route       = str_c(country_num, state_num, route, sep = "-"),
         observer_id = as.factor(obs_n)) %>%
  filter(rpid == 101) %>%
  select(year, route, route_name, observer_id)

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
# Saved so 03_bird_community.R does not need to re-download BBS
saveRDS(organized_data_bcr, p_clean("BBS_organized_counts.rds"))

# ---- 3. Route start-points (sf) ---------------------------------------------
routes_sf_bcr <- dat$routes %>%
  mutate(route = str_c(country_num, state_num, route, sep = "-")) %>%
  filter(route %in% unique(organized_data_bcr$route)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  dplyr::select(route, bcr) %>%
  distinct()
write_sf(routes_sf_bcr, file.path(route_out, "bbs_startpoints.shp"))

# ---- 4. Match BBS line routes to BCR routes -> lookups `a` and `b` ----------
bbsrtsl020 <- st_read(file.path(route_dir, "bbsrtsl020.shp")) %>%
  dplyr::select(RTENO, RTENAME) %>%
  st_transform(crs = TARGET_CRS)
n_distinct(st_drop_geometry(bbsrtsl020)[, 1:2]) == nrow(bbsrtsl020)   # QC

routes_sf_bcr <- st_read(file.path(route_out, "bbs_startpoints.shp")) %>%
  st_transform(crs = TARGET_CRS)

# For each line route, match whichever endpoint is nearer to a start point
start_pts <- st_startpoint(bbsrtsl020)
end_pts   <- st_endpoint(bbsrtsl020)

nearest_start <- st_nearest_feature(start_pts, routes_sf_bcr)
nearest_end   <- st_nearest_feature(end_pts,   routes_sf_bcr)

pt_start <- routes_sf_bcr[nearest_start, ]
pt_end   <- routes_sf_bcr[nearest_end, ]

dist_start <- st_distance(start_pts, pt_start, by_element = TRUE)
dist_end   <- st_distance(end_pts,   pt_end,   by_element = TRUE)

use_start <- dist_start <= dist_end
nearest_route <- ifelse(use_start, routes_sf_bcr$route[nearest_start],
                        routes_sf_bcr$route[nearest_end])
nearest_dist  <- ifelse(use_start, dist_start, dist_end)   # metres (units dropped)

result <- bbsrtsl020 %>%
  mutate(route = nearest_route, dist = nearest_dist)

# Closest match (< 1 km) per line, then one line per route
b <- result %>%
  st_drop_geometry() %>%
  group_by(RTENO, RTENAME) %>% arrange(dist) %>% slice(1) %>% ungroup() %>%
  filter(dist < 1000) %>%
  group_by(route) %>% arrange(dist) %>% slice(1) %>% ungroup()

# `a`: BBSRTSL020 <-> RTENO/RTENAME crosswalk (Earth Engine export)
# TODO(reviewer C3): `a` is later attached to the EH layers by ROW ORDER
# (cbind). Replace with a key join on BBSRTSL020 once the EH exports carry it.
export_table <- read.csv(p_gee("export_table.csv"))
n_distinct(export_table$BBSRTSL020)
a <- export_table %>% dplyr::select(RTENO, RTENAME, BBSRTSL020)

write_rds(a, p_clean("a.rds"))
write_rds(b, p_clean("b.rds"))

# ---- 5. Filter the route shapefiles to matched routes -----------------------
bbsrtsl020_filtered <- bbsrtsl020 %>%
  filter(paste(RTENO, RTENAME) %in% paste(b$RTENO, b$RTENAME))
routes_sf_bcr_filtered <- routes_sf_bcr %>%
  filter(route %in% b$route)

write_sf(bbsrtsl020_filtered,    file.path(route_out, "bbsrtsl020_filterd.shp"))
write_sf(routes_sf_bcr_filtered, file.path(route_out, "bbs_startpoints_filterd.shp"))
