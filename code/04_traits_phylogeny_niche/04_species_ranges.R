# =============================================================================
# 04_species_ranges.R
# Species range polygons clipped to the contiguous US (input to range size and
# range-map niche breadth)
#
# Source   : trees   - EmbeddingBio/analysis/FIA tree traits.R, lines 45-107
#            birds   - EmbeddingBio/analysis/BBS_traits_phylogeny.R, lines 417-456
#            mammals - EmbeddingBio/analysis/mammal_traits_phylogeny.R, lines 265-300
#            Interactive map checks were dropped; the steps are otherwise unchanged.
# Inputs   : data/clean/fia_sp.csv
#            HDR_EXT_DIR/SpeciesRange/BIEN_species_range/ (downloaded below)
#            HDR_EXT_DIR/SpeciesRange/IUCN_Trees{1,2}_Range.shp
#            HDR_RANGE_DIR/range_polygon_aves/BOTW_2024_2.gpkg   (BirdLife)
#            HDR_RANGE_DIR/range_polygon_mammalia/MAMMALS_TERRESTRIAL_ONLY.shp (IUCN)
#            HDR_EXT_DIR/tl_2024_us_state/tl_2024_us_state.shp
# Outputs  : data/clean/fia_species_range.shp, fia_species_range_us.shp
#            data/clean/bird_species_range_us.shp
#            data/clean/mammal_species_range_us.shp
# =============================================================================

source(here::here("code", "00_config.R"))

library(sf)
library(dplyr)
library(stringr)
library(BIEN)

us_boundary <- st_read(p_ext("tl_2024_us_state", "tl_2024_us_state.shp")) %>%
  st_union() %>%
  st_transform(4326) %>%
  st_make_valid()

# =============================================================================
# 1. Trees: BIEN ranges, gaps filled with IUCN
# =============================================================================
fia_sp <- read.csv(p_clean("fia_sp.csv"))

bien_dir <- p_ext("SpeciesRange", "BIEN_species_range")
BIEN_ranges_species_bulk(species = unique(fia_sp$SCIENTIFIC_NAME),
                         directory = bien_dir,
                         use_parallel = TRUE)
files <- list.files(path = file.path(bien_dir, "1"), pattern = ".shp", full.names = TRUE)

fia_sp_range <- lapply(files, st_read) %>%
  bind_rows() %>%
  select(species) %>%
  mutate(species = str_replace(species, "_", " "))   # 550 species

IUCN_range <- rbind(st_read(p_ext("SpeciesRange", "IUCN_Trees1_Range.shp")),
                    st_read(p_ext("SpeciesRange", "IUCN_Trees2_Range.shp"))) %>%
  select(species = sci_name) %>%
  filter(species %in% unique(fia_sp$SCIENTIFIC_NAME) &
           !(species %in% unique(fia_sp_range$species)))   # 47 species

fia_sp_range <- bind_rows(fia_sp_range, IUCN_range)   # 597 / 868 species
fia_sp_range <- st_make_valid(fia_sp_range)

inter_list <- st_intersects(fia_sp_range, us_boundary, sparse = TRUE)
fia_sp_range$intersects_us <- lengths(inter_list) > 0
st_write(fia_sp_range, p_clean("fia_species_range.shp"), append = FALSE)
fia_sp_range <- st_read(p_clean("fia_species_range.shp"))   # field becomes `intrsc_`

us_union <- st_union(us_boundary) %>% st_as_sf()
fia_sp_range_us <- fia_sp_range[fia_sp_range$intrsc_ == TRUE, ]
fia_sp_range_us <- st_intersection(fia_sp_range_us, us_union)
st_write(fia_sp_range_us, p_clean("fia_species_range_us.shp"), append = FALSE)

# =============================================================================
# 2. Birds: BirdLife BOTW 2024 (resident/breeding/non-breeding, native)
# =============================================================================
bird_sp_range <- st_read(
  p_range("range_polygon_aves", "BOTW_2024_2.gpkg"),
  layer = "all_species",
  wkt_filter = st_as_text(us_boundary)
)
bird_sp_range <- bird_sp_range %>%
  filter(presence == 1 & origin %in% 1:3 & seasonal %in% 1:3) %>%
  dplyr::select(species = sci_name)

bird_sp_range <- st_cast(bird_sp_range, "MULTIPOLYGON")
bird_sp_range <- st_make_valid(bird_sp_range)

sf::sf_use_s2(FALSE)
inter_list <- st_intersects(bird_sp_range, us_boundary, sparse = TRUE)
bird_sp_range$intersects_us <- lengths(inter_list) > 0

bird_sp_range_us <- bird_sp_range[bird_sp_range$intersects_us == TRUE, ] %>% st_transform(TARGET_CRS)
bird_sp_range_us <- st_intersection(bird_sp_range_us, st_transform(us_boundary, TARGET_CRS))
bird_sp_range_us <- bird_sp_range_us %>%
  group_by(species) %>%
  summarise(geom = st_union(geom), .groups = "drop")
st_write(bird_sp_range_us, p_clean("bird_species_range_us.shp"), append = FALSE)

# =============================================================================
# 3. Mammals: IUCN terrestrial mammals
# =============================================================================
mammal_sp_range <- st_read(p_range("range_polygon_mammalia", "MAMMALS_TERRESTRIAL_ONLY.shp")) %>%
  filter(presence == 1 & origin %in% 1:3 & seasonal %in% 1:3) %>%
  dplyr::select(species = sci_name, order = order_)
mammal_sp_range <- st_make_valid(mammal_sp_range)

sf::sf_use_s2(FALSE)
inter_list <- st_intersects(mammal_sp_range, us_boundary, sparse = TRUE)
mammal_sp_range$intersects_us <- lengths(inter_list) > 0

mammal_sp_range_us <- mammal_sp_range[mammal_sp_range$intersects_us == TRUE, ]
mammal_sp_range_us <- st_intersection(mammal_sp_range_us, us_boundary)
mammal_sp_range_us <- st_cast(mammal_sp_range_us, "MULTIPOLYGON")
mammal_sp_range_us <- mammal_sp_range_us %>%
  group_by(species, order) %>%
  summarise(geom = sf::st_union(geom), .groups = "drop")
st_write(mammal_sp_range_us, p_clean("mammal_species_range_us.shp"), append = FALSE)
