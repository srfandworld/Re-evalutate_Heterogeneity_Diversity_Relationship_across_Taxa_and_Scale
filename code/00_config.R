# =============================================================================
# 00_config.R - project paths and shared constants
# Re-evaluating the Environmental Heterogeneity-Diversity Relationship
#
# Every script starts with `source(here::here("code", "00_config.R"))`.
# All file paths in the pipeline are built with the helpers below, so the
# project can be moved without editing individual scripts.
#
# Folder roles
#   data/raw/     inputs that are downloaded or exported, never written by R
#                 (FIA state CSVs, SNAPSHOT USA, BBS route shapefiles, Google
#                 Earth Engine exports of EH and climate/NPP buffers)
#   data/clean/   intermediate files written by the pipeline
#   output/       figures and tables reported in the manuscript and SI
#
# External sources (too large to ship; point these at your local copies)
#   HDR_EXT_DIR    root holding environmental_data/, Traits_data/,
#                  tl_2024_us_state/, na_cec_eco_l1/, na_cec_eco_l2/,
#                  SpeciesRange/
#   HDR_RANGE_DIR  BirdLife / IUCN range polygons (BOTW gpkg, MAMMALS shp)
# Set them in .Renviron, e.g.  HDR_EXT_DIR=E:/BiodiversityEmbedding
# =============================================================================

library(here)

# ---- Project folders --------------------------------------------------------
DIR_RAW   <- here("data", "raw")
DIR_CLEAN <- here("data", "clean")
DIR_OUT   <- here("output")
DIR_FIG   <- here("output", "figures")
DIR_TAB   <- here("output", "tables")

# ---- External sources (EDIT via .Renviron, not here) ------------------------
DIR_EXT   <- Sys.getenv("HDR_EXT_DIR",   "E:/BiodiversityEmbedding")
DIR_RANGE <- Sys.getenv("HDR_RANGE_DIR", "D:/NaturalRangePolygon")

# ---- Path helpers -----------------------------------------------------------
p_raw   <- function(...) file.path(DIR_RAW, ...)
p_gee   <- function(...) file.path(DIR_RAW, "gee_exports", ...)  # Earth Engine outputs
p_clean <- function(...) file.path(DIR_CLEAN, ...)
p_model <- function(...) file.path(DIR_CLEAN, "models", ...)     # fitted model objects
p_fig   <- function(...) file.path(DIR_FIG, ...)
p_tab   <- function(...) file.path(DIR_TAB, ...)
p_ext   <- function(...) file.path(DIR_EXT, ...)
p_env   <- function(...) file.path(DIR_EXT, "environmental_data", ...)
p_trait <- function(...) file.path(DIR_EXT, "Traits_data", ...)
p_range <- function(...) file.path(DIR_RANGE, ...)

for (d in c(DIR_RAW, p_gee(), DIR_CLEAN, p_model(), p_clean("FIA_gamma_abundance"),
            p_clean("niche_results"), DIR_FIG, DIR_TAB)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ---- Analysis constants -----------------------------------------------------
TARGET_CRS    <- 5070   # EPSG:5070, CONUS Albers Equal Area
SCALES_TREE   <- c("1km", "5km", "10km", "20km", "50km")
SCALES_BIRD   <- c("400m", "1km", "5km", "10km", "20km")
SCALES_MAMMAL <- c("1km", "5km", "10km", "20km", "50km")

# Final tree-model set (reviewer item C2): coverage-based richness (iNEXT qD).
TREE_MODEL_SET <- "iNEXT"

# NLCD dominant land-cover codes -> labels (used by all three community tables)
LC_LEVELS <- c(1, 2, 3, 4, 5, 7, 9, 81, 82)
LC_LABELS <- c("Open water", "Developed", "Barren land", "Forest",
               "Shrub/Scrub", "Grassland/Herbaceous", "Wetlands",
               "Pasture/Hay", "Cultivated crops")
