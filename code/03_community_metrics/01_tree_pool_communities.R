# =============================================================================
# 01_tree_pool_communities.R
# Pool neighbouring FIA plots into buffer communities; richness/MSA per plot;
# assemble FIA_table and the scaled per-scale tables
#
# Source   : Code/Tree/01_FIA_data_download_cleanning.R, sections 2-5.
#            build_comm_abundance() now matches the version in
#            EmbeddingBio/analysis/FIA gamma diversity.R, which also stores
#            n_plots_species. The Code/ copy had dropped that column, which is
#            why the reviewer's iNEXT step returned an empty table (C5).
# Inputs   : data/clean/tpa_all5.rds, data/clean/FIA_env.csv
# Outputs  : data/clean/plot_pts.shp
#            data/clean/FIA_gamma_abundance/comm_abundance_r{1,5,10,20,50,110}km.rds
#            data/clean/FIA_table.csv
#            data/clean/FIA_table_{1,5,10,20,50}km.rds   (scaled covariates)
# Next     : 02_tree_richness_iNEXT.R
#
# TODO(reviewer C4): covariates are standardised here over ALL plots, but the
#   iNEXT step keeps only a subset per scale and Figure 2 back-transforms
#   turning points with the subset mean/SD. Pick one and use it consistently.
# TODO(reviewer C5): `NA_L1` is selected below, but FIA_env.csv carries
#   NA_L1CODE / NA_L1NAME (CEC level-1 ecoregion). Check which column name the
#   saved FIA_table_*.rds actually contain.
# NOTE: pooling at 50/110 km is memory- and time-intensive (hours).
# =============================================================================

source(here::here("code", "00_config.R"))

library(sf)
library(data.table)   # build_comm_abundance() relies on it
library(tidyverse)

tpa_all5 <- readRDS(p_clean("tpa_all5.rds"))
fia_sp   <- read.csv(p_clean("fia_sp.csv"))
FIA_env  <- read.csv(p_clean("FIA_env.csv"))   # 02_environment/02_environmental_covariates.R

# =============================================================================
# 2. Pool neighbouring plots into communities (radii 1, 5, 10, 20, 50, 110 km)
# =============================================================================
dt0 <- tpa_all5 %>%
  dplyr::select(pltID, SCIENTIFIC_NAME, BAA, geometry) %>%
  filter(!is.na(pltID), !is.na(SCIENTIFIC_NAME), !is.na(BAA))

plot_pts <- dt0 %>%
  distinct(pltID, geometry) %>%
  st_as_sf()
st_write(plot_pts, p_clean("plot_pts.shp"), append = FALSE)
plot_pts <- st_read(p_clean("plot_pts.shp"))

plot_sp <- dt0 %>%
  st_drop_geometry() %>%
  group_by(pltID, SCIENTIFIC_NAME) %>%
  summarise(BAA = sum(BAA, na.rm = TRUE), .groups = "drop") %>%
  as.data.table()

# Add an integer index to plot_pts so neighbour row numbers map back to pltID
plot_pts$plot_idx <- seq_len(nrow(plot_pts))
idx_to_pltID <- plot_pts$pltID

# Core function: find neighbours within a radius and pool them into a community
# abundance table (processed in batches to control memory).
build_comm_abundance <- function(plot_pts, plot_sp, radii_km,
                                 batch_size = 2000L,
                                 out_dir = NULL) {
  stopifnot(inherits(plot_pts, "sf"))
  plot_sp <- as.data.table(plot_sp)
  
  all_pts <- plot_pts
  n_all   <- nrow(all_pts)
  
  # Precompute pltID -> plot_idx mapping (for joins)
  map_dt <- data.table(pltID = all_pts$pltID, plot_idx = all_pts$plot_idx)
  setkey(map_dt, pltID)
  setkey(plot_sp, pltID)
  
  results_list <- list()
  
  for (r_km in radii_km) {
    message("==> radius: ", r_km, " km")
    dist_m <- r_km * 1000
    
    # Process focal plots in batches to avoid one huge neighbour list
    batch_starts <- seq(1L, n_all, by = batch_size)
    
    # One accumulator per radius (written to disk to avoid memory blow-up)
    out_accum <- list()
    kk <- 0L
    
    for (b0 in batch_starts) {
      b1    <- min(b0 + batch_size - 1L, n_all)
      focal <- all_pts[b0:b1, ]
      
      # Neighbour search: for each point in the focal batch, find every point
      # in all_pts within dist_m
      nb_list <- st_is_within_distance(focal, all_pts, dist = dist_m)
      
      # Neighbour count per focal plot (sampling effort)
      effort_dt <- data.table(
        focal_idx = focal$plot_idx,
        effort    = sapply(nb_list, length)
      )
      effort_dt[, focal_pltID := idx_to_pltID[focal_idx]]
      
      # Convert the neighbour list to an edge table: focal_idx -> neigh_idx
      # (nb_list[[i]] holds row numbers into all_pts)
      edges <- rbindlist(lapply(seq_along(nb_list), function(i) {
        ni <- nb_list[[i]]
        if (length(ni) == 0) return(NULL)
        data.table(
          focal_idx = focal$plot_idx[i],
          neigh_idx = all_pts$plot_idx[ni]
        )
      }))
      
      if (nrow(edges) == 0) next
      
      edges[, neigh_pltID := idx_to_pltID[neigh_idx]]
      edges[, focal_pltID := idx_to_pltID[focal_idx]]
      
      # Pull in each neighbour's per-species BAA, then pool by focal plot:
      # join edges (neigh_pltID) to plot_sp (pltID)
      setkey(edges, neigh_pltID)
      tmp <- plot_sp[
        edges,
        on = .(pltID = neigh_pltID),
        .(
          focal_pltID     = i.focal_pltID,
          neigh_pltID     = i.neigh_pltID,
          SCIENTIFIC_NAME = x.SCIENTIFIC_NAME,
          BAA             = x.BAA
        ),
        allow.cartesian = TRUE
      ]
      
      # Pool community abundance: for each focal plot and species, sum BAA
      # across all neighbours and count the plots the species occurs in
      # (n_plots_species = incidence frequency used by iNEXT in the next script)
      comm <- tmp[
        ,
        .(
          comm_BAA        = sum(BAA, na.rm = TRUE),
          n_plots_species = uniqueN(neigh_pltID)
        ),
        by = .(focal_pltID, SCIENTIFIC_NAME)
      ]
      
      # Attach effort
      setkey(comm, focal_pltID)
      setkey(effort_dt, focal_pltID)
      comm <- effort_dt[comm, on = .(focal_pltID)]
      comm[, c("focal_idx") := NULL]
      
      comm[, radius_km := r_km]
      
      kk <- kk + 1L
      out_accum[[kk]] <- comm
      
      message("  batch ", b0, "-", b1, " done. edges=", nrow(edges),
              " rows_out=", nrow(comm))
    }
    
    comm_r <- rbindlist(out_accum, use.names = TRUE, fill = TRUE)
    
    # Write to disk (strongly recommended for large radii)
    if (!is.null(out_dir)) {
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
      fout <- file.path(out_dir, paste0("comm_abundance_r", r_km, "km.rds"))
      saveRDS(comm_r, fout)   # or fwrite(comm_r, gsub(".rds", ".csv", fout))
      message("  wrote: ", fout)
    } else {
      results_list[[as.character(r_km)]] <- comm_r
    }
  }
  
  if (is.null(out_dir)) {
    return(rbindlist(results_list, use.names = TRUE, fill = TRUE))
  } else {
    return(invisible(TRUE))
  }
}

# Run for each radius
for (i in c(1, 5, 10, 20, 50, 110)) {
  build_comm_abundance(plot_pts, plot_sp, radii_km = c(i),
                       batch_size = 1500L,
                       out_dir = p_clean("FIA_gamma_abundance"))
}

# =============================================================================
# 3. Compute gamma diversity per focal plot and radius
# =============================================================================
# Read one per-radius community file and compute diversity metrics.
calc_gamma_from_rds <- function(rds_path) {   # reads .rds (not parquet)
  ds <- readRDS(rds_path)
  message("Read community file successfully; computing gamma diversity.")
  
  ds <- ds %>% filter(SCIENTIFIC_NAME %in% fia_sp$SCIENTIFIC_NAME)
  message("Filtered to identified species successfully.")
  
  # Normalise each species by its maximum pooled abundance, then summarise per
  # focal plot: richness, mean relative abundance, and mean species abundance
  # (MSA), each with and without an effort correction.
  ds %>%
    group_by(SCIENTIFIC_NAME) %>%
    mutate(max_abundance_species           = max(comm_BAA),
           max_abundance_species_corEffort = max(comm_BAA / effort)) %>%
    ungroup() %>%
    group_by(focal_pltID) %>%
    summarise(
      richness            = n_distinct(SCIENTIFIC_NAME),
      abundance           = sum(comm_BAA / max_abundance_species) / n_distinct(SCIENTIFIC_NAME),
      abundance_corEffort = sum((comm_BAA / effort) / max_abundance_species_corEffort) / n_distinct(SCIENTIFIC_NAME),
      MSA                 = sum(comm_BAA) / n_distinct(SCIENTIFIC_NAME),
      MSA_new             = sum(comm_BAA / effort) / n_distinct(SCIENTIFIC_NAME),
      effort              = first(effort),   # effort is constant within a focal_pltID
      .groups = "drop"
    )
}

buffer_names <- c("1km", "5km", "10km", "20km", "50km", "110km")
radii        <- c(1, 5, 10, 20, 50, 110)

FIA_gamma_list <- lapply(seq_along(radii), function(i) {
  r <- radii[i]
  f <- file.path(p_clean("FIA_gamma_abundance"),
                 paste0("comm_abundance_r", r, "km.rds"))
  calc_gamma_from_rds(f) %>%
    rename(
      pltID = focal_pltID,
      !!paste0("richness_",      buffer_names[i]) := richness,
      !!paste0("abundance_",     buffer_names[i]) := abundance,
      !!paste0("abundance_new_", buffer_names[i]) := abundance_corEffort,
      !!paste0("MSA_",           buffer_names[i]) := MSA,
      !!paste0("MSA_new_",       buffer_names[i]) := MSA_new,
      !!paste0("effort_",        buffer_names[i]) := effort
    )
})

# =============================================================================
# 4. Assemble FIA_table (gamma metrics + environment + coordinates)
# =============================================================================
coords <- st_coordinates(plot_pts)
plot_pts_xy <- plot_pts %>%
  st_drop_geometry() %>%
  mutate(X = coords[, 1], Y = coords[, 2]) %>%
  dplyr::select(pltID, X, Y)

FIA_table <- FIA_gamma_list %>%
  reduce(left_join, by = "pltID") %>%
  left_join(FIA_env,     by = "pltID") %>%
  left_join(plot_pts_xy, by = "pltID")

# Drop plots with negative mean elevation (artefacts)
hist(FIA_table$elev_mean_10km)   # quick sanity check
FIA_table <- FIA_table %>%
  filter(if_all(contains("elev_mean"), ~ . >= 0))


write.csv(FIA_table, p_clean("FIA_table.csv"), row.names = FALSE)
FIA_table <- read.csv(p_clean("FIA_table.csv"))

# Dominant land-cover codes -> labels
FIA_table <- FIA_table %>%
  mutate(across(contains("dominant_lc"), ~ factor(
    .x,
    levels = c(1, 2, 3, 4, 5, 7, 9, 81, 82),
    labels = c("Open water", "Developed", "Barren land", "Forest",
               "Shrub/Scrub", "Grassland/Herbaceous", "Wetlands",
               "Pasture/Hay", "Cultivated crops")
  )))

# =============================================================================
# 5. Per-scale scaled model tables (1, 5, 10, 20, 50 km)
# =============================================================================
FIA_table_50km <- FIA_table %>%
  dplyr::select(pltID, contains("_50km"), X, Y, NA_L1) %>%
  mutate(EH_50km_t        = scale(EH_50km)[, 1],
         soilph_50km_t    = scale(log(soilph_50km))[, 1],
         soilcec_50km_t   = scale(log(soilcec_50km))[, 1],
         elev_mean_50km_t = scale(sqrt(elev_mean_50km))[, 1],
         MAP_50km_t       = scale(sqrt(MAP_50km))[, 1],
         MAT_50km_t       = scale(MAT_50km)[, 1],
         X_t              = scale(X)[, 1],
         Y_t              = scale(Y)[, 1])

FIA_table_20km <- FIA_table %>%
  dplyr::select(pltID, contains("_20km"), X, Y, NA_L1) %>%
  mutate(EH_20km_t        = scale(EH_20km)[, 1],
         soilph_20km_t    = scale(log(soilph_20km))[, 1],
         soilcec_20km_t   = scale(log(soilcec_20km))[, 1],
         elev_mean_20km_t = scale(sqrt(elev_mean_20km))[, 1],
         MAP_20km_t       = scale(sqrt(MAP_20km))[, 1],
         MAT_20km_t       = scale(MAT_20km)[, 1],
         X_t              = scale(X)[, 1],
         Y_t              = scale(Y)[, 1])

FIA_table_10km <- FIA_table %>%
  dplyr::select(pltID, contains("_10km"), X, Y, NA_L1) %>%
  mutate(EH_10km_t        = scale(EH_10km)[, 1],
         soilph_10km_t    = scale(log(soilph_10km))[, 1],
         soilcec_10km_t   = scale(log(soilcec_10km))[, 1],
         elev_mean_10km_t = scale(sqrt(elev_mean_10km))[, 1],
         MAP_10km_t       = scale(sqrt(MAP_10km))[, 1],
         MAT_10km_t       = scale(MAT_10km)[, 1],
         X_t              = scale(X)[, 1],
         Y_t              = scale(Y)[, 1])

FIA_table_5km <- FIA_table %>%
  dplyr::select(pltID, contains("_5km"), X, Y, NA_L1) %>%
  mutate(EH_5km_t        = scale(EH_5km)[, 1],
         soilph_5km_t    = scale(log(soilph_5km))[, 1],
         soilcec_5km_t   = scale(log(soilcec_5km))[, 1],
         elev_mean_5km_t = scale(sqrt(elev_mean_5km))[, 1],
         MAP_5km_t       = scale(sqrt(MAP_5km))[, 1],
         MAT_5km_t       = scale(MAT_5km)[, 1],
         X_t             = scale(X)[, 1],
         Y_t             = scale(Y)[, 1])

FIA_table_1km <- FIA_table %>%
  dplyr::select(pltID, contains("_1km"), X, Y, NA_L1) %>%
  mutate(EH_1km_t        = scale(EH_1km)[, 1],
         soilph_1km_t    = scale(log(soilph_1km))[, 1],
         soilcec_1km_t   = scale(log(soilcec_1km))[, 1],
         elev_mean_1km_t = scale(sqrt(elev_mean_1km))[, 1],
         MAP_1km_t       = scale(sqrt(MAP_1km))[, 1],
         MAT_1km_t       = scale(MAT_1km)[, 1],
         X_t             = scale(X)[, 1],
         Y_t             = scale(Y)[, 1])

# ---- Write out the five per-scale tables ------------------------------------
FIA_tables_by_scale <- list(
  "1km"  = FIA_table_1km,
  "5km"  = FIA_table_5km,
  "10km" = FIA_table_10km,
  "20km" = FIA_table_20km,
  "50km" = FIA_table_50km
)

iwalk(FIA_tables_by_scale, function(tbl, scale) {
  # FIX(reorg): dropped `row.names = FALSE`, which saveRDS() rejects (reviewer C5)
  saveRDS(tbl, p_clean(paste0("FIA_table_", scale, ".rds")))
})
