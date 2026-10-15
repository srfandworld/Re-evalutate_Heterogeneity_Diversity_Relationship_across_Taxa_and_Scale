# =============================================================================
# fun_niche_breadth.R - 3-D KDE niche breadth and range size helpers
#
# Source   : calculate_niche_3d()       EmbeddingBio/analysis/niche_breadth_usingRangeMap.R, lines 112-260
#            calculate_niche_3d_occ()   EmbeddingBio/analysis/niche_breadth_usingOccurence.R, lines 46-194
#            calculate_range_convexhull() same file, lines 272-305
#            Copied verbatim. Niche space = first 3 PCs of the 14 CONUS
#            environmental layers (05_niche_space_PCA.R); volumes are the 95%
#            (total) and 50% (core) KDE contours with a global bandwidth H.
# =============================================================================

library(dplyr)
library(sf)
library(terra)
library(ks)
library(units)

vars <- c("MAT","MAP","Tseas","Tmax","Tmin","Pmax","Pmin","Pseas",
          "SoilPH","SOC","SoilNitro","SoilMoist","Elev","Slope")

# ---- Range-map niche breadth (birds, mammals) -------------------------------# ── 4. Core niche breadth function ────────────────────────────────────────────
calculate_niche_3d <- function(species_name,
                               range_sf,          # tree bird or mammal sp_range_us sf object
                               env_stack_coarse,  # coarsened env raster
                               pca_model,         # saved prcomp object
                               sc_center,         # scaling center vector
                               sc_sd,             # scaling SD vector
                               H,                 # global bandwidth matrix
                               grid_size  = 50,   # 50^3 grid
                               prob_total = 0.95,
                               prob_core  = 0.50,
                               max_points = 10000) {
  
  cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("Processing:", species_name, "\n")
  
  # ── 4a. Extract species range polygon ───────────────────────────────────────
  sp_range <- range_sf %>% filter(species == species_name)
  
  if (nrow(sp_range) == 0) {
    cat("No range found for species\n")
    return(NULL)
  }
  
  # ── 4b. Crop and mask env_stack to species range ────────────────────────────
  # This is the key step: extract env values from ALL cells within range
  # directly mimicking "extracting values from all cells within breeding range"
  sp_vect    <- terra::vect(sp_range)
  env_crop   <- terra::crop(env_stack_coarse, sp_vect)
  env_masked <- terra::mask(env_crop, sp_vect)
  
  # ── 4c. Convert to dataframe, drop NAs ──────────────────────────────────────
  sp_env_df <- as.data.frame(env_masked, na.rm = TRUE)
  colnames(sp_env_df) <- vars
  
  n_cells <- nrow(sp_env_df)
  cat("Grid cells in range:", n_cells, "\n")
  
  if (n_cells < 50) {
    cat("Too few cells — skipping\n")
    return(list(species = species_name, vol_total = NA,
                vol_core = NA, vol_periphery = NA,
                n_cells = n_cells, status = "too_few_cells"))
  }
  
  # ── 4d. Scale using GLOBAL parameters ───────────────────────────────────────
  # Critical: use same center/SD as global PCA scaling
  # so species PC scores are in the same space as the global PCA
  sp_scaled <- scale(sp_env_df, center = sc_center, scale = sc_sd)
  
  # Remove any rows with NA introduced by scaling
  sp_scaled <- sp_scaled[complete.cases(sp_scaled), ]
  
  # ── 4e. Project onto global PC axes ─────────────────────────────────────────
  sp_pc <- predict(pca_model, newdata = sp_scaled)[, 1:3]
  sp_pc <- as.data.frame(sp_pc)
  colnames(sp_pc) <- c("PC1","PC2","PC3")
  
  # ── 4f. Downsample if needed ─────────────────────────────────────────────────
  n_obs <- nrow(sp_pc)
  if (n_obs > max_points) {
    set.seed(42)
    sp_pc <- sp_pc[sample(n_obs, max_points), ]
    cat("Downsampled to:", max_points, "\n")
  } else {
    cat("PC scores extracted:", n_obs, "\n")
  }
  
  tryCatch({
    
    # ── 4g. Set grid limits to species' own extreme PC values ────────────────
    # Directly mimics: "X/Y limits set at most extreme PC values available"
    # Using species-specific limits focuses grid resolution on occupied space
    xmin <- min(sp_pc$PC1); xmax <- max(sp_pc$PC1)
    ymin <- min(sp_pc$PC2); ymax <- max(sp_pc$PC2)
    zmin <- min(sp_pc$PC3); zmax <- max(sp_pc$PC3)
    
    # 5% buffer so edge points are not cut off by grid boundary
    buf  <- 0.05
    xbuf <- (xmax - xmin) * buf; ybuf <- (ymax - ymin) * buf
    zbuf <- (zmax - zmin) * buf
    
    gridx <- seq(xmin - xbuf, xmax + xbuf, length.out = grid_size)
    gridy <- seq(ymin - ybuf, ymax + ybuf, length.out = grid_size)
    gridz <- seq(zmin - zbuf, zmax + zbuf, length.out = grid_size)
    
    # ── 4h. Fit 3D KDE ──────────────────────────────────────────────────────
    eval_grid <- as.matrix(expand.grid(gridx, gridy, gridz))
    
    fit <- ks::kde(
      x           = as.matrix(sp_pc),
      H           = H,               # global bandwidth matrix
      eval.points = eval_grid
    )
    
    # Reshape into 3D array
    density_3d <- array(fit$estimate,
                        dim = c(grid_size, grid_size, grid_size))
    
    # ── 4i. Voxel volume in PC space ─────────────────────────────────────────
    dx        <- gridx[2] - gridx[1]
    dy        <- gridy[2] - gridy[1]
    dz        <- gridz[2] - gridz[1]
    voxel_vol <- dx * dy * dz
    
    # ── 4j. Find probability contour thresholds ──────────────────────────────
    # Sort densities, find cutoff density that captures X% of total mass
    # This directly mirrors the 95%/50% contour logic from the paper
    dens_vec     <- as.vector(density_3d)
    total_mass   <- sum(dens_vec) * voxel_vol
    sorted_dens  <- sort(dens_vec, decreasing = TRUE)
    cum_prob     <- cumsum(sorted_dens) * voxel_vol / total_mass
    
    thresh_total <- sorted_dens[which(cum_prob >= prob_total)[1]]
    thresh_core  <- sorted_dens[which(cum_prob >= prob_core)[1]]
    
    # ── 4k. Calculate volumes ─────────────────────────────────────────────────
    vol_total     <- sum(density_3d >= thresh_total) * voxel_vol
    vol_core      <- sum(density_3d >= thresh_core)  * voxel_vol
    vol_periphery <- vol_total - vol_core
    
    cat("Total niche (95%):", round(vol_total,     4), "\n")
    cat("Core niche  (50%):", round(vol_core,      4), "\n")
    cat("Periphery        :", round(vol_periphery, 4), "\n")
    
    return(list(
      species       = species_name,
      vol_total     = vol_total,
      vol_core      = vol_core,
      vol_periphery = vol_periphery,
      n_cells       = n_cells,
      n_obs         = nrow(sp_pc),
      status        = "success"
    ))
    
  }, error = function(e) {
    cat("ERROR:", e$message, "\n")
    return(list(
      species       = species_name,
      vol_total     = NA,
      vol_core      = NA,
      vol_periphery = NA,
      n_cells       = n_cells,
      n_obs         = NA,
      error         = e$message,
      status        = "error"
    ))
  })
}

# ---- Occurrence-based niche breadth (trees, FIA plots) ---------------------
calculate_niche_3d_occ <- function(species_name,
                                   occ_data,            # dt0_filtered
                                   env_stack_coarse,    # coarsened env raster
                                   pca_model,           # saved prcomp object
                                   sc_center,           # scaling center vector
                                   sc_sd,               # scaling SD vector
                                   H,                   # global bandwidth matrix
                                   grid_size  = 50,
                                   prob_total = 0.95,
                                   prob_core  = 0.50,
                                   max_points = 10000) {
  
  cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("Processing:", species_name, "\n")
  
  # ── 3a. Filter occurrence points for this species ───────────────────────────
  sp_occ <- occ_data %>%
    filter(species == species_name)
  
  n_occ <- nrow(sp_occ)
  cat("Occurrence points:", n_occ, "\n")
  
  # ── 3b. Extract env values at occurrence points ─────────────────────────────
  # Convert to SpatVector for terra extraction
  sp_vect <- terra::vect(
    sp_occ,
    geom  = c("lon", "lat"),
    crs   = "EPSG:4326"
  )
  
  # Extract env values at each point — this directly replaces the
  # "extract all cells within range" step from the range-based method
  env_vals <- terra::extract(env_stack_coarse, sp_vect, ID = FALSE)
  colnames(env_vals) <- vars
  
  # Attach species and coordinates
  sp_env_df <- bind_cols(
    sp_occ %>% dplyr::select(species, lon, lat),
    env_vals
  ) %>%
    filter(complete.cases(.))                            # drop points outside raster
  
  n_cells <- nrow(sp_env_df)
  cat("Points with env data:", n_cells, "\n")
  
  if (n_cells < 50) {
    cat("Too few points with env data — skipping\n")
    return(list(
      species   = species_name,
      vol_total = NA, vol_core = NA, vol_periphery = NA,
      n_occ     = n_occ, n_cells = n_cells,
      status    = "too_few_cells"
    ))
  }
  
  # ── 3c. Scale using GLOBAL parameters ───────────────────────────────────────
  sp_scaled <- scale(
    sp_env_df %>% dplyr::select(all_of(vars)),
    center = sc_center,
    scale  = sc_sd
  )
  sp_scaled <- sp_scaled[complete.cases(sp_scaled), ]
  
  # ── 3d. Project onto global PC axes ─────────────────────────────────────────
  sp_pc <- predict(pca_model, newdata = sp_scaled)[, 1:3]
  sp_pc <- as.data.frame(sp_pc)
  colnames(sp_pc) <- c("PC1", "PC2", "PC3")
  
  # ── 3e. Downsample if needed ─────────────────────────────────────────────────
  n_obs <- nrow(sp_pc)
  if (n_obs > max_points) {
    set.seed(42)
    sp_pc <- sp_pc[sample(n_obs, max_points), ]
    cat("Downsampled to:", max_points, "\n")
  } else {
    cat("PC scores used:", n_obs, "\n")
  }
  
  tryCatch({
    
    # ── 3f. Species-specific grid limits ────────────────────────────────────
    xmin <- min(sp_pc$PC1); xmax <- max(sp_pc$PC1)
    ymin <- min(sp_pc$PC2); ymax <- max(sp_pc$PC2)
    zmin <- min(sp_pc$PC3); zmax <- max(sp_pc$PC3)
    
    buf  <- 0.05
    xbuf <- (xmax - xmin) * buf
    ybuf <- (ymax - ymin) * buf
    zbuf <- (zmax - zmin) * buf
    
    gridx <- seq(xmin - xbuf, xmax + xbuf, length.out = grid_size)
    gridy <- seq(ymin - ybuf, ymax + ybuf, length.out = grid_size)
    gridz <- seq(zmin - zbuf, zmax + zbuf, length.out = grid_size)
    
    # ── 3g. Fit 3D KDE ───────────────────────────────────────────────────────
    eval_grid <- as.matrix(expand.grid(gridx, gridy, gridz))
    
    fit <- ks::kde(
      x           = as.matrix(sp_pc),
      H           = H,
      eval.points = eval_grid
    )
    
    density_3d <- array(fit$estimate,
                        dim = c(grid_size, grid_size, grid_size))
    
    # ── 3h. Voxel volume ─────────────────────────────────────────────────────
    voxel_vol <- (gridx[2]-gridx[1]) * (gridy[2]-gridy[1]) * (gridz[2]-gridz[1])
    
    # ── 3i. Probability contour thresholds ───────────────────────────────────
    dens_vec    <- as.vector(density_3d)
    total_mass  <- sum(dens_vec) * voxel_vol
    sorted_dens <- sort(dens_vec, decreasing = TRUE)
    cum_prob    <- cumsum(sorted_dens) * voxel_vol / total_mass
    
    thresh_total <- sorted_dens[which(cum_prob >= prob_total)[1]]
    thresh_core  <- sorted_dens[which(cum_prob >= prob_core)[1]]
    
    # ── 3j. Calculate volumes ─────────────────────────────────────────────────
    vol_total     <- sum(density_3d >= thresh_total) * voxel_vol
    vol_core      <- sum(density_3d >= thresh_core)  * voxel_vol
    vol_periphery <- vol_total - vol_core
    
    cat("Total niche (95%):", round(vol_total,     4), "\n")
    cat("Core niche  (50%):", round(vol_core,      4), "\n")
    cat("Periphery        :", round(vol_periphery, 4), "\n")
    
    return(list(
      species       = species_name,
      vol_total     = vol_total,
      vol_core      = vol_core,
      vol_periphery = vol_periphery,
      n_occ         = n_occ,
      n_cells       = n_cells,
      status        = "success"
    ))
    
  }, error = function(e) {
    cat("ERROR:", e$message, "\n")
    return(list(
      species       = species_name,
      vol_total     = NA, vol_core = NA, vol_periphery = NA,
      n_occ         = n_occ, n_cells = n_cells,
      error         = e$message,
      status        = "error"
    ))
  })
}

# ---- Range size from occurrences: convex hull (trees) -----------------------
calculate_range_convexhull <- function(species_name, occ_data) {
  
  sp_occ <- occ_data %>% 
    filter(species == species_name) %>%
    dplyr::select(lon, lat) %>% 
    distinct()
  
  if (nrow(sp_occ) < 3) {
    return(data.frame(species = species_name, range_km2 = NA,
                      status = "too_few_points"))
  }
  
  tryCatch({
    # Build convex hull directly in sf
    range_sf <- sp_occ %>%
      st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
      st_union() %>%
      st_convex_hull() %>%
      st_transform(5070)                        # project for area
    
    range_km2 <- st_area(range_sf) %>% 
      set_units(km^2) %>% 
      as.numeric()
    
    return(data.frame(species = species_name, range_km2 = range_km2,
                      n_points = nrow(sp_occ), method = "convexhull",
                      status = "success"))
    
  }, error = function(e) {
    return(data.frame(species = species_name, range_km2 = NA,
                      n_points = nrow(sp_occ), method = "convexhull",
                      status = paste("error:", e$message)))
  })
}
