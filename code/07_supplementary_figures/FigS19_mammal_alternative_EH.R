# =============================================================================
# FigS19_mammal_alternative_EH.R
# Fig. S19 (mammals): same layout as Fig. S17 for camera-trap arrays (main
# scale 5 km). SAR models (spautolm, 8-NN weights on array centroids); response
# richness_Estimator / log(MSA) - same specification as 05_HDR_models/01.
#
# Source   : EmbeddingBio/analysis/sfigure17_physicallyVar_mammal.R (paths changed). Added so
#            it runs in a fresh session (the original relied on the workspace):
#              - libraries loaded up front; mammal_table read from data/clean
#              - QC hist() lines commented out (hist(dat$EH_10km_t) errors:
#                mammals use 5 km columns)
#              - map extent / graticule block (copied from FigS17) - the
#                original only had a placeholder comment
#              - usage block switched on (was if (FALSE))
# Inputs   : data/clean/mammal_table.csv, camera_array_centroids.shp
#            data/raw/gee_exports/mammal_EH_delta/, mammal_A57_EH/
# Outputs  : output/figures/FigS19_mammal_alternative_EH.png
# NOTE     : array-centroid weights use lon/lat coordinates, as in
#            05_HDR_models/01 (only the tree weights were moved to EPSG:5070).
# =============================================================================

source(here::here("code", "00_config.R"))

library(sf)
library(dplyr)
library(tidyr)
library(stringr)
library(glue)
library(ggplot2)
library(spatialreg)   # spautolm()
library(spdep)        # knearneigh(), knn2nb(), nb2listw()
mammal_table <- read.csv(p_clean("mammal_table.csv"))

#using physically interepretable variable

#mean-independence heterogeneity measures
mammal_5km_eledelta <- read.csv(p_gee("mammal_EH_delta/mammal_Elevation_EHdelta_5km.csv"))
mammal_5km_eledelta <- mammal_5km_eledelta %>% 
  dplyr::select(Camera_Trap_Array,elev_range,EH_delta)

#A57 EH
mammal_5km_A57 <- st_read(p_gee("mammal_A57_EH/mammal_A57_EH_5km.shp"))
mammal_5km_A57 <- mammal_5km_A57 %>% 
  st_drop_geometry() %>% 
  dplyr::select(Camera_Trap_Array = Camera_Tra,A57_sd)

#merge two table together
merged_table <- left_join(mammal_5km_eledelta,mammal_5km_A57,by="Camera_Trap_Array")

#data
dat <- mammal_table %>%
  dplyr::select(year_CTA,Year,Camera_Trap_Array,
                richness_Estimator,MSA,
                contains("_5km")) %>%  
  left_join(merged_table, by = c("Camera_Trap_Array")) %>%
  mutate(EH_5km_t = scale(EH_5km)[,1],
         EH_elevrange_t = scale(log(elev_range))[,1],
         EH_delta_t = scale(log(EH_delta))[,1],
         A57_sd_t = scale(A57_sd)[,1],
         elev_5km_t = scale(sqrt(elev_5km))[,1],
         PopDen_5km_t = scale(log1p(PopDen_5km))[,1],
         NPP_5km_t = scale(NPP_5km)[,1],
         Prec_5km_t = scale(sqrt(Prec_5km))[,1],
         Temp_5km_t = scale(Temp_5km)[,1]) %>%
  group_by(Camera_Trap_Array) %>% 
  summarise(across(c(richness_Estimator,MSA, contains("_t")), 
                   ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") %>% 
  drop_na(Camera_Trap_Array)
names(dat)


#histograms of response variables
# hist(dat$EH_5km_t)      # QC only (original referenced EH_10km_t, absent for mammals)
# hist(dat$EH_elevrange_t)
# hist(dat$EH_delta_t)
# hist(dat$A57_sd_t)

#lw
camera_array_centroids <- read_sf(p_clean("camera_array_centroids.shp")) %>% 
  rename(Camera_Trap_Array=Cmr_T_A)

resid_plots <-  dat %>% 
  dplyr::select(Camera_Trap_Array) %>% 
  left_join(camera_array_centroids,by="Camera_Trap_Array") %>% 
  st_as_sf()


#k nearest neighbors
coords <- st_coordinates(resid_plots)
nb <- knn2nb(knearneigh(coords, k = 8))
lw <- nb2listw(nb, style = "W")


scale <- "5km"
#richness-EH
form <- as.formula(glue("richness_Estimator~ EH_{scale}_t +I(EH_{scale}_t^2) +elev_{scale}_t+ PopDen_{scale}_t + NPP_{scale}_t + Prec_{scale}_t + Temp_{scale}_t  "))
#MSA-EH
form <- as.formula(glue("log(MSA)~ EH_{scale}_t +I(EH_{scale}_t^2) +elev_{scale}_t+ PopDen_{scale}_t + NPP_{scale}_t + Prec_{scale}_t + Temp_{scale}_t  "))

final_model <- spautolm(formula = form,
                      data = dat,
                      listw = lw,   #go to original files to get this
                      family = "SAR")


summary(final_model)
AIC(final_model)



#Draw the elev_range, mean-independent elev_range, EH_64d, EH_1d
states <- sf::st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))

a <- camera_array_centroids%>% 
  left_join(dat %>% 
              group_by(Camera_Trap_Array) %>% 
              summarise(richness = mean(richness_Estimator,na.rm=T), 
                        MSA = mean(MSA, na.rm = T)),
            by="Camera_Trap_Array")
conus_bbox <- st_bbox(c(
  xmin = -125,
  xmax = -66,
  ymin = 24,
  ymax = 50
), crs = st_crs(a))
a <- st_crop(a, conus_bbox) %>% st_transform(5070)

# Map extent / graticule / axis labels (added; copied from FigS17)
###map extent and graticule###
{
  # ── 1. Map extent ─────────────────────────────────────────────────────────────
  xlim_alb <- c(-2400000, 2300000)
  ylim_alb <- c(200000,   3200000)
  
  extent_poly <- st_bbox(
    c(xmin = xlim_alb[1], xmax = xlim_alb[2],
      ymin = ylim_alb[1], ymax = ylim_alb[2]),
    crs = st_crs(a)
  ) |> st_as_sfc()
  
  # ── 2. Graticule clipped inside box ──────────────────────────────────────────
  graticule <- st_graticule(
    lon = seq(-130, -60, by = 10),
    lat = seq(25, 55, by = 5),
    crs = st_crs(4326)
  ) |>
    st_transform(st_crs(a)) |>
    st_intersection(extent_poly)
  
  # ── 3. Build box border lines for intersection ────────────────────────────────
  # Bottom edge and left edge as sf linestrings in Albers
  bottom_edge <- st_linestring(matrix(
    c(xlim_alb[1], ylim_alb[1],
      xlim_alb[2], ylim_alb[1]), ncol = 2, byrow = TRUE
  )) |> st_sfc(crs = st_crs(a))
  
  left_edge <- st_linestring(matrix(
    c(xlim_alb[1], ylim_alb[1],
      xlim_alb[1], ylim_alb[2]), ncol = 2, byrow = TRUE
  )) |> st_sfc(crs = st_crs(a))
  
  # ── 4. Find exact intersection points of graticule with box edges ─────────────
  # Full graticule (unclipped) for intersection — need lines that reach the border
  graticule_full <- st_graticule(
    lon = seq(-130, -60, by = 10),
    lat = seq(25, 55, by = 5),
    crs = st_crs(4326)
  ) |> st_transform(st_crs(a))
  
  # Meridians (lon lines) intersect the bottom edge → lon labels
  meridians <- graticule_full |>
    filter(type == "E", degree %in% seq(-120, -70, by = 10))
  
  lon_pts <- st_intersection(meridians, bottom_edge) |>
    mutate(label = paste0(abs(degree), "°W"))
  
  # Parallels (lat lines) intersect the left edge → lat labels  
  parallels <- graticule_full |>
    filter(type == "N", degree %in% seq(30, 50, by = 5))
  
  lat_pts <- st_intersection(parallels, left_edge) |>
    mutate(label = paste0(degree, "°N"))
}


# =============================================================================
# EH figures (mammal): maps of the four EH measures + linear-vs-quadratic
# partial-residual plots from SAR models + inflection-point-vs-scale curves
#   (rows = measures, cols = map | richness partial | MSA partial | inflection~scale).
#
# Mammal-specific vs BBS (EH_figures_2.R):
#   - data is AGGREGATED to one row per Camera_Trap_Array (mean over years)
#   - covariates: elev, PopDen, NPP, Prec, Temp (all *_{scale}_t)
#   - response: raw richness_Estimator (no log) and log(MSA)
#   - map columns: EH_{scale}, elev_range, EH_delta, A57_sd  (elev_range has NO
#     scale suffix), so plot_eh_map() takes an explicit label + palette
#   - listw is built per scale from the camera-array centroids
#
# The four EH measures:
#   map column     predictor (_t)     legend / label                       palette
#   ------------   ----------------   ----------------------------------   -------
#   elev_range     EH_elevrange_t     Elevation_range                      rocket
#   EH_delta       EH_delta_t         Elevation range (mean-independent)   mako
#   EH_{scale}     EH_{scale}_t       EH (64-D)                            viridis
#   A57_sd         A57_sd_t           EH (1-D)                             plasma
# =============================================================================

library(ggplot2)
library(sf)
library(scales)        # squish
library(glue)
library(dplyr)
library(tidyr)         # drop_na
library(performance)   # r2
library(patchwork)
library(spatialreg)    # spautolm
library(spdep)         # knearneigh, knn2nb, nb2listw


# -----------------------------------------------------------------------------
# 1. Map one EH measure  (generalized: explicit column / label / palette)
# -----------------------------------------------------------------------------
# Colorbar clipped to `trim` quantiles with extremes squished into the end
# colors; NODATA (-9999) -> NA, drawn grey. Map "furniture" (states, graticule,
# lon_pts, lat_pts, xlim_alb, ylim_alb) is read from the workspace.
plot_eh_map <- function(a, measure,
                        label         = measure,
                        palette       = "viridis",
                        trim          = c(0.02, 0.98),
                        point_size    = 1.2,
                        drop_extremes = FALSE,
                        na_flag       = -9999) {
  
  if (!measure %in% names(a)) {
    stop(sprintf("plot_eh_map(): column '%s' is not in `a`. Columns present: %s.",
                 measure, paste(utils::head(names(a), 40), collapse = ", ")))
  }
  
  a[[measure]][a[[measure]] == na_flag] <- NA_real_
  lims <- stats::quantile(a[[measure]], probs = trim, na.rm = TRUE)
  
  if (drop_extremes) {
    keep <- !is.na(a[[measure]]) &
      a[[measure]] >= lims[1] & a[[measure]] <= lims[2]
    a <- a[keep, ]
  }
  
  ggplot() +
    geom_sf(data = states, fill = "grey80", color = "white", linewidth = 0.3) +
    geom_sf(data = graticule, color = "grey10", linewidth = 0.25, linetype = "dashed") +
    geom_sf(data = a, aes(color = .data[[measure]]), size = point_size) +
    scale_color_viridis_c(option = palette, name = label,
                          limits = lims, oob = scales::squish, na.value = "grey70") +
    geom_sf_text(data = lon_pts, aes(label = label),
                 size = 4, color = "black", vjust = 1.8, fontface = "plain") +
    geom_sf_text(data = lat_pts, aes(label = label),
                 size = 4, color = "black", hjust = 1.2, fontface = "plain") +
    coord_sf(crs = st_crs(a), xlim = xlim_alb, ylim = ylim_alb,
             expand = FALSE, clip = "off") +
    theme_minimal() +
    theme(
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid       = element_blank(),
      panel.background = element_rect(fill = "white"),
      axis.title       = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.margin    = margin(t = 10, b = 2),
      legend.title     = element_text(size = 18),
      legend.text      = element_text(size = 14),
      plot.margin      = margin(t = 5, r = 15, b = 20, l = 50)
    )
}


# -----------------------------------------------------------------------------
# 2. Partial-residual plot (linear vs quadratic), lm or SAR
# -----------------------------------------------------------------------------
# `data`/`lin_model` optional; for SAR pass lin_model (a pre-fitted no-quadratic
# model) so Delta AIC / R^2_lin work. Predictors are read from model$X.
sar_partial_plot <- function(model, var, scale, taxa,
                             quadratic    = FALSE,
                             data         = NULL,
                             lin_model    = NULL,
                             nsim         = 500,
                             n_seq        = 200,
                             seed         = 123,
                             point_size   = 0.4,
                             point_alpha  = 0.5,
                             line_color   = "steelblue",
                             linear_color = "firebrick",
                             ann_size     = 5,
                             x_label      = NULL,
                             y_label      = "Partial residual (log MSA)") {
  
  if (taxa %in% c("bird", "mammal")) {
    sar_coef_table <- summary(model)$Coef
    R2 <- cor(model$Y, model$fit$fitted.values)^2
  } else {
    sar_coef_table <- summary(model)$coefficients
    R2 <- unname(r2(model)[[2]])
  }
  
  coefs  <- coef(model)
  resids <- residuals(model)
  
  quad_name  <- paste0("I(", var, "^2)")
  has_quad   <- quadratic && (quad_name %in% names(coefs))
  coef_names <- if (has_quad) c(var, quad_name) else var
  
  if (!is.null(data)) {
    xvar <- data[[var]]
  } else if (taxa %in% c("bird", "mammal")) {
    xvar <- as.data.frame(model$X, check.names = FALSE)[[var]]
  } else {
    xvar <- stats::model.frame(model)[[var]]
  }
  
  x_seq <- seq(min(xvar, na.rm = TRUE), max(xvar, na.rm = TRUE), length.out = n_seq)
  
  if (has_quad) {
    fitted_line <- coefs[var] * x_seq + coefs[quad_name] * x_seq^2
    partial_res <- resids + coefs[var] * xvar + coefs[quad_name] * xvar^2
  } else {
    fitted_line <- coefs[var] * x_seq
    partial_res <- resids + coefs[var] * xvar
  }
  
  se_vec <- sar_coef_table[coef_names, "Std. Error"]
  set.seed(seed)
  coef_sims <- mapply(function(mu, se) rnorm(nsim, mean = mu, sd = se),
                      mu = coefs[coef_names], se = se_vec)
  if (!is.matrix(coef_sims)) coef_sims <- matrix(coef_sims, ncol = 1)
  X_mat <- if (has_quad) cbind(x_seq, x_seq^2) else matrix(x_seq, ncol = 1)
  sims  <- X_mat %*% t(coef_sims)
  
  fit_df <- data.frame(x = x_seq, fit = fitted_line,
                       lwr = apply(sims, 1, quantile, 0.025),
                       upr = apply(sims, 1, quantile, 0.975))
  res_df <- data.frame(x = xvar, res = partial_res)
  
  # linear fit: dashed line + Delta AIC + R^2_lin
  lin_ok <- FALSE; lin_df <- NULL
  stats_ok <- FALSE; dAIC <- NA_real_; R2_lin <- NA_real_
  if (has_quad) {
    lin_line <- lm(res ~ x, data = res_df)
    lin_df <- data.frame(x = x_seq,
                         fit = predict(lin_line, newdata = data.frame(x = x_seq)))
    lin_ok <- TRUE
    
    if (is.null(lin_model)) {
      lin_model <- tryCatch({
        if (taxa %in% c("bird", "mammal"))
          update(model, stats::as.formula(paste(". ~ . -", quad_name)))
        else {
          mf <- stats::model.frame(model); resp <- names(mf)[1]
          preds <- setdiff(names(mf)[-1], quad_name)
          lm(stats::as.formula(paste0("`", resp, "` ~ ",
                                      paste0("`", preds, "`", collapse = " + "))), data = mf)
        }
      }, error = function(e) NULL)
    }
    if (!is.null(lin_model)) {
      dAIC   <- AIC(lin_model) - AIC(model)
      R2_lin <- tryCatch(
        if (taxa %in% c("bird", "mammal"))
          cor(lin_model$Y, lin_model$fit$fitted.values)^2
        else unname(r2(lin_model)[[2]]),
        error = function(e) NA_real_)
      stats_ok <- TRUE
    }
  }
  
  inflection <- if (has_quad) -coefs[var] / (2 * coefs[quad_name]) else NA_real_
  infl_in_range <- has_quad && is.finite(inflection) &&
    inflection > min(xvar, na.rm = TRUE) && inflection < max(xvar, na.rm = TRUE)
  
  cat(sprintf("[%s] R2_quad=%.3f  R2_lin=%.3f  dAIC(lin-quad)=%.2f  inflection=%.3f%s\n",
              var, R2, R2_lin, dAIC, inflection,
              if (infl_in_range) "" else " (out of range)"))
  
  if (is.null(x_label)) x_label <- var
  
  p <- ggplot() +
    geom_point(data = res_df, aes(x = x, y = res),
               size = point_size, alpha = point_alpha, color = "grey50") +
    geom_ribbon(data = fit_df, aes(x = x, ymin = lwr, ymax = upr),
                fill = line_color, alpha = 0.25) +
    geom_line(data = fit_df, aes(x = x, y = fit), color = line_color, linewidth = 1)
  
  if (lin_ok) {
    p <- p + geom_line(data = lin_df, aes(x = x, y = fit),
                       color = linear_color, linetype = "dashed", linewidth = 0.8)
  }
  
  p <- p +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    labs(x = x_label, y = y_label) +
    theme_classic(base_size = 13) +
    theme(axis.title = element_text(size = 16),
          axis.text  = element_text(size = 13))
  
  if (infl_in_range) {
    p <- p + geom_vline(xintercept = inflection,
                        linetype = "dotted", color = "black", linewidth = 0.6)
  }
  
  fmt <- function(z, d = 2) formatC(z, format = "f", digits = d)
  lines <- character(0)
  if (stats_ok) {
    lines <- c(lines, paste0("Delta*AIC == ", fmt(dAIC, 1)))
    lines <- c(lines, paste0("R[quad]^2 == ", fmt(R2), '*","~~R[lin]^2 == ', fmt(R2_lin)))
  } else {
    lines <- c(lines, paste0("R[quad]^2 == ", fmt(R2)))
  }
  if (infl_in_range) lines <- c(lines, paste0("Inflection == ", fmt(inflection)))
  
  xr <- range(res_df$x, na.rm = TRUE)
  yr <- range(res_df$res, na.rm = TRUE)
  dy <- 0.10 * diff(yr)
  for (i in seq_along(lines)) {
    p <- p + annotate("text", x = xr[1], y = yr[2] - (i - 1) * dy,
                      label = lines[i], hjust = 0, vjust = 1,
                      size = ann_size, color = "black", parse = TRUE)
  }
  
  p
}


# -----------------------------------------------------------------------------
# 3. Fit the mammal SAR model for one measure
# -----------------------------------------------------------------------------
fit_eh_model_mammal <- function(eh_var,
                                response  = c("richness", "MSA"),
                                scale     = "10km",
                                quadratic = TRUE,
                                data      = dat,
                                listw     = lw) {
  response <- match.arg(response)
  covars <- glue("elev_{scale}_t + PopDen_{scale}_t + NPP_{scale}_t + ",
                 "Prec_{scale}_t + Temp_{scale}_t")
  lhs  <- if (response == "richness") "richness_Estimator" else "log(MSA)"
  quad <- if (quadratic) glue(" + I({eh_var}^2)") else ""
  form <- stats::as.formula(glue("{lhs} ~ {eh_var}{quad} + {covars}"))
  spatialreg::spautolm(formula = form, data = data, listw = listw, family = "SAR")
}


# -----------------------------------------------------------------------------
# 3b. Inflection point vs spatial scale  (the 4th column)
# -----------------------------------------------------------------------------
scale_to_km <- function(s) {
  num <- as.numeric(gsub("[^0-9.]", "", s))
  if (grepl("km$", s)) num else num / 1000
}
.fill <- function(tmpl, scale) as.character(glue::glue(tmpl))

#' Bootstrap the inflection point z = -b1 / (2 b2) of a quadratic EH term.
inflection_boot <- function(model, var, taxa = "mammal",
                            center = 0, scale_sd = 1, inv = identity,
                            var_range = NULL,
                            nsim = 500, seed = 123) {
  set.seed(seed)
  quad_name <- paste0("I(", var, "^2)")
  tab <- if (taxa %in% c("bird", "mammal")) summary(model)$Coef
  else summary(model)$coefficients
  
  b1 <- tab[var, "Estimate"];       s1 <- tab[var, "Std. Error"]
  b2 <- tab[quad_name, "Estimate"]; s2 <- tab[quad_name, "Std. Error"]
  
  z_obs  <- -b1 / (2 * b2)
  z_sims <- replicate(nsim, -rnorm(1, b1, s1) / (2 * rnorm(1, b2, s2)))
  z_sims <- z_sims[is.finite(z_sims)]
  
  in_range <- if (is.null(var_range)) TRUE else
    is.finite(z_obs) && z_obs >= var_range[1] && z_obs <= var_range[2]
  
  # clamp the CI to the observed predictor range so the back-transform (esp. exp
  # for log measures) can't blow up when b2~0 leaves the inflection ill-defined
  z_lo <- unname(quantile(z_sims, 0.025, na.rm = TRUE))
  z_hi <- unname(quantile(z_sims, 0.975, na.rm = TRUE))
  if (!is.null(var_range)) {
    clamp <- function(z) min(max(z, var_range[1]), var_range[2])
    z_lo <- clamp(z_lo); z_hi <- clamp(z_hi)
  }
  
  bt <- function(z) inv(center + scale_sd * z)
  data.frame(
    z        = bt(z_obs),
    lwr_95   = bt(z_lo),
    upr_95   = bt(z_hi),
    in_range = in_range
  )
}

#' Inflection point across scales for one measure & response (mammal SAR).
#' NOTE: mammal predictors are standardized at the year level then averaged per
#' array, so the raw-unit back-transform (original_scale = TRUE) is approximate;
#' use original_scale = FALSE for standardized (SD) units.
inflection_across_scales <- function(measure_cfg, response, scales,
                                     dat_by_scale, lw_by_scale,
                                     nsim = 500, seed = 123,
                                     response_label = response,
                                     original_scale = TRUE) {
  g    <- if (identical(measure_cfg$transform, "log")) log else identity
  ginv <- if (identical(measure_cfg$transform, "log")) exp else identity
  
  out <- lapply(scales, function(scale) {
    d     <- dat_by_scale[[scale]]
    var_t <- .fill(measure_cfg$var_tmpl, scale)
    raw   <- .fill(measure_cfg$raw_tmpl, scale)
    
    mdl <- fit_eh_model_mammal(var_t, response = response, scale = scale,
                               quadratic = TRUE, data = d, listw = lw_by_scale[[scale]])
    xr  <- range(as.data.frame(mdl$X, check.names = FALSE)[[var_t]], na.rm = TRUE)
    
    if (original_scale) {
      gv <- g(d[[raw]])
      ib <- inflection_boot(mdl, var = var_t, taxa = "mammal",
                            center = mean(gv, na.rm = TRUE),
                            scale_sd = sd(gv, na.rm = TRUE),
                            inv = ginv, var_range = xr, nsim = nsim, seed = seed)
    } else {
      ib <- inflection_boot(mdl, var = var_t, taxa = "mammal",
                            var_range = xr, nsim = nsim, seed = seed)
    }
    ib$scale <- scale; ib$scale_num <- scale_to_km(scale)
    ib
  })
  res <- do.call(rbind, out)
  res$response <- response_label
  res[res$in_range, , drop = FALSE]
}

#' Inflection-point-vs-scale panel (the 4th column), both responses overlaid.
plot_inflection_scale <- function(results_all, scales,
                                  y_lab = "Inflection point",
                                  x_lab = "Spatial scale",
                                  show_legend = TRUE) {
  scales_num <- vapply(scales, scale_to_km, numeric(1))
  cols <- c("Species richness"              = "steelblue",
            "log (Mean species abundance)"  = "firebrick")
  
  ggplot(results_all, aes(x = scale_num, y = z, color = response, fill = response)) +
    geom_ribbon(aes(ymin = lwr_95, ymax = upr_95), alpha = 0.15, color = NA) +
    geom_errorbar(aes(ymin = lwr_95, ymax = upr_95), width = 0.05, linewidth = 0.6) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.8) +
    scale_x_log10(breaks = scales_num, labels = scales) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    labs(x = x_lab, y = y_lab) +
    theme_classic(base_size = 13) +
    theme(
      axis.title       = element_text(size = 16),
      axis.text        = element_text(size = 13),
      legend.position  = if (show_legend) c(0.5, 1) else "none",
      legend.background = element_rect(fill = NA, color = NA),
      legend.direction = "vertical",
      legend.title     = element_blank(),
      legend.text      = element_text(size = 14)
    )
}

#' Build one mammal `dat` for a given scale (aggregated to one row per array).
#' ADJUST file paths / column names to your setup.
build_dat_mammal <- function(scale, mammal_table,
                             delta_csv = glue(p_gee("mammal_EH_delta/mammal_Elevation_EHdelta_{scale}.csv")),
                             a57_shp   = glue(p_gee("mammal_A57_EH/mammal_A57_EH_{scale}.shp"))) {
  eledelta <- read.csv(delta_csv) %>%
    dplyr::select(Camera_Trap_Array, elev_range, EH_delta)
  a57 <- sf::st_read(a57_shp, quiet = TRUE) %>%
    sf::st_drop_geometry() %>%
    dplyr::select(Camera_Trap_Array = Camera_Tra, A57_sd)
  merged <- dplyr::left_join(eledelta, a57, by = "Camera_Trap_Array")
  
  raw_cols <- c(glue("EH_{scale}"), "elev_range", "EH_delta", "A57_sd")
  
  mammal_table %>%
    dplyr::select(year_CTA, Year, Camera_Trap_Array,
                  richness_Estimator, MSA, dplyr::contains(glue("_{scale}"))) %>%
    dplyr::left_join(merged, by = "Camera_Trap_Array") %>%
    dplyr::filter(EH_delta > 0) %>%
    dplyr::mutate(
      "EH_{scale}_t"     := scale(.data[[glue("EH_{scale}")]])[, 1],
      EH_elevrange_t     = scale(log(elev_range))[, 1],
      EH_delta_t         = scale(log(EH_delta))[, 1],
      A57_sd_t           = scale(A57_sd)[, 1],
      "elev_{scale}_t"   := scale(sqrt(.data[[glue("elev_{scale}")]]))[, 1],
      "PopDen_{scale}_t" := scale(log1p(.data[[glue("PopDen_{scale}")]]))[, 1],
      "NPP_{scale}_t"    := scale(.data[[glue("NPP_{scale}")]])[, 1],
      "Prec_{scale}_t"   := scale(sqrt(.data[[glue("Prec_{scale}")]]))[, 1],
      "Temp_{scale}_t"   := scale(.data[[glue("Temp_{scale}")]])[, 1]
    ) %>%
    # aggregate to one row per array (keep raw EH cols for the maps)
    dplyr::group_by(Camera_Trap_Array) %>%
    dplyr::summarise(dplyr::across(c(richness_Estimator, MSA,
                                     dplyr::ends_with("_t"),
                                     dplyr::all_of(raw_cols)),
                                   ~ mean(.x, na.rm = TRUE)),
                     .groups = "drop") %>%
    tidyr::drop_na(richness_Estimator, MSA, dplyr::ends_with("_t"))
}

#' Build the SAR neighbour weights (listw) matching `dat`'s arrays.
build_lw_mammal <- function(dat, centroids, k = 8) {
  resid_plots <- dat %>%
    dplyr::select(Camera_Trap_Array) %>%
    dplyr::left_join(centroids, by = "Camera_Trap_Array") %>%
    sf::st_as_sf()
  coords <- sf::st_coordinates(resid_plots)
  nb <- spdep::knn2nb(spdep::knearneigh(coords, k = k))
  spdep::nb2listw(nb, style = "W")
}

#' One-call data prep for a scale: returns the aggregated `dat` + matching `lw`.
prepare_mammal <- function(scale, mammal_table, centroids, k = 8) {
  d <- build_dat_mammal(scale, mammal_table)
  list(dat = d, lw = build_lw_mammal(d, centroids, k = k))
}


# -----------------------------------------------------------------------------
# 4. Couple each EH measure's MAP with partial plots + inflection-vs-scale
# -----------------------------------------------------------------------------
# map_col / var_tmpl / raw_tmpl are glue templates resolved per scale.
eh_measures <- list(
  list(map_col = "elev_range", var_tmpl = "EH_elevrange_t", raw_tmpl = "elev_range",
       transform = "log",      palette = "rocket",  x_label = "Elevation_range"),
  list(map_col = "EH_delta",   var_tmpl = "EH_delta_t",     raw_tmpl = "EH_delta",
       transform = "log",      palette = "mako",    x_label = "Elevation range (mean-independent)"),
  list(map_col = "EH_{scale}", var_tmpl = "EH_{scale}_t",   raw_tmpl = "EH_{scale}",
       transform = "identity", palette = "viridis", x_label = "EH (64-D)"),
  list(map_col = "A57_sd",     var_tmpl = "A57_sd_t",       raw_tmpl = "A57_sd",
       transform = "identity", palette = "plasma",  x_label = "EH (1-D)")
)

#' Build the mammal figure: one row per EH measure, columns
#'   map | richness partial | MSA partial | inflection-point-vs-scale
#'
#' @param a            sf points object (camera-array centroids carrying the four
#'                     map columns at `main_scale`: EH_{main_scale}, elev_range,
#'                     EH_delta, A57_sd).
#' @param dat_by_scale named list of aggregated modelling data frames, per scale.
#' @param lw_by_scale  named list of listw objects, one per scale.
#' @param scales       scales for the inflection~scale column.
#' @param main_scale   scale used for the map + partial-plot columns (1-3).
#' @param responses    which partial-plot column(s) after the map.
#' @param ann_size     partial-plot annotation text size.
#' @param nsim,seed    inflection bootstrap settings.
#' @param inflection_original_scale  TRUE -> raw EH units (approx, see note);
#'                     FALSE -> standardized (SD) units.
#' @param ...          forwarded to plot_eh_map() (e.g. trim, point_size).
#' @return a patchwork object (rows = measures, cols = map + partials + inflection).
eh_map_partial_mammal <- function(a, dat_by_scale, lw_by_scale,
                                  scales     = c("1km", "5km", "10km", "20km", "50km"),
                                  main_scale = "10km",
                                  responses  = c("richness", "MSA"),
                                  ann_size   = 5,
                                  nsim       = 500,
                                  seed       = 123,
                                  inflection_original_scale = TRUE,
                                  ...) {
  responses <- match.arg(responses, c("richness", "MSA"), several.ok = TRUE)
  dat <- dat_by_scale[[main_scale]]
  lw  <- lw_by_scale[[main_scale]]
  y_lab <- function(resp) if (resp == "richness") "Partial residual (richness)"
  else "Partial residual (log MSA)"
  resp_label <- c(richness = "Species richness", MSA = "log (Mean species abundance)")
  
  rows <- lapply(seq_along(eh_measures), function(i) {
    m <- eh_measures[[i]]
    var_t_main <- .fill(m$var_tmpl, main_scale)
    map_col    <- .fill(m$map_col,  main_scale)
    
    # col 1: map (main scale)
    p_map <- plot_eh_map(a, map_col, label = m$x_label, palette = m$palette, ...)
    
    # cols 2..k: partial-residual plots (main scale)
    parts <- lapply(responses, function(resp) {
      mdl <- fit_eh_model_mammal(var_t_main, response = resp, scale = main_scale,
                                 quadratic = TRUE,  data = dat, listw = lw)
      lin <- fit_eh_model_mammal(var_t_main, response = resp, scale = main_scale,
                                 quadratic = FALSE, data = dat, listw = lw)
      sar_partial_plot(mdl, var = var_t_main, scale = main_scale, taxa = "mammal",
                       quadratic = TRUE, lin_model = lin, ann_size = ann_size,
                       x_label = m$x_label, y_label = y_lab(resp))
    })
    
    # last col: inflection point vs spatial scale (both responses)
    res <- do.call(rbind, lapply(responses, function(resp) {
      inflection_across_scales(m, resp, scales, dat_by_scale, lw_by_scale,
                               nsim = nsim, seed = seed,
                               response_label = resp_label[[resp]],
                               original_scale = inflection_original_scale)
    }))
    p_infl <- plot_inflection_scale(res, scales, show_legend = (i == 1))
    
    patchwork::wrap_plots(c(list(p_map), parts, list(p_infl)), nrow = 1)
  })
  
  patchwork::wrap_plots(rows, ncol = 1)
}


# =============================================================================
# 5. Usage
# =============================================================================
if (TRUE) {
  
  scales     <- c("1km", "5km", "10km", "20km", "50km")
  main_scale <- "5km"
  
  camera_array_centroids <- read_sf(p_clean("camera_array_centroids.shp")) %>%
    dplyr::rename(Camera_Trap_Array = Cmr_T_A)
  
  # 1. aggregated `dat` + matching `lw` per scale
  prep <- setNames(
    lapply(scales, prepare_mammal, mammal_table = mammal_table,
           centroids = camera_array_centroids),
    scales)
  dat_by_scale <- lapply(prep, `[[`, "dat")
  lw_by_scale  <- lapply(prep, `[[`, "lw")
  
  # 2. map object `a`: centroids + the main-scale EH columns
  states <- sf::st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))
  a <- camera_array_centroids %>%
    dplyr::left_join(dat_by_scale[[main_scale]], by = "Camera_Trap_Array")
  conus_bbox <- st_bbox(c(xmin = -125, xmax = -66, ymin = 24, ymax = 50),
                        crs = st_crs(a))
  a <- st_crop(a, conus_bbox) %>% st_transform(5070)
  # map extent / graticule / axis labels: built above (same extent, CRS 5070)
  
  # full 4-column panel: map | richness | MSA | inflection~scale ---------------
  fig <- eh_map_partial_mammal(a, dat_by_scale, lw_by_scale,
                               scales = scales, main_scale = main_scale)
  ggsave(p_fig("FigS19_mammal_alternative_EH.png"), fig, width = 24, height = 20, dpi = 300)
  
  # standardized (SD) inflection units instead of raw EH units:
  # eh_map_partial_mammal(a, dat_by_scale, lw_by_scale, scales = scales,
  #                       inflection_original_scale = FALSE)
  
}
