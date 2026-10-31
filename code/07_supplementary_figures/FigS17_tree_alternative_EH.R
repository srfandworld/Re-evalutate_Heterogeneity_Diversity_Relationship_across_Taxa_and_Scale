# =============================================================================
# FigS17_tree_alternative_EH.R
# Fig. S17 (trees): HDR under four EH measures (elevation range, mean-independent
# elevation range, AEF EH 64-D, A57 EH 1-D); rows = measures,
# cols = map | richness partial | MSA partial | inflection point vs scale.
# SAR models (spautolm, method = "LU", 8-NN weights on FIA plots in EPSG:5070),
# response log(qD) / log(MSA_new) - same specification as 05_HDR_models/01.
#
# Source   : EmbeddingBio/analysis/sfigure19_physicallyVar_FIA.R (2026-10-29; paths only changed)
# Inputs   : data/clean/FIA_table_{scale}_iNEXT.csv, FIA_plots.shp, plot_pts.shp
#            data/raw/gee_exports/FIA_embeddings/FIA_EH_delta/FIA_Elevation_EHdelta_{scale}.csv
#            data/raw/gee_exports/FIA_A57_EH/FIA_A57_EH_{scale}.shp
# Outputs  : output/figures/FigS17_tree_alternative_EH.png
# NOTE     : fits 16 SAR models at 10 km plus 8 per scale for the inflection
#            column (18,577-129,143 plots per scale) - slow and memory-heavy.
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

#using physically interepretable variable
library(sf)
library(dplyr)
library(glue)
library(spatialreg)   # spautolm()
library(spdep)        # knearneigh(), knn2nb(), nb2listw()

# FIA plot locations for the SAR spatial weights (as in 01_3taxa_modeling.R)
plots_sf_tree <- st_read(p_clean("FIA_plots.shp"), quiet = TRUE) %>%
  st_transform(crs = 5070)

#mean-independence heterogeneity measures
FIA_10km_eledelta <- read.csv(p_gee("FIA_embeddings/FIA_EH_delta/FIA_Elevation_EHdelta_10km.csv"))
FIA_10km_eledelta <- FIA_10km_eledelta %>% 
  dplyr::select(pltID,EH_delta)

#A57 EH
FIA_10km_A57 <- st_read(p_gee("FIA_A57_EH/FIA_A57_EH_10km.shp"))
FIA_10km_A57 <- FIA_10km_A57 %>% 
  st_drop_geometry() %>% 
  dplyr::select(pltID,A57_sd)
# 
# #FIA_table_10km
FIA_table_10km <- read.csv(p_clean("FIA_table_10km_iNEXT.csv"))
dat <- FIA_table_10km %>%
  rename(pltID = focal_pltID) %>%
  dplyr::select(pltID, qD, contains("_10km")) %>%
  left_join(FIA_10km_eledelta, by = "pltID") %>%   #mean-independence heterogeneity measures
  left_join(FIA_10km_A57, by = 'pltID') %>%        #A57 EH
  filter(EH_delta > 0) %>%
  mutate(EH_10km_t = scale(EH_10km)[,1],
         EH_elevrange_t = scale(log(elev_range_10km))[,1],
         EH_delta_t = scale(log(EH_delta))[,1],
         A57_sd_t = scale(A57_sd)[,1],
         soilph_10km_t = scale(log(soilph_10km))[,1],
         soilcec_10km_t = scale(log(soilcec_10km))[,1],
         elev_mean_10km_t = scale(sqrt(elev_mean_10km))[,1],
         MAP_10km_t = scale(sqrt(MAP_10km))[,1],
         MAT_10km_t = scale(MAT_10km)[,1]) %>%
  # complete cases for every modelled column, so spautolm drops no rows and the
  # spatial weights stay aligned with `dat`
  tidyr::drop_na(qD, MSA_new_10km, ends_with("_t"))
names(dat)

# SAR spatial weights: 8 nearest FIA plots, row-standardised (as build_lw_tree()
# in 01_3taxa_modeling.R)
# geom_dat <- dat %>%
#   dplyr::select(pltID) %>%
#   left_join(plots_sf_tree, by = "pltID") %>%
#   st_as_sf() %>%
#   st_transform(crs = 5070) 
#   
# stopifnot(nrow(geom_dat) == nrow(dat))
# lw <- nb2listw(knn2nb(knearneigh(st_coordinates(geom_dat), k = 8)), style = "W")

#histograms of response variables
# hist(dat$EH_10km_t)
# hist(dat$EH_elevrange_t)
# hist(dat$EH_delta_t)
# hist(dat$A57_sd_t)




#Draw the elev_range, mean-independent elev_range, EH_64d, EH_1d
plot_pts <- st_read(p_clean("plot_pts.shp"))
states <- sf::st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))

a <- plot_pts %>% 
  left_join(dat,by="pltID")

conus_bbox <- st_bbox(c(
  xmin = -125,
  xmax = -66,
  ymin = 24,
  ymax = 50
), crs = st_crs(a))
a <- st_crop(a, conus_bbox) %>% st_transform(5070)

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
# EH figures (FIA): maps of the four EH measures + linear-vs-quadratic partial-
# residual plots from SAR models + inflection-point-vs-scale curves, and a driver
# that couples them per measure
#   (rows = measures, cols = map | richness partial | MSA partial | inflection~scale).
#
# Models follow the main tree models in Code/Figure2/01_3taxa_modeling.R:
#   spatialreg::spautolm(family = "SAR", method = "LU") with row-standardised
#   8-nearest-neighbour weights on FIA plot locations (data_clean/FIA_plots.shp);
#   response log(qD) (iNEXT, 90% coverage) or log(MSA_new_<scale>).
#   As in sfigure_physicallyVar_BBS.R, the linear (no-quadratic) SAR model is
#   fitted in the driver and passed to sar_partial_plot() for Delta AIC / R^2_lin,
#   and the weights are rebuilt per scale (each scale keeps different plots).
#
# Contents
#   1. plot_eh_map()             - map one EH measure, outlier-robust colorbar
#   2. sar_partial_plot()        - partial-residual plot (quadratic vs dashed
#                                  linear; Delta AIC, R^2_quad, R^2_lin, inflection)
#   3. fit_eh_model()            - fit the quadratic EH SAR model for one measure
#   3b. inflection helpers       - inflection_boot / inflection_across_scales /
#                                  plot_inflection_scale / build_dat_fia /
#                                  build_lw_fia / prepare_fia
#   4. eh_map_partial()          - the combined 4-column panel for all measures
#   5. usage
#
# The four EH measures:
#   column            predictor (_t)     legend / label                       palette
#   ---------------   ----------------   ----------------------------------   -------
#   elev_range_10km   EH_elevrange_t     Elevation_range                      rocket
#   EH_delta          EH_delta_t         Elevation range (mean-independent)   mako
#   EH_10km           EH_10km_t          EH (64-D)                            viridis
#   A57_sd            A57_sd_t           EH (1-D)                             plasma
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
# 1. Map one EH measure
# -----------------------------------------------------------------------------
# Colorbar clipped to `trim` quantiles with extremes squished into the end
# colors; NODATA (-9999, from GEE) -> NA, drawn grey. Map "furniture" (states,
# graticule, lon_pts, lat_pts, xlim_alb, ylim_alb) is read from the workspace.
plot_eh_map <- function(a,
                        measure = c("EH_10km", "elev_range_10km",
                                    "EH_delta", "A57_sd"),
                        trim          = c(0.02, 0.98),
                        point_size    = 0.3,
                        drop_extremes = FALSE,
                        na_flag       = -9999) {
  
  measure <- match.arg(measure)
  
  meta <- switch(measure,
                 EH_10km         = list(label = "EH (64-D)",             option = "viridis"),
                 elev_range_10km = list(label = "Elevation_range",       option = "rocket"),
                 EH_delta        = list(label = "Elevation range (mean-independent)", option = "mako"),
                 A57_sd          = list(label = "EH (1-D)",              option = "plasma")
  )
  
  if (!measure %in% names(a)) {
    have <- intersect(c("EH_10km", "elev_range_10km", "EH_delta", "A57_sd"), names(a))
    stop(sprintf("plot_eh_map(): column '%s' is not in `a`. EH columns present: %s.",
                 measure, if (length(have)) paste(have, collapse = ", ") else "none"))
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
    scale_color_viridis_c(option = meta$option, name = meta$label,
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
# 2. Partial-residual plot (linear vs quadratic)
# -----------------------------------------------------------------------------
# `data`      : data frame holding the predictor `var`, aligned row-for-row with
#               residuals(model) (pass the `dat` the SAR model was fitted on).
#               If NULL, the predictor is taken from the model's design matrix.
# `lin_model` : pre-fitted SAR model WITHOUT the quadratic term, used for
#               Delta AIC and R^2_lin (a SAR built inside a helper can't be
#               update()-d reliably, so it is fitted in the driver and passed in).
sar_partial_plot <- function(model,        # spautolm (SAR)
                             var,           # predictor name e.g. "EH_delta_t"
                             scale,
                             taxa,
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
  
  # ── coef table + R^2 (quadratic / original): SAR for all taxa ──────────────
  sar_coef_table <- summary(model)$Coef
  R2 <- cor(model$Y, model$fit$fitted.values)^2

  # ── model components ───────────────────────────────────────────────────────
  coefs  <- coef(model)
  resids <- residuals(model)

  quad_name <- paste0("I(", var, "^2)")
  has_quad  <- quadratic && (quad_name %in% names(coefs))
  coef_names <- if (has_quad) c(var, quad_name) else var

  # ── predictor values, aligned row-for-row with residuals(model) ────────────
  xvar <- if (!is.null(data)) data[[var]] else
    as.data.frame(model$X, check.names = FALSE)[[var]]
  stopifnot(length(xvar) == length(resids))

  # ── predictor sequence ─────────────────────────────────────────────────────
  x_seq <- seq(min(xvar, na.rm = TRUE),
               max(xvar, na.rm = TRUE),
               length.out = n_seq)

  # ── quadratic fitted line + partial residuals ──────────────────────────────
  if (has_quad) {
    fitted_line <- coefs[var] * x_seq + coefs[quad_name] * x_seq^2
    partial_res <- resids + coefs[var] * xvar + coefs[quad_name] * xvar^2
  } else {
    fitted_line <- coefs[var] * x_seq
    partial_res <- resids + coefs[var] * xvar
  }
  
  # ── confidence band via coefficient simulation ─────────────────────────────
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

  # ── LINEAR fit (no quadratic term): dashed line + Delta AIC + linear R^2 ────
  lin_ok <- FALSE; lin_df <- NULL
  stats_ok <- FALSE; dAIC <- NA_real_; R2_lin <- NA_real_
  if (has_quad) {
    # dashed line = linear regression through the partial residuals
    lin_line <- lm(res ~ x, data = res_df)
    lin_df <- data.frame(x = x_seq,
                         fit = predict(lin_line, newdata = data.frame(x = x_seq)))
    lin_ok <- TRUE

    # full SAR model without the quadratic term (fitted in the driver)
    if (!is.null(lin_model)) {
      dAIC   <- AIC(lin_model) - AIC(model)        # >0 => quadratic preferred
      R2_lin <- tryCatch(cor(lin_model$Y, lin_model$fit$fitted.values)^2,
                         error = function(e) NA_real_)
      stats_ok <- TRUE
    }
  }

  # ── inflection point (quadratic only) ──────────────────────────────────────
  inflection <- if (has_quad) -coefs[var] / (2 * coefs[quad_name]) else NA_real_
  infl_in_range <- has_quad && is.finite(inflection) &&
    inflection > min(xvar, na.rm = TRUE) &&
    inflection < max(xvar, na.rm = TRUE)
  
  cat(sprintf("[%s] R2_quad=%.3f  R2_lin=%.3f  dAIC(lin-quad)=%.2f  inflection=%.3f%s\n",
              var, R2, R2_lin, dAIC, inflection,
              if (infl_in_range) "" else " (out of range)"))
  
  if (is.null(x_label)) x_label <- var
  
  # ── plot ───────────────────────────────────────────────────────────────────
  p <- ggplot() +
    geom_point(data = res_df, aes(x = x, y = res),
               size = point_size, alpha = point_alpha, color = "grey50") +
    geom_ribbon(data = fit_df, aes(x = x, ymin = lwr, ymax = upr),
                fill = line_color, alpha = 0.25) +
    geom_line(data = fit_df, aes(x = x, y = fit),
              color = line_color, linewidth = 1)
  
  if (lin_ok) {  # dashed linear fit
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
  
  # ── annotation: 1) Delta AIC  2) R^2 (quad & lin)  3) inflection ────────────
  fmt <- function(z, d = 2) formatC(z, format = "f", digits = d)
  lines <- character(0)
  if (stats_ok) {
    lines <- c(lines, paste0("Delta*AIC == ", fmt(dAIC, 1)))
    lines <- c(lines, paste0("R[quad]^2 == ", fmt(R2), '*","~~R[lin]^2 == ', fmt(R2_lin)))
  } else {
    lines <- c(lines, paste0("R[quad]^2 == ", fmt(R2)))
  }
  if (infl_in_range) {
    lines <- c(lines, paste0("Inflection == ", fmt(inflection)))
  }
  
  xr <- range(res_df$x,   na.rm = TRUE)
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
# 3. Fit the EH SAR model for one measure
# -----------------------------------------------------------------------------
# quadratic = TRUE  -> EH + I(EH^2) + covariates   (main model)
# quadratic = FALSE -> EH         + covariates     (linear, for Delta AIC / R^2_lin)
# Same specification as sar_modeling_tree() in Code/Figure2/01_3taxa_modeling.R.
fit_eh_model <- function(eh_var,
                         response  = c("richness", "MSA"),
                         scale     = "10km",
                         quadratic = TRUE,
                         data      = dat,
                         listw     = lw) {
  response <- match.arg(response)
  covars <- glue("MAT_{scale}_t + MAP_{scale}_t + soilph_{scale}_t + ",
                 "soilcec_{scale}_t + elev_mean_{scale}_t")
  lhs <- if (response == "richness")
    glue("log(qD)") else glue("log(MSA_new_{scale})")
  quad <- if (quadratic) glue(" + I({eh_var}^2)") else ""
  form <- stats::as.formula(glue("{lhs} ~ {eh_var}{quad} + {covars}"))
  spatialreg::spautolm(formula = form, data = data, listw = listw,
                       method = "LU", family = "SAR")
}


# -----------------------------------------------------------------------------
# 3b. Inflection point vs spatial scale  (the 4th column)
# -----------------------------------------------------------------------------

# "1km" -> 1, "50km" -> 50  (for the log x-axis)
scale_to_km <- function(s) {
  num <- as.numeric(gsub("[^0-9.]", "", s))
  if (grepl("km$", s)) num else num / 1000
}

# resolve a glue template like "EH_{scale}_t" for a given scale
.fill <- function(tmpl, scale) as.character(glue::glue(tmpl))

#' Bootstrap the inflection point z = -b1 / (2 b2) of a quadratic EH term.
#' Coefficients are simulated from their sampling distribution (no refit).
#' By default returns the STANDARDIZED inflection; supply center/scale_sd/inv to
#' map back to raw units as inv(center + scale_sd*z). `var_range` (the observed
#' standardized predictor range) flags whether the vertex is a real inflection.
inflection_boot <- function(model, var, taxa = "tree",
                            center = 0, scale_sd = 1, inv = identity,
                            var_range = NULL,
                            nsim = 500, seed = 123) {
  set.seed(seed)
  quad_name <- paste0("I(", var, "^2)")
  tab <- if (inherits(model, "Spautolm")) summary(model)$Coef
  else summary(model)$coefficients
  
  b1 <- tab[var, "Estimate"];       s1 <- tab[var, "Std. Error"]
  b2 <- tab[quad_name, "Estimate"]; s2 <- tab[quad_name, "Std. Error"]
  
  z_obs  <- -b1 / (2 * b2)
  z_sims <- replicate(nsim, -rnorm(1, b1, s1) / (2 * rnorm(1, b2, s2)))
  z_sims <- z_sims[is.finite(z_sims)]            # drop b2~0 blow-ups
  
  # a vertex is only a real inflection if it sits inside the observed range
  in_range <- if (is.null(var_range)) TRUE else
    is.finite(z_obs) && z_obs >= var_range[1] && z_obs <= var_range[2]
  
  bt <- function(z) inv(center + scale_sd * z)
  data.frame(
    z        = bt(z_obs),
    lwr_95   = bt(unname(quantile(z_sims, 0.025, na.rm = TRUE))),
    upr_95   = bt(unname(quantile(z_sims, 0.975, na.rm = TRUE))),
    in_range = in_range
  )
}

#' Inflection point across spatial scales for one EH measure & response.
#' Fits the quadratic SAR model at each scale (dat_by_scale[[scale]] with its
#' matching lw_by_scale[[scale]]) and bootstraps the inflection;
#' original_scale = TRUE maps back to the measure's raw units.
#' Out-of-range scales (no real inflection) are dropped.
inflection_across_scales <- function(measure_cfg, response, scales,
                                     dat_by_scale, lw_by_scale, taxa = "tree",
                                     nsim = 500, seed = 123,
                                     response_label = response,
                                     original_scale = TRUE) {
  g    <- if (identical(measure_cfg$transform, "log")) log else identity
  ginv <- if (identical(measure_cfg$transform, "log")) exp else identity

  out <- lapply(scales, function(scale) {
    d     <- dat_by_scale[[scale]]
    var_t <- .fill(measure_cfg$var_tmpl, scale)
    raw   <- .fill(measure_cfg$raw_tmpl, scale)

    mdl <- fit_eh_model(var_t, response = response, scale = scale,
                        quadratic = TRUE, data = d, listw = lw_by_scale[[scale]])
    # standardized predictor range (SAR design matrix) -> in-range check
    xr  <- range(as.data.frame(mdl$X, check.names = FALSE)[[var_t]], na.rm = TRUE)
    
    if (original_scale) {
      gv <- g(d[[raw]])
      ib <- inflection_boot(mdl, var = var_t, taxa = taxa,
                            center = mean(gv, na.rm = TRUE),
                            scale_sd = sd(gv, na.rm = TRUE),
                            inv = ginv, var_range = xr, nsim = nsim, seed = seed)
    } else {
      ib <- inflection_boot(mdl, var = var_t, taxa = taxa,
                            var_range = xr, nsim = nsim, seed = seed)
    }
    ib$scale <- scale; ib$scale_num <- scale_to_km(scale)
    ib
  })
  res <- do.call(rbind, out)
  res$response <- response_label
  res[res$in_range, , drop = FALSE]   # drop scales without a real inflection
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

#' Build one FIA `dat` for a given scale (mirrors your 10km pipeline).
#' ADJUST file paths / column names to your setup. Produces the *_t predictors
#' the models expect.
#' `FIA_table` is the per-scale iNEXT table (FIA_table_<scale>_iNEXT.csv); its
#' `focal_pltID` is renamed to `pltID`. Rows with any missing modelled column are
#' dropped so spautolm keeps every row and the spatial weights stay aligned.
build_dat_fia <- function(scale, FIA_table,
                          delta_csv = glue(p_gee("FIA_embeddings/FIA_EH_delta/FIA_Elevation_EHdelta_{scale}.csv")),
                          a57_shp   = glue(p_gee("FIA_A57_EH/FIA_A57_EH_{scale}.shp"))) {
  eledelta <- read.csv(delta_csv) %>% dplyr::select(pltID, EH_delta)
  a57 <- sf::st_read(a57_shp, quiet = TRUE) %>%
    sf::st_drop_geometry() %>% dplyr::select(pltID, A57_sd)

  if ("focal_pltID" %in% names(FIA_table))
    FIA_table <- dplyr::rename(FIA_table, pltID = focal_pltID)

  FIA_table %>%
    dplyr::select(pltID, qD, dplyr::contains(glue("_{scale}"))) %>%
    dplyr::left_join(eledelta, by = "pltID") %>%
    dplyr::left_join(a57,      by = "pltID") %>%
    dplyr::filter(EH_delta > 0) %>%
    dplyr::mutate(
      "EH_{scale}_t"        := scale(.data[[glue("EH_{scale}")]])[, 1],
      EH_elevrange_t        = scale(log(.data[[glue("elev_range_{scale}")]]))[, 1],
      EH_delta_t            = scale(log(EH_delta))[, 1],
      A57_sd_t              = scale(A57_sd)[, 1],
      "soilph_{scale}_t"    := scale(log(.data[[glue("soilph_{scale}")]]))[, 1],
      "soilcec_{scale}_t"   := scale(log(.data[[glue("soilcec_{scale}")]]))[, 1],
      "elev_mean_{scale}_t" := scale(sqrt(.data[[glue("elev_mean_{scale}")]]))[, 1],
      "MAP_{scale}_t"       := scale(sqrt(.data[[glue("MAP_{scale}")]]))[, 1],
      "MAT_{scale}_t"       := scale(.data[[glue("MAT_{scale}")]])[, 1]
    ) %>%
    tidyr::drop_na(qD, dplyr::all_of(glue("MSA_new_{scale}")), dplyr::ends_with("_t"))
}

#' SAR neighbour weights (listw) matching `dat`'s plots: 8 nearest FIA plots,
#' row-standardised. Same construction as build_lw_tree() in
#' Code/Figure2/01_3taxa_modeling.R. MUST be rebuilt per scale because each
#' scale's `dat` keeps a different set of plots.
build_lw_fia <- function(dat, plots_sf, k = 8) {
  geom <- dat %>%
    dplyr::select(pltID) %>%
    dplyr::left_join(plots_sf, by = "pltID") %>%
    sf::st_as_sf()
  stopifnot(nrow(geom) == nrow(dat))   # listw must align with data row-for-row
  coords <- sf::st_coordinates(geom)
  spdep::nb2listw(spdep::knn2nb(spdep::knearneigh(coords, k = k)), style = "W")
}

#' One-call data prep for a scale: returns the modelling `dat` and its matching
#' spatial weights `lw` together (use to build dat_by_scale / lw_by_scale).
prepare_fia <- function(scale, FIA_table, plots_sf, k = 8) {
  d <- build_dat_fia(scale, FIA_table)
  list(dat = d, lw = build_lw_fia(d, plots_sf, k = k))
}


# -----------------------------------------------------------------------------
# 4. Couple each EH measure's MAP with partial plots + inflection-vs-scale
# -----------------------------------------------------------------------------
# var_tmpl / raw_tmpl are glue templates resolved per scale; transform is the
# function applied before scale() (drives the inflection back-transform).
eh_measures <- list(
  list(map_col = "elev_range_10km", var_tmpl = "EH_elevrange_t", raw_tmpl = "elev_range_{scale}",
       transform = "log",      x_label = "Elevation_range"),
  list(map_col = "EH_delta",        var_tmpl = "EH_delta_t",     raw_tmpl = "EH_delta",
       transform = "log",      x_label = "Elevation range (mean-independent)"),
  list(map_col = "EH_10km",         var_tmpl = "EH_{scale}_t",   raw_tmpl = "EH_{scale}",
       transform = "identity", x_label = "EH (64-D)"),
  list(map_col = "A57_sd",          var_tmpl = "A57_sd_t",       raw_tmpl = "A57_sd",
       transform = "identity", x_label = "EH (1-D)")
)

#' Build the FIA figure: one row per EH measure, columns
#'   map | richness partial | MSA partial | inflection-point-vs-scale
#'
#' @param a            sf points object for the maps (carries the four map columns).
#' @param dat_by_scale named list of modelling data frames, one per scale.
#' @param scales       scales for the inflection~scale column (4th col).
#' @param main_scale   scale used for the map + partial-plot columns (1-3).
#' @param responses    which partial-plot column(s) after the map.
#' @param lw_by_scale  named list of SAR spatial weights (listw), one per scale,
#'                     each matching the rows of dat_by_scale[[scale]].
#' @param taxa         "tree" (FIA SAR).
#' @param ann_size     partial-plot annotation text size.
#' @param nsim,seed    inflection bootstrap settings.
#' @param inflection_original_scale  TRUE -> raw EH units; FALSE -> SD units.
#' @param ...          forwarded to plot_eh_map() (e.g. trim, drop_extremes).
#' @return a patchwork object (rows = measures, cols = map + partials + inflection).
eh_map_partial <- function(a, dat_by_scale, lw_by_scale,
                           scales     = c("1km", "5km", "10km", "20km", "50km"),
                           main_scale = "10km",
                           responses  = c("richness", "MSA"),
                           taxa       = "tree",
                           ann_size   = 5,
                           nsim       = 500,
                           seed       = 123,
                           inflection_original_scale = TRUE,
                           ...) {
  responses <- match.arg(responses, c("richness", "MSA"), several.ok = TRUE)
  dat <- dat_by_scale[[main_scale]]
  lw  <- lw_by_scale[[main_scale]]
  y_lab <- function(resp) if (resp == "richness") "Partial residual (log richness)"
  else "Partial residual (log MSA)"
  resp_label <- c(richness = "Species richness", MSA = "log (Mean species abundance)")

  rows <- lapply(seq_along(eh_measures), function(i) {
    m <- eh_measures[[i]]
    var_t_main <- .fill(m$var_tmpl, main_scale)

    # col 1: map (main scale)
    p_map <- plot_eh_map(a, m$map_col, ...)

    # cols 2..k: partial-residual plots (main scale). The linear SAR model is
    # fitted here and passed in for Delta AIC / R^2_lin.
    parts <- lapply(responses, function(resp) {
      mdl <- fit_eh_model(var_t_main, response = resp, scale = main_scale,
                          quadratic = TRUE,  data = dat, listw = lw)
      lin <- fit_eh_model(var_t_main, response = resp, scale = main_scale,
                          quadratic = FALSE, data = dat, listw = lw)
      sar_partial_plot(mdl, var = var_t_main, scale = main_scale, taxa = taxa,
                       quadratic = TRUE, data = dat, lin_model = lin,
                       ann_size = ann_size,
                       x_label = m$x_label, y_label = y_lab(resp))
    })

    # last col: inflection point vs spatial scale (both responses)
    res <- do.call(rbind, lapply(responses, function(resp) {
      inflection_across_scales(m, resp, scales, dat_by_scale, lw_by_scale, taxa = taxa,
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
# You need: a `dat` + matching SAR weights `lw` per scale (dat_by_scale,
# lw_by_scale), the map object `a`, and the map furniture (states, graticule,
# lon_pts, lat_pts, xlim_alb, ylim_alb) built at the top of this script.
# NOTE: SAR (method = "LU") on 18,577-129,143 plots per scale is slow and
# memory-heavy; the full panel fits 16 models at the main scale plus 8 per
# extra scale for the inflection column.
if (TRUE) {

  scales <- c("1km", "5km", "10km", "20km",'50km')

  # 1. one `dat` + matching `lw` per scale. Each scale has its own iNEXT table
  #    and keeps a different set of plots, so the weights MUST be rebuilt per scale.
  prep <- setNames(lapply(scales, function(s) {
    prepare_fia(s,
                FIA_table = read.csv(glue(p_clean("FIA_table_{s}_iNEXT.csv"))),
                plots_sf  = plots_sf_tree)
  }), scales)
  dat_by_scale <- lapply(prep, `[[`, "dat")
  lw_by_scale  <- lapply(prep, `[[`, "lw")

  # 2. the map object `a` and the map furniture are built at the top of this script.

  # single map / single partial plot ------------------------------------------
  plot_eh_map(a, "EH_delta")

  mq <- fit_eh_model("EH_delta_t", response = "richness", scale = "10km", quadratic = TRUE,
                     data = dat_by_scale[["10km"]], listw = lw_by_scale[["10km"]])
  ml <- fit_eh_model("EH_delta_t", response = "richness", scale = "10km", quadratic = FALSE,
                     data = dat_by_scale[["10km"]], listw = lw_by_scale[["10km"]])
  sar_partial_plot(mq, "EH_delta_t", scale = "10km", taxa = "tree",
                   quadratic = TRUE, data = dat_by_scale[["10km"]], lin_model = ml,
                   x_label = "Elevation range (mean-independent)",
                   y_label = "Partial residual (log richness)")

  # full 4-column panel: map | richness | MSA | inflection~scale ---------------
  fig <- eh_map_partial(a, dat_by_scale, lw_by_scale, scales = scales,
                        main_scale = "10km", taxa = "tree")
  ggsave(p_fig("FigS17_tree_alternative_EH.png"), fig, width = 24, height = 20, dpi = 300)

  # standardized (SD) inflection units instead of raw EH units:
  # eh_map_partial(a, dat_by_scale, lw_by_scale, scales = scales,
  #                inflection_original_scale = FALSE)

}
