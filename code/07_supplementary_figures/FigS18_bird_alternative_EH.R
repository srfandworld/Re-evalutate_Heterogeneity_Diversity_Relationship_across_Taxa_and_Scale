# =============================================================================
# FigS18_bird_alternative_EH.R
# Fig. S18 (birds): same layout as Fig. S17 for BBS routes. SAR models
# (spautolm, 8-NN weights on route start points, EPSG:5070); response richness /
# log(MSA) - same specification as 05_HDR_models/01.
#
# Source   : EmbeddingBio/analysis/sfigure18_physicallyVar_BBS.R (2026-10-29; paths only
#            changed; libraries loaded up front because the original used them
#            before its library() calls)
# Inputs   : data/clean/BBS_table.csv, bbsroutes_shp/bbsrtsl020_filterd.shp
#            data/raw/gee_exports/BBS_EH_delta/BBS_Elevation_EHdelta_{scale}.csv
#            data/raw/gee_exports/BBS_A57_EH/BBS_A57_EH_{scale}.shp
# Outputs  : output/figures/FigS18_bird_alternative_EH.png
# NOTE (reviewer C3): the EH-delta CSV and A57 shapefile are joined by ROW ORDER
#   (cbind), as in the original; confirm both exports share row order.
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

#mean-independence heterogeneity measures
BBS_10km_eledelta <- read.csv(p_gee("BBS_EH_delta/BBS_Elevation_EHdelta_10km.csv"))
BBS_10km_eledelta <- BBS_10km_eledelta %>%
  dplyr::select(BBSRTSL020,elev_range,EH_delta)

#A57 EH
BBS_10km_A57 <- st_read(p_gee("BBS_A57_EH/BBS_A57_EH_10km.shp"))
BBS_10km_A57 <- BBS_10km_A57 %>%
  st_drop_geometry() %>%
  dplyr::select(RTENAME,RTENO,A57_sd)

#merge two table together
merged_table <- cbind(BBS_10km_eledelta,BBS_10km_A57) %>%
  dplyr::select(RTENAME,RTENO,elev_range,EH_delta,A57_sd) %>%
  group_by(RTENAME,RTENO) %>%
  summarise(elev_range_10km = mean(elev_range),
            EH_delta = mean(EH_delta),
            A57_sd = mean(A57_sd)) %>%
  ungroup()

# BBS_table
BBS_table <- read.csv(p_clean("BBS_table.csv"))

dat <- BBS_table %>%
  dplyr::select(route,RTENAME,RTENO,richness,MSA, contains("_10km")) %>%
  left_join(merged_table, by = c("RTENAME",'RTENO')) %>%
  filter(EH_delta > 0) %>%
  mutate(EH_10km_t = scale(EH_10km)[,1],
         EH_elevrange_t = scale(log(elev_range_10km))[,1],
         EH_delta_t = scale(log(EH_delta))[,1],
         A57_sd_t = scale(A57_sd)[,1],
         MAT_10km_t = scale(MAT_10km)[,1],
         MAP_10km_t = scale(MAP_10km)[,1],
         Tsea_10km_t = scale(Tsea_10km)[,1],
         Psea_10km_t = scale(Psea_10km)[,1],
         elev_10km_t = scale(sqrt(elev_10km))[,1],
         PopDen_10km_t = scale(log1p(PopDen_10km))[,1])
names(dat)
# 
# 
# #histograms of response variables
# hist(dat$EH_10km_t)
# hist(dat$EH_elevrange_t)
# hist(dat$EH_delta_t)
# hist(dat$A57_sd_t)
# 
# scale <- '10km'
# 
# #neightbour list
# bbsrtsl020_filtered <- st_read(p_clean("bbsroutes_shp/bbsrtsl020_filterd.shp"))
# resid_plots <- bbsrtsl020_filtered %>% 
#   #filter to routes in dat
#   filter(str_c(RTENO,RTENAME,sep = "_") %in% str_c(dat$RTENO,dat$RTENAME,sep = "_")) %>% 
#   #join geometry when RTE is same one
#   mutate(RTE = str_c(RTENO, RTENAME, sep = "_")) %>%
#   group_by(RTE) %>%
#   summarise(geometry = st_union(geometry), .groups = "drop") %>% 
#   #reorder the RTE by dat
#   arrange(match(RTE,str_c(dat$RTENO,dat$RTENAME,sep = "_")))
# 
# #k nearest neighbors
# coords <- st_coordinates(resid_plots %>% lwgeom::st_startpoint())
# nb <- knn2nb(knearneigh(coords, k = 8))
# lw <- nb2listw(nb, style = "W")


# #richness-EH
# form <- as.formula(glue("richness ~ EH_delta_t +I(EH_delta_t^2) +MAT_{scale}_t+ MAP_{scale}_t + Tsea_{scale}_t + Psea_{scale}_t + elev_{scale}_t + PopDen_{scale}_t  "))
# #MSA-EH
# form <- as.formula(glue("log(MSA) ~ EH_{scale}_t +I(EH_{scale}_t^2) +MAT_{scale}_t+ MAP_{scale}_t + Tsea_{scale}_t + Psea_{scale}_t + elev_{scale}_t + PopDen_{scale}_t  "))
# 
# final_model <- spautolm(formula = form,
#                       data = dat,
#                       listw = lw,   #go to original files to get this
#                       family = "SAR")
# 
# 
# summary(final_model)
# AIC(final_model)

#Draw the elev_range, mean-independent elev_range, EH_64d, EH_1d
states <- sf::st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))
bbsrtsl020_filtered <- st_read(p_clean("bbsroutes_shp/bbsrtsl020_filterd.shp"))

a <- bbsrtsl020_filtered %>%
  left_join(dat, by=c("RTENO","RTENAME")) %>%
  drop_na()


conus_bbox <- st_bbox(c(
  xmin = -125,
  xmax = -66,
  ymin = 24,
  ymax = 50
), crs = st_crs(a))

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
# EH figures (BBS): maps of the four EH measures + linear-vs-quadratic partial-
# residual plots from SAR models + inflection-point-vs-scale curves, and a driver
# that couples them per measure
#   (rows = measures, cols = map | richness partial | MSA partial | inflection~scale).
#
# BBS-specific vs the FIA version (EH_figures.R):
#   - models are spatial SAR (spautolm, family = "SAR") fitted with a listw
#   - covariates: MAT, MAP, Tsea, Psea, elev, PopDen (all *_10km_t)
#   - response: raw richness (no log/effort) and log(MSA)
#   - linear (no-quadratic) model is fitted in the driver and passed in, so
#     Delta AIC / R^2_lin work for SAR (update() can't refit a model built
#     inside a helper)
#
# The four EH measures (same columns/labels as the FIA figure):
#   column            predictor (_t)     legend / label                       palette
#   ---------------   ----------------   ----------------------------------   -------
#   elev_range_10km   EH_elevrange_t     Elevation_range                      rocket
#   EH_delta          EH_delta_t         Elevation range (mean-independent)   mako
#   EH_10km           EH_10km_t          EH (64-D)                            viridis
#   A57_sd            A57_sd_t           EH (1-D)                             plasma
#
# Contents: 1) plot_eh_map  2) sar_partial_plot  3) fit_eh_model_bbs
#           4) eh_map_partial_bbs  5) usage
# =============================================================================

library(ggplot2)
library(sf)
library(scales)        # squish
library(glue)
library(dplyr)
library(tidyr)         # drop_na
library(stringr)       # str_c
library(performance)   # r2
library(patchwork)
library(spatialreg)    # spautolm
library(spdep)         # knearneigh, knn2nb, nb2listw
library(lwgeom)        # st_startpoint


# -----------------------------------------------------------------------------
# 1. Map one EH measure  (identical to the FIA version)
# -----------------------------------------------------------------------------
# Colorbar clipped to `trim` quantiles with extremes squished into the end
# colors; NODATA (-9999, from GEE) -> NA, drawn grey. Map "furniture" (states,
# graticule, lon_pts, lat_pts, xlim_alb, ylim_alb) is read from the workspace.
plot_eh_map <- function(a,
                        measure = c("EH_10km", "elev_range_10km",
                                    "EH_delta", "A57_sd"),
                        trim          = c(0.02, 0.98),
                        point_size    = 0.3,
                        linewidth     = 0.6,    # route (line) thickness
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
    geom_sf(data = a, aes(color = .data[[measure]]),
            size = point_size, linewidth = linewidth) +
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
# 2. Partial-residual plot (linear vs quadratic), lm or SAR
# -----------------------------------------------------------------------------
# `data`      : optional data frame holding the predictor `var` aligned with the
#               model's residuals (pass `dat` for SAR so partial residuals line
#               up exactly; otherwise derived from the model).
# `lin_model` : optional pre-fitted model WITHOUT the quadratic term, used for
#               Delta AIC and R^2_lin (needed for SAR, whose call can't be
#               update()-d when built inside a helper).
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
  
  # ── coef table + R^2 (quadratic / original) ────────────────────────────────
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
  
  # ── predictor values, aligned row-for-row with residuals(model) ────────────
  if (!is.null(data)) {
    xvar <- data[[var]]
  } else if (taxa %in% c("bird", "mammal")) {
    xvar <- as.data.frame(model$X, check.names = FALSE)[[var]]
  } else {
    xvar <- stats::model.frame(model)[[var]]
  }
  
  x_seq <- seq(min(xvar, na.rm = TRUE), max(xvar, na.rm = TRUE), length.out = n_seq)
  
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
  
  # ── linear fit: dashed line + Delta AIC + R^2_lin ──────────────────────────
  lin_ok <- FALSE; lin_df <- NULL
  stats_ok <- FALSE; dAIC <- NA_real_; R2_lin <- NA_real_
  if (has_quad) {
    # dashed line = linear regression through the partial residuals
    lin_line <- lm(res ~ x, data = res_df)
    lin_df <- data.frame(x = x_seq,
                         fit = predict(lin_line, newdata = data.frame(x = x_seq)))
    lin_ok <- TRUE
    
    # full model without the quadratic term (passed in for SAR; else refit)
    if (is.null(lin_model)) {
      lin_model <- tryCatch({
        if (taxa %in% c("bird", "mammal")) {
          update(model, stats::as.formula(paste(". ~ . -", quad_name)))
        } else {
          mf    <- stats::model.frame(model)
          resp  <- names(mf)[1]
          preds <- setdiff(names(mf)[-1], quad_name)
          lf <- stats::as.formula(
            paste0("`", resp, "` ~ ", paste0("`", preds, "`", collapse = " + ")))
          lm(lf, data = mf)
        }
      }, error = function(e) NULL)
    }
    if (!is.null(lin_model)) {
      dAIC   <- AIC(lin_model) - AIC(model)        # >0 => quadratic preferred
      R2_lin <- tryCatch(
        if (taxa %in% c("bird", "mammal"))
          cor(lin_model$Y, lin_model$fit$fitted.values)^2
        else unname(r2(lin_model)[[2]]),
        error = function(e) NA_real_)
      stats_ok <- TRUE
    }
  }
  
  # ── inflection point (quadratic only) ──────────────────────────────────────
  inflection <- if (has_quad) -coefs[var] / (2 * coefs[quad_name]) else NA_real_
  infl_in_range <- has_quad && is.finite(inflection) &&
    inflection > min(xvar, na.rm = TRUE) && inflection < max(xvar, na.rm = TRUE)
  
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
  
  # ── annotation (top-left, evenly spaced, all the SAME size) ────────────────
  # Row order: 1) Delta AIC   2) R^2 (quadratic & linear)   3) inflection.
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
# 3. Fit the BBS SAR model for one measure
# -----------------------------------------------------------------------------
# quadratic = TRUE  -> EH + I(EH^2) + covariates   (main model)
# quadratic = FALSE -> EH         + covariates     (linear, for Delta AIC / R^2_lin)
fit_eh_model_bbs <- function(eh_var,
                             response  = c("richness", "MSA"),
                             scale     = "10km",
                             quadratic = TRUE,
                             data      = dat,
                             listw     = lw) {
  response <- match.arg(response)
  covars <- glue("MAT_{scale}_t + MAP_{scale}_t + Tsea_{scale}_t + ",
                 "Psea_{scale}_t + elev_{scale}_t + PopDen_{scale}_t")
  lhs  <- if (response == "richness") "richness" else "log(MSA)"
  quad <- if (quadratic) glue(" + I({eh_var}^2)") else ""
  form <- stats::as.formula(glue("{lhs} ~ {eh_var}{quad} + {covars}"))
  spatialreg::spautolm(formula = form, data = data, listw = listw, family = "SAR")
}


# -----------------------------------------------------------------------------
# 3b. Inflection point vs spatial scale  (the 4th column)
# -----------------------------------------------------------------------------

# "400m" -> 0.4, "1km" -> 1, "20km" -> 20  (for the log x-axis)
scale_to_km <- function(s) {
  num <- as.numeric(gsub("[^0-9.]", "", s))
  if (grepl("km$", s)) num else num / 1000
}

# resolve a glue template like "EH_{scale}_t" for a given scale
.fill <- function(tmpl, scale) as.character(glue::glue(tmpl))

#' Bootstrap the inflection point z = -b1 / (2 b2) of a quadratic EH term.
#'
#' Coefficients are simulated from their sampling distribution (Estimate + Std.
#' Error in the model's coefficient table); the model is NOT refit. By default
#' the inflection is returned on the STANDARDIZED EH scale. Supply center /
#' scale_sd / inv to map it back to original units as inv(center + scale_sd*z):
#' for a predictor built as scale(g(raw)), use center = mean(g(raw)),
#' scale_sd = sd(g(raw)), and inv = the inverse of g (identity, or exp if
#' g = log). Coefficients are read by NAME, so works for lm and spautolm alike.
inflection_boot <- function(model, var, taxa = "bird",
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
  z_sims <- z_sims[is.finite(z_sims)]            # drop b2~0 blow-ups
  
  # a quadratic vertex is only a real inflection if it sits inside the observed
  # (standardized) predictor range; otherwise the relationship is monotonic here
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
#' Fits the quadratic SAR model at each scale (dat_by_scale[[scale]] + the shared
#' listw) and bootstraps the inflection. original_scale = TRUE maps back to the
#' measure's raw units using its transform (identity / log).
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
    
    mdl <- fit_eh_model_bbs(var_t, response = response, scale = scale,
                            quadratic = TRUE, data = d, listw = lw_by_scale[[scale]])
    
    # standardized predictor range -> drives the in-range check
    xr <- range(as.data.frame(mdl$X, check.names = FALSE)[[var_t]], na.rm = TRUE)
    
    if (original_scale) {
      gv <- g(d[[raw]])
      ib <- inflection_boot(mdl, var = var_t, taxa = "bird",
                            center = mean(gv, na.rm = TRUE),
                            scale_sd = sd(gv, na.rm = TRUE),
                            inv = ginv, var_range = xr, nsim = nsim, seed = seed)
    } else {
      ib <- inflection_boot(mdl, var = var_t, taxa = "bird",
                            var_range = xr, nsim = nsim, seed = seed)
    }
    ib$scale <- scale; ib$scale_num <- scale_to_km(scale)
    ib
  })
  res <- do.call(rbind, out)
  res$response <- response_label
  # drop scales where the inflection falls outside the data range (no inflection)
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

#' Build one BBS `dat` for a given scale (mirrors your 10km pipeline).
#' ADJUST the file paths / column names to your setup. Produces the *_t
#' predictors the models expect. The positional cbind() assumes the delta CSV
#' and A57 shapefile share row order (as in your original script).
build_dat_bbs <- function(scale, BBS_table,
                          delta_csv = glue(p_gee("BBS_EH_delta/BBS_Elevation_EHdelta_{scale}.csv")),
                          a57_shp   = glue(p_gee("BBS_A57_EH/BBS_A57_EH_{scale}.shp"))) {
  eledelta <- read.csv(delta_csv) %>%
    dplyr::select(BBSRTSL020, elev_range, EH_delta)
  a57 <- sf::st_read(a57_shp, quiet = TRUE) %>%
    sf::st_drop_geometry() %>%
    dplyr::select(RTENAME, RTENO, A57_sd)
  
  merged <- cbind(eledelta, a57) %>%
    dplyr::select(RTENAME, RTENO, elev_range, EH_delta, A57_sd) %>%
    dplyr::group_by(RTENAME, RTENO) %>%
    dplyr::summarise("elev_range_{scale}" := mean(elev_range),
                     EH_delta = mean(EH_delta),
                     A57_sd   = mean(A57_sd), .groups = "drop")
  
  BBS_table %>%
    dplyr::select(route, RTENAME, RTENO, richness, MSA,
                  dplyr::contains(glue("_{scale}"))) %>%
    dplyr::left_join(merged, by = c("RTENAME", "RTENO")) %>%
    dplyr::filter(EH_delta > 0) %>%
    dplyr::mutate(
      "EH_{scale}_t"     := scale(.data[[glue("EH_{scale}")]])[, 1],
      EH_elevrange_t     = scale(log(.data[[glue("elev_range_{scale}")]]))[, 1],
      EH_delta_t         = scale(log(EH_delta))[, 1],
      A57_sd_t           = scale(A57_sd)[, 1],
      "MAT_{scale}_t"    := scale(.data[[glue("MAT_{scale}")]])[, 1],
      "MAP_{scale}_t"    := scale(.data[[glue("MAP_{scale}")]])[, 1],
      "Tsea_{scale}_t"   := scale(.data[[glue("Tsea_{scale}")]])[, 1],
      "Psea_{scale}_t"   := scale(.data[[glue("Psea_{scale}")]])[, 1],
      "elev_{scale}_t"   := scale(sqrt(.data[[glue("elev_{scale}")]]))[, 1],
      "PopDen_{scale}_t" := scale(log1p(.data[[glue("PopDen_{scale}")]]))[, 1]
    ) %>%
    # complete cases for every modeled column, so the SAR doesn't drop rows and
    # the per-scale listw stays aligned (all four measures share these rows)
    tidyr::drop_na(richness, MSA, dplyr::ends_with("_t"))
}

#' Build the SAR neighbour weights (listw) matching `dat`'s routes.
#' Mirrors your resid_plots -> knn -> nb2listw construction. MUST be rebuilt per
#' scale because each scale's `dat` keeps a different set of routes.
build_lw_bbs <- function(dat, routes, k = 8) {
  key <- function(df) stringr::str_c(df$RTENO, df$RTENAME, sep = "_")
  resid_plots <- routes %>%
    dplyr::filter(stringr::str_c(RTENO, RTENAME, sep = "_") %in% key(dat)) %>%
    dplyr::mutate(RTE = stringr::str_c(RTENO, RTENAME, sep = "_")) %>%
    dplyr::group_by(RTE) %>%
    dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") %>%
    dplyr::arrange(match(RTE, key(dat)))
  coords <- sf::st_coordinates(lwgeom::st_startpoint(resid_plots))
  nb <- spdep::knn2nb(spdep::knearneigh(coords, k = k))
  spdep::nb2listw(nb, style = "W")
}

#' One-call data prep for a scale: returns the modelling `dat` and its matching
#' spatial weights `lw` together (use to build dat_by_scale / lw_by_scale).
prepare_bbs <- function(scale, BBS_table, routes, k = 8) {
  d  <- build_dat_bbs(scale, BBS_table)
  list(dat = d, lw = build_lw_bbs(d, routes, k = k))
}


# -----------------------------------------------------------------------------
# 4. Couple each EH measure's MAP with its richness & MSA partial plots
# -----------------------------------------------------------------------------
# plotting order: elevation range, mean-independent elevation range, EH 64-D, EH 1-D
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

#' Build the BBS figure: one row per EH measure, columns
#'   map | richness partial | MSA partial | inflection-point-vs-scale
#'
#' @param a            sf routes object for the maps (carries the four map columns).
#' @param dat_by_scale named list of modelling data frames, one per scale
#'                     (names must match `scales`, e.g. dat_by_scale[["10km"]]).
#' @param lw_by_scale  named list of spatial weights (listw), one per scale,
#'                     each matching the rows of dat_by_scale[[scale]].
#' @param scales       scales for the inflection~scale column (4th col).
#' @param main_scale   scale used for the map + partial-plot columns (1-3).
#' @param responses    which partial-plot column(s) after the map.
#' @param ann_size     partial-plot annotation text size.
#' @param nsim,seed    inflection bootstrap settings.
#' @param inflection_original_scale  TRUE -> 4th-col y-axis in raw EH units;
#'                     FALSE -> standardized (SD) units.
#' @param ...          forwarded to plot_eh_map() (e.g. trim, drop_extremes).
#' @return a patchwork object (rows = measures, cols = map + partials + inflection).
eh_map_partial_bbs <- function(a, dat_by_scale, lw_by_scale,
                               scales     = c("400m", "1km", "5km", "10km", "20km"),
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
    
    # col 1: map (main scale)
    p_map <- plot_eh_map(a, m$map_col, ...)
    
    # cols 2..k: partial-residual plots (main scale)
    parts <- lapply(responses, function(resp) {
      mdl <- fit_eh_model_bbs(var_t_main, response = resp, scale = main_scale,
                              quadratic = TRUE,  data = dat, listw = lw)
      lin <- fit_eh_model_bbs(var_t_main, response = resp, scale = main_scale,
                              quadratic = FALSE, data = dat, listw = lw)
      sar_partial_plot(mdl, var = var_t_main, scale = main_scale, taxa = "bird",
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
# You need: a `dat` per scale (dat_by_scale), the spatial weights `lw`, the map
# object `a`, and the map furniture (states, graticule, lon_pts, lat_pts,
# xlim_alb, ylim_alb) from your BBS script.
if (T) {
  
  scales <- c("400m", "1km", "5km", "10km", "20km")
  bbsrtsl020_filtered <- st_read(p_clean("bbsroutes_shp/bbsrtsl020_filterd.shp"))
  
  # 1. one `dat` + matching `lw` per scale (each scale keeps a different set of
  #    routes, so the neighbour list MUST be rebuilt per scale)
  prep <- setNames(
    lapply(scales, prepare_bbs, BBS_table = BBS_table, routes = bbsrtsl020_filtered),
    scales)
  dat_by_scale <- lapply(prep, `[[`, "dat")
  lw_by_scale  <- lapply(prep, `[[`, "lw")
  
  # 2. build the map object `a` and the map furniture as in your script.
  
  # single map / single partial plot ------------------------------------------
  plot_eh_map(a, "EH_delta")
  
  mq <- fit_eh_model_bbs("EH_delta_t", "richness", "10km", TRUE,
                         data = dat_by_scale[["10km"]], listw = lw_by_scale[["10km"]])
  ml <- fit_eh_model_bbs("EH_delta_t", "richness", "10km", FALSE,
                         data = dat_by_scale[["10km"]], listw = lw_by_scale[["10km"]])
  sar_partial_plot(mq, "EH_delta_t", scale = "10km", taxa = "bird",
                   quadratic = TRUE, lin_model = ml,
                   x_label = "Elevation range (mean-independent)",
                   y_label = "Partial residual (richness)")
  
  # full 4-column panel: map | richness | MSA | inflection~scale ---------------
  fig <- eh_map_partial_bbs(a, dat_by_scale, lw_by_scale,
                            scales = scales, main_scale = "10km")
  ggsave(p_fig("FigS18_bird_alternative_EH.png"), fig, width = 24, height = 20, dpi = 300)
  
  # standardized (SD) inflection units instead of raw EH units:
  # eh_map_partial_bbs(a, dat_by_scale, lw_by_scale, scales = scales,
  #                    inflection_original_scale = FALSE)
  
}
