# =============================================================================
# Figure 2 - Environmental heterogeneity effects on species richness and
# mean species abundance (MSA) across trees, birds, mammals.
#
# Outputs (in OUT_DIR):
#   Figure2A_maps_3taxa.png             3 rows x 2 cols (maps)        a,b/f,g/k,l
#   Figure2B_partial_scaling_3taxa.png  3 rows x 3 cols (analytics)   c,d,e/h,i,j/m,n,o
#
# Edit a single config to retune:
#   THEME_CONFIG  - text sizes, line widths, legend bar, partial-plot styling
#   COLOR_CONFIG  - per-panel palettes, dot size & alpha, fit-line colors
#   LAYOUT_CONFIG - map extent, output sizes, file paths, preview/cache toggles
#


## 1. Libraries --------------------------------------------------------------
suppressPackageStartupMessages({
  library(scales); library(sf); library(terra); library(patchwork)
  library(ggplot2); library(dplyr); library(tidyr); library(glue)
  library(maps); library(stringr); library(performance); library(paletteer)
})
# spatialreg loaded on-demand inside get_coef_table()
# ragg loaded on-demand inside save_figure()
# ggrastr loaded on-demand inside make_map_panel()


## 2. Paths ------------------------------------------------------------------
DATA_ROOT  <- "D:/BiodiversityEmbedding/data"
OUT_DIR    <- "D:/BiodiversityEmbedding/figures"
CACHE_DIR  <- file.path(OUT_DIR, "cache")
dir.create(OUT_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

PATHS <- list(
  tree   = list(plots = "plot_pts.shp",               table = "FIA_table.csv",
                sar_R = "sar_tree_richness.rds",      sar_M = "sar_tree_MSA.rds"),
  bird   = list(plots = "bbsrtsl020_filterd.shp",     table = "BBS_table.csv",
                sar_R = "sar_bird_richness.rds",      sar_M = "sar_bird_MSA.rds"),
  mammal = list(plots = "camera_array_centroids.shp", table = "mammal_table.csv",
                sar_R = "sar_mammal_richness.rds",    sar_M = "sar_mammal_MSA.rds")
)
PATHS <- lapply(PATHS, function(p) lapply(p, function(f) file.path(DATA_ROOT, f)))

# Per-scale RDS tables. Loaded into .GlobalEnv as <prefix>_table_<scale>
# so helpers can resolve get(glue("..._{scale}")).
SCALE_TABLE_PATHS <- list(
  tree   = setNames(file.path(DATA_ROOT, paste0("FIA_table_",    c("1km","5km","10km","20km","50km"), ".rds")),
                    c("1km","5km","10km","20km","50km")),
  bird   = setNames(file.path(DATA_ROOT, paste0("BBS_table_",    c("400m","1km","5km","10km","20km"), ".rds")),
                    c("400m","1km","5km","10km","20km")),
  mammal = setNames(file.path(DATA_ROOT, paste0("mammal_table_", c("1km","5km","10km","20km","50km"), ".rds")),
                    c("1km","5km","10km","20km","50km"))
)


## 3. Configuration ----------------------------------------------------------

# 3a. THEME_CONFIG - shared text/line/legend styling
#     v1.7: added scaling_legend_text_size and scaling_legend_labels for the
#     in-panel legend on e/j/o panels, and panel_spacing_B for figure_B
#     composite gutter width. These are exposed for easy tuning.
THEME_CONFIG <- list(
  base_size         = 13,
  axis_title_size   = 15, axis_text_size = 13, axis_color = "black",
  panel_bg_col      = "white",
  
  state_fill        = "grey80", state_border = "white", state_border_lwd = 0.3,
  latlon_label_size = 4,        latlon_label_col = "black",
  
  legend_title_size = 15, legend_text_size = 13,
  legend_bar_width  = unit(8, "cm"),
  legend_bar_height = unit(0.6, "cm"),
  legend_position   = "bottom", legend_direction = "horizontal",
  legend_ticks_col  = "black",  legend_frame_col = "black",
  
  partial_point_size    = 0.4, partial_point_alpha   = 0.5, partial_point_col = "grey50",
  partial_line_lwd      = 1.0, partial_ribbon_alpha  = 0.25,
  partial_zero_line_col = "grey40", partial_zero_line_lty = "dashed", partial_zero_line_lwd = 0.4,
  partial_annot_size    = 5,
  inflection_line_lty   = "dotted", inflection_line_col = "black", inflection_line_lwd = 0.6,
  
  scaling_point_size    = 3,    scaling_errorbar_lwd = 0.6, scaling_errorbar_wdth = 0.05,
  scaling_ribbon_alpha  = 0.15, scaling_trend_lwd    = 0.8, scaling_trend_lty     = "dashed",
  
  # v1.7 additions, panel_spacing_B widened in v1.8:
  # scaling_legend_text_size  Size of the in-panel legend text on e/j/o.
  #                           Reduced from 13 to 11 so the legend block fits
  #                           inside the column even when narrowed.
  # scaling_legend_labels     Abbreviated labels for the scaling-plot legend.
  #                           Shorter strings prevent right-edge overflow.
  # panel_spacing_B           Gutter applied symmetrically between rows AND
  #                           columns in figure_B. v1.8: 1.2 cm (was 0.7).
  #                           Increase to 1.5+ for looser spacing,
  #                           decrease to 0.5 for tighter packing.
  scaling_legend_text_size = 11,
  scaling_legend_labels    = c("Species richness"             = "Richness",
                               "log (Mean species abundance)" = "log (MSA)"),
  panel_spacing_B          = unit(2, "cm")
)

# 3b. COLOR_CONFIG - per-panel palette / dot styling and shared fit colors
.PAL_RICH <- rev(as.character(paletteer::paletteer_c("grDevices::Green-Yellow", 30)))
.PAL_ABUN <- as.character(paletteer::paletteer_c("grDevices::Earth", 6))

.MAP_PANEL <- function(palette, title, log_t, cex, alpha) {
  list(palette = palette, legend_title = title, log_transform = log_t,
       limits_probs = c(0.01, 0.99), point_cex = cex, point_alpha = alpha)
}

COLOR_CONFIG <- list(
  tree_richness    = .MAP_PANEL(.PAL_RICH, "Species richness",             FALSE, 0.1, 0.90),
  tree_abundance   = .MAP_PANEL(.PAL_ABUN, "log (Mean species abundance)", TRUE,  0.1, 0.90),
  bird_richness    = .MAP_PANEL(.PAL_RICH, "Species richness",             FALSE, 1.5, 0.70),
  bird_abundance   = .MAP_PANEL(.PAL_ABUN, "log (Mean species abundance)", TRUE,  1.5, 0.70),
  mammal_richness  = .MAP_PANEL(.PAL_RICH, "Species richness",             FALSE, 2.5, 0.90),
  mammal_abundance = .MAP_PANEL(.PAL_ABUN, "log (Mean species abundance)", TRUE,  2.5, 0.90),
  
  partial_richness_col  = "#006400",
  partial_abundance_col = "#8B4513",
  scaling_richness_col  = "#006400",
  scaling_abundance_col = "#8B4513"
)

# 3c. LAYOUT_CONFIG - map projection, output sizes, file paths, perf toggles
.GEO_ASPECT <- 3000000 / 4700000   # ~0.638

LAYOUT_CONFIG <- list(
  albers_crs   = 5070,
  xlim_alb     = c(-2400000, 2300000),
  ylim_alb     = c(  200000, 3200000),
  conus_bbox   = c(xmin = -125, xmax = -66, ymin = 24, ymax = 50),
  
  map_plot_margin = margin(t = 0, r = 0, b = 0, l = 0),
  
  out_width_A_in  = 9,
  out_width_B_in  = 9,
  out_height_A_in = 10,
  out_height_B_in = 10,
  out_dpi         = 300, draft_dpi = 100,
  
  output_file_A = file.path(OUT_DIR, "Figure2A_maps_3taxa.png"),
  output_file_B = file.path(OUT_DIR, "Figure2B_partial_scaling_3taxa.png"),
  draft_file_A  = file.path(OUT_DIR, "Figure2A_maps_3taxa_DRAFT.png"),
  draft_file_B  = file.path(OUT_DIR, "Figure2B_partial_scaling_3taxa_DRAFT.png"),
  
  preview_mode    = TRUE,
  preview_hex_km  = 5,
  rasterise_dpi   = 150,
  use_cache       = TRUE
)

# 3d. Per-taxon spatial scales used by the scaling plots
SCALES <- list(
  tree   = list(num = c(1, 5, 10, 20, 50),  lbl = c("1km","5km","10km","20km","50km")),
  bird   = list(num = c(0.4, 1, 5, 10, 20), lbl = c("400m","1km","5km","10km","20km")),
  mammal = list(num = c(1, 5, 10, 20, 50),  lbl = c("1km","5km","10km","20km","50km"))
)


## 4. Helpers ----------------------------------------------------------------

# 4a. Map scaffolding: states + Albers CRS.
build_map_scaffolding <- function() {
  list(states     = sf::st_as_sf(maps::map("state", plot = FALSE, fill = TRUE)),
       target_crs = sf::st_crs(LAYOUT_CONFIG$albers_crs))
}

# 4b. sf_cache(): RDS-backed cache keyed on source-file mtimes.
sf_cache <- function(cache_name, source_files, build_fn,
                     use_cache = LAYOUT_CONFIG$use_cache,
                     cache_dir = CACHE_DIR) {
  for (f in source_files)
    if (!file.exists(f)) stop(sprintf("sf_cache: missing source file: %s", f))
  
  cache_path <- file.path(cache_dir, paste0(cache_name, ".rds"))
  key_path   <- file.path(cache_dir, paste0(cache_name, ".key"))
  current_key <- digest_key(source_files, build_fn)
  
  if (use_cache && file.exists(cache_path) && file.exists(key_path)) {
    saved_key <- readLines(key_path, warn = FALSE)
    if (length(saved_key) && identical(saved_key, current_key)) {
      cat(sprintf("  [cache HIT ] %s\n", cache_name))
      return(readRDS(cache_path))
    }
  }
  cat(sprintf("  [cache MISS] %s -> building\n", cache_name))
  obj <- build_fn()
  saveRDS(obj, cache_path)
  writeLines(current_key, key_path)
  obj
}

digest_key <- function(source_files, build_fn) {
  mtimes <- vapply(source_files, function(f) as.character(file.info(f)$mtime),
                   character(1))
  body_str <- paste(deparse(build_fn), collapse = "\n")
  paste(c(mtimes, nchar(body_str), substr(body_str, 1, 200)), collapse = "|")
}

# 4c. thin_points_hex(): one representative point per hex cell.
thin_points_hex <- function(x, hex_km = LAYOUT_CONFIG$preview_hex_km,
                            min_n = 5000) {
  if (nrow(x) < min_n) return(x)
  geom_types <- unique(as.character(sf::st_geometry_type(x, by_geometry = FALSE)))
  if (any(grepl("LINE|POLYGON", geom_types))) return(x)
  
  cellsize_m <- hex_km * 1000
  hex <- sf::st_make_grid(x, cellsize = cellsize_m, square = FALSE)
  hex <- sf::st_sf(hex_id = seq_along(hex), geometry = hex)
  joined <- suppressWarnings(sf::st_join(x, hex, join = sf::st_intersects))
  thinned <- joined |>
    dplyr::group_by(hex_id) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::select(-hex_id)
  cat(sprintf("  thinned %d -> %d features (hex %.1f km)\n",
              nrow(x), nrow(thinned), hex_km))
  thinned
}

# 4d. Generic map panel.
make_map_panel <- function(data_sf, value_col, color_cfg, scaffold,
                           rasterise = TRUE,
                           rasterise_dpi = LAYOUT_CONFIG$rasterise_dpi) {
  th <- THEME_CONFIG; lc <- LAYOUT_CONFIG
  if (isTRUE(color_cfg$log_transform)) data_sf[[value_col]] <- log(data_sf[[value_col]])
  clim <- quantile(data_sf[[value_col]], color_cfg$limits_probs, na.rm = TRUE)
  
  color_scale <- scale_color_gradientn(
    colours = color_cfg$palette, limits = clim, oob = scales::squish,
    name = color_cfg$legend_title,
    guide = guide_colorbar(
      title.position = "top", title.hjust = 0,
      barwidth = th$legend_bar_width, barheight = th$legend_bar_height,
      ticks.colour = th$legend_ticks_col, frame.colour = th$legend_frame_col)
  )
  
  is_line <- any(grepl("LINE", as.character(sf::st_geometry_type(data_sf, by_geometry = FALSE))))
  data_layer <- if (is_line) {
    geom_sf(data = data_sf, aes(color = .data[[value_col]]),
            linewidth = color_cfg$point_cex, alpha = color_cfg$point_alpha)
  } else {
    geom_sf(data = data_sf, aes(color = .data[[value_col]]),
            size = color_cfg$point_cex,      alpha = color_cfg$point_alpha)
  }
  
  if (isTRUE(rasterise)) {
    if (!requireNamespace("ggrastr", quietly = TRUE)) {
      warning("make_map_panel(): 'ggrastr' not installed; falling back to vector layer.")
    } else {
      data_layer <- ggrastr::rasterise(data_layer, dpi = rasterise_dpi)
    }
  }
  
  ggplot() +
    geom_sf(data = scaffold$states, fill = th$state_fill,
            color = th$state_border, linewidth = th$state_border_lwd) +
    data_layer + color_scale +
    coord_sf(crs = scaffold$target_crs, xlim = lc$xlim_alb, ylim = lc$ylim_alb,
             expand = FALSE, clip = "off", datum = NA) +
    theme_classic(base_size = th$base_size) +
    theme(
      panel.border       = element_blank(),
      panel.grid         = element_blank(),
      panel.background   = element_rect(fill = th$panel_bg_col),
      axis.line          = element_blank(),
      axis.title         = element_blank(),
      axis.text          = element_blank(),
      axis.ticks         = element_blank(),
      legend.position    = th$legend_position,
      legend.direction   = th$legend_direction,
      legend.margin      = margin(t = 0, b = 0, l = 0, r = 0),
      legend.box.spacing = unit(0, "cm"),
      legend.title       = element_text(size = th$legend_title_size),
      legend.text        = element_text(size = th$legend_text_size),
      plot.margin        = lc$map_plot_margin
    )
}

# 4e. Coefficient-table extractor.
get_coef_table <- function(model) {
  if (inherits(model, "Spautolm")) {
    if (!requireNamespace("spatialreg", quietly = TRUE))
      stop("get_coef_table(): install 'spatialreg' for Spautolm models.")
    s <- tryCatch(summary(model), error = function(e) NULL)
    if (!is.list(s)) {
      s_fn <- getS3method("summary", "Spautolm", optional = TRUE)
      if (is.null(s_fn)) stop("get_coef_table(): summary.Spautolm not found.")
      s <- s_fn(model)
    }
    if (is.matrix(s$Coef))         return(s$Coef)
    if (is.matrix(s$coefficients)) return(s$coefficients)
    stop("get_coef_table(): no coefficient matrix in Spautolm summary.")
  }
  model_like <- c("Spautolm","Sarlm","sarlm","errorsarlm","lagsarlm","sacsarlm","lm","gls","glm")
  if (is.list(model) && !inherits(model, model_like)) {
    for (slot in c("model","fit","sarlm","sar")) {
      if (!is.null(model[[slot]]) && inherits(model[[slot]], model_like))
        return(get_coef_table(model[[slot]]))
    }
  }
  s <- summary(model)
  if (!is.list(s)) stop("get_coef_table(): summary() returned ", class(s)[1])
  if (is.matrix(s$Coef))         return(s$Coef)
  if (is.matrix(s$coefficients)) return(s$coefficients)
  tab <- try(stats::coef(s), silent = TRUE)
  if (inherits(tab, "try-error") || !is.matrix(tab))
    stop("get_coef_table(): no coefficient matrix found.")
  tab
}

# 4f. Partial residual plot.
sar_partial_plot <- function(model, var, scale, taxa,
                             quadratic = FALSE, nsim = 500, n_seq = 200, seed = 123,
                             line_color = "steelblue",
                             x_label = NULL, y_label = "Partial residual") {
  th <- THEME_CONFIG
  
  if (taxa == "bird") {
    data <- get(glue("BBS_table_{scale}")) %>% tidyr::drop_na()
    R2   <- cor(model$Y, model$fit$fitted.values)^2
  } else if (taxa == "mammal") {
    data <- get(glue("mammal_table_{scale}")) %>% tidyr::drop_na() %>%
      dplyr::group_by(Camera_Trap_Array) %>%
      dplyr::summarise(dplyr::across(c(richness_Estimator, abundance, MSA, contains("_t")),
                                     ~ mean(.x, na.rm = TRUE)), .groups = "drop")
    R2 <- cor(model$Y, model$fit$fitted.values)^2
  } else {
    data <- get(glue("FIA_table_{scale}")) %>% tidyr::drop_na()
    R2   <- unname(performance::r2(model)[[2]])
  }
  sar_coef_table <- get_coef_table(model)
  
  coefs   <- coef(model); resids <- residuals(model)
  x_seq   <- seq(min(data[[var]], na.rm = TRUE), max(data[[var]], na.rm = TRUE),
                 length.out = n_seq)
  quad_nm <- paste0("I(", var, "^2)")
  has_q   <- quadratic && (quad_nm %in% names(coefs))
  cnames  <- if (has_q) c(var, quad_nm) else var
  
  fit_line <- if (has_q) coefs[var] * x_seq + coefs[quad_nm] * x_seq^2 else coefs[var] * x_seq
  presid   <- if (has_q) resids + coefs[var] * data[[var]] + coefs[quad_nm] * data[[var]]^2
  else        resids + coefs[var] * data[[var]]
  
  set.seed(seed)
  se_col <- intersect(c("Std. Error", "Std.Error", "SE", "se"),
                      colnames(sar_coef_table))[1]
  if (is.na(se_col)) se_col <- colnames(sar_coef_table)[2]
  missing_cn <- setdiff(cnames, rownames(sar_coef_table))
  if (length(missing_cn) > 0)
    stop(sprintf("sar_partial_plot(): coefficient(s) '%s' not found in table rows.\nAvailable: %s",
                 paste(missing_cn, collapse = ", "),
                 paste(rownames(sar_coef_table), collapse = ", ")))
  coef_sims <- mapply(function(mu, se) rnorm(nsim, mu, se),
                      mu = coefs[cnames], se = sar_coef_table[cnames, se_col])
  if (!is.matrix(coef_sims)) coef_sims <- matrix(coef_sims, ncol = 1)
  X_mat <- if (has_q) cbind(x_seq, x_seq^2) else matrix(x_seq, ncol = 1)
  sims  <- X_mat %*% t(coef_sims)
  
  fit_df <- data.frame(x = x_seq, fit = fit_line,
                       lwr = apply(sims, 1, quantile, 0.025),
                       upr = apply(sims, 1, quantile, 0.975))
  res_df <- data.frame(x = data[[var]], res = presid)
  inflection <- if (has_q) -coefs[var] / (2 * coefs[quad_nm]) else NA_real_
  if (is.null(x_label)) x_label <- var
  
  p <- ggplot() +
    geom_point(data = res_df, aes(x, res),
               size = th$partial_point_size, alpha = th$partial_point_alpha,
               color = th$partial_point_col) +
    geom_ribbon(data = fit_df, aes(x = x, ymin = lwr, ymax = upr),
                fill = line_color, alpha = th$partial_ribbon_alpha) +
    geom_line(data = fit_df, aes(x, fit),
              color = line_color, linewidth = th$partial_line_lwd) +
    geom_hline(yintercept = 0, linetype = th$partial_zero_line_lty,
               color = th$partial_zero_line_col, linewidth = th$partial_zero_line_lwd) +
    labs(x = x_label, y = y_label) +
    theme_classic(base_size = th$base_size) +
    theme(axis.title = element_text(size = th$axis_title_size),
          axis.text  = element_text(size = th$axis_text_size),
          axis.line  = element_line(linewidth = 0.1),
          axis.ticks = element_line(linewidth = 0.1))
  
  if (has_q) {
    has_infl <- !is.na(inflection) && inflection < max(data[[var]]) && inflection > min(data[[var]])
    if (has_infl)
      p <- p + geom_vline(xintercept = inflection, linetype = th$inflection_line_lty,
                          color = th$inflection_line_col, linewidth = th$inflection_line_lwd)
    
    annot_hjust <- -0.08
    p <- p + annotate("text", x = -Inf, y = Inf,
                      label = paste0("R^2 == ", round(R2, 2)),
                      hjust = annot_hjust, vjust = 1.2,
                      size = th$partial_annot_size, color = "black", parse = TRUE)
    if (has_infl)
      p <- p + annotate("text", x = -Inf, y = Inf,
                        label = paste0("Turning~point == ", round(inflection, 2)),
                        hjust = annot_hjust, vjust = 4.0,
                        size = th$partial_annot_size, color = "black", parse = TRUE)
  }
  p
}

# 4g. Inflection bootstrap.
inflection_boot <- function(model, taxa, nsim = 500, seed = 123) {
  set.seed(seed)
  if (taxa == "bird") {
    scale  <- stringr::str_extract(colnames(model$X)[2], "\\d+(km|m)")
    EH_var <- glue("EH_{scale}")
    data   <- get(glue("BBS_table_{scale}")) %>% tidyr::drop_na()
  } else if (taxa == "mammal") {
    scale  <- stringr::str_extract(colnames(model$X)[2], "\\d+(km|m)")
    EH_var <- glue("EH_{scale}")
    data   <- get(glue("mammal_table_{scale}")) %>% tidyr::drop_na() %>%
      dplyr::group_by(Camera_Trap_Array) %>%
      dplyr::summarise(dplyr::across(where(is.numeric),
                                     ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  } else {
    scale  <- stringr::str_extract(colnames(model$model)[2], "\\d+(km|m)")
    EH_var <- glue("EH_{scale}")
    data   <- get(glue("FIA_table_{scale}")) %>% tidyr::drop_na()
  }
  tab <- get_coef_table(model)
  m   <- tab[, "Estimate"]; s <- tab[, "Std. Error"]
  z_sims <- replicate(nsim, { b1 <- rnorm(1, m[2], s[2]); b2 <- rnorm(1, m[3], s[3]); -b1/(2*b2) })
  z_obs  <- -m[2] / (2 * m[3])
  mu  <- mean(data[[EH_var]], na.rm = TRUE); sdv <- sd(data[[EH_var]], na.rm = TRUE)
  data.frame(z = mu + sdv*z_obs, se_z = mu + sdv*sd(z_sims),
             lwr_95 = mu + sdv*quantile(z_sims, 0.025),
             upr_95 = mu + sdv*quantile(z_sims, 0.975))
}

# 4h. Scaling plot: inflection point vs spatial scale, by response.
#     v1.7: shortened legend labels and reduced legend text size so the
#     in-panel legend doesn't overflow narrowed columns.
build_scaling_plot <- function(sar_R, sar_M, taxa, scales_info) {
  th <- THEME_CONFIG; cc <- COLOR_CONFIG
  res_R <- lapply(sar_R, inflection_boot, taxa = taxa) %>% dplyr::bind_rows() %>%
    dplyr::mutate(scale_num = scales_info$num, response = "Species richness")
  res_M <- lapply(sar_M, inflection_boot, taxa = taxa) %>% dplyr::bind_rows() %>%
    dplyr::mutate(scale_num = scales_info$num, response = "log (Mean species abundance)")
  res <- dplyr::bind_rows(res_R, res_M)
  
  # Use abbreviated legend labels (defined in THEME_CONFIG$scaling_legend_labels)
  # while keeping the underlying factor values intact for color matching.
  legend_labels <- th$scaling_legend_labels
  
  ggplot(res, aes(scale_num, z, color = response, fill = response)) +
    geom_ribbon(aes(ymin = lwr_95, ymax = upr_95),
                alpha = th$scaling_ribbon_alpha, color = NA) +
    geom_errorbar(aes(ymin = lwr_95, ymax = upr_95),
                  width = th$scaling_errorbar_wdth, linewidth = th$scaling_errorbar_lwd) +
    geom_point(size = th$scaling_point_size) +
    geom_smooth(method = "lm", se = FALSE,
                linetype = th$scaling_trend_lty, linewidth = th$scaling_trend_lwd) +
    scale_x_log10(breaks = scales_info$num, labels = scales_info$num) +
    scale_color_manual(values = c("Species richness"             = cc$scaling_richness_col,
                                  "log (Mean species abundance)" = cc$scaling_abundance_col),
                       labels = legend_labels,
                       name = "Response") +
    scale_fill_manual( values = c("Species richness"             = cc$scaling_richness_col,
                                  "log (Mean species abundance)" = cc$scaling_abundance_col),
                       labels = legend_labels,
                       name = "Response") +
    labs(x = "Spatial scale (km)", y = "Turning point") +
    theme_classic(base_size = th$base_size) +
    theme(axis.title           = element_text(size = th$axis_title_size),
          axis.text            = element_text(size = th$axis_text_size),
          axis.line            = element_line(linewidth = 0.1),
          axis.ticks           = element_line(linewidth = 0.1),
          legend.position      = c(0, 1),                              # top-left anchor
          legend.justification = c(0, 1),                              # legend's own top-left
          legend.background    = element_rect(fill = NA, color = NA),
          legend.direction     = "vertical",
          legend.margin        = margin(t = 4, r = 4, b = 2, l = 4),
          legend.title         = element_blank(),
          legend.text          = element_text(size = th$scaling_legend_text_size))
}

# 4i. Load per-scale RDS tables into .GlobalEnv so get(glue(...)) resolves.
load_scale_tables <- function(paths_by_taxa,
                              prefix_map = c(tree = "FIA", bird = "BBS", mammal = "mammal"),
                              envir = .GlobalEnv, verbose = TRUE) {
  for (taxa in names(paths_by_taxa)) {
    prefix <- prefix_map[[taxa]]
    for (scale in names(paths_by_taxa[[taxa]])) {
      fp <- paths_by_taxa[[taxa]][[scale]]
      if (!file.exists(fp)) stop(sprintf("Missing RDS for %s @ %s: %s", taxa, scale, fp))
      assign(paste0(prefix, "_table_", scale), readRDS(fp), envir = envir)
      if (verbose) cat(sprintf("  loaded %s_table_%s\n", prefix, scale))
    }
  }
}

# 4j. Save plot via ragg.
save_figure <- function(plot, filename, width_in, height_in, dpi, bg = "white") {
  if (!requireNamespace("ragg", quietly = TRUE)) stop("save_figure(): install 'ragg'.")
  setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  t0 <- Sys.time()
  ragg::agg_png(filename = filename, width = width_in, height = height_in,
                units = "in", res = dpi, bg = bg)
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)
  print(plot); dev.off(); on.exit()
  cat(sprintf("  -> %s (%.1fs, %d x %d px)\n", filename,
              as.numeric(difftime(Sys.time(), t0, units = "secs")),
              as.integer(width_in * dpi), as.integer(height_in * dpi)))
  invisible(filename)
}


## 5. Data loading (cached) --------------------------------------------------
cat("\n[1] Map scaffolding...\n")
scaffold <- build_map_scaffolding()

cat("[2] Trees...\n")
tree_sf <- sf_cache(
  cache_name   = "tree_sf",
  source_files = c(PATHS$tree$plots, PATHS$tree$table),
  build_fn     = function() {
    sf::st_read(PATHS$tree$plots, quiet = TRUE) %>%
      dplyr::left_join(read.csv(PATHS$tree$table), by = "pltID") %>%
      dplyr::select(pltID, dplyr::contains("richness"), dplyr::contains("abundance"),
                    dplyr::contains("MSA"), dplyr::contains("effort")) %>%
      sf::st_crop(sf::st_bbox(LAYOUT_CONFIG$conus_bbox, crs = sf::st_crs(.))) %>%
      sf::st_transform(LAYOUT_CONFIG$albers_crs)
  }
)
sar_tree_R <- readRDS(PATHS$tree$sar_R); sar_tree_M <- readRDS(PATHS$tree$sar_M)

cat("[3] Birds...\n")
bird_sf <- sf_cache(
  cache_name   = "bird_sf",
  source_files = c(PATHS$bird$plots, PATHS$bird$table),
  build_fn     = function() {
    sf::st_read(PATHS$bird$plots, quiet = TRUE) %>%
      dplyr::left_join(read.csv(PATHS$bird$table)[, 1:7], by = c("RTENO","RTENAME")) %>%
      tidyr::drop_na() %>%
      dplyr::mutate(total_niche = richness * MSA)
  }
)
sar_bird_R <- readRDS(PATHS$bird$sar_R); sar_bird_M <- readRDS(PATHS$bird$sar_M)

cat("[4] Mammals...\n")
mammal_sf <- sf_cache(
  cache_name   = "mammal_sf",
  source_files = c(PATHS$mammal$plots, PATHS$mammal$table),
  build_fn     = function() {
    sf::st_read(PATHS$mammal$plots, quiet = TRUE) %>%
      dplyr::rename(Camera_Trap_Array = Cmr_T_A) %>%
      dplyr::left_join(
        read.csv(PATHS$mammal$table) %>% dplyr::group_by(Camera_Trap_Array) %>%
          dplyr::summarise(richness  = mean(richness_Estimator, na.rm = TRUE),
                           shannon   = mean(shannon_Estimator,  na.rm = TRUE),
                           simpson   = mean(simpson_Estimator,  na.rm = TRUE),
                           abundance = mean(abundance,          na.rm = TRUE),
                           MSA       = mean(MSA,                na.rm = TRUE),
                           .groups = "drop"),
        by = "Camera_Trap_Array") %>%
      sf::st_crop(sf::st_bbox(LAYOUT_CONFIG$conus_bbox, crs = sf::st_crs(.))) %>%
      sf::st_transform(LAYOUT_CONFIG$albers_crs)
  }
)
sar_mammal_R <- readRDS(PATHS$mammal$sar_R); sar_mammal_M <- readRDS(PATHS$mammal$sar_M)

cat("\n[4b] Per-scale tables...\n")
load_scale_tables(SCALE_TABLE_PATHS)

# Optional preview thinning - applies only to dense point datasets.
if (isTRUE(LAYOUT_CONFIG$preview_mode)) {
  cat("\n[4c] Preview mode: thinning dense point layers...\n")
  tree_sf_plot   <- thin_points_hex(tree_sf,   LAYOUT_CONFIG$preview_hex_km)
  bird_sf_plot   <- bird_sf
  mammal_sf_plot <- mammal_sf
} else {
  tree_sf_plot   <- tree_sf
  bird_sf_plot   <- bird_sf
  mammal_sf_plot <- mammal_sf
}


## 6. Panel construction -----------------------------------------------------
# v1.8: Y_M shortened to use the "MSA" abbreviation now that the scaling
#       legend on e/j/o uses the same form. Y_R kept as "Residual richness"
#       since "richness" is already the short form. The newline (\n) is
#       no longer needed because the new label fits on a single line.
EH_LBL <- "Environmental heterogeneity \n(standardized)"
Y_R    <- "Residual richness"
Y_M    <- "Residual log (MSA)"

.NO_X_TITLE <- theme(axis.title.x = element_blank())

cat("\n[5] Tree panels...\n")
p_tree_R <- make_map_panel(tree_sf_plot, "richness_10km", COLOR_CONFIG$tree_richness,  scaffold)
p_tree_M <- make_map_panel(tree_sf_plot, "MSA_new_10km",  COLOR_CONFIG$tree_abundance, scaffold)
p1_tree  <- sar_partial_plot(sar_tree_R[[3]], "EH_10km_t", "10km", "tree",
                             quadratic = TRUE, line_color = COLOR_CONFIG$partial_richness_col,
                             x_label = EH_LBL, y_label = Y_R)
p2_tree  <- sar_partial_plot(sar_tree_M[[3]], "EH_10km_t", "10km", "tree",
                             quadratic = TRUE, line_color = COLOR_CONFIG$partial_abundance_col,
                             x_label = EH_LBL, y_label = Y_M)
p3_tree  <- build_scaling_plot(sar_tree_R, sar_tree_M, "tree", SCALES$tree)

cat("[6] Bird panels...\n")
p_bird_R <- make_map_panel(bird_sf_plot, "richness", COLOR_CONFIG$bird_richness,  scaffold)
p_bird_M <- make_map_panel(bird_sf_plot, "MSA",      COLOR_CONFIG$bird_abundance, scaffold)
p1_bird  <- sar_partial_plot(sar_bird_R[[4]], "EH_10km_t", "10km", "bird",
                             quadratic = TRUE, line_color = COLOR_CONFIG$partial_richness_col,
                             x_label = EH_LBL, y_label = Y_R)
p2_bird  <- sar_partial_plot(sar_bird_M[[4]], "EH_10km_t", "10km", "bird",
                             quadratic = TRUE, line_color = COLOR_CONFIG$partial_abundance_col,
                             x_label = EH_LBL, y_label = Y_M)
p3_bird  <- build_scaling_plot(sar_bird_R, sar_bird_M, "bird", SCALES$bird)

cat("[7] Mammal panels...\n")
p_mammal_R <- make_map_panel(mammal_sf_plot, "richness", COLOR_CONFIG$mammal_richness,  scaffold)
p_mammal_M <- make_map_panel(mammal_sf_plot, "MSA",      COLOR_CONFIG$mammal_abundance, scaffold)
p1_mammal  <- sar_partial_plot(sar_mammal_R[[2]], "EH_5km_t", "5km", "mammal",
                               quadratic = TRUE, line_color = COLOR_CONFIG$partial_richness_col,
                               x_label = EH_LBL, y_label = Y_R)
p2_mammal  <- sar_partial_plot(sar_mammal_M[[2]], "EH_5km_t", "5km", "mammal",
                               quadratic = TRUE, line_color = COLOR_CONFIG$partial_abundance_col,
                               x_label = EH_LBL, y_label = Y_M)
p3_mammal  <- build_scaling_plot(sar_mammal_R, sar_mammal_M, "mammal", SCALES$mammal)


## 7. Composite + save -------------------------------------------------------
# v1.7/v1.8: figure_B gets panel.spacing from THEME_CONFIG$panel_spacing_B
# which adds gutters between BOTH columns and rows of the 3x3 grid (patchwork
# applies it symmetrically). Canvas stays at 11 x 10 in but each panel's
# plot area shrinks slightly so axis labels, tags, and the in-panel legend
# on e/j/o no longer crowd against borders. v1.8 widened to 1.2 cm.
cat("\n[8] Assembling figures...\n")

figure_A <- (
  (p_tree_R   | p_tree_M)   /
    (p_bird_R   | p_bird_M)   /
    (p_mammal_R | p_mammal_M)
) +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(tag_levels = list(c("a","b","f","g","k","l"))) &
  theme(plot.tag    = element_text(size = 18, face = "bold"),
        plot.margin = margin(0, 0, 0, 0))

figure_B <- (
  (p1_tree   | p2_tree   | p3_tree)   /
    (p1_bird   | p2_bird   | p3_bird)   /
    (p1_mammal | p2_mammal | p3_mammal)
) +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(tag_levels = list(c("c","d","e","h","i","j","m","n","o"))) &
  theme(plot.tag      = element_text(size = 18, face = "bold"),
        plot.margin   = margin(0, 0, 0, 0),
        panel.spacing = THEME_CONFIG$panel_spacing_B)   # v1.8: row + col gutter

# Drafts (low DPI sanity check)
# cat("[9a] Figure A draft...\n")
# save_figure(figure_A, LAYOUT_CONFIG$draft_file_A,
#             LAYOUT_CONFIG$out_width_A_in, LAYOUT_CONFIG$out_height_A_in,
#             LAYOUT_CONFIG$draft_dpi)
# cat("[9b] Figure B draft...\n")
# save_figure(figure_B, LAYOUT_CONFIG$draft_file_B,
#             LAYOUT_CONFIG$out_width_B_in, LAYOUT_CONFIG$out_height_B_in,
#             LAYOUT_CONFIG$draft_dpi)

# Final renders.
cat("[10a] Figure A final...\n")
save_figure(figure_A, LAYOUT_CONFIG$output_file_A,
            LAYOUT_CONFIG$out_width_A_in, LAYOUT_CONFIG$out_height_A_in,
            LAYOUT_CONFIG$out_dpi)
cat("[10b] Figure B final...\n")
save_figure(figure_B, LAYOUT_CONFIG$output_file_B,
            LAYOUT_CONFIG$out_width_B_in, LAYOUT_CONFIG$out_height_B_in,
            LAYOUT_CONFIG$out_dpi)

cat("Done.\n")


## 8. Test code (small isolated checks) --------------------------------------
# These are intentionally separated from the main body. Run them ad-hoc.

# if (FALSE) {
#   # Test 8 (v1.7): verify scaling legend labels contain the abbreviated forms.
#   stopifnot("log (MSA)" %in% THEME_CONFIG$scaling_legend_labels)
#   stopifnot("Richness"  %in% THEME_CONFIG$scaling_legend_labels)
#   cat("v1.7 abbreviated-legend test OK\n")
#
#   # Test 9 (v1.8): verify Y_M uses MSA and panel_spacing_B was widened.
#   stopifnot(grepl("MSA", Y_M))
#   stopifnot(!grepl("mean species abundance", Y_M, ignore.case = TRUE))
#   stopifnot(as.numeric(THEME_CONFIG$panel_spacing_B) >= 1.0)  # cm
#   cat("v1.8 MSA-label + wider-spacing test OK\n")
# } .png