# =============================================================================
# FigS16_total_community_abundance.R
# Fig. S16: standardized effect of EH on log total community abundance
# (= log(richness x MSA)) across spatial scales for trees, birds and mammals.
#
# Source   : EmbeddingBio/analysis/sfigure16_highEH_highTotalAbundance.R
#            (2026-10-30; paths only changed)
# SAR models (spautolm) on the same rows, covariates and spatial weights as the
# richness / MSA models (05_HDR_models/01_fit_HDR_models.R); EH enters linearly.
# Inputs   : data/clean/model_data/{FIA,BBS,mammal}_model_<scale>.csv
#            data/clean/FIA_plots.shp, bbsroutes_shp/bbsrtsl020_filterd.shp,
#              camera_array_centroids.shp
# Outputs  : output/figures/FigureS16_EH_total_abundance.png
#            output/tables/FigureS16_EH_total_abundance.csv (estimates, CI, P, n)
#            data/clean/models/sar_models_figS16.rds (cache; tree fits ~15 min)
# Results (2026-10-30): positive for birds at every scale; trees positive at
#   5-50 km (NS at 1 km); mammals positive at 1 and 5 km only (NS at 10-50 km).
# TODO(reviewer, caption): "Across all three taxa, EH had a positive effect" is
#   too general - qualify for mammals (10-50 km) and trees (1 km).
# =============================================================================

source(here::here("code", "00_config.R"))
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(glue)
library(sf)
library(lwgeom)       # st_startpoint() (bird route start points)
library(spdep)        # knearneigh(), knn2nb(), nb2listw()
library(spatialreg)   # spautolm()
library(ggplot2)

data_dir   <- p_clean("model_data")
out_dir    <- DIR_FIG
cache_file <- p_model("sar_models_figS16.rds")

scales_tree   <- c("1km", "5km", "10km", "20km", "50km")
scales_bird   <- c("400m", "1km", "5km", "10km", "20km")
scales_mammal <- c("1km", "5km", "10km", "20km", "50km")


##########################1. spatial weights (same builders as 01_3taxa_modeling.R)###########
plots_sf_tree       <- st_read(p_clean("FIA_plots.shp"), quiet = TRUE) %>%
  st_transform(crs = 5070)
routes_sf_bird      <- st_read(p_clean("bbsroutes_shp", "bbsrtsl020_filterd.shp"), quiet = TRUE)
centroids_sf_mammal <- read_sf(p_clean("camera_array_centroids.shp")) %>%
  rename(Camera_Trap_Array = Cmr_T_A)

# trees: FIA plot locations
build_lw_tree <- function(dat, plots_sf, k = 8) {
  geom <- dat %>%
    dplyr::select(focal_pltID) %>%
    left_join(plots_sf, by = join_by(focal_pltID == pltID)) %>%
    st_as_sf()
  stopifnot(nrow(geom) == nrow(dat))   # listw must align with data row-for-row
  coords <- st_coordinates(geom)
  nb2listw(knn2nb(knearneigh(coords, k = k)), style = "W")
}

# birds: route lines united per route, start point as coordinate
build_lw_bird <- function(dat, routes_sf, k = 8) {
  dat_key <- str_c(dat$RTENO, dat$RTENAME, sep = "_")
  geom <- routes_sf %>%
    mutate(RTE = str_c(RTENO, RTENAME, sep = "_")) %>%
    filter(RTE %in% dat_key) %>%
    group_by(RTE) %>%
    summarise(geometry = st_union(geometry), .groups = "drop") %>%
    arrange(match(RTE, dat_key))
  stopifnot(nrow(geom) == nrow(dat))
  coords <- st_coordinates(lwgeom::st_startpoint(geom))
  nb2listw(knn2nb(knearneigh(coords, k = k)), style = "W")
}

# mammals: camera-trap-array centroids (one row per array)
build_lw_mammal <- function(dat, centroids_sf, k = 8) {
  stopifnot(all(dat$Camera_Trap_Array %in% centroids_sf$Camera_Trap_Array))
  geom <- dat %>%
    dplyr::select(Camera_Trap_Array) %>%
    left_join(centroids_sf, by = "Camera_Trap_Array") %>%
    st_as_sf()
  stopifnot(nrow(geom) == nrow(dat))
  coords <- st_coordinates(geom)
  nb2listw(knn2nb(knearneigh(coords, k = k)), style = "W")
}


##########################2. SAR models of log total abundance################################
f_tree <- function(scale){
  dat <- read.csv(file.path(data_dir, glue("FIA_model_{scale}.csv")))
  lw  <- build_lw_tree(dat, plots_sf_tree)
  form <- as.formula(glue("log(qD * MSA_new_{scale}) ~ EH_{scale}_t + MAT_{scale}_t + MAP_{scale}_t + ",
                          "soilph_{scale}_t + soilcec_{scale}_t + elev_mean_{scale}_t"))
  spautolm(formula = form, data = dat, listw = lw, method = "LU", family = "SAR")
}

f_bird <- function(scale){
  dat <- read.csv(file.path(data_dir, glue("BBS_model_{scale}.csv")))
  lw  <- build_lw_bird(dat, routes_sf_bird)
  form <- as.formula(glue("log(richness * MSA) ~ EH_{scale}_t + MAT_{scale}_t + MAP_{scale}_t + ",
                          "Tsea_{scale}_t + Psea_{scale}_t + elev_{scale}_t + PopDen_{scale}_t"))
  spautolm(formula = form, data = dat, listw = lw, family = "SAR")
}

f_mammal <- function(scale){
  dat <- read.csv(file.path(data_dir, glue("mammal_model_{scale}.csv")))
  lw  <- build_lw_mammal(dat, centroids_sf_mammal)
  form <- as.formula(glue("log(richness_Estimator * MSA) ~ EH_{scale}_t + elev_{scale}_t + PopDen_{scale}_t + ",
                          "NPP_{scale}_t + Prec_{scale}_t + Temp_{scale}_t"))
  spautolm(formula = form, data = dat, listw = lw, family = "SAR")
}

if (file.exists(cache_file)) {
  sar_models <- readRDS(cache_file)
} else {
  sar_models <- list(
    tree   = setNames(lapply(scales_tree,   f_tree),   scales_tree),
    bird   = setNames(lapply(scales_bird,   f_bird),   scales_bird),
    mammal = setNames(lapply(scales_mammal, f_mammal), scales_mammal)
  )
  saveRDS(sar_models, cache_file)
}


##########################3. EH coefficient per taxon and scale###############################
# Extract the EH coefficient, SE, and p-value from a fitted SAR model
extract_eh <- function(model, taxon, scale) {
  ct <- summary(model)$Coef
  eh_row <- grep("^EH", rownames(ct), value = TRUE)[1]
  est <- ct[eh_row, "Estimate"]
  se  <- ct[eh_row, "Std. Error"]
  data.frame(taxon, scale,
             estimate = est, se = se,
             lwr = est - 1.96 * se, upr = est + 1.96 * se,
             p = ct[eh_row, "Pr(>|z|)"],
             significant = ct[eh_row, "Pr(>|z|)"] < 0.05,
             n = length(model$Y))
}

eh_df <- bind_rows(
  map2_dfr(sar_models$tree,   scales_tree,   ~ extract_eh(.x, "Tree",   .y)),
  map2_dfr(sar_models$bird,   scales_bird,   ~ extract_eh(.x, "Bird",   .y)),
  map2_dfr(sar_models$mammal, scales_mammal, ~ extract_eh(.x, "Mammal", .y))
) %>%
  mutate(
    scale = factor(scale, levels = c("400m","1km","5km","10km","20km","50km")),
    taxon = factor(taxon, levels = c("Tree","Bird","Mammal"))
  )
write.csv(eh_df, p_tab("FigureS16_EH_total_abundance.csv"), row.names = FALSE)


##########################4. Figure S16########################################################
p <- ggplot(eh_df, aes(x = estimate, y = scale, color = taxon)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(xmin = lwr, xmax = upr), orientation = "y", width = 0.25, linewidth = 0.9) +
  geom_point(aes(shape = significant), size = 3.5, fill = "white", stroke = 1) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                     labels = c(`TRUE` = "p < 0.05", `FALSE` = "non-significant"),
                     name = NULL) +
  scale_color_manual(values = c(Tree = "#4a7c3f", Bird = "#3b7ba8", Mammal = "#d98841"),
                     guide = "none") +
  facet_wrap(~ taxon, ncol = 1, scales = "free_y") +
  labs(x = "Standardized effect of EH on log(total community abundance)",
       y = "Spatial scale") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "top")

ggsave(file.path(out_dir, "FigureS16_EH_total_abundance.png"), p, width = 6, height = 8, dpi = 300)
p
