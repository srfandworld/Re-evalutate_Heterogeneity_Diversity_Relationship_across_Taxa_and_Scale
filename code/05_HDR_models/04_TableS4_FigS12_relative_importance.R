# =============================================================================
# 04_TableS4_FigS12_relative_importance.R
# LMG relative importance of predictors (relaimpo) -> donut charts (Fig. S12)
#
# Source   : EmbeddingBio/analysis/sfigure_donut_relaimpor.R (paths only changed)
# Inputs   : data/clean/models/sar_*.rds
# Outputs  : donut plots printed to the device (no files written)
#
# TODO(reviewer C2 / Table S4): this script is NOT yet a full reproduction of
#   Table S4 / Fig. S12:
#   1. relaimpo::boot.relimp() needs an lm object. The tree models are lm, but
#      the bird/mammal objects are spautolm - refit the same formula with lm()
#      for the LMG decomposition (and say so in the Methods).
#   2. Bird "MSA" block bootstraps sar_tree_MSA[[3]] (copy-paste slip) - should
#      be sar_bird_MSA[[...]].
#   3. Mammal block never computes `ci` before using ci@lmg.
#   4. Table S4 needs every scale x taxon x response; add a loop that writes
#      output/tables/TableS4_LMG.csv.
#   5. Reported tree values (22.0-23.8% EH) match the OLD non-iNEXT models;
#      the _iNEXT models give ~6.4% at 10 km. Update text after re-running.
# =============================================================================

source(here::here("code", "00_config.R"))
library(tidyverse)
library(relaimpo)

library(ggplot2)
library(dplyr)


tree_donut_relative_importance <- function(ri_vals,scale="10km"){
  names(ri_vals) <- c(
    paste0("EH_", scale, "_t"),
    paste0("I(EH_", scale, "_t^2)"),
    paste0("MAT_", scale, "_t"),
    paste0("MAP_", scale, "_t"),
    paste0("soilph_", scale, "_t"),
    paste0("soilcec_", scale, "_t"),
    paste0("elev_mean_", scale, "_t")
  )
  
  # --- 2. Merge EH terms and build data frame ---
  ri_df <- tibble(
    predictor = names(ri_vals),
    importance = as.numeric(ri_vals)
  ) %>%
    mutate(predictor = case_when(
      predictor %in% c(
        paste0("EH_", scale, "_t"),
        paste0("I(EH_", scale, "_t^2)")
      ) ~ "EH",
      predictor == paste0("MAT_", scale, "_t") ~ "MAT",
      predictor == paste0("MAP_", scale, "_t") ~ "MAP",
      predictor == paste0("soilph_", scale, "_t") ~ "SoilPH",
      predictor == paste0("soilcec_", scale, "_t") ~ "SoilCEC",
      predictor == paste0("elev_mean_", scale, "_t") ~ "Elevation",
      TRUE ~ predictor
    )) %>%
    group_by(predictor) %>%
    summarise(importance = sum(importance), .groups = "drop") %>%
    mutate(
      pct   = importance / sum(importance),
      label = paste0(predictor, "\n", scales::percent(pct, accuracy = 0.1)),
      ymax  = cumsum(pct),
      ymin  = lag(ymax, default = 0),
      ymid  = (ymin + ymax) / 2
    )
  
  # --- 3. Colour palette ---
  my_colors <- c(
    "EH" = "#c52a20",
    "MAT" = "#508ab2",
    "MAP" = "#a1d0c7",
    "SoilPH" = "#f0bb41",
    "SoilCEC" = "#d5b282",
    "Elevation" = "#b36a6f"
  )
  
  # --- 4. Draw the donut ---
  ggplot(ri_df) +
    geom_rect(
      aes(xmin = 3, xmax = 4.5, ymin = ymin, ymax = ymax, fill = predictor),
      color = "white", linewidth = 0.5
    ) +
    geom_text(
      aes(x = 5, y = ymid, label = label),
      size = 3, lineheight = 0.9, hjust = 0
    ) +
    scale_fill_manual(values = my_colors) +
    coord_polar(theta = "y", start = 0) +
    xlim(2, 6.5) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10)
    )
}

sar_tree_richness <- readRDS(p_model("sar_tree_richness_iNEXT.rds"))
sar_tree_MSA <- readRDS(p_model("sar_tree_MSA_iNEXT.rds"))

#richness
#relative importance
bootresults <- relaimpo::boot.relimp(sar_tree_richness[[3]],rela=T) 
ci <- relaimpo::booteval.relimp(bootresults, norank=T)
ci
ri_vals <- ci@lmg  # named vector of LMG proportions
tree_donut_relative_importance(ri_vals,scale="10km")

#MSA
bootresults <- relaimpo::boot.relimp(sar_tree_MSA[[3]],rela=T) 
ci <- relaimpo::booteval.relimp(bootresults, norank=T)
ci
ri_vals <- ci@lmg  # named vector of LMG proportions
tree_donut_relative_importance(ri_vals,scale="10km")





bird_donut_relative_importance <- function(ri_vals, scale = "10km") {
  
  names(ri_vals) <- c(
    paste0("EH_", scale, "_t"),
    paste0("I(EH_", scale, "_t^2)"),
    paste0("MAT_", scale, "_t"),
    paste0("MAP_", scale, "_t"),
    paste0("Tsea_", scale, "_t"),
    paste0("Psea_", scale, "_t"),
    paste0("elev_", scale, "_t"),
    paste0("PopDen_", scale, "_t")
  )
  
  # --- 2. Merge EH terms and build data frame ---
  ri_df <- tibble(
    predictor = names(ri_vals),
    importance = as.numeric(ri_vals)
  ) %>%
    mutate(
      predictor = case_when(
        predictor %in% c(
          paste0("EH_", scale, "_t"),
          paste0("I(EH_", scale, "_t^2)")
        ) ~ "EH",
        predictor == paste0("MAT_", scale, "_t") ~ "MAT",
        predictor == paste0("MAP_", scale, "_t") ~ "MAP",
        predictor == paste0("Tsea_", scale, "_t") ~ "Tsea",
        predictor == paste0("Psea_", scale, "_t") ~ "Psea",
        predictor == paste0("elev_", scale, "_t") ~ "Elevation",
        predictor == paste0("PopDen_", scale, "_t") ~ "PopDen",
        TRUE ~ predictor
      )
    ) %>%
    group_by(predictor) %>%
    summarise(importance = sum(importance), .groups = "drop") %>%
    mutate(
      pct   = importance / sum(importance),
      label = paste0(
        predictor,
        "\n",
        scales::percent(pct, accuracy = 0.1)
      ),
      ymax = cumsum(pct),
      ymin = lag(ymax, default = 0),
      ymid = (ymin + ymax) / 2
    )
  
  # --- 3. Colour palette ---
  my_colors <- c(
    "EH"        = "#c52a20",
    "MAT"       = "#508ab2",
    "MAP"       = "#a1d0c7",
    "Tsea"      = "#f0bb41",
    "Psea"      = "#d5b282",
    "Elevation" = "#b36a6f",
    "PopDen"    = "#639922"
  )
  
  # --- 4. Draw the donut ---
  ggplot(ri_df) +
    geom_rect(
      aes(
        xmin = 3, xmax = 4.5,
        ymin = ymin, ymax = ymax,
        fill = predictor
      ),
      color = "white",
      linewidth = 0.5
    ) +
    geom_text(
      aes(x = 5, y = ymid, label = label),
      size = 3,
      lineheight = 0.9,
      hjust = 0
    ) +
    scale_fill_manual(values = my_colors) +
    coord_polar(theta = "y", start = 0) +
    xlim(2, 6.5) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10)
    )
}


sar_bird_richness <- readRDS(p_model("sar_bird_richness.rds"))
sar_bird_MSA <- readRDS(p_model("sar_bird_MSA.rds"))

#richness
bootresults <- relaimpo::boot.relimp(sar_bird_richness[[3]],rela=T) 
ci <- relaimpo::booteval.relimp(bootresults, norank=T)
ci
ri_vals <- ci@lmg  # named vector of LMG proportions
bird_donut_relative_importance(ri_vals,scale="10km")
#MSA
bootresults <- relaimpo::boot.relimp(sar_tree_MSA[[3]],rela=T) 
ci <- relaimpo::booteval.relimp(bootresults, norank=T)
ci
ri_vals <- ci@lmg  # named vector of LMG proportions
bird_donut_relative_importance(ri_vals,scale="10km")





mammal_donut_relative_importance <- function(ri_vals, scale = "1km") {
  
  names(ri_vals) <- c(
    paste0("EH_", scale, "_t"),
    paste0("I(EH_", scale, "_t^2)"),
    paste0("elev_", scale, "_t"),
    paste0("PopDen_", scale, "_t"),
    paste0("NPP_", scale, "_t"),
    paste0("MAP_", scale, "_t"),
    paste0("MAT_", scale, "_t")
  )
  
  # --- 2. Merge EH terms and build data frame ---
  ri_df <- tibble(
    predictor = names(ri_vals),
    importance = as.numeric(ri_vals)
  ) %>%
    mutate(
      predictor = case_when(
        predictor %in% c(
          paste0("EH_", scale, "_t"),
          paste0("I(EH_", scale, "_t^2)")
        ) ~ "EH",
        predictor == paste0("MAT_", scale, "_t") ~ "MAT",
        predictor == paste0("MAP_", scale, "_t") ~ "MAP",
        predictor == paste0("NPP_", scale, "_t") ~ "NPP",
        predictor == paste0("PopDen_", scale, "_t") ~ "PopDen",
        predictor == paste0("elev_", scale, "_t") ~ "Elevation",
        TRUE ~ predictor
      )
    ) %>%
    group_by(predictor) %>%
    summarise(importance = sum(importance), .groups = "drop") %>%
    mutate(
      pct   = importance / sum(importance),
      label = paste0(
        predictor,
        "\n",
        scales::percent(pct, accuracy = 0.1)
      ),
      ymax = cumsum(pct),
      ymin = lag(ymax, default = 0),
      ymid = (ymin + ymax) / 2
    )
  
  # --- 3. Colour palette ---
  my_colors <- c(
    "EH"        = "#c52a20",
    "MAT"       = "#508ab2",
    "MAP"       = "#a1d0c7",
    "NPP"       = "#f0bb41",
    "Elevation" = "#b36a6f",
    "PopDen"    = "#639922"
  )
  
  # --- 4. Draw the donut ---
  ggplot(ri_df) +
    geom_rect(
      aes(
        xmin = 3, xmax = 4.5,
        ymin = ymin, ymax = ymax,
        fill = predictor
      ),
      color = "white",
      linewidth = 0.5
    ) +
    geom_text(
      aes(x = 5, y = ymid, label = label),
      size = 3,
      lineheight = 0.9,
      hjust = 0
    ) +
    scale_fill_manual(values = my_colors) +
    coord_polar(theta = "y", start = 0) +
    xlim(2, 6.5) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10)
    )
}
#richness
ri_vals <- ci@lmg  # named vector of LMG proportions
mammal_donut_relative_importance(ri_vals,scale="50km")
#MSA
ri_vals <- ci@lmg  # named vector of LMG proportions
mammal_donut_relative_importance(ri_vals,scale="20km")
