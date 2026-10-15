# =============================================================================
# FigS05_niche_env_layers.R  -  STUB (reviewer: no code located)
# Fig. S5: maps of the 14 environmental layers that define niche space
#
# Input    : data/clean/envstack_for_FIAniche.tif (04_traits_phylogeny_niche/05)
# Suggested: plot each layer of the stack (aggregated for display) in a 4 x 4
#            panel with terra::plot() or tidyterra + ggplot2, and save to
#            output/figures/FigS5_niche_env_layers.png.
# NOTE     : the caption says maps are global "for visualization purposes";
#            check that this matches the stack extent (see the CRS TODO in
#            05_niche_space_PCA.R).
# =============================================================================

source(here::here("code", "00_config.R"))
library(terra)

env_stack <- rast(p_clean("envstack_for_FIAniche.tif"))
# TODO: plotting code
