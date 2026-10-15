# =============================================================================
# FigS01_phylogenies.R  -  STUB (reviewer: plotting script not located)
# Fig. S1: phylogenies of (a) trees, (b) birds, (c) mammals
#
# Only the tree panel had code (commented out in
# EmbeddingBio/analysis/FIA_speciesabundance-EH.R, lines 165-172); it is
# reproduced below. Add the bird and mammal panels and save the figure.
# TODO(reviewer): the caption says the trees were used "for trait imputation
#   and species-specific analysis": the imputation tree has 864 tips, while the
#   species-level tree has 263. State which one is shown.
# =============================================================================

source(here::here("code", "00_config.R"))
library(tidyverse)
library(ape)
library(ggtree)

# (a) Trees: 263 species, coloured by family
FIA_phylo <- readRDS(p_clean("FIA_phylo_263species.rds"))
# sp.list.sp_ge_fa (species, genus, family) is built in
# 06_species_trait_responses/01_species_EH_trait_models.R via TNRS
# sp.list.sp_ge_fa <- sp.list.sp_ge_fa %>% mutate(tip_label = str_replace(species, " ", "_"))
# groupInfo <- split(sp.list.sp_ge_fa$tip_label, sp.list.sp_ge_fa$family)
# viz_FIA_phylo <- groupOTU(FIA_phylo, groupInfo, group_name = "Family")
# ggtree(viz_FIA_phylo, aes(color = Family), layout = "circular") + geom_tiplab(size = 2)

# (b) Birds: data/clean/BBS_phylo_520species.rds      TODO
# (c) Mammals: mammal_4999species_phylotree-pruner/output.nex (pruned to the
#     modelled species)                                TODO
