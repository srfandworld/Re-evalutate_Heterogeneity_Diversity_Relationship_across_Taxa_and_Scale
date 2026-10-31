# =============================================================================
# FigS01_phylogenies.R
# Fig. S1: phylogenies of (a) tree, (b) bird and (c) mammal species used in the
# species-specific analysis (Fig. 3, Table S5); circular layout, branches
# coloured by family (trees) or order (birds, mammals); panels stacked vertically.
#
# Source   : EmbeddingBio/analysis/sfigure1_phylogenies.R (2026-10-30; paths only
#            changed)
# Inputs   : data/clean/FIA_phylo_263species.rds (V.PhyloMaker2 tree, built in
#              06_species_trait_responses/01_species_EH_trait_models.R)
#            data/clean/FIA_phylo_263species_family.csv  genus -> family lookup
#              (from V.PhyloMaker2's tips.info.TPL; Gymnocladus = Fabaceae added)
#            data/clean/BBS_phylo_520species.rds, BBS_species2.csv, BBS_data.csv
#            data/clean/mammal_RA_weighted.csv, data/raw/Mammals/ssusa_finalsequences.csv
#            HDR_EXT_DIR/Traits_data/mammal_4999species_phylotree-pruner/output.nex,
#              mammal_species_122.csv
# Outputs  : output/figures/FigureS1_phylogenies.png (a/b/c stacked) and
#            FigureS1a_tree_phylogeny.png, FigureS1b_bird_phylogeny.png,
#            FigureS1c_mammal_phylogeny.png
# Species  : 263 trees, 267 birds (modelled species present in the McTavish
#            tree), 17 mammals
# Packages : ape (CRAN), ggtree (Bioconductor: renv::install("bioc::ggtree"))
# =============================================================================

source(here::here("code", "00_config.R"))
library(dplyr)
library(stringr)
library(ape)
library(ggtree)
library(ggplot2)
library(patchwork)

out_dir <- DIR_FIG
min_obs <- 100   # minimum observations per species in the species-level models


##########################1. trees (263 species)###############################################
FIA_phylo <- readRDS(p_clean("FIA_phylo_263species.rds"))
#family of each tip: genus -> family from V.PhyloMaker2's tips.info.TPL (the
#taxonomy that built FIA_phylo), saved once as a small lookup table;
#Gymnocladus (absent from tips.info.TPL) was added by hand as Fabaceae
tree_tax <- read.csv(p_clean("FIA_phylo_263species_family.csv")) %>%
  filter(tip_label %in% FIA_phylo$tip.label)
if (any(is.na(tree_tax$family))) warning("No family for: ", paste(tree_tax$tip_label[is.na(tree_tax$family)], collapse = ", "))
groupInfo_tree <- split(tree_tax$tip_label, tree_tax$family)
viz_FIA_phylo  <- groupOTU(FIA_phylo, groupInfo_tree, group_name = "Family")


##########################2. birds (modelled species)##########################################
BBS_phylo    <- readRDS(p_clean("BBS_phylo_520species.rds"))   #tip names = eBird taxa
BBS_species2 <- read.csv(p_clean("BBS_species2.csv"))

bird_modelled <- read.csv(p_clean("BBS_data.csv")) %>%
  distinct(route, aou) %>%
  count(aou, name = "n_obs") %>%
  filter(n_obs > min_obs)

bird_tax <- bird_modelled %>%
  left_join(BBS_species2 %>% dplyr::select(aou, ebird_species, ORDER), by = "aou") %>%
  mutate(tip_label = str_replace_all(ebird_species, " ", "_")) %>%
  filter(tip_label %in% BBS_phylo$tip.label) %>%
  distinct(tip_label, .keep_all = TRUE)
BBS_phylo_pruned <- keep.tip(BBS_phylo, bird_tax$tip_label)
groupInfo_bird   <- split(bird_tax$tip_label, bird_tax$ORDER)
viz_BBS_phylo    <- groupOTU(BBS_phylo_pruned, groupInfo_bird, group_name = "Order")


##########################3. mammals (modelled species)########################################
mammal_trees <- ape::read.nexus(p_trait("mammal_4999species_phylotree-pruner", "output.nex"))
species_122  <- read.csv(p_trait("mammal_species_122.csv"))   #SNAPSHOT -> VertLife names

mammal_modelled <- read.csv(p_clean("mammal_RA_weighted.csv")) %>%
  mutate(year_CTA = str_c(Year, Camera_Trap_Array, sep = "_")) %>%
  group_by(SCIENTIFIC_NAME) %>%
  summarise(n_obs = n_distinct(year_CTA), .groups = "drop") %>%
  filter(n_obs > min_obs)

mammal_order <- read.csv(p_raw("Mammals", "ssusa_finalsequences.csv")) %>%
  filter(Class == "Mammalia") %>%
  transmute(SCIENTIFIC_NAME = paste(Genus, Species), Order) %>%
  distinct(SCIENTIFIC_NAME, .keep_all = TRUE)

mammal_tax <- mammal_modelled %>%
  left_join(species_122 %>% dplyr::select(SCIENTIFIC_NAME, Vertlife_NAME), by = "SCIENTIFIC_NAME") %>%
  left_join(mammal_order, by = "SCIENTIFIC_NAME") %>%
  mutate(tip_label = str_replace_all(coalesce(Vertlife_NAME, SCIENTIFIC_NAME), " ", "_")) %>%
  filter(tip_label %in% mammal_trees[[1]]$tip.label)
mammal_phylo    <- keep.tip(mammal_trees[[1]], mammal_tax$tip_label)
groupInfo_mamm  <- split(mammal_tax$tip_label, mammal_tax$Order)
viz_mammal_phylo <- groupOTU(mammal_phylo, groupInfo_mamm, group_name = "Order")


##########################4. Figure S1#########################################################
#circular tree coloured by group; the backbone branches that groupOTU() cannot
#assign to one group ("0") are drawn grey and left out of the legend
plot_tree <- function(phy, group, tiplab_size = 2, title = NULL) {
  grp <- sort(setdiff(unique(as.character(attr(phy, group))), "0"))
  ggtree(phy, aes(color = .data[[group]]), layout = "circular") +
    geom_tiplab(size = tiplab_size, show.legend = FALSE) +
    scale_color_discrete(breaks = grp, name = group, na.value = "grey60") +
    labs(title = title) +
    theme(legend.position = "right",
          legend.text  = element_text(size = 8),
          legend.title = element_text(size = 10),
          plot.title   = element_text(face = "bold", size = 16)) +
    guides(color = guide_legend(ncol = 2, override.aes = list(linewidth = 1)))
}

p_tree   <- plot_tree(viz_FIA_phylo,    "Family", tiplab_size = 1.6, title = "a")
p_bird   <- plot_tree(viz_BBS_phylo,    "Order",  tiplab_size = 1.6, title = "b")
p_mammal <- plot_tree(viz_mammal_phylo, "Order",  tiplab_size = 3.0, title = "c")

ggsave(file.path(out_dir, "FigureS1a_tree_phylogeny.png"),   p_tree,   width = 11, height = 9, dpi = 300)
ggsave(file.path(out_dir, "FigureS1b_bird_phylogeny.png"),   p_bird,   width = 11, height = 9, dpi = 300)
ggsave(file.path(out_dir, "FigureS1c_mammal_phylogeny.png"), p_mammal, width = 8,  height = 6, dpi = 300)
#panels stacked vertically: a (trees) / b (birds) / c (mammals)
ggsave(file.path(out_dir, "FigureS1_phylogenies.png"),
       p_tree / p_bird / p_mammal + plot_layout(ncol = 1, heights = c(1, 1, 0.7)),
       width = 11, height = 26, dpi = 300)

cat(sprintf("Fig. S1: %d tree, %d bird, %d mammal species\n",
            Ntip(FIA_phylo), Ntip(BBS_phylo_pruned), Ntip(mammal_phylo)))
