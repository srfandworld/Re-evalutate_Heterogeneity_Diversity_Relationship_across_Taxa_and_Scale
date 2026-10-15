# =============================================================================
# 03_mammal_traits_imputation.R
# Mammal traits (PanTHERIA + EltonTraits), VertLife names, phylogenetic imputation
#
# Source   : Code/Mammal/02_Mammal_traits_phylogeny.R, sections 1-3.
#            Section 4 (trait PCoA -> mammal_traits_space.csv, Fig. S4) moved to
#            06_species_trait_responses/01_species_EH_trait_models.R, because it
#            needs the 17-species trait table built there. This removes the
#            circular dependency flagged by the reviewer (C5).
# Inputs   : data/clean/mammal_RA_weighted.csv   (03_community_metrics/04)
#            HDR_EXT_DIR/Traits_data/{PanTHERIA, EltonTraits 1.0, vertlife_taxonomies.csv}
#            HDR_EXT_DIR/Traits_data/mammal_4999species_phylotree-pruner/output.nex
#              (VertLife phylotree pruner, run on mammal_species_5416.csv)
# Outputs  : HDR_EXT_DIR/Traits_data/mammal_species_122.csv, mammal_species_5416.csv
#            data/clean/mammal_imputed_traits.csv
#
# TODO(reviewer M4): Methods describe imputation across 4,999 species and
#   10 trees, but mammal_imputed_traits.csv is not used by the Fig. 3 / Table S5
#   models (they use raw PanTHERIA/EltonTraits). Either use it or revise Methods.
# =============================================================================

source(here::here("code", "00_config.R"))

# ---- Libraries --------------------------------------------------------------
library(tidyverse)   # dplyr, tidyr, stringr, ggplot2, purrr, readr, tibble
library(ape)         # read.nexus()
library(Rphylopars)  # phylopars()
library(abind)       # abind()
library(FD)          # gowdis()
library(grid)        # unit(), arrow()  (used in the PCoA biplot)


traits_dir <- p_trait()   # raw trait sources (HDR_EXT_DIR/Traits_data)

# =============================================================================
# 1. PanTHERIA + EltonTraits -> sp_trait
# =============================================================================
Pan_mammal_traits <- read.delim(file.path(traits_dir, "PanTHERIA", "PanTHERIA_1-0_WR05_Aug2008.txt")) %>%
  dplyr::select(-MSW05_Order, -MSW05_Family, -MSW05_Genus, -MSW05_Species, -References) %>%
  mutate(across(where(is.numeric), ~ na_if(., -999)))

# EltonTraits 1.0: species-level foraging attributes of birds and mammals
Elton_mammal_trait <- read.delim(file.path(traits_dir, "EltonTraits 1.0", "MamFuncDat.txt")) %>%
  dplyr::select(-MSW3_ID, -MSWFamilyLatin, -Diet.Source, -Diet.Certainty, -ForStrat.Certainty,
                -ForStrat.Comment, -Activity.Source, -Activity.Certainty, -BodyMass.Source, -BodyMass.SpecLevel)

# Merge; encode categorical/ordinal traits numerically for imputation
sp_trait <- Pan_mammal_traits %>%
  left_join(Elton_mammal_trait, join_by(MSW05_Binomial == Scientific)) %>%
  rename(SCIENTIFIC_NAME = MSW05_Binomial) %>%
  mutate(
    activity             = as.numeric(factor(X1.1_ActivityCycle)),
    terrestriality       = as.numeric(factor(X12.2_Terrestriality)),
    ForStrat             = as.numeric(factor(ForStrat.Value)),
    X6.1_DietBreadth     = as.numeric(factor(X6.1_DietBreadth, ordered = TRUE)),
    X12.1_HabitatBreadth = as.numeric(factor(X12.1_HabitatBreadth, ordered = TRUE)),
    X6.2_TrophicLevel    = as.numeric(factor(X6.2_TrophicLevel, ordered = TRUE)),
    DietGuild = case_when(
      Diet.Inv + Diet.Vend + Diet.Vect + Diet.Vfish + Diet.Vunk > 50 ~ "Carnivore",
      Diet.Scav > 50                                                 ~ "Scavenger",
      Diet.Fruit + Diet.Nect + Diet.Seed + Diet.PlantO > 50          ~ "Herbivore",
      TRUE                                                           ~ "Omnivore"),
    DietGuild = as.numeric(factor(DietGuild))) %>%
  dplyr::select(SCIENTIFIC_NAME, activity, terrestriality, ForStrat, DietGuild, BodyMass.Value,
                X13.1_AdultHeadBodyLen_mm, X6.1_DietBreadth, X9.1_GestationLen_d, X12.1_HabitatBreadth,
                X22.1_HomeRange_km2, X15.1_LitterSize, X6.2_TrophicLevel, X27.2_HuPopDen_Mean_n.km2,
                X28.1_Precip_Mean_mm, X28.2_Temp_Mean_01degC, X30.1_AET_Mean_mm, X30.2_PET_Mean_mm)

# =============================================================================
# 2. Harmonise names to VertLife taxonomy
#    Mammal phylogeny: https://doi.org/10.1371/journal.pbio.3000494 (VertLife.org)
# =============================================================================
# NOTE: mammal_RA_weighted is produced by the mammal prep script; read it in.
mammal_RA_weighted <- read.csv(p_clean("mammal_RA_weighted.csv"))

VertLife <- read.csv(file.path(traits_dir, "vertlife_taxonomies.csv")) %>%
  filter(group == "Mammals")

# 122 focal species
species_122 <- mammal_RA_weighted %>%
  dplyr::select(SCIENTIFIC_NAME) %>%
  distinct() %>%
  left_join(VertLife, join_by(SCIENTIFIC_NAME == scientificname)) %>%
  mutate(Vertlife_NAME = case_when(
    str_detect(SCIENTIFIC_NAME, "Neotamias ") ~ str_replace(SCIENTIFIC_NAME, "Neotamias", "Tamias"),
    SCIENTIFIC_NAME == "Urva javanica"     ~ "Herpestes javanicus",
    SCIENTIFIC_NAME == "Pekania pennanti"  ~ "Martes pennanti",
    SCIENTIFIC_NAME == "Cervus canadensis" ~ "Cervus elaphus",
    .default = SCIENTIFIC_NAME))
write.csv(species_122, file.path(traits_dir, "mammal_species_122.csv"), row.names = FALSE)

# 4999 species with a VertLife match (for imputation)
species_5416 <- sp_trait %>%
  dplyr::select(SCIENTIFIC_NAME) %>%
  distinct() %>%
  left_join(VertLife, join_by(SCIENTIFIC_NAME == scientificname)) %>%
  drop_na()
n_distinct(species_5416$SCIENTIFIC_NAME)   # 4999
write.csv(species_5416, file.path(traits_dir, "mammal_species_5416.csv"), row.names = FALSE)

# NOTE: mammal_species_5416.csv is fed to the VertLife phylotree pruner (an
#       external step) which produces output.nex, read back in Section 3.
sp_trait_4999 <- sp_trait %>%
  filter(SCIENTIFIC_NAME %in% unique(species_5416$SCIENTIFIC_NAME))
colMeans(is.na(sp_trait_4999))   # QC: missingness per trait

# =============================================================================
# 3. Phylogeny + phylopars imputation (averaged over 10 trees)
# =============================================================================
mammal_trees <- ape::read.nexus(file.path(traits_dir, "mammal_4999species_phylotree-pruner", "output.nex"))
# plot(mammal_trees[[1]])   # QC

# NOTE: `1:4999` hard-codes the tip count (keep in sync with the species total).
imputed_list <- lapply(mammal_trees[1:10], function(tr) {
  fit <- Rphylopars::phylopars(
    trait_data = sp_trait_4999 %>%
      rename(species = SCIENTIFIC_NAME) %>%
      mutate(species = str_replace_all(species, " ", "_")) %>%
      # drop traits highly correlated with precipitation
      dplyr::select(-X28.2_Temp_Mean_01degC, -X30.1_AET_Mean_mm, -X30.2_PET_Mean_mm, -DietGuild),
    tree  = tr,
    model = "lambda")
  fit$anc_recon[1:4999, ]
})

imp_array     <- abind::abind(imputed_list, along = 3)
trait_imputed <- apply(imp_array, c(1, 2), mean) %>% as.data.frame()

trait_imputed$SCIENTIFIC_NAME <- str_replace(rownames(trait_imputed), "_", " ")
trait_imputed <- trait_imputed %>% relocate(SCIENTIFIC_NAME, .before = activity)
write.csv(trait_imputed, p_clean("mammal_imputed_traits.csv"), row.names = FALSE)
mammal_trait_imputed <- read.csv(p_clean("mammal_imputed_traits.csv"))
