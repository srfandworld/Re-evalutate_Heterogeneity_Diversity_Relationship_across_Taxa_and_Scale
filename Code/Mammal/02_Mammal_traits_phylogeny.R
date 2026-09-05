# =============================================================================
# Mammal traits: PanTHERIA + EltonTraits, phylogenetic imputation, and PCoA
# Re-evaluating the Environmental Heterogeneity-Diversity Relationship
#
# Sections:
#   1. Load PanTHERIA + EltonTraits -> sp_trait
#   2. Harmonise names to VertLife taxonomy (122 focal + 4999 for imputation)
#   3. Phylogeny + phylopars imputation across trees -> mammal_imputed_traits.csv
#   4. Trait PCoA -> mammal_traits_space.csv
# =============================================================================

# ---- Libraries --------------------------------------------------------------
library(tidyverse)   # dplyr, tidyr, stringr, ggplot2, purrr, readr, tibble
library(ape)         # read.nexus()
library(Rphylopars)  # phylopars()
library(abind)       # abind()
library(FD)          # gowdis()
library(grid)        # unit(), arrow()  (used in the PCoA biplot)


traits_dir <- "D:/BiodiversityEmbedding/Traits_data"   # raw trait sources (EDIT ME)

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
mammal_RA_weighted <- read.csv("data_clean/mammal_RA_weighted.csv")

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
write.csv(trait_imputed, "data_clean/mammal_imputed_traits.csv", row.names = FALSE)
mammal_trait_imputed <- read.csv("data_clean/mammal_imputed_traits.csv")

# =============================================================================
# 4. Trait PCoA -> mammal_traits_space.csv
# =============================================================================
# NOTE: EH_trait_plot is NOT created here - it comes from the mammal abundance-EH
#       script (the trait-response script), in its pre-PCA form (raw trait
#       columns, before that script joins mammal_traits_space back on). Run that
#       script up to the EH_trait_plot construction, then this section, which
#       produces mammal_traits_space.csv for it to consume. (Circular by design.)
trait_df <- EH_trait_plot %>%
  select(ActivityCycle, AdultBodyMass_g, AdultHeadBodyLen_mm, DietBreadth, GestationLen_d,
         HabitatBreadth, InterbirthInterval_d, LitterSize, MaxLongevity_m, NeonateBodyMass_g,
         PopulationDensity_n.km2, SexualMaturityAge_d, TeatNumber, Terrestriality,
         TrophicLevel, WeaningAge_d, HuPopDen_Mean_n.km2, Precip_Mean_mm, Temp_Mean_01degC,
         AET_Mean_mm, PET_Mean_mm, ForStrat.Value, logbodymass, loghomerange, DietGuild)

# --- QC correlation panels (interactive) -------------------------------------
panel.cor <- function(x, y, digits = 2, cex.cor = 1.2) {
  usr <- par("usr"); on.exit(par(usr)); par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use = "complete.obs")
  text(0.5, 0.5, formatC(abs(r), format = "f", digits = digits), cex = cex.cor)
}
panel.hist <- function(x, ...) {
  usr <- par("usr"); par(usr = c(usr[1:2], 0, 1.5))
  h <- hist(x, plot = FALSE); breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts / max(h$counts)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
# pairs(trait_df, diag.panel = panel.hist, upper.panel = panel.cor)

# Keep traits with |r| < 0.8 (drop the correlated ones)
trait_df <- trait_df %>%
  dplyr::select(-AdultBodyMass_g, -AdultHeadBodyLen_mm, -GestationLen_d,
                -SexualMaturityAge_d, -Temp_Mean_01degC, -AET_Mean_mm, -PET_Mean_mm)

# PCoA on Gower distances of scaled traits
trait_scaled <- trait_df %>% mutate(across(where(is.numeric), ~ as.numeric(scale(.))))
gower_dist   <- FD::gowdis(trait_scaled)
pcoa_res     <- cmdscale(gower_dist, k = 3, eig = TRUE)   # classical (metric) MDS

traits_space <- as.data.frame(pcoa_res$points)
colnames(traits_space) <- c("PCo1", "PCo2", "PCo3")
traits_space <- cbind(traits_space, trait_df[, sapply(trait_df, is.factor)])

# Export the 3-axis trait space (one row per species)
traits_space_output <- EH_trait_plot %>% select(SCIENTIFIC_NAME) %>% cbind(traits_space[, 1:3])
write.csv(traits_space_output, "data_clean/mammal_traits_space.csv", row.names = FALSE)

# --- QC biplot (interactive) -------------------------------------------------
var_exp <- pcoa_res$eig / sum(pcoa_res$eig)
xlab <- paste0("PCo1 (", round(var_exp[1] * 100, 1), "%)")
ylab <- paste0("PCo2 (", round(var_exp[2] * 100, 1), "%)")

traits_num    <- trait_df[, sapply(trait_df, is.numeric)]
trait_vectors <- as.data.frame(cor(traits_num, traits_space[, 1:2], use = "pairwise.complete.obs"))
trait_vectors$trait <- rownames(trait_vectors)
arrow_scale <- 0.4   # adjustable, ~1 to 5
trait_vectors <- trait_vectors %>% mutate(PCo1 = PCo1 * arrow_scale, PCo2 = PCo2 * arrow_scale)

ggplot() +
  geom_point(data = traits_space, aes(PCo1, PCo2, color = DietGuild), size = 3, alpha = 0.7) +
  stat_ellipse(data = traits_space, aes(PCo1, PCo2, color = DietGuild)) +
  geom_segment(data = trait_vectors, aes(0, 0, xend = PCo1, yend = PCo2),
               arrow = arrow(length = unit(0.25, "cm")), linewidth = 0.8, color = "black") +
  geom_text(data = trait_vectors, aes(PCo1, PCo2, label = trait), size = 4, hjust = 0.5, vjust = -0.5) +
  geom_hline(yintercept = 0, col = "grey", linetype = "dashed") +
  geom_vline(xintercept = 0, col = "grey", linetype = "dashed") +
  labs(x = xlab, y = ylab) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(), axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black")) +
  coord_equal()