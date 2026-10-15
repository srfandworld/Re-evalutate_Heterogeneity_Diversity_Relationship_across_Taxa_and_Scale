# =============================================================================
# 01_species_EH_trait_models.R
# Species-level abundance ~ EH slopes (step 1) and coef_EH ~ trait
# meta-regression (step 2): trees 10 km, birds 10 km, mammals 5 km
#
# Source   : Code/Figure3/01_figure3_analysis.R, sections 2-4. Reorganised:
#            - meta-regression functions -> code/functions/fun_beta_trait.R
#            - tree phylogeny build restored (was commented out elsewhere)
#            - mammal trait PCoA moved in from Code/Mammal/02 (breaks the cycle)
# Inputs   : data/clean/FIA_gamma_abundance/comm_abundance_r10km.rds,
#            FIA_table_10km.rds, fia_sp.csv, FIA_864sp_imputed_traits.csv,
#            niche_results/*_niche_breadth_3d.csv, BBS_data.csv, BBS_env.csv,
#            BBS_traits.csv, BBS_traits_space.csv, BBS_phylo_520species.rds,
#            BBS_species2.csv, mammal_RA_weighted.csv, Camera_Trap_Array_env.csv
#            HDR_EXT_DIR/Traits_data/{PanTHERIA, EltonTraits 1.0,
#              mammal_4999species_phylotree-pruner/output.nex}
# Outputs  : data/clean/FIA_phylo_263species.rds, mammal_traits_space.csv,
#            tree_data_10km.csv, bird_data_10km.csv, mammal_data_5km.csv
#            (-> 02_Figure3_forestplot.R and Table S5)
#
# TODO(reviewer M1): bird filter (n_obs > 100) gives 281 species; Methods and
#   Table S5 say 295.
# TODO(reviewer M4): tree response is log(comm_BAA) summed within 10 km with no
#   plot-count offset; bird/mammal predictors are standardised within each
#   species' sites, trees globally. Document both choices.
# TODO(reviewer M4): manual PanTHERIA fills for Sciurus niger, Tamias striatus,
#   Odocoileus hemionus, Sus scrofa need a cited source (Tamias striatus coded
#   above-ground vs EltonTraits ground foraging).
# TODO(reviewer M4): Table S5 tree N varies by trait (178-263) although the
#   imputed trait table is complete - document where N drops.
# =============================================================================

source(here::here("code", "00_config.R"))
source(here::here("code", "functions", "fun_beta_trait.R"))

# ---- 0. Libraries -----------------------------------------------------------
library(tidyverse)   # dplyr, tidyr, stringr, purrr, readr, ggplot2, tibble
library(broom)       # tidy()
library(nlme)        # gls(), varFixed()
library(ape)         # keep.tip(), corPagel(), read.nexus()
library(phytools)    # phylosig()

library(glue)
library(FD)          # gowdis() for the mammal PCoA
library(grid)        # unit(), arrow() in the PCoA biplot

# =============================================================================
# 2. Trees (10 km)
# =============================================================================
FIA_gamma_10km   <- readRDS(p_clean("FIA_gamma_abundance/comm_abundance_r10km.rds"))
sp_trait         <- read.csv(p_clean("FIA_864sp_raw_trait.csv"))
imputed_sp_trait <- read.csv(p_clean("FIA_864sp_imputed_traits.csv"))
fia_niche        <- read.csv(p_clean("niche_results/FIA_niche_breadth_3d.csv"))
# ADD: used below but produced by earlier scripts - read them in
# TODO(reviewer C5): the original read FIA_table_10km.csv, which no script writes.
# The scaled all-plot table (03_community_metrics/01) is used here; confirm this
# (vs. the iNEXT subset FIA_table_10km_iNEXT.csv) reproduces the 263 species.
FIA_table_10km   <- readRDS(p_clean("FIA_table_10km.rds"))
fia_sp           <- read.csv(p_clean("fia_sp.csv"))

FIA_abundance_10km <- FIA_gamma_10km %>%
  rename(pltID = focal_pltID) %>%
  filter(pltID %in% unique(FIA_table_10km$pltID)) %>%
  left_join(FIA_table_10km, by = "pltID")

FIA_abundance_10km_filter <- FIA_abundance_10km %>%
  drop_na(comm_BAA, EH_10km_t, MAT_10km_t, MAP_10km_t, soilph_10km_t, soilcec_10km_t, elev_mean_10km_t) %>%
  filter(SCIENTIFIC_NAME %in% unique(fia_sp$SCIENTIFIC_NAME)) %>%
  group_by(SCIENTIFIC_NAME) %>%
  mutate(n_obs = n_distinct(pltID)) %>%
  ungroup() %>%
  filter(n_obs > 100)                       # minimum observations for modelling
n_distinct(FIA_abundance_10km_filter$SCIENTIFIC_NAME)   # 263 tree species

# Per-species abundance ~ EH slopes
model_results <- FIA_abundance_10km_filter %>%
  group_by(SCIENTIFIC_NAME) %>%
  nest() %>%
  mutate(
    model = purrr::map(data, ~ lm(log(comm_BAA) ~ EH_10km_t + MAT_10km_t + MAP_10km_t +
                                    soilph_10km_t + soilcec_10km_t + elev_mean_10km_t, data = .x)),
    tidy  = purrr::map(model, tidy)) %>%
  unnest(tidy) %>% ungroup() %>%
  dplyr::select(-data, -model)

EH_coef <- model_results %>%
  filter(term == "EH_10km_t") %>%
  dplyr::select(SCIENTIFIC_NAME, estimate, std.error, statistic, p.value) %>%
  rename(coef_EH = estimate)

EH_trait <- EH_coef %>%
  left_join(imputed_sp_trait, by = "SCIENTIFIC_NAME") %>%              # imputed 864-sp traits
  left_join(fia_niche %>% dplyr::select(-n_points, -n_occ), join_by(SCIENTIFIC_NAME == species))

EH_trait_plot <- EH_trait %>%
  mutate(PFT_binary = case_when(
    PFT %in% c("broad_deciduous", "needle_deciduous")                       ~ "Deciduous",
    PFT %in% c("broad_evergreen", "needle_evergreen", "scalelike_evergreen") ~ "Evergreen",
    .default = PFT))


# Tree phylogeny for the modelled species (V.PhyloMaker2, scenario 3).
# Source: EmbeddingBio/analysis/FIA_speciesabundance-EH.R lines 153-163, where
# it was commented out after the first run; restored so FIA_phylo_263species.rds
# can be regenerated. Fig. S1a tree: same file, lines 165-172 (ggtree).
phylo_file <- p_clean("FIA_phylo_263species.rds")
if (!file.exists(phylo_file)) {
  taxonomic_name <- EH_trait_plot %>%
    rownames_to_column() %>%
    dplyr::select(rowname, SCIENTIFIC_NAME) %>%
    ungroup()
  sp.list.sp_ge_fa <- TNRS::TNRS(taxonomic_names = taxonomic_name) %>%
    dplyr::select(species = Name_submitted, genus = Genus_matched, family = Accepted_family)
  phylo <- V.PhyloMaker2::phylo.maker(sp.list = sp.list.sp_ge_fa)
  saveRDS(phylo$scenario.3, phylo_file)
}
FIA_phylo <- readRDS(phylo_file)   # 263 species
tree_data <- tree_beta_trait(data = EH_trait_plot, tree_phylo = FIA_phylo)
write.csv(tree_data, p_clean("tree_data_10km.csv"), row.names = FALSE)

# =============================================================================
# 3. Birds (10 km)
# =============================================================================
BBS_data          <- read.csv(p_clean("BBS_data.csv"))
BBS_env           <- read.csv(p_clean("BBS_env.csv"))
BBS_traits        <- read.csv(p_clean("BBS_traits.csv"))
BBS_traits_space  <- read.csv(p_clean("BBS_traits_space.csv"))
BBS_niche_breadth <- read.csv(p_clean("niche_results/BBSBird_niche_breadth_3d.csv"))
BBS_phylo         <- readRDS(p_clean("BBS_phylo_520species.rds"))   # tips = eBird taxa
BBS_species2      <- read.csv(p_clean("BBS_species2.csv"))

BBS_abundance <- BBS_data %>%
  group_by(route, aou) %>%
  summarise(annual_mean_count = unique(annual_mean_count), .groups = "drop") %>%
  left_join(BBS_env, by = "route")

BBS_abundance_filter <- BBS_abundance %>%
  group_by(aou) %>%
  mutate(n_obs = n_distinct(route)) %>%
  ungroup() %>%
  filter(n_obs > 100)
n_distinct(BBS_abundance_filter$aou)   # 281 species

model_results <- BBS_abundance_filter %>%
  group_by(aou) %>%
  nest() %>%
  mutate(
    model = purrr::map(data, ~ lm(log(annual_mean_count) ~ scale(EH_10km)[, 1] + scale(MAT_10km)[, 1] +
                                    scale(MAP_10km)[, 1] + scale(Tsea_10km)[, 1] + scale(Psea_10km)[, 1] +
                                    scale(log1p(PopDen_10km))[, 1] + scale(sqrt(elev_10km))[, 1], data = .x)),
    tidy  = purrr::map(model, tidy)) %>%
  unnest(tidy) %>% ungroup() %>%
  dplyr::select(-data, -model)

EH_coef <- model_results %>%
  filter(term == "scale(EH_10km)[, 1]") %>%
  dplyr::select(aou, estimate, std.error, statistic, p.value) %>%
  rename(coef_EH = estimate)

EH_trait <- EH_coef %>%
  left_join(BBS_traits, by = "aou") %>%
  relocate(ScientificName, .after = aou) %>%
  left_join(BBS_traits_space, by = "aou") %>%
  left_join(BBS_niche_breadth %>% dplyr::select(aou, vol_total, vol_core, vol_periphery, range_km2),
            by = "aou", multiple = "any")

EH_trait_plot <- EH_trait %>%
  mutate(
    Habitat           = factor(Habitat),
    Habitat.Density   = factor(Habitat.Density, levels = c(1, 2, 3), labels = c("Dense", "Semi-open", "Open")),
    Migration         = factor(Migration, levels = c(1, 2, 3), labels = c("Sedentary", "Partially migratory", "Migratory")),
    Trophic.Level     = factor(Trophic.Level, levels = c("Herbivore", "Omnivore", "Carnivore", "Scavenger")),
    Trophic.Niche     = factor(Trophic.Niche),
    Primary.Lifestyle = factor(Primary.Lifestyle))

EH_trait_plot_PGLS <- EH_trait_plot %>%
  left_join(BBS_species2 %>% dplyr::select(aou, ebird_species), by = "aou") %>%
  relocate(ebird_species, .after = ScientificName) %>%
  mutate(ebird_species = str_replace(ebird_species, " ", "_"))

bird_data <- bird_beta_trait(data = EH_trait_plot_PGLS, bird_phylo = BBS_phylo)
write.csv(bird_data, p_clean("bird_data_10km.csv"), row.names = FALSE)

# =============================================================================
# 4. Mammals (5 km)
# =============================================================================
mammal_RA_weighted   <- read.csv(p_clean("mammal_RA_weighted.csv"))
CTA_env              <- read.csv(p_clean("Camera_Trap_Array_env.csv"))
mammal_niche_breadth <- read.csv(p_clean("niche_results/Mammal_niche_breadth_3d.csv"))
mammal_trees         <- ape::read.nexus(p_trait("mammal_4999species_phylotree-pruner", "output.nex"))
#mammal traits
Pan_mammal_traits <- read.delim(p_trait("PanTHERIA/PanTHERIA_1-0_WR05_Aug2008.txt"))
Pan_mammal_traits <- Pan_mammal_traits%>%
  dplyr::select(-MSW05_Order,-MSW05_Family,-MSW05_Genus,-MSW05_Species,-References) %>% 
  mutate(across(where(is.numeric), ~ na_if(., -999)))

#EltonTraits 1.0: Species-level foraging attributes of the world's birds and mammals
Elton_mammal_trait <- read.delim(p_trait("EltonTraits 1.0/MamFuncDat.txt"))
Elton_mammal_trait <- Elton_mammal_trait %>% 
  dplyr::select(-MSW3_ID,-MSWFamilyLatin,-Diet.Source,-Diet.Certainty,-ForStrat.Certainty,-ForStrat.Comment,-Activity.Source,-Activity.Certainty,-BodyMass.Source,-BodyMass.SpecLevel)

mammal_abundance <- mammal_RA_weighted %>%
  mutate(year_CTA = str_c(Year, Camera_Trap_Array, sep = "_")) %>%
  relocate(year_CTA, .before = Year) %>%
  left_join(CTA_env, by = c("year_CTA", "Year", "Camera_Trap_Array")) %>%
  mutate(NA_L1CODE = factor(as.integer(NA_L1CODE), levels = sort(unique(as.integer(NA_L1CODE)))))

mammal_abundance_filter <- mammal_abundance %>%
  group_by(SCIENTIFIC_NAME) %>%
  mutate(n_obs = n_distinct(year_CTA)) %>%
  ungroup() %>%
  filter(n_obs > 100)
n_distinct(mammal_abundance_filter$SCIENTIFIC_NAME)   # 17 species (>100 obs)

model_results <- mammal_abundance_filter %>%
  group_by(SCIENTIFIC_NAME) %>%
  nest() %>%
  mutate(
    model = purrr::map(data, ~ lm(log(RA_weighted) ~ scale(EH_5km)[, 1] + scale(Temp_5km)[, 1] +
                                    scale(sqrt(Prec_5km))[, 1] + scale(NPP_5km)[, 1] +
                                    scale(log1p(PopDen_5km))[, 1] + scale(sqrt(elev_5km))[, 1], data = .x)),
    tidy  = purrr::map(model, tidy)) %>%
  unnest(tidy) %>% ungroup() %>%
  dplyr::select(-data, -model)

EH_coef <- model_results %>%
  filter(term == "scale(EH_5km)[, 1]") %>%
  dplyr::select(SCIENTIFIC_NAME, estimate, std.error, statistic, p.value) %>%
  rename(coef_EH = estimate)

EH_trait <- EH_coef %>%
  left_join(Pan_mammal_traits,  join_by(SCIENTIFIC_NAME == MSW05_Binomial)) %>%   # unimputed
  left_join(Elton_mammal_trait, join_by(SCIENTIFIC_NAME == Scientific)) %>%       # unimputed
  mutate(
    X1.1_ActivityCycle = case_when(
      SCIENTIFIC_NAME == "Sciurus niger"       ~ 3,
      SCIENTIFIC_NAME == "Tamias striatus"     ~ 3,
      SCIENTIFIC_NAME == "Odocoileus hemionus" ~ 2,
      .default = X1.1_ActivityCycle),
    X12.2_Terrestriality = case_when(
      SCIENTIFIC_NAME == "Sciurus niger"   ~ 2,
      SCIENTIFIC_NAME == "Sus scrofa"      ~ 1,
      SCIENTIFIC_NAME == "Tamias striatus" ~ 2,
      .default = X12.2_Terrestriality),
    X6.2_TrophicLevel = case_when(
      SCIENTIFIC_NAME == "Sciurus niger"   ~ 2,
      SCIENTIFIC_NAME == "Tamias striatus" ~ 2,
      .default = X6.2_TrophicLevel)) %>%
  left_join(mammal_niche_breadth %>% dplyr::select(species, vol_total, vol_core, vol_periphery, range_km2),
            join_by(SCIENTIFIC_NAME == species))

colnames(EH_trait) <- sub("^X\\d+\\.\\d+_", "", colnames(EH_trait))

EH_trait_plot <- EH_trait %>%
  mutate(
    ActivityCycle  = factor(ActivityCycle,  levels = c(1, 2, 3), labels = c("Nocturnal only", "Mixed", "Diurnal only")),
    Terrestriality = factor(Terrestriality, levels = c(1, 2),    labels = c("Ground dwelling", "Above-ground dwelling")),
    ForStrat.Value = factor(ForStrat.Value, levels = c("G", "S", "Ar"), labels = c("Ground", "Scansorial", "Arboreal")),
    HabitatBreadth = factor(HabitatBreadth),
    TrophicLevel   = factor(TrophicLevel,   levels = c(1, 2, 3), labels = c("Herbivore", "Omnivore", "Carnivore")),
    logbodymass    = log(BodyMass.Value),
    loghomerange   = log(HomeRange_km2),
    DietGuild = case_when(
      Diet.Inv + Diet.Vend + Diet.Vect + Diet.Vfish + Diet.Vunk > 50 ~ "Carnivore",
      Diet.Scav > 50                                                 ~ "Scavenger",
      Diet.Fruit + Diet.Nect + Diet.Seed + Diet.PlantO > 50          ~ "Herbivore",
      TRUE                                                           ~ "Omnivore"),
    DietGuild = factor(DietGuild, levels = c("Herbivore", "Omnivore", "Carnivore", "Scavenger"))) %>%
  dplyr::select(-contains("Diet."))


# ---- Mammal trait PCoA -> mammal_traits_space.csv (Fig. S4) ---------------
# Moved here from Code/Mammal/02_Mammal_traits_phylogeny.R section 4 (lines
# 122-192): it needs the 17-species EH_trait_plot built just above, which
# removes the circular dependency between the two scripts (reviewer C5).
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
write.csv(traits_space_output, p_clean("mammal_traits_space.csv"), row.names = FALSE)
mammal_traits_space <- traits_space_output

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

# Add PCA components; drop columns that are >50% NA
# NOTE: this drop can remove trait columns the model expects. If a trait listed
#       in mammal_beta_trait() gets dropped here, select(all_of(trait)) errors -
#       confirm the retained columns cover the trait list.
EH_trait_plot <- EH_trait_plot %>% left_join(mammal_traits_space, by = "SCIENTIFIC_NAME")
EH_trait_plot <- EH_trait_plot[, colMeans(is.na(EH_trait_plot)) <= 0.5]

# FIX: build mammal_phylo HERE (after the mammal EH_trait_plot exists). The
#      original built it near the top of this section, where EH_trait_plot still
#      held the BIRD table, so it pruned to the wrong (or missing) species.
mammal_phylo <- ape::keep.tip(
  phy = mammal_trees[[1]],
  tip = str_replace(EH_trait_plot$SCIENTIFIC_NAME, " ", "_"))

mammal_data <- mammal_beta_trait(data = EH_trait_plot, mammal_phylo = mammal_phylo)
write.csv(mammal_data, p_clean("mammal_data_5km.csv"), row.names = FALSE)
