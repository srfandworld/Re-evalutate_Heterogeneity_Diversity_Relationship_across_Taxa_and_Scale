# =============================================================================
# Trait modulation of the abundance-EH relationship
#   step 1: per-species slope of abundance ~ EH  (coef_EH)
#   step 2: coef_EH ~ trait, via phylogenetic meta-regression (lambda test ->
#           weighted lm if no phylogenetic signal, else PGLS)
# Re-evaluating the Environmental Heterogeneity-Diversity Relationship
#
# Sections:
#   0. Libraries
#   1. Meta-regression functions (tree / bird / mammal)
#   2. Trees   (10 km)
#   3. Birds   (10 km)
#   4. Mammals (5 km)
# =============================================================================

# ---- 0. Libraries -----------------------------------------------------------
library(tidyverse)   # dplyr, tidyr, stringr, purrr, readr, ggplot2, tibble
library(broom)       # tidy()
library(nlme)        # gls(), varFixed()
library(ape)         # keep.tip(), corPagel(), read.nexus()
library(phytools)    # phylosig()


# =============================================================================
# 1. Meta-regression functions
#    Each: scale trait covariates, loop over traits, test phylogenetic signal
#    in coef_EH, then fit a variance-weighted lm (no signal) or PGLS (signal).
#    NOTE: the three are ~95% identical; they differ only in the species-ID
#          column and whether Pagel's lambda is clamped to [0, 1]. They can be
#          merged into a single beta_trait(species_col, clamp_lambda) helper -
#          kept separate here to preserve each taxon's exact behaviour.
# =============================================================================

## Trees: coef_EH ~ trait -----------------------------------------------------
tree_beta_trait <- function(data = EH_trait_plot, tree_phylo = FIA_phylo) {
  traits       <- c("swd", "maxHt", "leafnitro", "SLA", "Seed.dry.mass",
                    "vol_total", "vol_core", "vol_periphery", "range_km2", "PFT_binary")
  scale_traits <- c("swd", "maxHt", "leafnitro", "SLA", "Seed.dry.mass",
                    "vol_total", "vol_core", "vol_periphery", "range_km2")
  
  data_scaled <- data %>% mutate(across(all_of(scale_traits), ~ scale(.x)[, 1]))
  result_list <- list()
  
  for (trait in traits) {
    dat <- data_scaled %>%
      dplyr::select(SCIENTIFIC_NAME, coef_EH, std.error, all_of(trait)) %>%
      mutate(species = str_replace_all(SCIENTIFIC_NAME, " ", "_")) %>%
      tidyr::drop_na()
    
    phylo_pruned <- keep.tip(tree_phylo, dat$species)
    dat <- dat %>% arrange(match(species, phylo_pruned$tip.label))
    
    signal_test <- phytools::phylosig(phylo_pruned,
                                      setNames(dat$coef_EH, dat$species),
                                      method = "lambda", test = TRUE)
    
    if (signal_test$P > 0.05) {
      meta_model <- lm(as.formula(paste("coef_EH ~", trait)),
                       data = dat, weights = 1 / std.error^2)
      coefs <- summary(meta_model)$coefficients[-1, , drop = FALSE]
      df    <- summary(meta_model)$df[2]
      result <- data.frame(
        model = "lm", lambda = 0, trait = trait, level = rownames(coefs),
        estimate = coefs[, "Estimate"], std.error = coefs[, "Std. Error"], p = coefs[, "Pr(>|t|)"],
        lwr = coefs[, "Estimate"] - qt(0.975, df) * coefs[, "Std. Error"],
        upr = coefs[, "Estimate"] + qt(0.975, df) * coefs[, "Std. Error"],
        n = nobs(meta_model), row.names = NULL)
    } else {
      cor_struct <- ape::corPagel(value = 1, phy = phylo_pruned, fixed = FALSE, form = ~ species)
      pgls_fit   <- nlme::gls(as.formula(paste("coef_EH ~", trait)), data = dat,
                              correlation = cor_struct, weights = varFixed(~ I(std.error^2)),
                              method = "REML")
      coefs      <- summary(pgls_fit)$tTable[-1, , drop = FALSE]
      df         <- summary(pgls_fit)$dims$N - summary(pgls_fit)$dims$p
      lambda_val <- as.numeric(coef(pgls_fit$modelStruct$corStruct, unconstrained = FALSE))
      result <- data.frame(
        model = "PGLS", lambda = lambda_val, trait = trait, level = rownames(coefs),
        estimate = coefs[, "Value"], std.error = coefs[, "Std.Error"], p = coefs[, "p-value"],
        lwr = coefs[, "Value"] - qt(0.975, df) * coefs[, "Std.Error"],
        upr = coefs[, "Value"] + qt(0.975, df) * coefs[, "Std.Error"],
        n = summary(pgls_fit)$dims$N, row.names = NULL)
    }
    result_list[[trait]] <- result
  }
  do.call(rbind, result_list)
}

## Birds ----------------------------------------------------------------------
bird_beta_trait <- function(data = EH_trait_plot_PGLS, bird_phylo = BBS_phylo) {
  traits       <- c("Beak.Length_Culmen", "Beak.Length_Nares", "Beak.Width", "Beak.Depth", "Tarsus.Length",
                    "Wing.Length", "Kipps.Distance", "Secondary1", "Hand.Wing.Index", "Tail.Length", "Mass",
                    "Habitat", "Habitat.Density", "Migration", "Trophic.Level", "Trophic.Niche", "Primary.Lifestyle",
                    "Clutch_Max", "Incu2", "PCo1", "PCo2", "PCo3", "vol_total", "vol_core", "vol_periphery", "range_km2")
  scale_traits <- c("Beak.Length_Culmen", "Beak.Length_Nares", "Beak.Width", "Beak.Depth", "Tarsus.Length",
                    "Wing.Length", "Kipps.Distance", "Secondary1", "Hand.Wing.Index", "Tail.Length", "Mass",
                    "Clutch_Max", "Incu2", "PCo1", "PCo2", "PCo3", "vol_total", "vol_core", "vol_periphery", "range_km2")
  
  data_scaled <- data %>% mutate(across(all_of(scale_traits), ~ scale(.x)[, 1]))
  result_list <- list()
  
  for (trait in traits) {
    dat <- data_scaled %>%
      dplyr::select(ebird_species, coef_EH, std.error, all_of(trait)) %>%
      tidyr::drop_na()
    
    phylo_pruned <- keep.tip(bird_phylo, dat$ebird_species)   # FIX: was global BBS_phylo
    dat <- dat %>% arrange(match(ebird_species, phylo_pruned$tip.label))
    
    signal_test <- phytools::phylosig(phylo_pruned,
                                      setNames(dat$coef_EH, dat$ebird_species),
                                      method = "lambda", test = TRUE)
    
    if (signal_test$P > 0.05) {
      meta_model <- lm(as.formula(paste("coef_EH ~", trait)),
                       data = dat, weights = 1 / std.error^2)
      coefs <- summary(meta_model)$coefficients[-1, , drop = FALSE]
      df    <- summary(meta_model)$df[2]
      result <- data.frame(
        model = "lm", lambda = 0, trait = trait, level = rownames(coefs),
        estimate = coefs[, "Estimate"], std.error = coefs[, "Std. Error"], p = coefs[, "Pr(>|t|)"],
        lwr = coefs[, "Estimate"] - qt(0.975, df) * coefs[, "Std. Error"],
        upr = coefs[, "Estimate"] + qt(0.975, df) * coefs[, "Std. Error"],
        n = nobs(meta_model), row.names = NULL)
    } else {
      cor_struct <- ape::corPagel(value = 1, phy = phylo_pruned, fixed = FALSE, form = ~ ebird_species)
      pgls_fit   <- nlme::gls(as.formula(paste("coef_EH ~", trait)), data = dat,
                              correlation = cor_struct, weights = varFixed(~ I(std.error^2)),
                              method = "REML")
      lambda <- coef(pgls_fit$modelStruct$corStruct, unconstrained = FALSE)
      lambda <- min(max(lambda, 0), 1)                          # clamp to [0, 1]
      coefs  <- summary(pgls_fit)$tTable[-1, , drop = FALSE]
      df     <- summary(pgls_fit)$dims$N - summary(pgls_fit)$dims$p
      result <- data.frame(
        model = "PGLS", lambda = lambda, trait = trait, level = rownames(coefs),
        estimate = coefs[, "Value"], std.error = coefs[, "Std.Error"], p = coefs[, "p-value"],
        lwr = coefs[, "Value"] - qt(0.975, df) * coefs[, "Std.Error"],
        upr = coefs[, "Value"] + qt(0.975, df) * coefs[, "Std.Error"],
        n = summary(pgls_fit)$dims$N, row.names = NULL)
    }
    result_list[[trait]] <- result
  }
  do.call(rbind, result_list)
}

## Mammals --------------------------------------------------------------------
mammal_beta_trait <- function(data = EH_trait_plot, mammal_phylo = mammal_phylo) {
  traits       <- c("ActivityCycle", "ForStrat.Value", "logbodymass", "AdultHeadBodyLen_mm", "DietBreadth", "DietGuild",
                    "HabitatBreadth", "loghomerange", "TrophicLevel",
                    "GestationLen_d", "InterbirthInterval_d", "LitterSize", "MaxLongevity_m", "NeonateBodyMass_g",
                    "SexualMaturityAge_d", "TeatNumber", "WeaningAge_d",
                    "Terrestriality", "HuPopDen_Mean_n.km2",
                    "Precip_Mean_mm", "Temp_Mean_01degC", "AET_Mean_mm", "PET_Mean_mm",
                    "PCo1", "PCo2", "PCo3", "vol_total", "vol_core", "vol_periphery", "range_km2")
  scale_traits <- c("logbodymass", "AdultHeadBodyLen_mm", "DietBreadth", "loghomerange",
                    "GestationLen_d", "InterbirthInterval_d", "LitterSize", "MaxLongevity_m", "NeonateBodyMass_g",
                    "SexualMaturityAge_d", "TeatNumber", "WeaningAge_d",
                    "HuPopDen_Mean_n.km2",
                    "Precip_Mean_mm", "Temp_Mean_01degC", "AET_Mean_mm", "PET_Mean_mm",
                    "PCo1", "PCo2", "PCo3", "vol_total", "vol_core", "vol_periphery", "range_km2")
  
  data_scaled <- data %>% mutate(across(all_of(scale_traits), ~ scale(.x)[, 1]))
  result_list <- list()
  
  for (trait in traits) {
    dat <- data_scaled %>%
      dplyr::select(SCIENTIFIC_NAME, coef_EH, std.error, all_of(trait)) %>%
      mutate(SCIENTIFIC_NAME = str_replace(SCIENTIFIC_NAME, " ", "_")) %>%
      tidyr::drop_na()
    
    phylo_pruned <- keep.tip(mammal_phylo, dat$SCIENTIFIC_NAME)
    dat <- dat %>% arrange(match(SCIENTIFIC_NAME, phylo_pruned$tip.label))
    
    signal_test <- phytools::phylosig(phylo_pruned,
                                      setNames(dat$coef_EH, dat$SCIENTIFIC_NAME),
                                      method = "lambda", test = TRUE)
    
    if (signal_test$P > 0.05) {
      meta_model <- lm(as.formula(paste("coef_EH ~", trait)),
                       data = dat, weights = 1 / std.error^2)
      coefs <- summary(meta_model)$coefficients[-1, , drop = FALSE]
      df    <- summary(meta_model)$df[2]
      result <- data.frame(
        model = "lm", lambda = 0, trait = trait, level = rownames(coefs),   # FIX: was "0" (string)
        estimate = coefs[, "Estimate"], std.error = coefs[, "Std. Error"], p = coefs[, "Pr(>|t|)"],
        lwr = coefs[, "Estimate"] - qt(0.975, df) * coefs[, "Std. Error"],
        upr = coefs[, "Estimate"] + qt(0.975, df) * coefs[, "Std. Error"],
        n = nobs(meta_model), row.names = NULL)
    } else {
      cor_struct <- ape::corPagel(value = 1, phy = phylo_pruned, fixed = FALSE, form = ~ SCIENTIFIC_NAME)
      pgls_fit   <- nlme::gls(as.formula(paste("coef_EH ~", trait)), data = dat,
                              correlation = cor_struct, weights = varFixed(~ I(std.error^2)),
                              method = "REML")
      lambda <- coef(pgls_fit$modelStruct$corStruct, unconstrained = FALSE)
      lambda <- min(max(lambda, 0), 1)
      coefs  <- summary(pgls_fit)$tTable[-1, , drop = FALSE]
      df     <- summary(pgls_fit)$dims$N - summary(pgls_fit)$dims$p
      result <- data.frame(
        model = "PGLS", lambda = lambda, trait = trait, level = rownames(coefs),
        estimate = coefs[, "Value"], std.error = coefs[, "Std.Error"], p = coefs[, "p-value"],
        lwr = coefs[, "Value"] - qt(0.975, df) * coefs[, "Std.Error"],
        upr = coefs[, "Value"] + qt(0.975, df) * coefs[, "Std.Error"],
        n = summary(pgls_fit)$dims$N, row.names = NULL)
    }
    result_list[[trait]] <- result
  }
  do.call(rbind, result_list)
}

# =============================================================================
# 2. Trees (10 km)
# =============================================================================
FIA_gamma_10km   <- readRDS("output/FIA_gamma_abundance/comm_abundance_r10km.rds")
sp_trait         <- read.csv("data_clean/FIA_864sp_raw_trait.csv")
imputed_sp_trait <- read.csv("data_clean/FIA_864sp_imputed_traits.csv")
fia_niche        <- read.csv("data_clean/niche_results/FIA_niche_breadth_3d.csv")
FIA_phylo        <- readRDS("data_clean/FIA_phylo_263species.rds")   # 263 species
# ADD: used below but produced by earlier scripts - read them in
FIA_table_10km   <- read.csv("data_clean/FIA_table_10km.csv")
fia_sp           <- read.csv("data_clean/fia_sp.csv")

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

tree_data <- tree_beta_trait(data = EH_trait_plot, tree_phylo = FIA_phylo)
write.csv(tree_data, "data_clean/tree_data_10km.csv", row.names = FALSE)

# =============================================================================
# 3. Birds (10 km)
# =============================================================================
BBS_data          <- read.csv("data_clean/BBS_data.csv")
BBS_env           <- read.csv("data_clean/BBS_env.csv")
BBS_traits        <- read.csv("data_clean/BBS_traits.csv")
BBS_traits_space  <- read.csv("data_clean/BBS_traits_space.csv")
BBS_niche_breadth <- read.csv("data_clean/niche_results/BBSBird_niche_breadth_3d.csv")
BBS_phylo         <- readRDS("data_clean/BBS_phylo_520species.rds")   # tips = eBird taxa
BBS_species2      <- read.csv("data_clean/BBS_species2.csv")

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
write.csv(bird_data, "data_clean/bird_data_10km.csv", row.names = FALSE)

# =============================================================================
# 4. Mammals (5 km)
# =============================================================================
mammal_RA_weighted   <- read.csv("data_clean/mammal_RA_weighted.csv")
CTA_env              <- read.csv("data_clean/Camera_Trap_Array_env.csv")
mammal_traits_space  <- read.csv("data_clean/mammal_traits_space.csv")
mammal_niche_breadth <- read.csv("data_clean/niche_results/Mammal_niche_breadth_3d.csv")
# NOTE: absolute path made relative; adjust to where the phylogeny lives.
mammal_trees         <- ape::read.nexus("data/Mammals/mammal_4999species_phylotree-pruner/output.nex")
#mammal traits
Pan_mammal_traits <- read.delim("D:/BiodiversityEmbedding/Traits_data/PanTHERIA/PanTHERIA_1-0_WR05_Aug2008.txt")
Pan_mammal_traits <- Pan_mammal_traits%>%
  dplyr::select(-MSW05_Order,-MSW05_Family,-MSW05_Genus,-MSW05_Species,-References) %>% 
  mutate(across(where(is.numeric), ~ na_if(., -999)))

#EltonTraits 1.0: Species-level foraging attributes of the world's birds and mammals
Elton_mammal_trait <- read.delim("D:/BiodiversityEmbedding/Traits_data/EltonTraits 1.0/MamFuncDat.txt")
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
write.csv(mammal_data, "data_clean/mammal_data_5km.csv", row.names = FALSE)