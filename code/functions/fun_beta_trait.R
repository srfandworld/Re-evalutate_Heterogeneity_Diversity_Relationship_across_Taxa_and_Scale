# =============================================================================
# fun_beta_trait.R - second-stage meta-regression: coef_EH ~ trait
#   Pagel's lambda test on coef_EH; weighted lm (1/SE^2) if no signal,
#   otherwise PGLS (gls + corPagel + varFixed(SE^2)).
#
# Source   : Code/Figure3/01_figure3_analysis.R, lines 24-207 (verbatim)
#
# TODO(reviewer M4): trees report the unconstrained lambda while birds/mammals
#   clamp it to [0, 1] for reporting only; coefficients always come from the
#   unconstrained fit. Decide on one rule and state it in the Methods.
# =============================================================================

library(dplyr)
library(stringr)
library(nlme)
library(ape)
library(phytools)
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
