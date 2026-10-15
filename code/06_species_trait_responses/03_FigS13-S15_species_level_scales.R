# =============================================================================
# 03_FigS13-S15_species_level_scales.R
# Fig. S13: pooled cross-species effect of EH on species-level abundance
#           (lmer, random EH slope by species) at every scale and taxon
#
# Source   : marginal_effect()  EmbeddingBio/analysis/fun_betaEH_trait.R, lines 295-529
#            panel layout       EmbeddingBio/analysis/figure3.R, lines 17-51
#            (paths only changed)
# Inputs   : data/clean/FIA_gamma_abundance/comm_abundance_r<scale>.rds,
#            FIA_table_<scale>.rds, fia_sp.csv, BBS_*, mammal_* (see function)
# Outputs  : patchwork printed to the device (see TODO at bottom to save)
#
# TODO(reviewer M1 / Fig. S13): figure labels 281 bird species; Methods and
#   Table S5 say 295.
# TODO(reviewer, "Abundance relationships"): bird CIs include zero at all scales
#   and mammals are positive at 5 km only; align the Results text and the SI
#   caption ("positive responses across 1-10 km").
# =============================================================================

source(here::here("code", "00_config.R"))
library(tidyverse)
library(glue)
library(lme4)
library(ggeffects)
library(broom.mixed)
library(performance)
library(patchwork)

fia_sp <- read.csv(p_clean("fia_sp.csv"))
# Tree per-scale covariate tables (all plots, scaled) used by marginal_effect().
# TODO(reviewer C2/C5): confirm whether these or the *_iNEXT subsets were used.
for (s in SCALES_TREE) assign(glue("FIA_table_{s}"), readRDS(p_clean(glue("FIA_table_{s}.rds"))))

#####Supplemantary figure####
marginal_effect <- function(taxa='tree',scale='10km'){
  if(taxa == "tree"){
    term = paste0('EH_',scale,"_t")
    
    FIA_gamma <- readRDS(paste0(p_clean("FIA_gamma_abundance/comm_abundance_r"),scale,".rds"))
    sp_trait <- read.csv(p_clean("FIA_864sp_raw_trait.csv"))
    fia_niche <- read.csv(p_clean("niche_results/FIA_niche_breadth_3d.csv"))
    
    
    FIA_abundance <- FIA_gamma %>% 
      rename(pltID = focal_pltID) %>% 
      filter(pltID %in% unique(get(glue('FIA_table_',scale))$pltID)) %>% 
      left_join(get(glue('FIA_table_',scale)), by="pltID")
  
    
    FIA_abundance_filter <- FIA_abundance %>%  
      drop_na(comm_BAA,contains('EH'), contains('MAP'),contains('MAT'),contains('soilph'), contains('soilcec'),contains('elev'))%>%
      filter(SCIENTIFIC_NAME %in% unique(fia_sp$SCIENTIFIC_NAME)) %>%  #keep effective species name
      group_by(SCIENTIFIC_NAME) %>% 
      mutate(n_obs = n_distinct(pltID)) %>% 
      ungroup() %>% 
      filter(n_obs > 100)         # minimum observations for modelling
    
    n_sp = n_distinct(FIA_abundance_filter$SCIENTIFIC_NAME)
    
    #########################################global analysis with species as random effects##############################################
    # to decrease the effect of baseline abundance of each species, we log-tranform abundance
    form <- as.formula(glue("log(comm_BAA)~ EH_{scale}_t +MAT_{scale}_t+ MAP_{scale}_t +soilph_{scale}_t + soilcec_{scale}_t + elev_mean_{scale}_t + (1 + EH_{scale}_t || SCIENTIFIC_NAME) "))
    
    global_model <- lme4::lmer(
      form,
      data = FIA_abundance_filter
    )
    
    conditional_R2 = r2(global_model)[[1]]
    
    fixed_ci <- broom.mixed::tidy(global_model,
                                  effects = "fixed",
                                  conf.int = TRUE,
                                  conf.method = "Wald")
    beta_EH = fixed_ci[2,'estimate']
    conf.low = fixed_ci[2,'conf.low']
    conf.high = fixed_ci[2,'conf.high']
    
    
    #plot marginal effects of EH
    eh_effect <- ggpredict(global_model, terms = term)
    ggplot(eh_effect, aes(x = x, y = predicted)) +
      geom_line(color = "#3B6D11", linewidth = 1.2) +
      geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
                  fill = "#3B6D11", alpha = 0.2) +
      labs(x = "Environmental heterogeneity (Standardized)",
           y = "Predicted species-level relative abundance") +
      annotate("text",
               x     = -Inf,
               y     = Inf,
               label = sprintf(
                 "No. of species (> 100 obs.) = %d\nConditional R^2 = %.2f\nβ = %.3f (95%% CI: %.3f, %.3f)",
                 n_sp, conditional_R2, beta_EH, conf.low, conf.high
               ),
               hjust = -0.1,    # small positive = slight right offset from left edge
               vjust = 1.5,     # positive = slight down offset from top edge
               size = 4, color = "black") +
      theme_classic(base_size = 13) +
      theme(
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 11)
      )
  }else if (taxa == 'bird'){
    
    term = paste0('EH_',scale,"_t")
    BBS_data <- read.csv(p_clean("BBS_data.csv"))
    BBS_env <- read.csv(p_clean("BBS_env.csv"))
    BBS_traits <- read.csv(p_clean("BBS_traits.csv"))
    BBS_traits_space <- read.csv(p_clean("BBS_traits_space.csv"))
    BBS_niche_breadth <- read.csv(p_clean("niche_results/BBSBird_niche_breadth_3d.csv"))
    
    
    BBS_abundance <- BBS_data %>% 
      group_by(route,aou) %>% 
      summarise(annual_mean_count = unique(annual_mean_count),
                .groups = 'drop') %>% 
      left_join(BBS_env,by="route")
    n_distinct(BBS_abundance$aou)  #584 species
    
    #choose species with enough observations
    BBS_abundance_filter <- BBS_abundance %>% 
      group_by(aou) %>% 
      mutate(n_obs = n_distinct(route)) %>% 
      ungroup() %>% 
      filter(n_obs > 100)   
    n_sp <- n_distinct(BBS_abundance_filter$aou)  
    
  
    ############################################global analysis with species as random effects####################################
    
    dat <- BBS_abundance_filter %>%
      mutate(
        "{glue('EH_{scale}_t')}" := base::scale(.data[[glue("EH_{scale}")]])[, 1],
        "{glue('MAT_{scale}_t')}" := base::scale(.data[[glue("MAT_{scale}")]])[, 1],
        "{glue('MAP_{scale}_t')}" := base::scale(.data[[glue("MAP_{scale}")]])[, 1],
        "{glue('Tsea_{scale}_t')}" := base::scale(.data[[glue("Tsea_{scale}")]])[, 1],
        "{glue('Psea_{scale}_t')}" := base::scale(.data[[glue("Psea_{scale}")]])[, 1],
        "{glue('PopDen_{scale}_t')}" := base::scale(log1p(.data[[glue("PopDen_{scale}")]]))[, 1],
        "{glue('elev_{scale}_t')}" := base::scale(sqrt(.data[[glue("elev_{scale}")]]))[, 1]
      )
    
    global_model <- lme4::lmer(
      as.formula(glue(
        "log(annual_mean_count) ~
     EH_{scale}_t + MAT_{scale}_t + MAP_{scale}_t + Tsea_{scale}_t +
     Psea_{scale}_t + PopDen_{scale}_t + elev_{scale}_t +
     (1 + EH_{scale}_t | aou)"
      )),
      data = dat
    )
    
    conditional_R2 = r2(global_model)[[1]]
    
    fixed_ci <- broom.mixed::tidy(global_model,
                                  effects = "fixed",
                                  conf.int = TRUE,
                                  conf.method = "Wald")
    beta_EH = fixed_ci[2,'estimate']
    conf.low = fixed_ci[2,'conf.low']
    conf.high = fixed_ci[2,'conf.high']
  
    
    eh_effect <- ggeffects::ggpredict(global_model, terms = glue("EH_{scale}_t"))
    ggplot(eh_effect, aes(x = x, y = predicted)) +
      geom_line(color = "#185FA5", linewidth = 1.2) +
      geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
                  fill = "#185FA5", alpha = 0.2) +
      labs(x = "Environmental heterogeneity (Standardized)",
           y = "Predicted species-level relative abundance") +
      annotate("text",
               x     = -Inf,
               y     = Inf,
               label = sprintf(
                 "No. of species (> 100 obs.) = %d\nConditional R^2 = %.2f\nβ = %.3f (95%% CI: %.3f, %.3f)",
                 n_sp, conditional_R2, beta_EH, conf.low, conf.high
               ),
               hjust = -0.1,    # small positive = slight right offset from left edge
               vjust = 1.5,     # positive = slight down offset from top edge
               size = 4, color = "black") +
      theme_classic(base_size = 13) +
      theme(
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 11)
      )
    
  }else{ #mammal
    term = paste0('EH_',scale,"_t")
    
    mammal_RA_weighted <- read.csv(p_clean("mammal_RA_weighted.csv"))
    CTA_env <- read.csv(p_clean("Camera_Trap_Array_env.csv"))
    mammal_traits_space <- read.csv(p_clean("mammal_traits_space.csv"))
    mammal_niche_breadth <- read.csv(p_clean("niche_results/Mammal_niche_breadth_3d.csv"))
    
    mammal_abundance <- mammal_RA_weighted %>% 
      mutate(year_CTA = str_c(Year,Camera_Trap_Array,sep = "_")) %>% 
      relocate(year_CTA,.before = Year) %>% 
      left_join(CTA_env,by=c("year_CTA","Year","Camera_Trap_Array")) %>% 
      mutate(NA_L1CODE = factor(as.integer(NA_L1CODE),
                                levels = sort(unique(as.integer(NA_L1CODE)))))
    
    #choose species with enough observations
    mammal_abundance_filter <- mammal_abundance %>% 
      group_by(SCIENTIFIC_NAME) %>% 
      mutate(n_obs = n_distinct(year_CTA)) %>% 
      ungroup() %>% 
      filter(n_obs > 100)
    
    n_sp = n_distinct(mammal_abundance_filter$SCIENTIFIC_NAME)
    
    
    ############################################global analysis with species as random effects####################################
    
    dat <- mammal_abundance_filter %>%
      mutate(
        "{glue('EH_{scale}_t')}" := base::scale(.data[[glue("EH_{scale}")]])[, 1],
        "{glue('Temp_{scale}_t')}" := base::scale(.data[[glue("Temp_{scale}")]])[, 1],
        "{glue('Prec_{scale}_t')}" := base::scale(sqrt(.data[[glue("Prec_{scale}")]]))[, 1],
        "{glue('NPP_{scale}_t')}" := base::scale(.data[[glue("NPP_{scale}")]])[, 1],
        "{glue('PopDen_{scale}_t')}" := base::scale(log1p(.data[[glue("PopDen_{scale}")]]))[, 1],
        "{glue('elev_{scale}_t')}" := base::scale(sqrt(.data[[glue("elev_{scale}")]]))[, 1]
      )
    
    global_model <- lme4::lmer(
      as.formula(glue(
        "log(RA_weighted) ~
     EH_{scale}_t + Temp_{scale}_t + Prec_{scale}_t + NPP_{scale}_t + PopDen_{scale}_t + elev_{scale}_t +
     (1 + EH_{scale}_t | SCIENTIFIC_NAME)"
      )),
      data = dat
    )
    
    conditional_R2 = r2(global_model)[[1]]
    
    fixed_ci <- broom.mixed::tidy(global_model,
                                  effects = "fixed",
                                  conf.int = TRUE,
                                  conf.method = "Wald")
    beta_EH = fixed_ci[2,'estimate']
    conf.low = fixed_ci[2,'conf.low']
    conf.high = fixed_ci[2,'conf.high']
    
    
    eh_effect <- ggeffects::ggpredict(global_model, terms = glue("EH_{scale}_t"))
    ggplot(eh_effect, aes(x = x, y = predicted)) +
      geom_line(color = "#993C1D", linewidth = 1.2) +
      geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
                  fill = "#993C1D", alpha = 0.2) +
      labs(x = "Environmental heterogeneity (Standardized)",
           y = "Predicted species-level relative abundance") +
      annotate("text",
               x     = -Inf,
               y     = Inf,
               label = sprintf(
                 "No. of species (> 100 obs.) = %d\nConditional R^2 = %.2f\nβ = %.3f (95%% CI: %.3f, %.3f)",
                 n_sp, conditional_R2, beta_EH, conf.low, conf.high
               ),
               hjust = -0.1,    # small positive = slight right offset from left edge
               vjust = 1.5,     # positive = slight down offset from top edge
               size = 4, color = "black") +
      theme_classic(base_size = 13) +
      theme(
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 11)
      )
  }
  
  
}

##########################supplementary figure: marginal effect of EH on species-level abundance #################
p_tree <- lapply(c("1km","5km","10km","20km","50km"),function(s) marginal_effect(taxa = "tree",scale=s) +
                   labs(title = s, size = 14) +
                   theme(plot.title = element_text(hjust = 0.5, size = 14)))
p_bird <- lapply(c('400m',"1km","5km","10km","20km"),function(s) marginal_effect(taxa = "bird",scale=s) +
                   labs(title = s) +
                   theme(plot.title = element_text(hjust = 0.5, size = 14)))
p_mammal <- lapply(c("1km","5km","10km","20km","50km"),function(s) marginal_effect(taxa = "mammal",scale=s) +
                     labs(title = s) +
                     theme(plot.title = element_text(hjust = 0.5, size = 14)))

#layout
p_tree <- lapply(seq_along(p_tree), function(i) {
  
  p_tree[[i]] + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())
})

p_bird <- lapply(seq_along(p_bird), function(i) {
  if (i == 1) {
    p_bird[[i]] + theme(axis.title.x = element_blank())
  } else {
    p_bird[[i]] + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())
  }
})

p_mammal <- lapply(seq_along(p_mammal), function(i) {
  if (i == 3) {
    p_mammal[[i]] + theme(axis.title.y = element_blank())
  } else {
    p_mammal[[i]] + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank())
  }
})

# combine
wrap_plots(c(p_tree, p_bird,p_mammal), nrow = 3)

# TODO: save the composite, e.g.
# ggsave(p_fig("FigS13_pooled_EH_abundance.png"), wrap_plots(c(p_tree, p_bird, p_mammal), nrow = 3),
#        width = 16, height = 10, dpi = 300)

# =============================================================================
# Fig. S14 - trait effects across ALL scales                (NOT YET SCRIPTED)
# Fig. S15 - raw vs imputed tree traits across scales       (NOT YET SCRIPTED)
# -----------------------------------------------------------------------------
# TODO(reviewer, Fig. S14/S15): the existing code fits the trait meta-regression
#   at one scale per taxon only (01_species_EH_trait_models.R). Fig. S14-S15
#   need a loop over scales that rebuilds EH_trait_plot per scale and calls
#   tree_beta_trait() / bird_beta_trait() / mammal_beta_trait()
#   (code/functions/fun_beta_trait.R), once with raw and once with imputed tree
#   traits. Per-scale results exist in EmbeddingBio/data_clean/
#   {tree,bird,mammal}_data_<scale>.csv, tree_data_<scale>_imputed.csv and
#   EmbeddingBio/output/supplementarymaterials/
#   betaEHwithTraits_unimputed_imputed.xlsx, but the script that wrote them was
#   not found. Recover it or write the loop, then plot S14/S15 from those CSVs.
# TODO(reviewer, Fig. S15 caption): imputation was across 864 tree species
#   (01_tree_traits_imputation.R); the caption's "4,999 species" is the mammal count.
# =============================================================================
