# =============================================================================
# 04_TableS4_LMG.R
# Table S4: relative importance (% of model-explained variance) of each predictor
# for species richness and log MSA, trees / birds / mammals, five scales (LMG).
#
# Source   : EmbeddingBio/analysis/tableS4_LMG.R (2026-10-30; paths only changed)
# LMG (relaimpo) is defined for ordinary least squares only, so each SAR model of
# 01_fit_HDR_models.R is refitted with lm() using the SAME formula on the SAME rows
# (data/clean/model_data/*_model_<scale>.csv) and the lm R^2 is decomposed.
# EH and EH^2 are summed into one "EH" share.
# Methods: "Because LMG is defined for ordinary least squares, relative importance
# was computed for OLS models with the same predictors and observations as the
# corresponding SAR models."
# Outputs  : output/tables/TableS4_LMG.csv
# =============================================================================

source(here::here("code", "00_config.R"))
library(dplyr)
library(tidyr)
library(glue)

data_dir <- p_clean("model_data")
out_dir  <- DIR_TAB

#response and covariates per taxon: same formulas as 01_3taxa_modeling.R
spec <- list(
  tree   = list(prefix = "FIA",    scales = c("1km", "5km", "10km", "20km", "50km"),
                lhs = c(richness = "log(qD)", MSA = "log(MSA_new_{scale})"),
                covars = "MAT_{scale}_t + MAP_{scale}_t + soilph_{scale}_t + soilcec_{scale}_t + elev_mean_{scale}_t"),
  bird   = list(prefix = "BBS",    scales = c("400m", "1km", "5km", "10km", "20km"),
                lhs = c(richness = "richness", MSA = "log(MSA)"),
                covars = "MAT_{scale}_t + MAP_{scale}_t + Tsea_{scale}_t + Psea_{scale}_t + elev_{scale}_t + PopDen_{scale}_t"),
  mammal = list(prefix = "mammal", scales = c("1km", "5km", "10km", "20km", "50km"),
                lhs = c(richness = "richness_Estimator", MSA = "log(MSA)"),
                covars = "elev_{scale}_t + PopDen_{scale}_t + NPP_{scale}_t + Prec_{scale}_t + Temp_{scale}_t")
)

#predictor term -> column name in Table S4
predictor_name <- function(term) {
  case_when(
    grepl("^(I\\()?EH_", term)            ~ "EH",
    grepl("^elev", term)                  ~ "Elevation",
    grepl("^(MAP|Prec)_", term)           ~ "MAP",
    grepl("^(MAT|Temp)_", term)           ~ "MAT",
    grepl("^soilcec_", term)              ~ "SoilCEC",
    grepl("^soilph_", term)               ~ "SoilPH",
    grepl("^Tsea_", term)                 ~ "Tsea",
    grepl("^Psea_", term)                 ~ "Psea",
    grepl("^PopDen_", term)               ~ "PopDen",
    grepl("^NPP_", term)                  ~ "NPP"
  )
}

#LMG shares (%) for one taxon, response and scale
lmg_row <- function(taxa, response, scale) {
  s    <- spec[[taxa]]
  dat  <- read.csv(file.path(data_dir, glue("{s$prefix}_model_{scale}.csv")))
  form <- as.formula(glue("{glue(s$lhs[[response]])} ~ EH_{scale}_t + I(EH_{scale}_t^2) + {glue(s$covars)}"))
  m    <- lm(form, data = dat)

  ri <- relaimpo::calc.relimp(m, type = "lmg", rela = TRUE)@lmg   # proportions, sum to 1
  data.frame(term = names(ri), share = 100 * unname(ri)) %>%
    mutate(predictor = predictor_name(term)) %>%
    group_by(predictor) %>%
    summarise(share = sum(share), .groups = "drop") %>%      # EH + EH^2 -> EH
    mutate(Taxa = taxa, Response = response, Scale = scale,
           R2_lm = 100 * summary(m)$r.squared, n = nrow(dat))
}

table_S4_long <- bind_rows(lapply(names(spec), function(taxa) {
  bind_rows(lapply(c("richness", "MSA"), function(resp) {
    bind_rows(lapply(spec[[taxa]]$scales, function(s) lmg_row(taxa, resp, s)))
  }))
}))

#wide format as in the SI (one row per taxon x response x scale)
table_S4 <- table_S4_long %>%
  mutate(share = round(share, 1)) %>%
  pivot_wider(names_from = predictor, values_from = share) %>%
  relocate(Taxa, Response, Scale, EH, Elevation, MAP, MAT,
           any_of(c("SoilCEC", "SoilPH", "Tsea", "Psea", "PopDen", "NPP")), R2_lm, n)

write.csv(table_S4, file.path(out_dir, "TableS4_LMG.csv"), row.names = FALSE)
table_S4
