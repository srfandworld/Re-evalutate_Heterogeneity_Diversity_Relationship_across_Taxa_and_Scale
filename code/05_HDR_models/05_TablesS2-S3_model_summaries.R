# =============================================================================
# 05_TablesS2-S3_model_summaries.R  -  STUB (reviewer: no table/ΔAIC/
# classification code located)
# Table S2 (richness) and Table S3 (log MSA): EH and EH^2 coefficients, SE, p,
# R2, AIC, ΔAIC (linear - quadratic), turning point and HDR shape, for every
# taxon x scale.
#
# The reviewer matched all 30 fits to the saved models, so the numbers come
# from data/clean/models/sar_*.rds. What's missing is the code that:
#   1. refits each model WITHOUT the EH^2 term and computes
#      ΔAIC = AIC(linear) - AIC(quadratic);
#   2. computes turning points -b1 / (2 b2) and back-transforms them to raw EH
#      (use the SAME mean/SD as the fitting data, reviewer C4);
#   3. applies the HDR classification rule (unimodal / U-shaped / linear / none).
#      TODO(reviewer): state the rule; ΔAIC of 1.8-2.0 (mammal richness 5 km,
#      mammal MSA 50 km) is weak support for "unimodal".
#   4. writes output/tables/TableS2_richness.csv and TableS3_MSA.csv.
# TODO(reviewer C2): Results say ΔAIC ranges 18-8,756.5, but Tables S2-S3 reach
#   20,104.6 (8,756.5 belongs to the old non-iNEXT tree MSA model). Update the
#   text from the regenerated tables.
# Helpers: inflection_boot() in code/functions/fun_sar_partial_plot.R
#   bootstraps turning points.
# =============================================================================

source(here::here("code", "00_config.R"))
source(here::here("code", "functions", "fun_sar_partial_plot.R"))
