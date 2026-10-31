# =============================================================================
# 01_fit_HDR_models.R
# Richness / MSA ~ EH + EH^2 + covariates, per taxon and spatial scale; SAR
# models (spatialreg::spautolm) with row-standardised 8-nearest-neighbour weights:
#   trees   - FIA plot locations (EPSG:5070), method = "LU"
#   birds   - BBS route start points (EPSG:5070)
#   mammals - camera-trap-array centroids (one row per array; lon/lat)
#
# Source   : Code/Figure2/01_3taxa_modeling.R (2026-10-30; paths only changed).
#            Data, weights and model functions: functions/fun_hdr_sar_models.R
# Inputs   : data/clean/FIA_table_<s>_iNEXT.csv, BBS_table.csv, mammal_table.csv
#            data/raw/gee_exports/{FIA_embeddings/FIA_EH_delta, FIA_A57_EH,
#              BBS_EH_delta, BBS_A57_EH, mammal_EH_delta, mammal_A57_EH}
#            data/clean/FIA_plots.shp, bbsroutes_shp/bbsrtsl020_filterd.shp,
#              camera_array_centroids.shp
# Outputs  : data/clean/model_data/{FIA,BBS,mammal}_model_<s>.csv
#            data/clean/models/sar_{tree,bird,mammal}_{richness,MSA}*.rds
# Feeds    : Fig. 2, Tables S2-S4, Figs. S9-S11
# NOTE     : tree fits (18,588-129,178 plots) take ~25 min; run in a fresh R
#            session (the LU factorisation needs several GB of free memory).
# =============================================================================

source(here::here("code", "00_config.R"))
source(here::here("code", "functions", "fun_hdr_sar_models.R"))   # data, weights, model functions

# =============================================================================
# 4. Run all models (named by scale) and save
# =============================================================================

# ---- Trees ----
sar_tree_richness <- setNames(lapply(scales_tree, sar_modeling_tree, y = "richness"),  scales_tree)
sar_tree_MSA      <- setNames(lapply(scales_tree, sar_modeling_tree, y = "log(MSA)"),  scales_tree)
saveRDS(sar_tree_richness, p_model("sar_tree_richness_iNEXT.rds"))
saveRDS(sar_tree_MSA,      p_model("sar_tree_MSA_iNEXT.rds"))

# ---- Birds ----
sar_bird_richness <- setNames(lapply(scales_bird, sar_modeling_bird, y = "richness"),  scales_bird)
sar_bird_MSA      <- setNames(lapply(scales_bird, sar_modeling_bird, y = "log(MSA)"),  scales_bird)
saveRDS(sar_bird_richness, p_model("sar_bird_richness.rds"))
saveRDS(sar_bird_MSA,      p_model("sar_bird_MSA.rds"))

# ---- Mammals ----
sar_mammal_richness <- setNames(lapply(scales_mammal, sar_modeling_mammal, y = "richness_Estimator"), scales_mammal)
sar_mammal_MSA      <- setNames(lapply(scales_mammal, sar_modeling_mammal, y = "log(MSA)"),            scales_mammal)
saveRDS(sar_mammal_richness, p_model("sar_mammal_richness.rds"))
saveRDS(sar_mammal_MSA,      p_model("sar_mammal_MSA.rds"))
