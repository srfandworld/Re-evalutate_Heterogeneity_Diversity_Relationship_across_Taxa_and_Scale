# =============================================================================
# run_all.R - run the pipeline in dependency order
#
# Each stage reads only files written by earlier stages (see README.md), so
# stages can also be run one at a time. Slow / online steps are off by default;
# with data/ staged by 00_stage_data.R, the analysis stages (05-07) run from
# the saved intermediates.
# =============================================================================

source(here::here("code", "00_config.R"))

RUN <- list(
  acquisition  = FALSE,  # FIA / BBS downloads (hours, needs internet)
  environment  = FALSE,  # raster buffer extraction (hours; needs HDR_EXT_DIR)
  community    = FALSE,  # FIA pooling at 1-110 km + iNEXT (many hours)
  traits_niche = FALSE,  # TRY/BIEN/TNRS online; KDE niche breadth (hours)
  hdr_models   = TRUE,   # Fig. 2, Tables S2-S4, Figs. S9-S12
  species      = TRUE,   # Fig. 3, Table S5, Figs. S13-S15
  supplement   = TRUE    # Figs. S1, S5-S8, S16-S19
)

run <- function(...) {
  f <- here::here("code", ...)
  message("\n==== ", file.path(...), " ====")
  source(f, local = new.env(), echo = FALSE)
}

if (RUN$acquisition) {
  run("01_data_acquisition", "01_tree_FIA_download.R")
  run("01_data_acquisition", "02_bird_BBS_routes.R")
  run("01_data_acquisition", "03_mammal_SNAPSHOT_sites.R")
}
if (RUN$environment) {
  # EH itself comes from Earth Engine: see 02_environment/README_GEE.md
  run("02_environment", "01_CTA_climate_NPP.R")
  run("02_environment", "02_environmental_covariates.R")
}
if (RUN$community) {
  run("03_community_metrics", "01_tree_pool_communities.R")
  run("03_community_metrics", "02_tree_richness_iNEXT.R")
  run("03_community_metrics", "03_bird_community.R")
  run("03_community_metrics", "04_mammal_community.R")
}
if (RUN$traits_niche) {
  run("04_traits_phylogeny_niche", "01_tree_traits_imputation.R")
  run("04_traits_phylogeny_niche", "02_bird_traits_phylogeny.R")
  run("04_traits_phylogeny_niche", "03_mammal_traits_imputation.R")
  run("04_traits_phylogeny_niche", "04_species_ranges.R")
  run("04_traits_phylogeny_niche", "05_niche_space_PCA.R")
  run("04_traits_phylogeny_niche", "06_niche_breadth_KDE.R")
}
if (RUN$hdr_models) {
  run("05_HDR_models", "01_fit_HDR_models.R")
  run("05_HDR_models", "02_Figure2.R")
  run("05_HDR_models", "03_FigS9-S11_partial_residuals.R")
  run("05_HDR_models", "04_TableS4_LMG.R")             # Fig. S12 donuts: 04b (TODO)
  run("05_HDR_models", "05_TablesS2-S3_model_summaries.R")   # refits linear SAR once (cached)
}
if (RUN$species) {
  run("06_species_trait_responses", "01_species_EH_trait_models.R")
  run("06_species_trait_responses", "02_Figure3_forestplot.R")
  run("06_species_trait_responses", "03_FigS13_species_level_scales.R")
  run("06_species_trait_responses", "04_FigS14_forestplot_multiscale.R")   # needs *_data_<scale>.csv for all scales
  run("06_species_trait_responses", "05_FigS15_forestplot_treeimputation.R")   # needs tree_data_<scale>(_imputed).csv
}
if (RUN$supplement) {
  run("07_supplementary_figures", "FigS01_phylogenies.R")
  run("07_supplementary_figures", "FigS05_niche_env_layers.R")
  run("07_supplementary_figures", "FigS06_EH_characterization.R")
  run("07_supplementary_figures", "FigS07_EH_vs_conventional_metrics.R")
  run("07_supplementary_figures", "FigS08_EH_across_scales.R")
  run("07_supplementary_figures", "FigS16_total_community_abundance.R")   # SAR, cached
  run("07_supplementary_figures", "FigS17_tree_alternative_EH.R")    # SAR, slow
  run("07_supplementary_figures", "FigS18_bird_alternative_EH.R")
  run("07_supplementary_figures", "FigS19_mammal_alternative_EH.R")
}

sessionInfo()
