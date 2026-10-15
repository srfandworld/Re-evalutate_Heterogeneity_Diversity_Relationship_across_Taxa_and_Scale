# =============================================================================
# FigS16_total_community_abundance.R  -  STUB (no code located)
# Fig. S16: standardised EH effect on log total community abundance
#           (= richness x MSA) across scales, same covariates as the HDR models
#
# Searched EmbeddingBio/analysis, Code/ and the user's Downloads/Documents
# (2026-10-15) - no script found. Partial clue: 05_HDR_models/02_Figure2.R
# creates `total_niche = richness * MSA` for birds.
#
# To implement: reuse the model functions in 05_HDR_models/01_fit_HDR_models.R
# with the response log(richness * MSA) (trees: log(qD * MSA_new_<s>)) and a
# linear EH term only, collect the EH coefficient + 95% CI per taxon x scale,
# and draw a forest plot -> output/figures/FigS16_total_abundance.png.
# TODO(reviewer): mammal estimates at 10/20/50 km overlap zero, so the caption's
#   "Across all three taxa, EH had a positive effect" is too general.
# =============================================================================

source(here::here("code", "00_config.R"))
