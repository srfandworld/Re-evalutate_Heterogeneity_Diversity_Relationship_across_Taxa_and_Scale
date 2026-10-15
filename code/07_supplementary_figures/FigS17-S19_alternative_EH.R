# =============================================================================
# FigS17-S19_alternative_EH.R  -  STUB (no code located)
# Figs. S17 (trees), S18 (birds), S19 (mammals): HDR with alternative EH
#   1. elevation range                      (FIA_table: elev_range_<s>)
#   2. mean-independent elevation range (delta, Eq. S1)
#        data/raw/gee_exports/{FIA,BBS,mammal}_EH_delta/*_Elevation_EHdelta_<s>.csv
#   3. A57 EH (alternative embedding dissimilarity)
#        data/raw/gee_exports/{FIA,BBS,mammal}_A57_EH/*_A57_EH_<s>.shp
#   4. AEF EH (main metric), for comparison
#
# The input layers exist in EmbeddingBio/data/, but no script reading them was
# found (searched 2026-10-15).
# TODO(reviewer): tree panels use log(richness/effort) and "multiple linear
#   regressions", which matches the OLD (non-iNEXT) tree models. Re-run with the
#   final model set (TREE_MODEL_SET = "iNEXT") so S17 matches Fig. 2 / Tables S2-S3.
# To implement: join each metric to the per-scale tables, refit the
#   01_fit_HDR_models.R specification with EH swapped for the metric, and plot
#   partial residuals with fun_sar_partial_plot.R.
# =============================================================================

source(here::here("code", "00_config.R"))
