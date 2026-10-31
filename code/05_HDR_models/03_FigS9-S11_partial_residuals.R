# =============================================================================
# 03_FigS9-S11_partial_residuals.R
# Figs. S9 (trees), S10 (birds), S11 (mammals): partial residual plots of species
# richness (top row) and log MSA (bottom row) against standardized EH across five
# spatial scales, from the SAR models of 01_fit_HDR_models.R
#
# Source   : EmbeddingBio/analysis/figureS9_S11.R, sections 2-3 (verbatim). That
#            script refits the same models via the physicallyVar functions; here
#            the saved models are read instead (identical data, weights and
#            specification since 2026-10-30).
# Inputs   : data/clean/models/sar_{tree,bird,mammal}_{richness,MSA}*.rds
# Outputs  : output/figures/FigureS9_tree_partial_residuals.png,
#            FigureS10_bird_partial_residuals.png, FigureS11_mammal_partial_residuals.png
# NOTE     : R^2 = cor(y, fitted)^2 includes the SAR spatial term (tree values
#            0.4-1.0); decide on the R^2 definition together with Table S2.
# =============================================================================

source(here::here("code", "00_config.R"))
source(here::here("code", "functions", "fun_sar_partial_plot.R"))
library(stringr)
library(patchwork)
library(spatialreg)   # summary() / coef() methods for spautolm objects

read_models <- function(file, scales) {
  m <- readRDS(p_model(file))
  if (is.null(names(m))) names(m) <- scales
  m
}
sar_models <- list(
  tree   = list(richness = read_models("sar_tree_richness_iNEXT.rds", SCALES_TREE),
                MSA      = read_models("sar_tree_MSA_iNEXT.rds",      SCALES_TREE)),
  bird   = list(richness = read_models("sar_bird_richness.rds",       SCALES_BIRD),
                MSA      = read_models("sar_bird_MSA.rds",            SCALES_BIRD)),
  mammal = list(richness = read_models("sar_mammal_richness.rds",     SCALES_MAMMAL),
                MSA      = read_models("sar_mammal_MSA.rds",          SCALES_MAMMAL))
)

##########################2. richness + MSA across scales for one taxon######################
#same code as the "#4. Supplementary figures" block of figure2_fit.R
plot_scales <- function(models_richness, models_MSA,
                        y_lab_richness = "Residual richness",
                        y_lab_MSA      = "Residual log (mean species abundance)") {
  n <- length(models_richness)

  p_richness <- lapply(models_richness, function(m){
    var <- colnames(m$X)[2]                       # e.g. "EH_10km_t"
    scale <- str_extract(var, "\\d+(km|m)")
    sar_partial_plot(m, var = var, quadratic = TRUE,
                     line_color = 'steelblue',
                     x_label   = "Environmental heterogeneity \n(standardized)",
                     y_label   = y_lab_richness) +
      labs(title = scale) +
      theme(plot.title = element_text(hjust = 0.5))
  })
  p_MSA <- lapply(models_MSA, function(m){
    var <- colnames(m$X)[2]
    sar_partial_plot(m, var = var, quadratic = TRUE,
                     line_color = 'firebrick',
                     x_label   = "Environmental heterogeneity \n(standardized)",
                     y_label   = y_lab_MSA)
  })

  #layout
  # remove y-axis titles except first column
  p_richness <- lapply(seq_along(p_richness), function(i) {
    if (i %% n != 1) p_richness[[i]] + theme(axis.title.y = element_blank()) else p_richness[[i]]
  })
  p_MSA <- lapply(seq_along(p_MSA), function(i) {
    if (i %% n != 1) p_MSA[[i]] + theme(axis.title.y = element_blank()) else p_MSA[[i]]
  })
  # remove x-axis titles for first row (richness row)
  p_richness <- lapply(p_richness, function(p) p + theme(axis.title.x = element_blank()))

  # combine
  wrap_plots(c(p_richness, p_MSA), nrow = 2)
}


##########################3. Figures S9-S11#####################################################
#tree richness is modelled as log(qD) (coverage-based richness, iNEXT)
fig_S9  <- plot_scales(sar_models$tree$richness,   sar_models$tree$MSA,
                       y_lab_richness = "Residual log (richness)")
fig_S10 <- plot_scales(sar_models$bird$richness,   sar_models$bird$MSA)
fig_S11 <- plot_scales(sar_models$mammal$richness, sar_models$mammal$MSA)

ggsave(p_fig("FigureS9_tree_partial_residuals.png"),    fig_S9,  width = 20, height = 8, dpi = 300)
ggsave(p_fig("FigureS10_bird_partial_residuals.png"),   fig_S10, width = 20, height = 8, dpi = 300)
ggsave(p_fig("FigureS11_mammal_partial_residuals.png"), fig_S11, width = 20, height = 8, dpi = 300)
