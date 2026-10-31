# =============================================================================
# FigS05_niche_env_layers.R
# Fig. S5: global maps of the 14 environmental layers defining the niche space
# for niche-breadth estimation (4 columns, viridis, one panel per layer).
#
# Source   : EmbeddingBio/analysis/sfigure5_niche_env_layers.R (2026-10-31;
#            paths only changed). Not yet run.
# Inputs   : data/clean/envstack_for_FIAniche.tif (~10 GB;
#            04_traits_phylogeny_niche/05_niche_space_PCA.R)
# Outputs  : output/figures/FigureS5_niche_env_layers.png
# Options  : rescale_soil (SoilGrids storage units -> caption units: pH, g/kg),
#            crop_antarctica (drop < 60 S); plotted from a regular sample of
#            ~1e6 cells per layer.
# =============================================================================

source(here::here("code", "00_config.R"))
library(terra)

out_png <- p_fig("FigureS5_niche_env_layers.png")

rescale_soil    <- TRUE
crop_antarctica <- TRUE
n_cells         <- 1e6   # cells per layer used for plotting

vars <- c("MAT","MAP","Tseas","Tmax","Tmin","Pmax","Pmin","Pseas",
          "SoilPH","SOC","SoilNitro","SoilMoist","Elev","Slope")

#panel titles with units (as in the Fig. S5 caption)
titles <- c(MAT       = "MAT (°C)",
            MAP       = "MAP (mm)",
            Tseas     = "Tseas (SD x 100)",
            Tmax      = "Tmax (°C)",
            Tmin      = "Tmin (°C)",
            Pmax      = "Pmax (mm)",
            Pmin      = "Pmin (mm)",
            Pseas     = "Pseas (CV)",
            SoilPH    = if (rescale_soil) "SoilPH (pH)"           else "SoilPH (pH x 10)",
            SOC       = if (rescale_soil) "SOC (g/kg)"            else "SOC (dg/kg)",
            SoilNitro = if (rescale_soil) "SoilNitro (g/kg)"      else "SoilNitro (cg/kg)",
            SoilMoist = "SoilMoist (m³/m³)",
            Elev      = "Elev (m)",
            Slope     = "Slope (°)")


##########################1. read and sample the stack#########################################
env_stack <- rast(p_clean("envstack_for_FIAniche.tif"))
names(env_stack) <- vars

env_plot <- spatSample(env_stack, size = n_cells, method = "regular", as.raster = TRUE)
if (crop_antarctica) env_plot <- crop(env_plot, ext(-180, 180, -60, 90))

if (rescale_soil) {
  env_plot[["SoilPH"]]    <- env_plot[["SoilPH"]]    / 10    # pH x 10 -> pH
  env_plot[["SOC"]]       <- env_plot[["SOC"]]       / 10    # dg/kg   -> g/kg
  env_plot[["SoilNitro"]] <- env_plot[["SoilNitro"]] / 100   # cg/kg   -> g/kg
}


##########################2. Figure S5#########################################################
png(out_png, width = 16, height = 11, units = "in", res = 300, bg = "white")
plot(env_plot, nc = 4, main = unname(titles[names(env_plot)]),
     col = hcl.colors(100, "viridis"), cex.main = 1.4, mar = c(2, 2, 2.5, 5))
dev.off()
