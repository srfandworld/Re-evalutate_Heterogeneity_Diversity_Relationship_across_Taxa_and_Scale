# =============================================================================
# 00_packages.R - install the packages used by the pipeline
#
# EmbeddingBio/renv.lock was never snapshotted (it only lists renv), so there
# are no pinned versions yet. After a successful run on the analysis machine:
#   renv::init(bare = TRUE); renv::snapshot()
# and commit the resulting renv.lock (reviewer: "required packages/versions").
# Analysis was run on R 4.5.1.
# =============================================================================

cran <- c(
  # data handling
  "here", "tidyverse", "data.table", "glue", "readxl", "abind", "units",
  # spatial
  "sf", "terra", "lwgeom", "spatialEco", "exactextractr", "maps", "fields",
  "spdep", "spatialreg",
  # biodiversity data
  "rFIA", "iNEXT", "BIEN", "TNRS", "rtry",
  # traits / phylogeny
  "ape", "phytools", "Rphylopars", "FD", "nlme",
  # niche breadth
  "ks",
  # models
  "lme4", "broom", "broom.mixed", "performance", "relaimpo", "ggeffects",
  "emmeans", "multcomp", "multcompView", "corrplot",
  # plotting
  "ggplot2", "patchwork", "cowplot", "ggpubr", "ggtext", "ggimage", "ggrastr",
  "gtable", "magick", "rsvg", "ragg", "paletteer", "scales", "viridis"
)
install.packages(setdiff(cran, rownames(installed.packages())))

# Non-CRAN
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

install.packages("bbsBayes2", repos = c(bbsbayes = "https://bbsbayes.r-universe.dev",
                                        CRAN = getOption("repos")))   # BBS data
remotes::install_github("jinyizju/V.PhyloMaker2")                    # tree phylogeny
remotes::install_github("OpenTreeOfLife/clootl")                     # bird phylogeny (McTavish et al. 2025)
BiocManager::install("ggtree")                                       # Fig. S1
