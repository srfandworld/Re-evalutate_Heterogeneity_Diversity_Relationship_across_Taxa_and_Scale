# =============================================================================
# 01_tree_FIA_download.R
# Download FIA, compute plot-level basal area per species, export plot points
#
# Source   : Code/Tree/01_FIA_data_download_cleanning.R, section 1
# Inputs   : FIA DataMart (downloaded by rFIA::getFIA into data/raw/FIA/)
# Outputs  : data/clean/tpa_all5.rds     plot x species basal area (DIA > 5 in)
#            data/clean/FIA_plots.shp    plot locations -> 02_environment
#            data/clean/fia_sp.csv       SPCD -> scientific name (868 species)
# Next     : 02_environment/02_environmental_covariates.R needs FIA_plots.shp
# =============================================================================

source(here::here("code", "00_config.R"))

library(rFIA)         # getFIA(), readFIA(), tpa(), clipFIA()
library(sf)
library(tidyverse)

options(timeout = 3600)   # long download timeout for the FIA pull

# ---- 1. Download / load FIA -------------------------------------------------
# getFIA() only needs to run once; afterwards readFIA() loads the local copy.
# TODO(reviewer M1): this list includes AK, HI and territories (GU, FM, MP, PW,
# AS, PR, VI). Non-CONUS plots are dropped later (models are CONUS only), but
# the species tally (868) in the Methods should state the CONUS filter used.
fia_dir <- p_raw("FIA")
dir.create(fia_dir, showWarnings = FALSE, recursive = TRUE)
if (length(list.files(fia_dir, pattern = "_TREE\\.csv$")) == 0) {
  db <- getFIA(
    states = c("KS","LA","MN","MO","MT","NE","NM","ND","OK","SD","TX","UT","WY",
               "AL","CT","DE","FL","GA",
               "AK","CA","HI","ID","NV","OR","WA","AZ","AR","CO","IA",
               "IL","IN","KY",
               "ME","MD","MA",
               "MI","MS","NH","NJ","NY","NC","OH","PA","RI","SC","TN","VT","VA","WV","WI",
               "GU","FM","MP","PW","AS","PR","VI"),
    dir = fia_dir)
}
db <- readFIA(dir = fia_dir)

# Most recent survey per plot; basal area per species for trees with DIA > 5
# TODO(reviewer M2): Methods say "at least 5 inches"; code uses DIA > 5.
tpa_all5 <- tpa(clipFIA(db), bySpecies = TRUE, byPlot = TRUE,
                treeDomain = DIA > 5, returnSpatial = TRUE)
saveRDS(tpa_all5, p_clean("tpa_all5.rds"))

# ---- 2. Plot locations (consumed by the environmental-covariate script) -----
FIA_plots <- tpa_all5 %>%
  dplyr::select(pltID, PLT_CN) %>%
  distinct()
st_write(FIA_plots, p_clean("FIA_plots.shp"), append = FALSE)

# ---- 3. Species lookup: drop genus-only / unidentified records --------------
fia_spcd_to_name <- tpa_all5 %>%
  st_drop_geometry() %>%
  filter(!is.na(SPCD), !is.na(SCIENTIFIC_NAME)) %>%
  distinct(SPCD, SCIENTIFIC_NAME)
fia_sp <- fia_spcd_to_name %>%
  filter(!str_detect(SCIENTIFIC_NAME, "spp\\.|Tree"))
n_distinct(fia_sp$SCIENTIFIC_NAME)   # 868
write.csv(fia_sp, p_clean("fia_sp.csv"), row.names = FALSE)
