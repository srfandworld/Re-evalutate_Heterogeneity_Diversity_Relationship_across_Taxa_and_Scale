# Environmental heterogeneity (EH) from AlphaEarth embeddings: Earth Engine step

EH (manuscript Eq. 2) is the mean pairwise cosine dissimilarity among 10 m
AlphaEarth Foundations embedding pixels within each buffer, averaged over
2017–2023. It is computed in Google Earth Engine, not in R. The EE scripts live
in the `/gee` repository linked from the main README.

> **TODO (reviewer C3):** The EE scripts are not in this folder and were not
> found on the analysis machine. Before submission:
> 1. Add the EE code (or a read-only repository link) to `/gee`.
> 2. Name the exact script that produced each export listed below.

## Exports the R pipeline expects (`data/raw/gee_exports/`)

| File(s) | Used by | Unit / key |
|---|---|---|
| `FIA_embeddings/FIA_Dissimilarity_{1,5,10,20}km_2017_2023_mean_TEST.shp`, `..._50km_..._TEST.csv` | `02_environmental_covariates.R` | FIA plot, `pltID`; field `mean_cosin` |
| `BBS_embeddings/BBS_Dissimilarity_{400m,1km,5km,10km,20km}_2017_2023_mean.shp` | `02_environmental_covariates.R` | BBS route segment; field `mean_cosin` (joined by **row order** to `export_table.csv`, see C3) |
| `export_table.csv` | `01_data_acquisition/02_bird_BBS_routes.R` | `BBSRTSL020`, `RTENO`, `RTENAME` |
| `Camera_Trap_Array_EH_2019-2023.csv` | `02_environmental_covariates.R` | `Year` × `Camera_Trap_Array` |
| `camera_array_centroids_NPP_Climate_1_5_10_20_50km/*.shp` | `01_CTA_climate_NPP.R` | `Cmr_T_A`; columns `<var>_<year>` |
| `EH_US_5km_4326.shp` | Fig. 2 and Fig. S6 maps | 5 km CONUS grid; field `EH` |
| `*_A57_EH/`, `*_EH_delta/` | Fig. S17–S19 (alternative EH metrics) | see `07_supplementary_figures/FigS17-S19_alternative_EH.R` |

Open questions to settle in the EE documentation:
- What does the `_TEST` suffix on the FIA files mean, and are they the final runs?
- Which embedding year range and pixel size were used for each export?
- Was the BBS export keyed by `BBSRTSL020`? If so, the row-order `cbind(a)` join
  can be replaced with a key join.
