# =============================================================================
# Mammals (SNAPSHOT USA): camera arrays, relative abundance, iNEXT richness,
# and per-scale tables
# Re-evaluating the Environmental Heterogeneity-Diversity Relationship
# SNAPSHOT USA reference: https://doi.org/10.1111/geb.13941
#
# Sections:
#   1. Load SNAPSHOT data; build camera sites and array centroids
#   2. Filter to wild mammal detections
#   3. Relative abundance (RA) + body mass (EltonTraits) + weighting
#   4. Richness via iNEXT (incidence-based Chao estimators)
#   5. Assemble mammal_table (diversity + abundance + environment)
#   6. Per-scale standardised tables (1, 5, 10, 20, 50 km) + export
#
# NOTE: Section 1 WRITES data_clean/camera_array_centroids.shp, which the
#       environmental-covariate script needs; that script in turn produces
#       Camera_Trap_Array_env.csv, which Section 5 reads. So the run order is:
#       Section 1 here -> covariate script -> Sections 2-6 here.
# =============================================================================

# ---- Libraries --------------------------------------------------------------
library(tidyverse)   # dplyr, tidyr, stringr, purrr, ggplot2, readr, tibble
library(sf)
library(iNEXT)       # iNEXT(), ggiNEXT()  (was loaded mid-script)


dir.create("data_clean", showWarnings = FALSE)

# =============================================================================
# 1. Load SNAPSHOT data; camera sites and array centroids
# =============================================================================
ssusa_seq <- read.csv("data/Mammals/ssusa_finalsequences.csv")
ssusa_dep <- read.csv("data/Mammals/ssusa_finaldeployments.csv")

# Camera site points
mammal_sites <- ssusa_dep %>%
  distinct(Camera_Trap_Array, Site_Name, Deployment_ID, Latitude, Longitude) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)

# Dominant habitat per camera-trap array
CTA_habitat <- ssusa_dep %>%
  count(Camera_Trap_Array, Habitat) %>%
  group_by(Camera_Trap_Array) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(Camera_Trap_Array, Habitat)

# Array centroids (union of member sites), with habitat
camera_array_centroids <- mammal_sites %>%
  group_by(Camera_Trap_Array) %>%
  summarise(geometry = st_centroid(st_union(geometry))) %>%
  left_join(CTA_habitat, by = "Camera_Trap_Array")
write_sf(camera_array_centroids, "data_clean/camera_array_centroids.shp")
camera_array_centroids <- read_sf("data_clean/camera_array_centroids.shp") %>%
  rename(Camera_Trap_Array = Cmr_T_A)   # shapefile truncates the field name

# =============================================================================
# 2. Filter to wild mammal detections
# =============================================================================
mammal_seq <- ssusa_seq %>%
  filter(
    Class == "Mammalia",
    if_all(c(Order, Family, Genus, Species), ~ !is.na(.x)),  # was complete.cases(across(...))
    Genus != " ", Species != " ", Genus != "Homo",
    !grepl("Domestic", Common_Name, ignore.case = TRUE)
  ) %>%
  mutate(SCIENTIFIC_NAME = paste0(Genus, " ", Species))

n_distinct(mammal_seq$SCIENTIFIC_NAME)   # 122 species (131 - 1 human - 8 domestic types)
mammal_scinam_comnam <- mammal_seq %>% distinct(SCIENTIFIC_NAME, Common_Name)

# Deployments retained for those detections
mammal_dep <- ssusa_dep %>%
  filter(Deployment_ID %in% unique(mammal_seq$Deployment_ID))

# =============================================================================
# 3. Relative abundance (RA), body mass, and weighting
#    RA / detection rate: see https://doi.org/10.1111/2041-210X.13370
# =============================================================================
#EltonTraits 1.0: Species-level foraging attributes of the world's birds and mammals
Elton_mammal_trait <- read.delim("D:/BiodiversityEmbedding/Traits_data/EltonTraits 1.0/MamFuncDat.txt")
Elton_mammal_trait <- Elton_mammal_trait %>% 
  dplyr::select(-MSW3_ID,-MSWFamilyLatin,-Diet.Source,-Diet.Certainty,-ForStrat.Certainty,-ForStrat.Comment,-Activity.Source,-Activity.Certainty,-BodyMass.Source,-BodyMass.SpecLevel)

# Independent detections per year x array x species
mammal_IndDect <- mammal_seq %>%
  group_by(Year, Camera_Trap_Array, SCIENTIFIC_NAME) %>%
  summarise(IndDect = n_distinct(Sequence_ID), .groups = "drop")

# Survey nights per year x array
mammal_surveynight <- mammal_dep %>%
  group_by(Year, Camera_Trap_Array) %>%
  summarise(SurNights = sum(Survey_Nights), .groups = "drop")

mammal_RA <- mammal_IndDect %>%
  left_join(mammal_surveynight, by = c("Year", "Camera_Trap_Array")) %>%
  drop_na() %>%
  mutate(RA = IndDect / SurNights) %>%
  left_join(Elton_mammal_trait %>% dplyr::select(Scientific, bodymass = BodyMass.Value),
            by = join_by(SCIENTIFIC_NAME == Scientific))

# Species missing body mass (QC), then fill manually
missing_bodymass <- mammal_RA[is.na(mammal_RA$bodymass), c("SCIENTIFIC_NAME", "bodymass")] %>% distinct()

mammal_RA <- mammal_RA %>%
  mutate(bodymass = case_when(
    SCIENTIFIC_NAME == "Pekania pennanti"              ~ 3500,
    SCIENTIFIC_NAME == "Cervus canadensis"             ~ 329000,
    SCIENTIFIC_NAME == "Otospermophilus variegatus"    ~ 680.389,
    SCIENTIFIC_NAME == "Neotamias minimus"             ~ 45.359,
    SCIENTIFIC_NAME == "Canis rufus"                   ~ 27000,
    SCIENTIFIC_NAME == "Neotamias senex"               ~ 90.718,
    SCIENTIFIC_NAME == "Neotamias dorsalis"            ~ 70,
    SCIENTIFIC_NAME == "Urva javanica"                 ~ 771.107,
    SCIENTIFIC_NAME == "Neotamias townsendii"          ~ 73.708,
    SCIENTIFIC_NAME == "Neotamias umbrinus"            ~ 67,
    SCIENTIFIC_NAME == "Otospermophilus beecheyi"      ~ 589.67,
    SCIENTIFIC_NAME == "Neotamias amoenus"             ~ 51.029,
    SCIENTIFIC_NAME == "Neotamias speciosus"           ~ 60,
    SCIENTIFIC_NAME == "Neotamias panamintinus"        ~ 90,
    SCIENTIFIC_NAME == "Callospermophilus lateralis"   ~ 256,
    SCIENTIFIC_NAME == "Ictidomys tridecemlineatus"    ~ 190,
    SCIENTIFIC_NAME == "Neotamias merriami"            ~ 73.708,
    SCIENTIFIC_NAME == "Neotamias sonomae"             ~ 75,
    SCIENTIFIC_NAME == "Urocitellus richardsonii"      ~ 317.5,
    SCIENTIFIC_NAME == "Xerospermophilus tereticaudus" ~ 125,
    SCIENTIFIC_NAME == "Neotamias obscurus"            ~ 73.708,
    SCIENTIFIC_NAME == "Neotamias ruficaudus"          ~ 56,
    SCIENTIFIC_NAME == "Urocitellus armatus"           ~ 210,
    .default = bodymass
  ))
write.csv(mammal_RA, "data_clean/mammal_RA.csv", row.names = FALSE)
mammal_RA <- read.csv("data_clean/mammal_RA.csv")

# Body-mass-weighted RA (divide by mass-adjusted home-range area)
mammal_RA_weighted <- mammal_RA %>%
  mutate(RA_weighted = RA / (1.65 * bodymass^(1/3)))
write.csv(mammal_RA_weighted, "data_clean/mammal_RA_weighted.csv", row.names = FALSE)
mammal_RA_weighted <- read.csv("data_clean/mammal_RA_weighted.csv")

# =============================================================================
# 4. Richness via iNEXT (incidence-based)
# =============================================================================
# Build per-detection dates relative to deployment/array start (for SAC checks)
SAC_data <- mammal_seq %>%
  dplyr::select(Year, Camera_Trap_Array, Deployment_ID, Sequence_ID,
                Start_Time, End_Time, SCIENTIFIC_NAME) %>%
  mutate(date = as.Date(substr(Start_Time, 1, 10), format = "%Y/%m/%d")) %>%  # was mammal_seq$Start_Time
  left_join(
    mammal_dep %>%
      mutate(Start_Date = as.Date(Start_Date), End_Date = as.Date(End_Date)) %>%
      dplyr::select(Year, Camera_Trap_Array, Deployment_ID,
                    Start_Date_Dep = Start_Date, End_Date_Dep = End_Date),
    by = c("Year", "Camera_Trap_Array", "Deployment_ID")) %>%
  mutate(Days_Dep = as.numeric(date - Start_Date_Dep)) %>%
  group_by(Year, Camera_Trap_Array) %>%
  mutate(Start_Date_CTA  = min(Start_Date_Dep),
         Days_CTA        = as.numeric(date - Start_Date_CTA),
         NumberOfCameras = n_distinct(Deployment_ID)) %>%
  ungroup()

# --- Optional QC: species-accumulation curve for one array-year --------------
dat_sub <- SAC_data %>% filter(Year == 2021, Camera_Trap_Array == "Rentz")
dep_sac <- dat_sub %>%
  arrange(Deployment_ID, Days_Dep) %>%
  group_by(Deployment_ID) %>%
  mutate(cum_rich = sapply(seq_along(Days_Dep), function(i) length(unique(SCIENTIFIC_NAME[1:i])))) %>%
  ungroup()
total_sac <- dat_sub %>%
  arrange(Days_CTA) %>%
  mutate(cum_rich = sapply(seq_along(Days_CTA), function(i) length(unique(SCIENTIFIC_NAME[1:i]))))
ggplot() +
  geom_line(data = dep_sac,  aes(Days_Dep, cum_rich, group = Deployment_ID), color = "grey70", alpha = 0.6) +
  geom_line(data = total_sac, aes(Days_CTA, cum_rich), color = "black", linewidth = 1.2) +
  theme_classic() +
  labs(x = "Days of deployment", y = "Species accumulation",
       title = "Species Accumulation Curve (2021-Rentz)")

# --- Incidence-raw matrices per array-year -----------------------------------
inc_data <- SAC_data %>%
  distinct(Year, Camera_Trap_Array, Deployment_ID, SCIENTIFIC_NAME) %>%
  mutate(year_CTA = str_c(Year, Camera_Trap_Array, sep = "_")) %>%
  dplyr::select(year_CTA, Deployment_ID, SCIENTIFIC_NAME)
all_sp <- sort(unique(inc_data$SCIENTIFIC_NAME))

inc_list <- inc_data %>%
  split(.$year_CTA) %>%
  map(~ .x %>%
        mutate(pres = 1) %>%
        pivot_wider(names_from = SCIENTIFIC_NAME, values_from = pres, values_fill = 0) %>%
        dplyr::select(-year_CTA, -Deployment_ID) %>%
        as.matrix() %>%
        t())

# Pad every matrix to the full species set (rows = species, cols = cameras)
inc_list_full <- map(inc_list, function(m) {
  missing <- setdiff(all_sp, rownames(m))
  if (length(missing) > 0) {
    add <- matrix(0L, nrow = length(missing), ncol = ncol(m),
                  dimnames = list(missing, colnames(m)))
    m <- rbind(m, add)
  }
  m <- m[all_sp, , drop = FALSE]
  storage.mode(m) <- "integer"
  m
})

# Identify array-years iNEXT can run on without error
results <- sapply(names(inc_list_full), function(nm) {
  n_species <- sum(rowSums(inc_list_full[[nm]]) > 0)
  n_samples <- ncol(inc_list_full[[nm]])
  err <- tryCatch({
    iNEXT(inc_list_full[nm], q = c(0, 1, 2), datatype = "incidence_raw")
    "OK"
  }, error = function(e) "ERROR")
  c(n_species = n_species, n_samples = n_samples, status = err)
})
results <- t(results)
valid_sites <- rownames(results)[which(results[, 3] == "OK")]

# Assumes survey nights are roughly equal across cameras within a CTA
out <- iNEXT(inc_list_full[valid_sites], q = c(0, 1, 2), datatype = "incidence_raw")
write_rds(out, "data_clean/mammal_chaoDiversity.rds")
out <- readRDS("data_clean/mammal_chaoDiversity.rds")

# --- Optional QC plots -------------------------------------------------------
# ggiNEXT(out, type = 1) + theme(legend.position = "none")  # size-based R/E
# ggiNEXT(out, type = 2) + theme(legend.position = "none")  # sample completeness
# ggiNEXT(out, type = 3) + theme(legend.position = "none")  # coverage-based R/E

# =============================================================================
# 5. Assemble mammal_table (diversity + abundance + environment)
# =============================================================================
mammal_diversity <- out$AsyEst %>%
  rename(year_CTA = Assemblage) %>%
  mutate(Diversity = case_when(
    Diversity == "Species richness"  ~ "richness",
    Diversity == "Shannon diversity" ~ "shannon",
    Diversity == "Simpson diversity" ~ "simpson"
  )) %>%
  pivot_wider(
    id_cols    = year_CTA,
    names_from = Diversity,
    values_from = c(Observed, Estimator, s.e., LCL, UCL),
    names_glue = "{Diversity}_{.value}"
  ) %>%
  dplyr::select(year_CTA, richness_Observed, shannon_Observed, simpson_Observed,
                richness_Estimator, shannon_Estimator, simpson_Estimator)

# Community mean species abundance
mammal_meanabundance <- mammal_RA_weighted %>%
  mutate(year_CTA = str_c(Year, Camera_Trap_Array, sep = "_")) %>%
  dplyr::select(year_CTA, SCIENTIFIC_NAME, RA_weighted) %>%
  group_by(SCIENTIFIC_NAME) %>%
  mutate(max_abundance_species = max(RA_weighted)) %>%
  ungroup() %>%
  group_by(year_CTA) %>%
  summarise(abundance = sum(RA_weighted / max_abundance_species) / n_distinct(SCIENTIFIC_NAME),
            MSA       = sum(RA_weighted) / n_distinct(SCIENTIFIC_NAME),
            .groups = "drop")

CTA_env <- read.csv("data_clean/Camera_Trap_Array_env.csv")

mammal_table <- mammal_diversity %>%
  left_join(mammal_meanabundance, by = "year_CTA") %>%
  left_join(CTA_env,             by = "year_CTA") %>%
  left_join(CTA_habitat,         by = "Camera_Trap_Array") %>%   # was st_drop_geometry() on a non-sf frame
  mutate(NA_L1CODE = factor(as.integer(NA_L1CODE),
                            levels = sort(unique(as.integer(NA_L1CODE))))) %>%
  relocate(Year, Camera_Trap_Array, .after = year_CTA)
write.csv(mammal_table, "data_clean/mammal_table.csv", row.names = FALSE)
mammal_table <- read.csv("data_clean/mammal_table.csv")

# Dominant land-cover codes -> labels
mammal_table <- mammal_table %>%
  mutate(across(contains("dominant_lc"), ~ factor(
    .x,
    levels = c(1, 2, 3, 4, 5, 7, 9, 81, 82),
    labels = c("Open water", "Developed", "Barren land", "Forest",
               "Shrub/Scrub", "Grassland/Herbaceous", "Wetlands",
               "Pasture/Hay", "Cultivated crops")
  )))

# =============================================================================
# 6. Per-scale standardised tables (1, 5, 10, 20, 50 km)
# =============================================================================
keep_cols <- c("year_CTA", "Year", "Camera_Trap_Array",
               "richness_Estimator", "shannon_Estimator", "simpson_Estimator",
               "abundance", "MSA")

mammal_table_50km <- mammal_table %>%
  dplyr::select(all_of(keep_cols), contains("_50km"), NA_L1CODE, Habitat) %>%
  mutate(EH_50km_t     = scale(EH_50km)[, 1],
         elev_50km_t   = scale(sqrt(elev_50km))[, 1],
         PopDen_50km_t = scale(log1p(PopDen_50km))[, 1],
         NPP_50km_t    = scale(NPP_50km)[, 1],
         Prec_50km_t   = scale(sqrt(Prec_50km))[, 1],
         Temp_50km_t   = scale(Temp_50km)[, 1])

mammal_table_20km <- mammal_table %>%
  dplyr::select(all_of(keep_cols), contains("_20km"), NA_L1CODE, Habitat) %>%
  mutate(EH_20km_t     = scale(EH_20km)[, 1],
         elev_20km_t   = scale(sqrt(elev_20km))[, 1],
         PopDen_20km_t = scale(log1p(PopDen_20km))[, 1],
         NPP_20km_t    = scale(NPP_20km)[, 1],
         Prec_20km_t   = scale(sqrt(Prec_20km))[, 1],
         Temp_20km_t   = scale(Temp_20km)[, 1])

mammal_table_10km <- mammal_table %>%
  dplyr::select(all_of(keep_cols), contains("_10km"), NA_L1CODE, Habitat) %>%
  mutate(EH_10km_t     = scale(EH_10km)[, 1],
         elev_10km_t   = scale(sqrt(elev_10km))[, 1],
         PopDen_10km_t = scale(log1p(PopDen_10km))[, 1],
         NPP_10km_t    = scale(NPP_10km)[, 1],
         Prec_10km_t   = scale(sqrt(Prec_10km))[, 1],
         Temp_10km_t   = scale(Temp_10km)[, 1])

mammal_table_5km <- mammal_table %>%
  dplyr::select(all_of(keep_cols), contains("_5km"), NA_L1CODE, Habitat) %>%
  mutate(EH_5km_t     = scale(EH_5km)[, 1],
         elev_5km_t   = scale(sqrt(elev_5km))[, 1],
         PopDen_5km_t = scale(log1p(PopDen_5km))[, 1],
         NPP_5km_t    = scale(NPP_5km)[, 1],
         Prec_5km_t   = scale(sqrt(Prec_5km))[, 1],
         Temp_5km_t   = scale(Temp_5km)[, 1])

mammal_table_1km <- mammal_table %>%
  dplyr::select(all_of(keep_cols), contains("_1km"), NA_L1CODE, Habitat) %>%
  mutate(EH_1km_t     = scale(EH_1km)[, 1],
         elev_1km_t   = scale(sqrt(elev_1km))[, 1],
         PopDen_1km_t = scale(log1p(PopDen_1km))[, 1],
         NPP_1km_t    = scale(NPP_1km)[, 1],
         Prec_1km_t   = scale(sqrt(Prec_1km))[, 1],
         Temp_1km_t   = scale(Temp_1km)[, 1])

# ---- Write out the five per-scale tables ------------------------------------
mammal_tables_by_scale <- list(
  "1km"  = mammal_table_1km,
  "5km"  = mammal_table_5km,
  "10km" = mammal_table_10km,
  "20km" = mammal_table_20km,
  "50km" = mammal_table_50km
)

iwalk(mammal_tables_by_scale, function(tbl, scale) {
  write.csv(tbl, file.path("data_clean", paste0("mammal_table_", scale, ".csv")),
            row.names = FALSE)
})