# =============================================================================
# FIA tree traits: species list,  TRY + BIEN traits, and phylogenetic imputation
# Re-evaluating the Environmental Heterogeneity-Diversity Relationship
#
# Sections:
#   0. Unit-conversion lookup (defined first - it is used in the TRY section)
#   1. FIA species list (fia_sp)
#   2. Traits from TRY
#   3. Traits from BIEN
#   4. Combine TRY + BIEN
#   5. Phylogenetic imputation
# =============================================================================

# ---- Libraries --------------------------------------------------------------
library(tidyverse)      # dplyr, tidyr, stringr, tibble, purrr, ggplot2, readr
library(sf)
library(terra)
library(maps)
library(rtry)           # rtry_import()
library(BIEN)           # BIEN ranges + traits
library(TNRS)           # taxonomic name resolution (was library()'d mid-script)
library(V.PhyloMaker2)  # phylo.maker()
library(Rphylopars)     # phylopars()

# =============================================================================
# 0. Unit-conversion lookup for TRY traits
# =============================================================================
conversion_rules_corrected <- tribble(
  ~TraitName, ~OrigUnitStr, ~commonUnit, ~scaleVal,
  
  # ---- SLA (specific leaf area) - petiole excluded; target m2/kg ----
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded", "cm2 / g",   "m2/kg", 0.1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded", "cm2 g-1",   "m2/kg", 0.1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded", "cm2/g",     "m2/kg", 0.1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded", "cm2\u00b7g-1", "m2/kg", 0.1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded", "g m-2",     "m2/kg", -1000,  # LMA: invert
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded", "g.m2",      "m2/kg", -1000,  # LMA: invert
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded", "g/m2",      "m2/kg", -1000,  # LMA: invert
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded", "m2 kg-1",   "m2/kg", 1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded", "m2/g",      "m2/kg", 1000,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded", "m2/kg",     "m2/kg", 1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded", "m^2/kg",    "m2/kg", 1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded", "mg/cm2",    "m2/kg", -100,   # LMA: invert
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded", "mm2 / mg",  "m2/kg", 1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded", "mm2 mg-1",  "m2/kg", 1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded", "mm2*mg-1",  "m2/kg", 1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded", "mm2/g",     "m2/kg", 0.001,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded", "mm2/mg",    "m2/kg", 1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded", "mm2\u00b7mg-1", "m2/kg", 1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded", "mm^2/mg",   "m2/kg", 1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded", "m\u00b2 kg-1", "m2/kg", 1,
  
  # ---- SLA - petiole undefined; target m2/kg ----
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "(cm2/g)/100", "m2/kg", 0.001,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "(g m-2)",     "m2/kg", -1,      # LMA
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "(g/cm2)",     "m2/kg", -10000,  # LMA
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "SLA cm2/g",   "m2/kg", 0.1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "cm2 g-1",     "m2/kg", 0.1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "cm2 mg-1",    "m2/kg", 100,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "cm2.g-1",     "m2/kg", 0.1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "cm2/ g",      "m2/kg", 0.1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "cm2/g",       "m2/kg", 0.1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "cm2\u00b7g-1", "m2/kg", 0.1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "cm^2/g",      "m2/kg", 0.1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "g m-2",       "m2/kg", -1,      # LMA
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "g/cm2",       "m2/kg", -10000,  # LMA
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "g/m2",        "m2/kg", -1,      # LMA
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "gm-2",        "m2/kg", -1,      # LMA
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "m2 / kg",     "m2/kg", 1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "m2 kg-1",     "m2/kg", 1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "m2/g",        "m2/kg", 1000,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "m2/kg",       "m2/kg", 1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "mg/cm2",      "m2/kg", -0.1,    # LMA
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "mm2 / mg",    "m2/kg", 1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "mm2 mg-1",    "m2/kg", 1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "mm2/g",       "m2/kg", 0.001,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "mm2/mg",      "m2/kg", 1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "mm2mg-1",     "m2/kg", 1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "mm^2/mg",     "m2/kg", 1,
  "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded", "mm\u00b2/mg", "m2/kg", 1,
  
  # ---- Leaf nitrogen; target mg/g ----
  "Leaf nitrogen (N) content per leaf dry mass", "%",                 "mg/g", 10,
  "Leaf nitrogen (N) content per leaf dry mass", "% (100 * g g-1 )",  "mg/g", 10,
  "Leaf nitrogen (N) content per leaf dry mass", "% dry mass",        "mg/g", 10,
  "Leaf nitrogen (N) content per leaf dry mass", "% mass/mass",       "mg/g", 10,
  "Leaf nitrogen (N) content per leaf dry mass", "%DM",               "mg/g", 10,
  "Leaf nitrogen (N) content per leaf dry mass", "(mg g-1)",          "mg/g", 1,
  "Leaf nitrogen (N) content per leaf dry mass", "g N g-1 DW",        "mg/g", 1000,
  "Leaf nitrogen (N) content per leaf dry mass", "g/100g",            "mg/g", 10,
  "Leaf nitrogen (N) content per leaf dry mass", "g/g",               "mg/g", 1000,
  "Leaf nitrogen (N) content per leaf dry mass", "g/kg",              "mg/g", 1,
  "Leaf nitrogen (N) content per leaf dry mass", "g\u00b7kg-1",       "mg/g", 1,
  "Leaf nitrogen (N) content per leaf dry mass", "kg/kg",             "mg/g", 1000000,
  "Leaf nitrogen (N) content per leaf dry mass", "mg / g",            "mg/g", 1,
  "Leaf nitrogen (N) content per leaf dry mass", "mg N g-1",          "mg/g", 1,
  "Leaf nitrogen (N) content per leaf dry mass", "mg g-1",            "mg/g", 1,
  "Leaf nitrogen (N) content per leaf dry mass", "mg/g",              "mg/g", 1,
  "Leaf nitrogen (N) content per leaf dry mass", "mg/g dry mass",     "mg/g", 1,
  "Leaf nitrogen (N) content per leaf dry mass", "mg/mg *100",        "mg/g", 100,
  "Leaf nitrogen (N) content per leaf dry mass", "mg_g-1",            "mg/g", 1,
  "Leaf nitrogen (N) content per leaf dry mass", "percent",           "mg/g", 10,
  
  # ---- Plant lifespan; target year ----
  "Plant lifespan (longevity)", "NN",     "year", 1,
  "Plant lifespan (longevity)", "year",   "year", 1,
  "Plant lifespan (longevity)", "years",  "year", 1,
  "Plant lifespan (longevity)", "years?", "year", 1,
  "Plant lifespan (longevity)", "yr",     "year", 1,
  "Plant lifespan (longevity)", "yrs",    "year", 1,
  
  # ---- Seed dry mass; target mg ----
  "Seed dry mass", "(mg)",                   "mg", 1,
  "Seed dry mass", "1/g",                    "mg", -1000,
  "Seed dry mass", "1/kg",                   "mg", -1000000,
  "Seed dry mass", "1/lb",                   "mg", -454000,
  "Seed dry mass", "1/lb (lb~454g)",         "mg", -454000,
  "Seed dry mass", "1/liter (liter~1,7kg)",  "mg", -1700000,
  "Seed dry mass", "1/pound",                "mg", -454000,
  "Seed dry mass", "g",                      "mg", 1000,
  "Seed dry mass", "g / 1000 seeds",         "mg", 1,
  "Seed dry mass", "g/1000",                 "mg", 1,
  "Seed dry mass", "g/1000 seeds",           "mg", 1,
  "Seed dry mass", "gr",                     "mg", 1000,
  "Seed dry mass", "mg",                     "mg", 1
)

# =============================================================================
# 1. FIA species list
# =============================================================================
tpa_all5 <- readRDS("data_clean/tpa_all5.rds")

# Which species do we have?
fia_spcd_to_name <- tpa_all5 %>%
  st_drop_geometry() %>%
  filter(!is.na(SPCD), !is.na(SCIENTIFIC_NAME)) %>%
  distinct(SPCD, SCIENTIFIC_NAME)

n_distinct(fia_spcd_to_name$SCIENTIFIC_NAME)  # 902

# Drop genus-only ("spp.") and unidentified ("Tree") entries.
fia_sp <- fia_spcd_to_name %>%
  filter(!str_detect(SCIENTIFIC_NAME, "spp\\.|Tree"))
write.csv(fia_sp, "data_clean/fia_sp.csv", row.names = FALSE)

n_distinct(fia_sp$SCIENTIFIC_NAME)  # 868 = 902 - 34
n_distinct(fia_sp$SPCD)             # 880

# Occurrence count per species (QC)
fia_sp_occurrence <- table(tpa_all5$SCIENTIFIC_NAME[tpa_all5$SCIENTIFIC_NAME %in% fia_sp$SCIENTIFIC_NAME])
quantile(fia_sp_occurrence, c(0, 0.25, 0.5, 0.75, 1))



# =============================================================================
# 2. Traits from TRY
#    TraitIDs: 59 lifespan, 14 leaf N, 3115/3117 SLA, 26 seed mass
#    (38 Plant woodiness is dropped - categorical)
# =============================================================================
TRYSpecies <- read.delim("D:/BiodiversityEmbedding/Traits_data/TryAccSpecies.txt") %>%
  dplyr::select(AccSpeciesID, AccSpeciesName)

# Species already present in TRY by exact name
id  <- which(is.na(match(fia_sp$SCIENTIFIC_NAME, TRYSpecies$AccSpeciesName)))  # 107 unmatched
sp1 <- left_join(fia_sp[-id, ], TRYSpecies, by = join_by(SCIENTIFIC_NAME == AccSpeciesName)) %>%
  mutate(AccSpeciesName = SCIENTIFIC_NAME)

# Resolve the remaining names with TNRS, then match to TRY
fia_sp_tomatchTRY <- fia_sp[id, ]
results    <- TNRS(taxonomic_names = fia_sp_tomatchTRY) %>%
  dplyr::select(ID, Name_submitted, Overall_score, Accepted_species)
resultsTRY <- TNRS(taxonomic_names = TRYSpecies) %>%
  dplyr::select(ID, Name_submitted, Overall_score, Accepted_species)
write.csv(resultsTRY, "data/TNRS_TRYSpecies.csv", row.names = FALSE)
resultsTRY <- read.csv("data/TNRS_TRYSpecies.csv")

id  <- which(is.na(match(results$Accepted_species, resultsTRY$Accepted_species)))
sp2 <- results[-id, ] %>%
  dplyr::select(SPCD = ID, SCIENTIFIC_NAME = Name_submitted) %>%
  cbind(AccSpeciesID   = resultsTRY$ID[na.omit(match(results$Accepted_species, resultsTRY$Accepted_species))],
        AccSpeciesName = resultsTRY$Name_submitted[na.omit(match(results$Accepted_species, resultsTRY$Accepted_species))])

fia_sp_matchedTRY <- rbind(sp1, sp2) %>%
  filter(AccSpeciesName != "Aa sp")
n_distinct(fia_sp_matchedTRY$SCIENTIFIC_NAME)  # 835
n_distinct(fia_sp_matchedTRY$AccSpeciesID)     # 821
cat(paste(unique(fia_sp_matchedTRY$AccSpeciesID), collapse = ",\n"))  # IDs to request from TRY

# Import the TRY release and standardise units
FIA_traits <- rtry_import("D:/BiodiversityEmbedding/Traits_data/47273_11022026043624/47273.txt")
FIA_traits <- FIA_traits %>%
  dplyr::select(AccSpeciesID, AccSpeciesName, TraitID, TraitName, DataID, DataName,
                OriglName, OrigValueStr, OrigUnitStr) %>%
  filter(!is.na(TraitID)) %>%
  filter(TraitName != "Plant woodiness") %>%
  mutate(OrigValueStr = str_trim(OrigValueStr)) %>%
  filter(
    !is.na(OrigValueStr),                       # drop NA
    OrigValueStr != "",                         # drop empty
    OrigValueStr != "0",
    !str_detect(OrigValueStr, "[A-Za-z]|>|-")   # keep numeric-only values
  ) %>%
  mutate(OrigValueNum = as.numeric(OrigValueStr))

# Inspect the units present per trait (QC)
trait_units <- FIA_traits %>%
  group_by(TraitName, OrigUnitStr) %>%
  summarise(num = n(), .groups = "drop")
# View(trait_units)   # FIX: original called view() (no such function); use View()

# Apply the conversion table -> StdValue
trait_units_converted <- trait_units %>%
  left_join(conversion_rules_corrected, by = c("TraitName", "OrigUnitStr"))

FIA_traits <- FIA_traits %>%
  left_join(trait_units_converted %>% dplyr::select(-num), by = c("TraitName", "OrigUnitStr")) %>%
  mutate(StdValue = case_when(
    scaleVal < 0                    ~ abs(scaleVal) / OrigValueNum,   # LMA -> SLA
    !is.na(scaleVal) & scaleVal > 0 ~ OrigValueNum * scaleVal,
    TRUE                            ~ NA_real_
  )) %>%
  mutate(TraitName = if_else(str_detect(TraitName, "Leaf area per leaf dry mass"), "SLA", TraitName))

# One mean value per species x trait, wide format
FIA_traits <- FIA_traits %>%
  group_by(AccSpeciesID, TraitName) %>%
  summarise(Trait_mean_value = mean(StdValue), .groups = "drop") %>%
  pivot_wider(names_from = TraitName, values_from = Trait_mean_value, values_fill = NA) %>%
  filter(AccSpeciesID != 1)

# Cap implausible values, then keep the four target traits
FIA_traits <- FIA_traits %>%
  mutate(
    `Leaf nitrogen (N) content per leaf dry mass` =
      if_else(`Leaf nitrogen (N) content per leaf dry mass` > 65, NA_real_,
              `Leaf nitrogen (N) content per leaf dry mass`),
    SLA = if_else(SLA > 86, NA_real_, SLA)
  ) %>%
  rename(leafnitro = `Leaf nitrogen (N) content per leaf dry mass`) %>%
  left_join(fia_sp_matchedTRY %>% dplyr::select(-AccSpeciesName) %>%
              mutate(AccSpeciesID = as.integer(AccSpeciesID)),
            by = "AccSpeciesID") %>%
  dplyr::select(SPCD, SCIENTIFIC_NAME, leafnitro, SLA, `Seed dry mass`, `Plant lifespan (longevity)`)

write.csv(FIA_traits, "data_clean/FIA_traits_fromTRY.csv", row.names = FALSE)
n_distinct(FIA_traits$SCIENTIFIC_NAME)  # 648
cor(FIA_traits[, 3:6], use = "pairwise.complete.obs")

# =============================================================================
# 3. Traits from BIEN (to supplement TRY)
# =============================================================================
BIEN_trait_list()   # lists available trait types

bien_sla <- BIEN_trait_trait(trait = "leaf area per leaf dry mass") %>%   # m2/kg
  dplyr::select(scrubbed_species_binomial, trait_value) %>%
  mutate(trait_value = as.numeric(trait_value)) %>%
  na.omit() %>%
  group_by(scrubbed_species_binomial) %>%
  summarise(trait_value = mean(trait_value), .groups = "drop")

bien_leafnitro <- BIEN_trait_trait(trait = "leaf nitrogen content per leaf dry mass") %>%
  dplyr::select(scrubbed_species_binomial, trait_value) %>%
  mutate(trait_value = as.numeric(trait_value)) %>%
  na.omit() %>%
  group_by(scrubbed_species_binomial) %>%
  summarise(trait_value = mean(trait_value), .groups = "drop")

bien_seedmass <- BIEN_trait_trait(trait = "seed mass") %>%
  dplyr::select(scrubbed_species_binomial, trait_value) %>%
  mutate(trait_value = as.numeric(trait_value)) %>%
  na.omit() %>%
  group_by(scrubbed_species_binomial) %>%
  summarise(trait_value = mean(trait_value), .groups = "drop")

bien_traits <- full_join(bien_sla       %>% rename(SLA = trait_value),
                         bien_leafnitro %>% rename(leafnitro = trait_value),
                         by = "scrubbed_species_binomial") %>%
  full_join(bien_seedmass %>% rename(`Seed dry mass` = trait_value),
            by = "scrubbed_species_binomial") %>%
  mutate(bien_ID = row_number())

# Match FIA species to BIEN (exact, then via TNRS)
id  <- which(is.na(match(fia_sp$SCIENTIFIC_NAME, bien_traits$scrubbed_species_binomial)))  # 361
sp1 <- left_join(fia_sp[-id, ], bien_traits %>% dplyr::select(bien_ID, scrubbed_species_binomial),
                 by = join_by(SCIENTIFIC_NAME == scrubbed_species_binomial)) %>%
  mutate(scrubbed_species_binomial = SCIENTIFIC_NAME)

fia_sp_tomatchBIEN <- fia_sp[id, ]
results <- TNRS(taxonomic_names = fia_sp_tomatchBIEN) %>%
  dplyr::select(ID, Name_submitted, Overall_score, Accepted_species) %>%
  filter(Accepted_species != "")
resultsBIEN <- TNRS(taxonomic_names = bien_traits[, c("bien_ID", "scrubbed_species_binomial")]) %>%
  dplyr::select(ID, Name_submitted, Overall_score, Accepted_species)
write.csv(resultsBIEN, "data/TNRS_BIENSpecies.csv", row.names = FALSE)
resultsBIEN <- read.csv("data/TNRS_BIENSpecies.csv")

id  <- which(is.na(match(results$Accepted_species, resultsBIEN$Accepted_species)))
sp2 <- results[-id, ] %>%
  dplyr::select(SPCD = ID, SCIENTIFIC_NAME = Name_submitted) %>%
  cbind(bien_ID                   = resultsBIEN$ID[na.omit(match(results$Accepted_species, resultsBIEN$Accepted_species))],
        scrubbed_species_binomial = resultsBIEN$Name_submitted[na.omit(match(results$Accepted_species, resultsBIEN$Accepted_species))])

fia_sp_matchedBIEN <- rbind(sp1, sp2)
BIEN_traits <- fia_sp_matchedBIEN %>%
  left_join(bien_traits %>% mutate(bien_ID = as.character(bien_ID)),
            by = c("scrubbed_species_binomial", "bien_ID")) %>%
  dplyr::select(-bien_ID, -scrubbed_species_binomial)
write.csv(BIEN_traits, "data_clean/FIA_traits_fromBIEN.csv", row.names = FALSE)

# =============================================================================
# 4. Combine TRY + BIEN
# =============================================================================
combined_traits <- bind_rows(
  BIEN_traits,
  FIA_traits[, -6]   # drop the extra lifespan column so the sets align
)
final_traits <- combined_traits %>%
  group_by(SCIENTIFIC_NAME) %>%
  summarise(
    SPCD            = first(SPCD),   # keep one code
    leafnitro       = mean(leafnitro, na.rm = TRUE),
    SLA             = mean(SLA, na.rm = TRUE),
    `Seed dry mass` = mean(`Seed dry mass`, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(final_traits, "data_clean/FIA_final_traits.csv", row.names = FALSE)

# =============================================================================
# 5. Phylogenetic trait imputation
# =============================================================================
# Build a phylogeny for the trait species
sp.list <- TNRS(taxonomic_names = final_traits %>% dplyr::select(SPCD, SCIENTIFIC_NAME))
sp.list.sp_ge_fa <- sp.list %>%
  dplyr::select(species = Name_submitted, genus = Genus_matched, family = Accepted_family)
phylo <- phylo.maker(sp.list = sp.list.sp_ge_fa)

# Impute with phylopars
imputed_trait <- Rphylopars::phylopars(
  trait_data = final_traits %>%
    dplyr::select(-SPCD, species = SCIENTIFIC_NAME) %>%
    mutate(species = str_replace(species, " ", "_")),
  tree = phylo$scenario.3,
  pheno_error = FALSE, phylo_correlated = TRUE, pheno_correlated = FALSE)

# NOTE: `1:663` hard-codes the number of tips (663 tips + 598 internal nodes).
#       If the species count changes, update this index.
final_imputed_traits <- as.data.frame(imputed_trait$anc_recon[1:663, ]) %>%
  rownames_to_column(var = "SCIENTIFIC_NAME") %>%
  mutate(SCIENTIFIC_NAME = str_replace(SCIENTIFIC_NAME, "_", " ")) %>%
  left_join(final_traits[, 1:2], by = "SCIENTIFIC_NAME") %>%
  relocate(SPCD, .after = SCIENTIFIC_NAME)

final_imputed_traits <- final_imputed_traits[
  complete.cases(final_imputed_traits) &
    final_imputed_traits$SLA > 0 &
    final_imputed_traits$Seed.dry.mass > 0 &
    final_imputed_traits$leafnitro > 0, ]
write.csv(final_imputed_traits, "data_clean/FIA_final_imputed_traits.csv", row.names = FALSE)
final_imputed_traits <- read.csv("data_clean/FIA_final_imputed_traits.csv")
corrplot::corrplot(cor(final_imputed_traits[, 3:5]), type = "upper")

# ---- Alternative: impute from the 864-species raw trait table ---------------
# NOTE: FIA_864sp_raw_trait.csv is an external input, not produced above.
sp_trait <- read.csv("data_clean/FIA_864sp_raw_trait.csv")
apply(sp_trait, 2, function(x) length(which(is.na(x))) / nrow(sp_trait))  # missingness

sp.list <- TNRS(taxonomic_names = sp_trait %>% dplyr::select(rownumber, SCIENTIFIC_NAME))
sp.list.sp_ge_fa <- sp.list %>%
  dplyr::select(species = Name_submitted, genus = Genus_matched, family = Accepted_family)
phylo <- V.PhyloMaker2::phylo.maker(sp.list = sp.list.sp_ge_fa)

imputed_fit <- Rphylopars::phylopars(   # renamed from `a` for clarity
  trait_data = sp_trait %>%
    dplyr::select(species = SCIENTIFIC_NAME, swd, SLA, maxHt, leafnitro, Seed.dry.mass) %>%
    mutate(species = str_replace_all(species, " ", "_")),
  tree = phylo$scenario.3,
  pheno_error = FALSE, phylo_correlated = TRUE, pheno_correlated = FALSE)

# NOTE: `1:864` hard-codes the tip count - keep in sync with the species total.
imputed_sp_trait <- as.data.frame(imputed_fit$anc_recon[1:864, ]) %>%
  rownames_to_column(var = "SCIENTIFIC_NAME") %>%
  mutate(SCIENTIFIC_NAME = str_replace_all(SCIENTIFIC_NAME, "_", " ")) %>%
  left_join(sp_trait[, 2:4], by = "SCIENTIFIC_NAME") %>%
  relocate(species_tnrs, .after = SCIENTIFIC_NAME)
write.csv(imputed_sp_trait, "data_clean/FIA_864sp_imputed_traits.csv", row.names = FALSE)
imputed_sp_trait <- read.csv("data_clean/FIA_864sp_imputed_traits.csv")

# QC: raw vs imputed correlations
par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
cols <- c("swd", "SLA", "maxHt", "leafnitro", "Seed.dry.mass")
corrplot::corrplot(cor(sp_trait[, cols], use = "pairwise.complete.obs")); title("Incomplete traits", line = -1)
corrplot::corrplot(cor(imputed_sp_trait[, cols]));                        title("Imputed traits",   line = -1)