# =============================================================================
# BBS bird traits: AVONET + BIRDBASE, trait PCoA, and phylogeny
# Re-evaluating the Environmental Heterogeneity-Diversity Relationship
#
# Sections:
#   1. BBS species list (AOS names) + scientific-name cleaning
#   2. Harmonise to eBird taxonomy -> BBS_species2
#   3. AVONET morphological / ecological traits -> BBS_AVONET
#   4. BIRDBASE life-history traits -> BBS_BIRDBASE
#   5. Combine -> BBS_traits
#   6. Trait PCoA -> BBS_traits_space
#   7. Bird phylogeny -> BBS_phylo
# =============================================================================

# ---- Libraries --------------------------------------------------------------
library(tidyverse)   # dplyr, tidyr, stringr, ggplot2, readr, tibble, purrr
library(readxl)      # read_excel(), excel_sheets()
library(grid)        # unit(), arrow()  (used in the PCoA biplot)
library(FD)          # gowdis()
library(clootl)      # clootl_data, extractTree()  (were loaded mid-script)


traits_dir <- "D:/BiodiversityEmbedding/Traits_data"   # raw trait sources (EDIT ME)

# =============================================================================
# 1. BBS species list + scientific-name cleaning
# =============================================================================
BBS_data     <- read.csv("data_clean/BBS_data.csv")
species_hash <- read.csv("data_clean/BBS_species_hash.csv")

BBS_species <- BBS_data %>%
  distinct(aou) %>%
  left_join(species_hash %>% distinct(aou, ScientificName), by = "aou")

# Names with hybrids (" x "), slashes ("/"), or "sp." need special handling
clean_species <- function(x) {
  # Leave "sp." untouched
  if (str_detect(x, "sp\\.")) return(x)
  
  has_x     <- str_detect(x, " x ")
  has_slash <- str_detect(x, "/")
  
  # Case A: exactly two "/" and no "x" -> pair words as 1+3 / 2+4
  if (!has_x && str_count(x, "/") == 2) {
    parts <- str_split(x, "/")[[1]] %>% str_trim()
    words <- unlist(str_split(parts, " "))
    words <- words[words != ""]
    if (length(words) == 4) {
      return(c(paste(words[1], words[3]), paste(words[2], words[4])))
    }
  }
  
  # Other cases (has "x" or multiple "/") -> original logic
  x_clean <- str_replace_all(x, " x ", "/")     # replace " x " with "/"
  parts   <- str_trim(str_split(x_clean, "/")[[1]])
  genus   <- word(parts[1], 1)
  
  parts_full <- sapply(parts, function(p) {
    if (str_count(p, " ") == 0) p <- paste(genus, p)   # species only -> prepend genus
    word(p, 1, 2)                                       # drop subspecies
  })
  parts_full
}

BBS_species <- BBS_species %>%
  mutate(species_list = lapply(ScientificName, clean_species)) %>%
  unnest(species_list) %>%
  rename(species = species_list) %>%
  mutate(species = case_when(
    species == "Dryocopus pileatus"          ~ "Hylatomus pileatus",
    species == "Butorides virescens"         ~ "Butorides striata",
    species == "Dryobates villosus"          ~ "Leuconotopicus villosus",
    species == "Nannopterum auritum"         ~ "Nannopterum auritus",
    species == "Dryobates borealis"          ~ "Leuconotopicus borealis",
    species == "Leucophaeus atricilla"       ~ "Larus atricilla",
    species == "Himantopus mexicanus"        ~ "Himantopus himantopus",
    species == "Aphelocoma woodhouseii"      ~ "Aphelocoma californica",
    species == "Icterus bullockii"           ~ "Icterus bullockiorum",
    species == "Picoides dorsalis"           ~ "Picoides tridactylus",
    species == "Anas diazi"                  ~ "Anas platyrhynchos",
    species == "Dryobates arizonae"          ~ "Leuconotopicus arizonae",
    species == "Accipiter atricapillus"      ~ "Accipiter gentilis",
    species == "Coccothraustes vespertinus"  ~ "Hesperiphona vespertina",
    species == "Corthylio calendula"         ~ "Regulus calendula",
    species == "Sturnella lilianae"          ~ "Sturnella magna",
    species == "Nannopterum brasilianum"     ~ "Nannopterum brasilianus",
    species == "Haematopus bachmani"         ~ "Haematopus ater",
    species == "Pica nuttalli"               ~ "Pica nutalli",
    species == "Dryobates albolarvatus"      ~ "Leuconotopicus albolarvatus",
    species == "Streptopelia chinensis"      ~ "Spilopelia chinensis",
    species == "Phalaropus tricolor"         ~ "Steganopus tricolor",
    species == "Leucophaeus pipixcan"        ~ "Larus pipixcan",
    species == "Carpodacus cassinii"         ~ "Haemorhous cassinii",
    species == "purpureus mexicanus"         ~ "Haemorhous purpureus",
    species == "Loxia sinesciuris"           ~ "Loxia curvirostra",
    species == "Centronyx henslowii"         ~ "Passerculus henslowii",
    species == "Canachites canadensis"       ~ "Falcipennis canadensis",
    species == "Centronyx bairdii"           ~ "Passerculus bairdii",
    species == "auratus cafer"               ~ "Colaptes auratus",
    .default = species
  ))

# =============================================================================
# 2. Harmonise to eBird taxonomy -> BBS_species2
# =============================================================================
eBird_taxa <- clootl_data$taxonomy.files$year2025

BBS_species2 <- BBS_species %>%
  mutate(ebird_species = coalesce(
    eBird_taxa$SCI_NAME[match(ScientificName, eBird_taxa$SCI_NAME)],
    eBird_taxa$SCI_NAME[match(species,        eBird_taxa$SCI_NAME)])) %>%
  mutate(ebird_species = case_when(
    species == "Bubulcus ibis"          ~ "Ardea ibis",
    species == "Accipiter cooperii"     ~ "Astur cooperii",
    species == "Nannopterum auritus"    ~ "Nannopterum auritum",
    species == "Ixobrychus exilis"      ~ "Botaurus exilis",
    species == "Accipiter gentilis"     ~ "Astur gentilis",
    species == "Nannopterum brasilianus"~ "Nannopterum brasilianum",
    species == "Charadrius montanus"    ~ "Anarhynchus montanus",
    species == "Porphyrio martinicus"   ~ "Porphyrio martinica",
    species == "Charadrius nivosus"     ~ "Anarhynchus nivosus",
    species == "Icterus bullockiorum"   ~ "Icterus bullockii",
    species == "Charadrius wilsonia"    ~ "Anarhynchus wilsonia",
    .default = ebird_species)) %>%
  left_join(eBird_taxa %>% select(SCI_NAME, ORDER), by = join_by(ebird_species == SCI_NAME))

write.csv(BBS_species2, "data_clean/BBS_species2.csv", row.names = FALSE)
length(na.omit(unique(BBS_species2$ebird_species)))   # 520 species

# =============================================================================
# 3. AVONET traits -> BBS_AVONET
#    BBS uses AOS taxonomy; AVONET uses BirdLife. Join on the cleaned `species`.
# =============================================================================
avonet_file <- file.path(traits_dir, "AVONET", "AVONET Supplementary dataset 1.xlsx")
AVONET <- read_excel(avonet_file, sheet = "AVONET1_BirdLife")

# Most frequent non-NA value (for categorical traits)
Mode1 <- function(x) {
  x <- na.omit(as.character(x))
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

BBS_AVONET <- BBS_species %>%
  left_join(AVONET, by = join_by(species == Species1)) %>%
  dplyr::select(-Sequence, -Family1, -Order1, -Avibase.ID1, -Total.individuals, -Female, -Male,
                -Unknown, -Complete.measures, -Mass.Source, -Mass.Refs.Other, -Inference,
                -Traits.inferred, -Reference.species) %>%
  mutate(Min.Latitude = as.numeric(Min.Latitude), Max.Latitude = as.numeric(Max.Latitude),
         Centroid.Latitude = as.numeric(Centroid.Latitude), Centroid.Longitude = as.numeric(Centroid.Longitude),
         Range.Size = as.numeric(Range.Size), Habitat.Density = as.character(Habitat.Density)) %>%
  group_by(aou, ScientificName) %>%
  summarise(
    # numeric traits (mean)
    Beak.Length_Culmen = mean(Beak.Length_Culmen, na.rm = TRUE),
    Beak.Length_Nares  = mean(Beak.Length_Nares,  na.rm = TRUE),
    Beak.Width         = mean(Beak.Width,         na.rm = TRUE),
    Beak.Depth         = mean(Beak.Depth,         na.rm = TRUE),
    Tarsus.Length      = mean(Tarsus.Length,      na.rm = TRUE),
    Wing.Length        = mean(Wing.Length,        na.rm = TRUE),
    Kipps.Distance     = mean(Kipps.Distance,     na.rm = TRUE),
    Secondary1         = mean(Secondary1,         na.rm = TRUE),
    `Hand-Wing.Index`  = mean(`Hand-Wing.Index`,  na.rm = TRUE),
    Tail.Length        = mean(Tail.Length,        na.rm = TRUE),
    Mass               = mean(Mass,               na.rm = TRUE),
    Min.Latitude       = mean(Min.Latitude,       na.rm = TRUE),
    Max.Latitude       = mean(Max.Latitude,       na.rm = TRUE),
    Centroid.Latitude  = mean(Centroid.Latitude,  na.rm = TRUE),
    Centroid.Longitude = mean(Centroid.Longitude, na.rm = TRUE),
    Range.Size         = mean(Range.Size,         na.rm = TRUE),
    # categorical traits (mode)
    Habitat           = Mode1(Habitat),
    Habitat.Density   = Mode1(Habitat.Density),
    Migration         = Mode1(Migration),
    Trophic.Level     = Mode1(Trophic.Level),
    Trophic.Niche     = Mode1(Trophic.Niche),
    Primary.Lifestyle = Mode1(Primary.Lifestyle),
    .groups = "drop"
  )
write.csv(BBS_AVONET, "data_clean/BBS_AVONET.csv", row.names = FALSE)

# =============================================================================
# 4. BIRDBASE life-history traits -> BBS_BIRDBASE
# =============================================================================
birdbase_file <- file.path(traits_dir, "BIRDBASE", "BIRDBASE v2025.1 Sekercioglu et al. Final.xlsx")
BIRDBASE <- read_excel(birdbase_file, sheet = "Data", skip = 1) %>%
  select(species = `HBW/BirdLife International (v9.1)`, Clutch_Max, Incu2, Fldg2)

BBS_BIRDBASE <- BBS_species2 %>%          # uses the in-memory table from Section 2
  left_join(BIRDBASE, by = "species")

# =============================================================================
# 5. Combine AVONET + BIRDBASE -> BBS_traits
# =============================================================================
# NOTE: BBS_species2 can have >1 row per aou (from the unnested hybrid/slash
#       names), so joining the BIRDBASE columns by aou may duplicate aou rows.
#       If that happens, summarise BBS_BIRDBASE to one row per aou first.
BBS_traits <- BBS_AVONET %>%
  left_join(BBS_BIRDBASE %>% select(aou, Clutch_Max, Incu2, Fldg2), by = "aou")
write.csv(BBS_traits, "data_clean/BBS_traits.csv", row.names = FALSE)
BBS_traits <- read.csv("data_clean/BBS_traits.csv")

# =============================================================================
# 6. Trait PCoA -> BBS_traits_space
# =============================================================================
trait_df <- BBS_traits %>%
  drop_na() %>%
  mutate(
    logmass  = log(Mass),
    logrange = log(Range.Size),
    Habitat           = factor(Habitat),
    Habitat.Density   = factor(Habitat.Density, levels = c(1, 2, 3), labels = c("Dense", "Semi-open", "Open")),
    Migration         = factor(Migration, levels = c(1, 2, 3), labels = c("Sedentary", "Partially migratory", "Migratory")),
    Trophic.Level     = factor(Trophic.Level, levels = c("Herbivore", "Omnivore", "Carnivore", "Scavenger")),
    Trophic.Niche     = factor(Trophic.Niche),
    Primary.Lifestyle = factor(Primary.Lifestyle)) %>%
  select(-Min.Latitude, -Max.Latitude, -Centroid.Latitude, -Centroid.Longitude, -Mass, -Range.Size)

# --- QC correlation panels (interactive) -------------------------------------
panel.cor <- function(x, y, digits = 2, cex.cor = 1.2) {
  usr <- par("usr"); on.exit(par(usr)); par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use = "complete.obs")
  text(0.5, 0.5, formatC(abs(r), format = "f", digits = digits), cex = cex.cor)
}
panel.hist <- function(x, ...) {
  usr <- par("usr"); par(usr = c(usr[1:2], 0, 1.5))
  h <- hist(x, plot = FALSE); breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts / max(h$counts)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
# pairs(trait_df[, -c(1:2, 13:18)], diag.panel = panel.hist, upper.panel = panel.cor)

# Keep traits with |r| < 0.8 and log-transform right-skewed ones
trait_df <- trait_df %>%
  dplyr::select(Beak.Length_Culmen, Beak.Width, Tarsus.Length, Hand.Wing.Index,
                Tail.Length, logmass, logrange, Clutch_Max, Incu2, Fldg2,
                Habitat, Habitat.Density, Migration, Trophic.Level, Trophic.Niche, Primary.Lifestyle) %>%
  mutate(
    logbeaklength = log(Beak.Length_Culmen),
    logbeakwidth  = log(Beak.Width),
    logtarsus     = log(Tarsus.Length),
    logtail       = log(Tail.Length),
    logclutch     = log(Clutch_Max),
    logincubation = log(Incu2),
    logfledging   = log(Fldg2)) %>%
  dplyr::select(logbeaklength, logbeakwidth, logtarsus, logtail, Hand.Wing.Index,
                logmass, logrange, logclutch, logincubation, logfledging,   # removed duplicate logrange
                Habitat, Habitat.Density, Migration, Trophic.Level, Trophic.Niche, Primary.Lifestyle)

# PCoA on Gower distances of scaled traits
trait_scaled <- trait_df %>% mutate(across(where(is.numeric), ~ as.numeric(scale(.))))
gower_dist   <- FD::gowdis(trait_scaled)
pcoa_res     <- cmdscale(gower_dist, k = 3, eig = TRUE)   # classical (metric) MDS

traits_space <- as.data.frame(pcoa_res$points)
colnames(traits_space) <- c("PCo1", "PCo2", "PCo3")
traits_space <- cbind(traits_space, trait_df[, sapply(trait_df, is.factor)])

# Export the 3-axis trait space (one row per aou)
traits_space_output <- BBS_traits %>% drop_na() %>% select(aou) %>% cbind(traits_space[, 1:3])
write.csv(traits_space_output, "data_clean/BBS_traits_space.csv", row.names = FALSE)

# --- QC biplot (interactive) -------------------------------------------------
var_exp <- pcoa_res$eig / sum(pcoa_res$eig)
ylab <- paste0("PCo2 (", round(var_exp[2] * 100, 2), "%)")
zlab <- paste0("PCo3 (", round(var_exp[3] * 100, 2), "%)")

traits_num    <- trait_df[, sapply(trait_df, is.numeric)]
trait_vectors <- as.data.frame(cor(traits_num, traits_space[, 2:3]))
trait_vectors$trait <- rownames(trait_vectors)
arrow_scale <- 0.5   # adjustable, ~1 to 5
trait_vectors <- trait_vectors %>% mutate(PCo2 = PCo2 * arrow_scale, PCo3 = PCo3 * arrow_scale)

ggplot() +
  geom_point(data = traits_space, aes(PCo2, PCo3, color = Migration), size = 3, alpha = 0.7) +
  stat_ellipse(data = traits_space, aes(PCo2, PCo3, color = Migration)) +
  geom_segment(data = trait_vectors, aes(0, 0, xend = PCo2, yend = PCo3),
               arrow = arrow(length = unit(0.25, "cm")), linewidth = 0.8, color = "black") +
  geom_text(data = trait_vectors, aes(PCo2, PCo3, label = trait), size = 4, hjust = 0.5, vjust = -0.5) +
  geom_hline(yintercept = 0, col = "grey", linetype = "dashed") +
  geom_vline(xintercept = 0, col = "grey", linetype = "dashed") +
  labs(x = ylab, y = zlab) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(), axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black")) +
  coord_equal()

# =============================================================================
# 7. Bird phylogeny -> BBS_phylo
#    Birds of the World phylogeny: https://www.pnas.org/doi/10.1073/pnas.2409658122
# =============================================================================
BBS_phylo <- clootl::extractTree(species = na.omit(unique(BBS_species2$ebird_species)))
saveRDS(BBS_phylo, "data_clean/BBS_phylo_520species.rds")