# =============================================================================
# FigS06_EH_characterization.R
# Fig. S6: EH map (5 km CONUS grid), correlation of EH with land-cover
# diversity / canopy structure / elevation / climate, LMG of EH drivers, and
# EH across dominant land-cover groups (ANOVA + emmeans letters)
#
# Source   : EmbeddingBio/analysis/EH_representation.R (paths only changed)
# Inputs   : data/raw/gee_exports/EH_US_5km_4326.shp (EE export, see README_GEE.md)
#            HDR_EXT_DIR/environmental_data/{Elevation_US_5070.tif, tcc_US_5070.tif,
#              NLCD_US_5070.tif, EH_elev_mean.csv, EH_tcc_mean.csv,
#              EH_US_5km_with_lc_frac.csv, EH_MAT_mean.csv, EH_MAP_mean.csv,
#              EH_US_5km_CHM_stats.csv, EH_soilph_std.csv}
# TODO(reviewer M5): the EH_*_mean/std CSVs are zonal statistics of the 5 km
#   grid made outside R (EE/ArcGIS?). Add the code or document the tool.
# TODO: panels are printed, not saved; add ggsave(p_fig("FigS6_*.png")).
# =============================================================================

source(here::here("code", "00_config.R"))

#to see how EH relates with tree cover, land cover, elevation, 
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(fields)
library(dplyr)
library(tidyr)

us_boundary <- st_read(p_ext("tl_2024_us_state/tl_2024_us_state.shp"))
us_boundary <- st_transform(us_boundary, "epsg:5070")

#Environmental heterogeneity
EH_US_5km <- st_read(p_gee("EH_US_5km_4326.shp"))

EH_vect <- vect(EH_US_5km)

#visualize EH
summary(EH_US_5km$mean_cosin)
cols <- rev(hcl.colors(100, "YlOrRd"))
vals <- pmax(EH_US_5km$mean_cosin, 0)
plot(EH_vect,
     col = cols[cut(vals, breaks = 100, include.lowest = TRUE)],
     border = NA)
states <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))
plot(st_geometry(states), add = TRUE, col = NA, border = "black", lwd = 0.5)


image.plot(
  legend.only = TRUE,
  zlim = range(vals, na.rm = TRUE),
  col = cols,
  legend.lab = "EH",
  smallplot = c(0.86, 0.89, 0.18, 0.78)
)


#Elevation (30m)
Elev_proj <- rast(p_env("Elevation_US_5070.tif"))
Elev_proj
plot(Elev_proj)


#USDA forest services Tree Canopy Cover (TCC) Datasets：https://data.fs.usda.gov/geodata/rastergateway/treecanopycover/ (30m)
tcc_proj <- rast(p_env("tcc_US_5070.tif"))
tcc_proj
plot(tcc_proj)


#NLCD land cover (30m) https://www.usgs.gov/node/279743
NLCD_proj <- rast(p_env("NLCD_US_5070.tif"))
NLCD_proj
plot(NLCD_proj)


#Temp, precp
MAT_proj <- rast(p_clean("MAT_usa_5070.tif"))
MAP_proj <- rast(p_clean("MAP_usa_5070.tif"))



#statistics
# elev_mean <- terra::extract(Elev_proj, EH_vect, fun = "mean", na.rm = TRUE)
# tcc_mean <- terra::extract(tcc_proj, EH_vect, fun = "mean", na.rm = TRUE)
# lc <- terra::extract(NLCD_proj, EH_vect)

elev_mean <- read.csv(p_env("EH_elev_mean.csv"))
tcc_mean <- read.csv(p_env("EH_tcc_mean.csv"))
lc_frac <- read.csv(p_env("EH_US_5km_with_lc_frac.csv")) %>% drop_na()
MAT_mean <- read.csv(p_env("EH_MAT_mean.csv"))
MAP_mean <- read.csv(p_env("EH_MAP_mean.csv"))
# TPsea_mean <- read.csv(p_env("EH_climate_seasonality_mean.csv"))
CHM_stat <- read.csv(p_env("EH_US_5km_CHM_stats.csv")) %>% dplyr::select(ID,EH,chm_mean,chm_std) %>% arrange(ID)
soilph_sd <- read.csv(p_env("EH_soilph_std.csv"))

#compute land cover shannon diversity and dominant lc type
lc_cols <- grep("^nlcd_", names(lc_frac), value = TRUE)
lc_frac$shannon <- apply(lc_frac[, lc_cols], 1, function(x) {
  x <- x[x > 0]
  -sum(x * log(x))
})
lc_frac$dom_lc <- lc_cols[max.col(lc_frac[, lc_cols], ties.method = "first")]
lc_frac$dom_lc <- as.numeric(sub("nlcd_", "", lc_frac$dom_lc))

#NLCD code
lc_frac$dom_lc_f <- factor(
  lc_frac$dom_lc,
  levels = c(11, 21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95),
  labels = c(
    "Open Water",
    "Developed, Open Space",
    "Developed, Low Intensity",
    "Developed, Medium Intensity",
    "Developed, High Intensity",
    "Barren Land",
    "Deciduous Forest",
    "Evergreen Forest",
    "Mixed Forest",
    "Shrub/Scrub",
    "Herbaceous",
    "Pasture/Hay",
    "Cultivated Crops",
    "Woody Wetlands",
    "Emergent Herbaceous Wetlands"
  )
)
#reclassification
lc_frac <- lc_frac %>%
  mutate(
    dom_lc_group = case_when(
      
      # Developed
      dom_lc_f %in% c(
        "Developed, Open Space",
        "Developed, Low Intensity",
        "Developed, Medium Intensity",
        "Developed, High Intensity"
      ) ~ "Developed",
      
      # Forest
      dom_lc_f %in% c(
        "Deciduous Forest",
        "Evergreen Forest",
        "Mixed Forest"
      ) ~ "Forest",
      
      # Wetlands
      dom_lc_f %in% c(
        "Woody Wetlands",
        "Emergent Herbaceous Wetlands"
      ) ~ "Wetlands",
      
     
      TRUE ~ as.character(dom_lc_f)
    )
  )
lc_frac$dom_lc_group <- factor(lc_frac$dom_lc_group)

str(lc_frac)

#join vars together
dat <- lc_frac %>%  
  left_join(elev_mean,by="ID") %>% 
  left_join(tcc_mean, by="ID") %>% 
  left_join(MAT_mean %>% dplyr::select(-EH),by="ID") %>% 
  left_join(MAP_mean %>% dplyr::select(-EH),by="ID") %>% 
  left_join(CHM_stat %>% dplyr::select(-EH),by="ID") %>% 
  left_join(soilph_sd %>% dplyr::select(-EH),by="ID") %>% 
  dplyr::select(-ID) %>% 
  filter(EH>=0) %>% 
  dplyr::select(-contains('nlcd'))
str(dat)

write.csv(dat,p_clean("EH_representation.csv"),row.names = FALSE)
dat <- read.csv(p_clean("EH_representation.csv"))
str(dat)

#analysis: How EH correlate with elev, tcc, and land cover shannon index; 
cor(dat[,c(1,2,6:12)],use = "pairwise.complete.obs")
corrplot::corrplot(cor(dat[,c(1,2,6:12)],use = "pairwise.complete.obs"),type="lower")

form <- as.formula(EH ~ elev_mean + tcc_mean + shannon + MAT_mean + MAP_mean + chm_mean+ chm_std + soilph_std)
model <- lm(form,data = dat)
summary(model)

#variance explanation
library(relaimpo)
a <- calc.relimp(model, type = "lmg", rela = TRUE)
imp <- data.frame(
  Variable = names(a$lmg),
  Importance = a$lmg
) %>%
  arrange(Importance) %>%
  mutate(Variable = factor(Variable, levels = Variable))
ggplot(imp, aes(x = Variable, y = Importance)) +
  geom_col(fill = "#4C9F70", width = 0.7) +
  geom_text(
    aes(label = round(Importance, 2)),
    hjust = 1,                  
    size = 4
  ) +
  coord_flip() +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 13),
    axis.line = element_line(color = "black")
  ) +
  labs(
    x = NULL,
    y = "Relative importance (LMG)",
    title = "Relative contribution of predictors"
  )+
  annotate(
    "text",
    x = 1.5,                  # bottom category after flip
    y = max(imp$Importance) * 0.98,            # near right edge
    label = "R^2 == 56.34*'%'",
    hjust = 1,
    vjust = 1.2,
    size = 5,
    color = "black",
    parse = TRUE
  )




#visualize
ggplot(dat, aes(x = dom_lc_group, y = EH, fill = dom_lc_group)) +
  geom_violin(trim = FALSE, color = NA, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(x = "Dominant land cover", y = "EH")

#significance letter
library(emmeans)
library(dplyr)

#overall ANOVA null hypothesis: EH means are equal across groups'
fit <- aov(EH ~ dom_lc_group, data = dat)
#post-hoc pairwise comparison: EH means are equal between two groups 
cld_res <- multcomp::cld(emmeans(fit, ~ dom_lc_group), Letters = letters)

label_df <- dat %>%
  left_join(
    as.data.frame(cld_res) %>% dplyr::select(dom_lc_group, .group),
    by = "dom_lc_group"
  )
label_df$.group <- trimws(label_df$.group)


ggplot(dat, aes(x = dom_lc_group, y = EH, fill = dom_lc_group)) +
  geom_violin(trim = FALSE, color = NA, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  geom_text(data = label_df,
            aes(x = dom_lc_group, y = 0.65, label = .group),
            inherit.aes = FALSE,
            size = 5) +
  scale_fill_manual(values = nlcd_colors) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = 'black'),
    legend.position = "none"
  ) +
  labs(x = "Dominant land cover", y = "Environmental heterogeneity")





nlcd_colors <- c(
  "Open Water" = "#466b9f",
  
  "Developed" = "#d99282",  
  
  "Barren Land" = "#d1b26f",
  
  "Forest" = "#4daf4a",     
  
  "Shrub/Scrub" = "#dcd939",
  "Herbaceous" = "#c8d29b",
  
  "Pasture/Hay" = "#e6c84f",
  "Cultivated Crops" = "#e6a141",
  
  "Wetlands" = "#66c2a5"   
)
