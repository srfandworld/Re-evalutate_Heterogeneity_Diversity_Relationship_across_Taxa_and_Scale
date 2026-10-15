# =============================================================================
# 03_FigS9-S11_partial_residuals.R
# Partial-residual plots of richness and log(MSA) vs EH at every scale
#   Fig. S9 trees, Fig. S10 birds, Fig. S11 mammals
#   (the 10 km / 5 km panels are the same as Fig. 2c-d, h-i, m-n)
#
# Source   : EmbeddingBio/analysis/figure2_fit.R (paths only changed)
# Inputs   : data/clean/models/sar_*.rds, per-scale tables
# Outputs  : patchwork objects printed to the device.
#            TODO: add ggsave(p_fig("FigS9_tree_partial.png"), ...) etc. so the
#            SI figures are written to output/figures.
# NOTE     : reviewer listed S9-S11 as "helper functions only / not tested".
# =============================================================================

source(here::here("code", "00_config.R"))
source(here::here("code", "functions", "fun_sar_partial_plot.R"))
library(tidyverse)
library(patchwork)
library(glue)

for (s in SCALES_TREE)   assign(glue("FIA_table_{s}"),    read.csv(p_clean(glue("FIA_table_{s}_iNEXT.csv"))))
for (s in SCALES_BIRD)   assign(glue("BBS_table_{s}"),    read.csv(p_clean(glue("BBS_table_{s}.csv"))))
for (s in SCALES_MAMMAL) assign(glue("mammal_table_{s}"), read.csv(p_clean(glue("mammal_table_{s}.csv"))))

#################tree###############################
# sar_tree_richness <- lapply(c("1km","5km","10km","20km","50km"),function(s){
#   sar_modeling_tree(scale = s,y="richness")
# })
# sar_tree_MSA <- lapply(c("1km","5km","10km","20km","50km"),function(s){
#   sar_modeling_tree(scale = s,y="log(MSA)")
# })
# write_rds(sar_tree_richness,p_model("sar_tree_richness.rds"))
# write_rds(sar_tree_MSA,p_model("sar_tree_MSA.rds"))


sar_tree_richness <- readRDS(p_model("sar_tree_richness_iNEXT.rds"))
sar_tree_MSA <- readRDS(p_model("sar_tree_MSA_iNEXT.rds"))

#1. richness ~ EH partial resiudal plot at 10km
p1_tree <- sar_partial_plot(
  model     = sar_tree_richness[[3]],
  taxa      = "tree",
  var       = "EH_10km_t",
  scale     = "10km",
  quadratic = TRUE,
  line_color = "steelblue",
  x_label   = "Environmental heterogeneity \n(standardized)",
  y_label   = "Residual richness"
)
p1_tree

#2. log(MSA) ~ EH partial residual plot at 10km

p2_tree <- sar_partial_plot(
  model     = sar_tree_MSA[[3]],
  taxa      = "tree",
  var       = "EH_10km_t",
  scale     = "10km",
  quadratic = TRUE,
  line_color = "firebrick",
  x_label   = "Environmental heterogeneity \n(standardized)",
  y_label   = "Residual log (mean species abundance)"
)

p2_tree


#3. positive scaling effects
scales_num <- c( 1, 5, 10, 20, 50)   # adjust — your last scale is 20km not 50km
results_richness <- lapply(sar_tree_richness, function(m){inflection_boot(m,'tree')}) %>%
  bind_rows() %>%
  mutate(scale_num = scales_num,
         response  = "Species richness")
results_MSA <- lapply(sar_tree_MSA,  function(m){inflection_boot(m,'tree')}) %>%
  bind_rows() %>%
  mutate(scale_num = scales_num,
         response  = "log (Mean species abundance)")
results_all <- bind_rows(results_richness, results_MSA)

# ── Plot ──────────────────────────────────────────────────────────────────────
p3_tree <- ggplot(results_all, aes(x = scale_num, y = z,
                                   color = response, fill = response)) +
  # CI band
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95),
              alpha = 0.15, color = NA) +
  # Error bars
  geom_errorbar(aes(ymin = lwr_95, ymax = upr_95),
                width = 0.05, linewidth = 0.6) +
  # Points
  geom_point(size = 3) +
  # Trend line per response
  geom_smooth(method = "lm", se = FALSE,
              linetype = "dashed", linewidth = 0.8) +
  
  scale_x_log10(
    breaks = scales_num,
    labels = c( "1km", "5km", "10km", "20km", "50km")
  ) +
  scale_color_manual(
    values = c("Species richness"  = "steelblue",
               "log (Mean species abundance)"  = "firebrick"),
    name   = "Response"
  ) +
  scale_fill_manual(
    values = c("Species richness"  = "steelblue",
               "log (Mean species abundance)"  = "firebrick"),
    name   = "Response"
  ) +
  labs(
    x = "Spatial scale",
    y = "Inflection point"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x    = element_text(size = 11),
    legend.position = c(0.5,1),
    legend.background = element_rect(fill=NA,color = NA),
    legend.direction = "vertical",
    legend.margin    = margin(t = 10, b = 2),
    legend.title     = element_blank(),
    legend.text      = element_text(size = 12),
  )

p3_tree

#layout 1 2 3
p1_tree + p2_tree + p3_tree+
  plot_layout(nrow = 1)


#4. Supplementary figures: richness across scales + MSA across scales
scales_num <- c( 1, 5, 10, 20, 50)   # adjust — your last scale is 20km not 50km
p_richness <- lapply(sar_tree_richness, function(m){
  var <- colnames(m$model)[2]
  scale <- str_extract(var, "\\d+(km|m)")
  sar_partial_plot(m, taxa='tree',  
                   var = var,  scale= scale, 
                   quadratic = TRUE,
                   line_color = 'steelblue',  
                   x_label   = "Environmental heterogeneity \n(standardized)",
                   y_label   = "Residual richness") + 
    labs(title = scale) +
    theme(plot.title = element_text(hjust = 0.5))
}) 
p_MSA <- lapply(sar_tree_MSA,  function(m){
  var <- colnames(m$model)[2]
  scale <- str_extract(var, "\\d+(km|m)")
  sar_partial_plot(m, taxa='tree',  
                   var = var,  scale= scale, 
                   quadratic = TRUE,
                   line_color = 'firebrick',  
                   x_label   = "Environmental heterogeneity \n(standardized)",
                   y_label   = "Residual log (mean species abundance)")}) 

#layout
# remove y-axis titles except first column (plots 1 and 6)
p_richness <- lapply(seq_along(p_richness), function(i) {
  if (i %% 5 != 1) {
    p_richness[[i]] + theme(axis.title.y = element_blank())
  } else {
    p_richness[[i]]
  }
})

p_MSA <- lapply(seq_along(p_MSA), function(i) {
  if (i %% 5 != 1) {
    p_MSA[[i]] + theme(axis.title.y = element_blank())
  } else {
    p_MSA[[i]]
  }
})

# remove x-axis titles for first row (richness row)
p_richness <- lapply(p_richness, function(p) {
  p + theme(axis.title.x = element_blank())
})

# combine
wrap_plots(c(p_richness, p_MSA), nrow = 2)












#################bird###############################
# sar_bird_richness <- lapply(c("400m","1km","5km","10km","20km"),function(s){
#   sar_modeling_bird(scale = s,y="richness")
# })
# sar_bird_MSA <- lapply(c("400m","1km","5km","10km","20km"),function(s){
#   sar_modeling_bird(scale = s,y="log(MSA)")
# })
# write_rds(sar_bird_richness,p_model("sar_bird_richness.rds"))
# write_rds(sar_bird_MSA,p_model("sar_bird_MSA.rds"))


sar_bird_richness <- readRDS(p_model("sar_bird_richness.rds"))
sar_bird_MSA <- readRDS(p_model("sar_bird_MSA.rds"))

#1. richness ~ EH partial resiudal plot at 10km
p1_bird <- sar_partial_plot(
  model     = sar_bird_richness[[4]],
  taxa      = "bird",
  var       = "EH_10km_t",
  scale     = "10km",
  quadratic = TRUE,
  line_color = "steelblue",
  x_label   = "Environmental heterogeneity \n(standardized)",
  y_label   = "Residual richness"
)
p1_bird
  
#2. log(MSA) ~ EH partial residual plot at 10km

p2_bird <- sar_partial_plot(
  model     = sar_bird_MSA[[4]],
  taxa      = "bird",
  var       = "EH_10km_t",
  scale     = "10km",
  quadratic = TRUE,
  line_color = "firebrick",
  x_label   = "Environmental heterogeneity \n(standardized)",
  y_label   = "Residual log (mean species abundance)"
)

p2_bird


#3. positive scaling effects
scales_num <- c(0.4, 1, 5, 10, 20)   # adjust — your last scale is 20km not 50km
results_richness <- lapply(sar_bird_richness, function(m){inflection_boot(m,'bird')}) %>%
  bind_rows() %>%
  mutate(scale_num = scales_num,
         response  = "Species richness")
results_MSA <- lapply(sar_bird_MSA, function(m){inflection_boot(m,'bird')}) %>%
  bind_rows() %>%
  mutate(scale_num = scales_num,
         response  = "log (Mean species abundance)")
results_all <- bind_rows(results_richness, results_MSA)

# ── Plot ──────────────────────────────────────────────────────────────────────
p3_bird <- ggplot(results_all, aes(x = scale_num, y = z,
                        color = response, fill = response)) +
  # CI band
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95),
              alpha = 0.15, color = NA) +
  # Error bars
  geom_errorbar(aes(ymin = lwr_95, ymax = upr_95),
                width = 0.05, linewidth = 0.6) +
  # Points
  geom_point(size = 3) +
  # Trend line per response
  geom_smooth(method = "lm", se = FALSE,
              linetype = "dashed", linewidth = 0.8) +
  
  scale_x_log10(
    breaks = scales_num,
    labels = c("400m", "1km", "5km", "10km", "20km")
  ) +
  scale_color_manual(
    values = c("Species richness"  = "steelblue",
               "log (Mean species abundance)"  = "firebrick"),
    name   = "Response"
  ) +
  scale_fill_manual(
    values = c("Species richness"  = "steelblue",
               "log (Mean species abundance)"  = "firebrick"),
    name   = "Response"
  ) +
  labs(
    x = "Spatial scale",
    y = "Inflection point"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x    = element_text(size = 11),
    legend.position = c(0.5,1),
    legend.background = element_rect(fill=NA,color = NA),
    legend.direction = "vertical",
    legend.margin    = margin(t = 10, b = 2),
    legend.title     = element_blank(),
    legend.text      = element_text(size = 12),
  )


 p1_bird + p2_bird + p3_bird+
  plot_layout(nrow = 1)



 #4. Supplementary figures: richness across scales + MSA across scales
 scales_num <- c(0.4, 1, 5, 10, 20) 
 p_richness <- lapply(sar_bird_richness, function(m){
   var <- colnames(m$X)[2]
   scale <- str_extract(var, "\\d+(km|m)")
   sar_partial_plot(m, taxa='bird',  
                    var = var,  scale= scale, 
                    quadratic = TRUE,
                    line_color = 'steelblue',  
                    x_label   = "Environmental heterogeneity \n(standardized)",
                    y_label   = "Residual richness") + 
     labs(title = scale) +
     theme(plot.title = element_text(hjust = 0.5))}) 
 
 p_MSA <- lapply(sar_bird_MSA,  function(m){
   var <- colnames(m$X)[2]
   scale <- str_extract(var, "\\d+(km|m)")
   sar_partial_plot(m, taxa='bird',  
                    var = var,  scale= scale, 
                    quadratic = TRUE,
                    line_color = 'firebrick',  
                    x_label   = "Environmental heterogeneity \n(standardized)",
                    y_label   = "Residual log (mean species abundance)")}) 
 
 
 #layout
 # remove y-axis titles except first column (plots 1 and 6)
 p_richness <- lapply(seq_along(p_richness), function(i) {
   if (i %% 5 != 1) {
     p_richness[[i]] + theme(axis.title.y = element_blank())
   } else {
     p_richness[[i]]
   }
 })
 
 p_MSA <- lapply(seq_along(p_MSA), function(i) {
   if (i %% 5 != 1) {
     p_MSA[[i]] + theme(axis.title.y = element_blank())
   } else {
     p_MSA[[i]]
   }
 })
 
 # remove x-axis titles for first row (richness row)
 p_richness <- lapply(p_richness, function(p) {
   p + theme(axis.title.x = element_blank())
 })
 
 # combine
 wrap_plots(c(p_richness, p_MSA), nrow = 2)
 
 










#################mammal###############################
# sar_mammal_richness <- lapply(c("1km","5km","10km","20km","50km"),function(s){
#   sar_modeling_mammal(scale = s,y="richness_Estimator")
# })
# sar_mammal_MSA <- lapply(c("1km","5km","10km","20km","50km"),function(s){
#   sar_modeling_mammal(scale = s,y="log(MSA)")
# })
# write_rds(sar_mammal_richness,p_model("sar_mammal_richness.rds"))
# write_rds(sar_mammal_MSA,p_model("sar_mammal_MSA.rds"))


sar_mammal_richness <- readRDS(p_model("sar_mammal_richness.rds"))
sar_mammal_MSA <- readRDS(p_model("sar_mammal_MSA.rds"))

#1. richness ~ EH partial resiudal plot at 10km
p1_mammal <- sar_partial_plot(
  model     = sar_mammal_richness[[1]],
  taxa      = "mammal",
  var       = "EH_5km_t",
  scale     = "1km",
  quadratic = TRUE,
  line_color = "steelblue",
  x_label   = "Environmental heterogeneity \n(standardized)",
  y_label   = "Residual richness"
)
p1_mammal

#2. log(MSA) ~ EH partial residual plot at 10km

p2_mammal <- sar_partial_plot(
  model     = sar_mammal_MSA[[1]],
  taxa      = "mammal",
  var       = "EH_5km_t",
  scale     = "1km",
  quadratic = TRUE,
  line_color = "firebrick",
  x_label   = "Environmental heterogeneity \n(standardized)",
  y_label   = "Residual log (mean species abundance)"
)

p2_mammal


#3. positive scaling effects
scales_num <- c(1, 5, 10, 20, 50) 
results_richness <- lapply(sar_mammal_richness,  function(m){inflection_boot(m,'mammal')}) %>%
  bind_rows() %>%
  mutate(scale_num = scales_num,
         response  = "Species richness")
results_MSA <- lapply(sar_mammal_MSA,  function(m){inflection_boot(m,'mammal')}) %>%
  bind_rows() %>%
  mutate(scale_num = scales_num,
         response  = "log (Mean species abundance)")
results_all <- bind_rows(results_richness, results_MSA)

# ── Plot ──────────────────────────────────────────────────────────────────────
p3_mammal <- ggplot(results_all, aes(x = scale_num, y = z,
                                   color = response, fill = response)) +
  # CI band
  geom_ribbon(aes(ymin = lwr_95, ymax = upr_95),
              alpha = 0.15, color = NA) +
  # Error bars
  geom_errorbar(aes(ymin = lwr_95, ymax = upr_95),
                width = 0.05, linewidth = 0.6) +
  # Points
  geom_point(size = 3) +
  # Trend line per response
  geom_smooth(method = "lm", se = FALSE,
              linetype = "dashed", linewidth = 0.8) +
  
  scale_x_log10(
    breaks = scales_num,
    labels = c( "1km", "5km", "10km", "20km", "50km")
  ) +
  scale_color_manual(
    values = c("Species richness"  = "steelblue",
               "log (Mean species abundance)"  = "firebrick"),
    name   = "Response"
  ) +
  scale_fill_manual(
    values = c("Species richness"  = "steelblue",
               "log (Mean species abundance)"  = "firebrick"),
    name   = "Response"
  ) +
  labs(
    x = "Spatial scale",
    y = "Inflection point"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x    = element_text(size = 11),
    legend.position = c(0.5,1),
    legend.background = element_rect(fill=NA,color = NA),
    legend.direction = "vertical",
    legend.margin    = margin(t = 10, b = 2),
    legend.title     = element_blank(),
    legend.text      = element_text(size = 12),
  )

p3_mammal

#p_mammal_richness + p_mammal_abundance +
p1_mammal + p2_mammal + p3_mammal+
  plot_layout(nrow = 1)






#4. Supplementary figures: richness across scales + MSA across scales
scales_num <- c(1, 5, 10, 20, 50)
p_richness <- lapply(sar_mammal_richness, function(m){
  var <- colnames(m$X)[2]
  scale <- str_extract(var, "\\d+(km|m)")
  
  sar_partial_plot(m, taxa='mammal',  
                   var = var,  scale= scale, 
                   quadratic = TRUE,
                   line_color = 'steelblue',  
                   x_label   = "Environmental heterogeneity \n(standardized)",
                   y_label   = "Residual richness") + 
    labs(title = scale) +
    theme(plot.title = element_text(hjust = 0.5))
  }) 

p_MSA <- lapply(sar_mammal_MSA,  function(m){
  var <- colnames(m$X)[2]
  scale <- str_extract(var, "\\d+(km|m)")
  
  sar_partial_plot(m, taxa='mammal',  
                   var = var,  scale= scale, 
                   quadratic = TRUE,
                   line_color = 'firebrick',  
                   x_label   = "Environmental heterogeneity \n(standardized)",
                   y_label   = "Residual log (mean species abundance)")}) 


#layout
# remove y-axis titles except first column (plots 1 and 6)
p_richness <- lapply(seq_along(p_richness), function(i) {
  if (i %% 5 != 1) {
    p_richness[[i]] + theme(axis.title.y = element_blank())
  } else {
    p_richness[[i]]
  }
})

p_MSA <- lapply(seq_along(p_MSA), function(i) {
  if (i %% 5 != 1) {
    p_MSA[[i]] + theme(axis.title.y = element_blank())
  } else {
    p_MSA[[i]]
  }
})

# remove x-axis titles for first row (richness row)
p_richness <- lapply(p_richness, function(p) {
  p + theme(axis.title.x = element_blank())
})

# combine
wrap_plots(c(p_richness, p_MSA), nrow = 2)

