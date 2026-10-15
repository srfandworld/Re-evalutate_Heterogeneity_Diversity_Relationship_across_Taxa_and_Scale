# =============================================================================
# FigS08_EH_across_scales.R
# Fig. S8: EH distribution by spatial scale (violin + box) with paired
# Wilcoxon letters (BH-adjusted), one panel per taxon
#
# Source   : EmbeddingBio/analysis/sfigure_violin_EH.R. The original overwrote
#            EH_long three times (tree, bird, mammal) and plotted whichever ran
#            last; the three tables are loaded below so each block can be run
#            in turn. Logic unchanged.
# Inputs   : data/clean/FIA_table.csv, BBS_table.csv, mammal_table.csv
# TODO(reviewer, Fig. S8 caption): the 400 m-50 km range is not available for
#   every taxon (trees lack 400 m; birds lack 50 km); x-axis label hard-codes
#   n = 132,262 (tree plots). Save each panel with ggsave(p_fig(...)).
# =============================================================================

source(here::here("code", "00_config.R"))
library(tidyverse)

FIA_table    <- read.csv(p_clean("FIA_table.csv"))
BBS_table    <- read.csv(p_clean("BBS_table.csv"))
mammal_table <- read.csv(p_clean("mammal_table.csv"))

####violin plot for EH comparison at different spatial scales

library(ggpubr)

EH_long <- FIA_table %>%
  dplyr::select(matches("^EH_\\d+km$")) %>%
  pivot_longer(
    cols = everything(),
    names_to = "scale",
    values_to = "EH"
  ) %>%
  mutate(
    scale = str_extract(scale, "\\d+km"),
    scale = factor(scale, levels = c("1km", "5km", "10km", "20km", "50km"))
  )
EH_long <- BBS_table %>%
  dplyr::select(matches("^EH_(400m|1km|5km|10km|20km)$")) %>%
  pivot_longer(
    cols = everything(),
    names_to = "scale",
    values_to = "EH"
  ) %>%
  mutate(
    scale = str_extract(scale, "\\d+(km|m)"),
    scale = factor(scale, levels = c("400m", "1km", "5km", "10km", "20km"))
  )
EH_long <- mammal_table %>%
  dplyr::select(matches("^EH_\\d+km$")) %>%
  pivot_longer(
    cols = everything(),
    names_to = "scale",
    values_to = "EH"
  ) %>%
  mutate(
    scale = str_extract(scale, "\\d+km"),
    scale = factor(scale, levels = c("1km", "5km", "10km", "20km", "50km"))
  )

library(multcompView)

pw <- pairwise.wilcox.test(
  EH_long$EH,
  EH_long$scale,
  paired = TRUE,
  p.adjust.method = "BH"
)

lvls <- levels(EH_long$scale)

# build full matrix
p_full <- matrix(1, length(lvls), length(lvls),
                 dimnames = list(lvls, lvls))

for (r in rownames(pw$p.value)) {
  for (c in colnames(pw$p.value)) {
    if (!is.na(pw$p.value[r, c])) {
      p_full[r, c] <- pw$p.value[r, c]
      p_full[c, r] <- pw$p.value[r, c]
    }
  }
}

diag(p_full) <- 1
p_full[is.na(p_full)] <- 1

letters <- multcompLetters(p_full)$Letters
letters

letter_df <- data.frame(
  scale = names(letters),
  label = letters
)

# position letters above each violin
y_pos <- EH_long %>%
  group_by(scale) %>%
  summarise(y = max(EH, na.rm = TRUE) + 0.06, .groups = "drop")

letter_df <- left_join(letter_df, y_pos, by = "scale")

c("Tree"   = "#62b187",
  "Bird"   = "#65afc1",
  "Mammal" = "#e8964c")

ggplot(EH_long, aes(x = scale, y = EH)) +
  geom_violin(fill = "#62b187", color = "black", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  
  geom_text(
    data = letter_df,
    aes(x = scale, y = y, label = label),
    inherit.aes = FALSE,
    size = 6
  ) +
  
  theme_classic(base_size = 13) +
  labs(
    x = "Spatial scale (n = 132,262)",
    y = "Environmental heterogeneity (EH)"
  ) +
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 12)
  )

