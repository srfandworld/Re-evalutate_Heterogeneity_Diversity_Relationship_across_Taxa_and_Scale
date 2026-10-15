# =============================================================================
# 02_tree_richness_iNEXT.R
# Coverage-based tree richness (qD at 90% sample coverage) per buffer
#
# Source   : EmbeddingBio/analysis/FIA_richness_recalc.R, lines 9-164 (the
#            version that wrote the *_iNEXT tables behind Tables S2-S3 / Fig. 2).
#            Code/Tree/01 held a copy that could not run (n_plots_species missing).
# Inputs   : data/clean/FIA_gamma_abundance/comm_abundance_r{5,10,20,50}km.rds
#            data/clean/FIA_table_{1,5,10,20,50}km.rds, data/clean/fia_sp.csv
# Outputs  : data/clean/FIA_table_{1,5,10,20,50}km_iNEXT.csv   (model input)
#
# Scale-specific rules (TODO reviewer M2: state these in the Methods):
#   1 km     : focal plots with no other plot within 1 km; observed richness of
#              that single plot (no pooling, no iNEXT)
#   5 km     : buffers with >= 6 plots   (n_target = 6)
#   10-50 km : buffers with >= 10 plots  (n_target = 10)
# These rules make the plot set differ by scale (reviewer C4: 18,577 plots at
# 5 km vs 129,143 at 50 km).
# =============================================================================

source(here::here("code", "00_config.R"))

library(data.table)
library(dplyr)
library(purrr)
library(iNEXT)
library(tidyr)
estimate_FIA_richness <- function(
    scale,
    fia_sp,
    gamma_file = NULL,
    table_file = NULL,
    output_file = NULL,
    n_target = 10,
    coverage = 0.90
) {
  
  scale <- as.character(scale)
  
  #-----------------------------
  # Special case: 1 km
  #-----------------------------
  if (scale == "1km") {
    
    FIA_table <- readRDS(table_file)
    
    result <- FIA_table %>%
      filter(.data[[paste0("effort_", scale)]] == 1) %>%
      mutate(
        qD = .data[[paste0("richness_", scale)]]
      ) %>%
      rename(focal_pltID = pltID) %>%
      relocate(qD, .before = all_of(paste0("richness_", scale)))
    
    if (!is.null(output_file)) {
      write.csv(result, output_file, row.names = FALSE)
    }
    
    return(result)
  }
  
  
  #-----------------------------
  # Other scales: iNEXT
  #-----------------------------
  gamma <- readRDS(gamma_file) %>%
    filter(SCIENTIFIC_NAME %in% fia_sp$SCIENTIFIC_NAME) %>%
    as.data.table()
  
  all_species <- sort(unique(gamma$SCIENTIFIC_NAME))
  
  incidence_list <- lapply(
    split(gamma, gamma$focal_pltID),
    function(d) {
      
      freq <- setNames(
        d$n_plots_species,
        d$SCIENTIFIC_NAME
      )
      
      freq <- freq[all_species]
      freq[is.na(freq)] <- 0
      
      c(
        nT = d$effort[1],
        freq
      )
    }
  )
  
  # keep enough buffers to estimate richness at 90% coverage
  valid_buffers <- incidence_list[
    sapply(incidence_list, function(x) x[["nT"]] >= n_target)
  ]
  
  richness_results <- purrr::imap_dfr(
    valid_buffers,
    function(one_buffer, focal_pltID) {
      
      tryCatch({
        
        out2 <- iNEXT::estimateD(
          one_buffer,
          q = 0,
          datatype = "incidence_freq",
          base = "coverage",
          level = coverage
        )
        
        out2 %>%
          mutate(
            focal_pltID = focal_pltID,
            n_plots = one_buffer[["nT"]]
          )
        
      }, error = function(e) {
        
        message("Skipped buffer: ", focal_pltID)
        message("Reason: ", e$message)
        
        NULL
      })
    }
  )
  
  
  #-----------------------------
  # Join with FIA table
  #-----------------------------
  FIA_table <- readRDS(table_file)
  
  result <- richness_results %>%
    select(focal_pltID, n_plots, qD) %>%
    left_join(
      FIA_table,
      by = join_by(focal_pltID == pltID)
    )
  
  if (!is.null(output_file)) {
    write.csv(result, output_file, row.names = FALSE)
  }
  
  return(result)
}

fia_sp <- read.csv(p_clean("fia_sp.csv"))
gamma_path <- function(scale) p_clean("FIA_gamma_abundance", paste0("comm_abundance_r", scale, ".rds"))

FIA_table_1km <- estimate_FIA_richness(
  scale = "1km",
  fia_sp = fia_sp,
  table_file = p_clean("FIA_table_1km.rds"),
  output_file = p_clean("FIA_table_1km_iNEXT.csv")
)
FIA_table_5km <- estimate_FIA_richness(
  scale = "5km",
  fia_sp = fia_sp,
  gamma_file = gamma_path("5km"),
  table_file = p_clean("FIA_table_5km.rds"),
  n_target = 6,
  output_file = p_clean("FIA_table_5km_iNEXT.csv")
)
FIA_table_10km <- estimate_FIA_richness(
  scale = "10km",
  fia_sp = fia_sp,
  gamma_file = gamma_path("10km"),
  table_file = p_clean("FIA_table_10km.rds"),
  n_target = 10,
  output_file = p_clean("FIA_table_10km_iNEXT.csv")
)
FIA_table_20km <- estimate_FIA_richness(
  scale = "20km",
  fia_sp = fia_sp,
  gamma_file = gamma_path("20km"),
  table_file = p_clean("FIA_table_20km.rds"),
  n_target = 10,
  output_file = p_clean("FIA_table_20km_iNEXT.csv")
)
FIA_table_50km <- estimate_FIA_richness(
  scale = "50km",
  fia_sp = fia_sp,
  gamma_file = gamma_path("50km"),
  table_file = p_clean("FIA_table_50km.rds"),
  n_target = 10,
  output_file = p_clean("FIA_table_50km_iNEXT.csv")
)