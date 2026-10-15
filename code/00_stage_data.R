# =============================================================================
# 00_stage_data.R - copy existing inputs/intermediates into data/
#
# Reads data/data_manifest.csv (dest, source, role, ...) and copies every file
# or folder that exists at `source` to data/<dest>. Use it to:
#   - run later stages (models, figures) from the saved intermediates without
#     re-running the slow upstream steps, and
#   - collect the files for the Figshare data deposit.
# Rows with exists = FALSE are listed at the end - those are the gaps.
#
# The FIA DataMart row (raw/FIA/) is skipped by default because it is the
# whole EmbeddingBio/data folder; set COPY_FIA_RAW <- TRUE to copy only the
# FIA state CSVs (<STATE>_<TABLE>.csv).
# =============================================================================

source(here::here("code", "00_config.R"))

COPY_FIA_RAW <- FALSE
manifest <- read.csv(here::here("data", "data_manifest.csv"), stringsAsFactors = FALSE)

copy_one <- function(src, dest) {
  dest <- here::here("data", dest)
  if (dir.exists(src)) {
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    file.copy(list.files(src, full.names = TRUE), dest, recursive = TRUE, overwrite = FALSE)
  } else {
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    file.copy(src, dest, overwrite = FALSE)
  }
}

for (i in seq_len(nrow(manifest))) {
  r <- manifest[i, ]
  if (!file.exists(r$source)) next
  if (r$dest == "raw/FIA/") {
    if (!COPY_FIA_RAW) next
    fia_csv <- list.files(r$source, pattern = "^[A-Z]{2}_[A-Z_0-9]+\\.csv$", full.names = TRUE)
    dir.create(p_raw("FIA"), showWarnings = FALSE)
    file.copy(fia_csv, p_raw("FIA"), overwrite = FALSE)
    next
  }
  message("copy: ", r$dest)
  copy_one(r$source, r$dest)
}

missing <- manifest[!file.exists(manifest$source), c("dest", "note")]
if (nrow(missing) > 0) {
  message("\nNot found (see note):")
  print(missing, row.names = FALSE)
}
