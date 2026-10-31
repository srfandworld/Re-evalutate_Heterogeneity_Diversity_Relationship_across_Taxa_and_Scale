# =============================================================================
# 05_TablesS2-S3_model_summaries.R
# Table S2 (species richness) and Table S3 (log MSA): per taxon x scale
#   dAIC (linear - quadratic), R2 (%), AIC (quadratic), beta_EH, beta_EH2 (+ P),
#   turning point (raw EH units) and HDR shape (SI Section D.1)
#
# Source   : EmbeddingBio/analysis/tableS2_S3_HDR.R (2026-10-30): same
#            calculations, classification and formatting. That script refits
#            both models via the physicallyVar functions; here the quadratic
#            models are read from 01_fit_HDR_models.R (identical data, weights
#            and specification) and only the linear models are refitted.
#   - Linear models (EH^2 dropped): same rows, weights and covariates; cached in
#     data/clean/models/sar_*_linear.rds (tree refits: method = "LU", slow).
#   - R2 = cor(y, fitted)^2 (includes the SAR spatial term).
#   - Turning point = -b1 / (2 b2), back-transformed as center + sd * z with the
#     center/sd used to standardize EH (recovered exactly from the model data).
# Inputs   : data/clean/models/sar_*.rds; model data + weights via
#            functions/fun_hdr_sar_models.R
# Outputs  : output/tables/TableS2_richness.csv, TableS3_MSA.csv (formatted;
#            "*" = P < 0.05, bold in the SI) and *_full.csv (full precision)
# Run in a fresh R session (several GB of free memory needed for tree SAR fits).
# =============================================================================

source(here::here("code", "00_config.R"))
source(here::here("code", "functions", "fun_hdr_sar_models.R"))   # model data, weights, model functions


##########################1. HDR classification (SI Section D.1)##############################
#position of the turning point relative to the observed EH range
#  dAIC = AIC(linear) - AIC(quadratic)
classify_hdr <- function(dAIC, b1, p1, b2, p2, tp_position, alpha = 0.05) {
  sig1 <- p1 < alpha
  sig2 <- p2 < alpha
  within <- tp_position == "within"
  left   <- tp_position == "left"
  right  <- tp_position == "right"

  if (dAIC > 0 && sig2 && b2 < 0 && within) return("Unimodal")
  if (dAIC > 0 && sig2 && b2 > 0 && within) return("U-shaped")
  if ((!sig2 && sig1 && b1 > 0) ||
      (sig2 && b2 > 0 && left)  ||           #only the ascending arm of a U
      (sig2 && b2 < 0 && right))             #only the ascending arm of a hump
    return("Positive")
  if ((!sig2 && sig1 && b1 < 0) ||
      (sig2 && b2 > 0 && right) ||           #only the descending arm of a U
      (sig2 && b2 < 0 && left))              #only the descending arm of a hump
    return("Negative")
  "NS"
}


##########################2. models, linear refits and one table row##########################
spec <- list(
  tree   = list(label = "Tree",   prefix = "FIA",    scales = scales_tree,
                files = c(richness = "sar_tree_richness_iNEXT.rds", MSA = "sar_tree_MSA_iNEXT.rds"),
                lhs   = c(richness = "log(qD)", MSA = "log(MSA_new_{scale})")),
  bird   = list(label = "Bird",   prefix = "BBS",    scales = scales_bird,
                files = c(richness = "sar_bird_richness.rds", MSA = "sar_bird_MSA.rds"),
                lhs   = c(richness = "richness", MSA = "log(MSA)")),
  mammal = list(label = "Mammal", prefix = "mammal", scales = scales_mammal,
                files = c(richness = "sar_mammal_richness.rds", MSA = "sar_mammal_MSA.rds"),
                lhs   = c(richness = "richness_Estimator", MSA = "log(MSA)"))
)

#linear SAR: same response, rows, weights and covariates as the quadratic model
fit_linear <- function(m, taxa, scale, lhs) {
  d  <- get(glue("{spec[[taxa]]$prefix}_table_{scale}"))      # model data (fun_hdr_sar_models.R)
  stopifnot(nrow(d) == length(m$Y))
  lw <- switch(taxa,
               tree   = build_lw_tree(d, plots_sf_tree),
               bird   = build_lw_bird(d, routes_sf_bird),
               mammal = build_lw_mammal(d, centroids_sf_mammal))
  quad  <- as.character(glue("I(EH_{scale}_t^2)"))
  preds <- setdiff(colnames(m$X)[-1], quad)
  form  <- as.formula(paste(glue(lhs), "~", paste(preds, collapse = " + ")))
  if (taxa == "tree") {
    spautolm(formula = form, data = d, listw = lw, method = "LU", family = "SAR")
  } else {
    spautolm(formula = form, data = d, listw = lw, family = "SAR")
  }
}

#one row of Table S2/S3 (same calculations as tableS2_S3_HDR.R)
hdr_row <- function(m_quad, m_lin, taxa, scale) {
  d  <- get(glue("{spec[[taxa]]$prefix}_table_{scale}"))
  eh <- as.character(glue("EH_{scale}_t"))

  coef_tab <- summary(m_quad)$Coef
  quad <- as.character(glue("I({eh}^2)"))
  b1 <- coef_tab[eh,   "Estimate"]; p1 <- coef_tab[eh,   "Pr(>|z|)"]
  b2 <- coef_tab[quad, "Estimate"]; p2 <- coef_tab[quad, "Pr(>|z|)"]
  dAIC <- AIC(m_lin) - AIC(m_quad)          # > 0: quadratic model fits better

  #turning point on the standardized scale, and where it falls in the data
  z  <- -b1 / (2 * b2)
  xr <- range(as.data.frame(m_quad$X, check.names = FALSE)[[eh]])
  tp_position <- if (!is.finite(z)) "none" else if (z < xr[1]) "left" else
    if (z > xr[2]) "right" else "within"

  #back-transform to raw EH: EH = center + sd * EH_t (exact)
  bt <- unname(coef(lm(d[[glue("EH_{scale}")]] ~ d[[eh]])))
  turning_point <- bt[1] + bt[2] * z

  data.frame(
    Taxa          = spec[[taxa]]$label,
    Scale         = scale,
    dAIC          = dAIC,
    R2            = 100 * cor(m_quad$Y, m_quad$fit$fitted.values)^2,
    AIC           = AIC(m_quad),
    beta_EH       = b1,
    p_EH          = p1,
    beta_EH2      = b2,
    p_EH2         = p2,
    turning_point = turning_point,
    tp_position   = tp_position,
    HDR           = classify_hdr(dAIC, b1, p1, b2, p2, tp_position),
    n             = length(m_quad$Y)
  )
}

#all taxa x scales for one response
hdr_table <- function(response) {
  bind_rows(lapply(names(spec), function(taxa) {
    s    <- spec[[taxa]]
    quad <- readRDS(p_model(s$files[[response]]))
    if (is.null(names(quad))) names(quad) <- s$scales

    lin_file <- p_model(sub("\\.rds$", "_linear.rds", s$files[[response]]))
    lin <- if (file.exists(lin_file)) readRDS(lin_file) else list()
    for (sc in setdiff(s$scales, names(lin))) {
      message("refitting linear SAR: ", taxa, " ", response, " ", sc)
      lin[[sc]] <- fit_linear(quad[[sc]], taxa, sc, s$lhs[[response]])
      saveRDS(lin, lin_file)                  # cache after every fit
    }
    bind_rows(lapply(s$scales, function(sc) hdr_row(quad[[sc]], lin[[sc]], taxa, sc)))
  }))
}

table_S2 <- hdr_table("richness")
table_S3 <- hdr_table("MSA")


##########################3. format as in the SI tables and save##############################
#significant coefficients (P < 0.05) are marked with "*" (bold in the Word table)
format_table <- function(tab) {
  sig <- function(x, p) paste0(formatC(x, format = "f", digits = 3), ifelse(p < 0.05, "*", ""))
  tab %>%
    transmute(Taxa, Scale,
              dAIC            = round(dAIC, 1),
              `R2 (%)`        = round(R2, 2),
              AIC             = round(AIC, 0),
              beta_EH         = sig(beta_EH, p_EH),
              beta_EH2        = sig(beta_EH2, p_EH2),
              `Turning point` = round(turning_point, 3),
              HDR)
}

write.csv(table_S2, p_tab("TableS2_richness_full.csv"), row.names = FALSE)
write.csv(table_S3, p_tab("TableS3_MSA_full.csv"),      row.names = FALSE)
write.csv(format_table(table_S2), p_tab("TableS2_richness.csv"), row.names = FALSE)
write.csv(format_table(table_S3), p_tab("TableS3_MSA.csv"),      row.names = FALSE)

format_table(table_S2)
format_table(table_S3)

#dAIC range reported in the Results
range(c(table_S2$dAIC, table_S3$dAIC))
