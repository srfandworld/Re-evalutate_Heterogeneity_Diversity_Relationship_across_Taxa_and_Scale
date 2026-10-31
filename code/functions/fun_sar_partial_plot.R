# =============================================================================
# fun_sar_partial_plot.R  (source: EmbeddingBio/analysis/fun_sar_partial_plot.R,
#                         2026-10-30; verbatim)
# Used by: 05_HDR_models/03_FigS9-S11_partial_residuals.R
# =============================================================================

#partial residual plot of richness / MSA ~ EH from a SAR model (spautolm)
#restored from the old fun_sar_partial_plot.R (lost); same plotting style:
#  grey partial residuals, fitted EH (+ EH^2) curve with 95% band from coefficient
#  simulation, dotted line at the turning point, R^2 and turning point annotated
#changes vs the old version:
#  - all three taxa are SAR now (trees too), so the same code path is used
#  - EH values are taken from the model's design matrix (model$X), so the points
#    line up row-for-row with residuals(model) whatever data the model was fit on
#    (no need for FIA_table_<s> / BBS_table_<s> / mammal_table_<s> in the workspace)
library(dplyr)
library(glue)
library(ggplot2)

sar_partial_plot <- function(model,        # spautolm object
                             var,          # predictor name e.g. "EH_400m_t"
                             scale      = NULL,  # kept for compatibility (not used)
                             taxa       = NULL,  # kept for compatibility (not used)
                             quadratic  = FALSE, # TRUE if var has I(var^2) term
                             nsim       = 500,
                             n_seq      = 200,
                             seed       = 123,
                             point_size = 0.4,
                             point_alpha= 0.5,
                             line_color = "steelblue",
                             x_label    = NULL,
                             y_label    = "Partial residual (log MSA)") {

  # ── SAR coefficient table and R2 ─────────────────────────────────────────────
  sar_coef_table <- summary(model)$Coef
  R2 <- cor(model$Y, model$fit$fitted.values)^2

  # ── Extract model components ────────────────────────────────────────────────
  coefs  <- coef(model)
  resids <- residuals(model)
  xvar   <- model$X[, var]                 # EH values used in the fit

  # ── Predictor sequence ──────────────────────────────────────────────────────
  x_seq <- seq(min(xvar, na.rm = TRUE),
               max(xvar, na.rm = TRUE),
               length.out = n_seq)

  # ── Coefficient and SE names ────────────────────────────────────────────────
  quad_name  <- paste0("I(", var, "^2)")   # e.g. "I(EH_400m_t^2)"
  has_quad   <- quadratic && (quad_name %in% names(coefs))

  if (has_quad) {
    coef_names <- c(var, quad_name)
  } else {
    coef_names <- var
  }

  # ── Fitted line (predictor contribution only, no intercept) ─────────────────
  if (has_quad) {
    fitted_line <- coefs[var] * x_seq + coefs[quad_name] * x_seq^2
  } else {
    fitted_line <- coefs[var] * x_seq
  }

  # ── Partial residuals ────────────────────────────────────────────────────────
  if (has_quad) {
    partial_res <- resids +
      coefs[var]       * xvar +
      coefs[quad_name] * xvar^2
  } else {
    partial_res <- resids + coefs[var] * xvar
  }

  # ── Confidence band via SE simulation ───────────────────────────────────────
  se_vec <- sar_coef_table[coef_names, "Std. Error"]

  set.seed(seed)
  coef_sims <- mapply(
    function(mu, se) rnorm(nsim, mean = mu, sd = se),
    mu = coefs[coef_names],
    se = se_vec
  )
  # Ensure matrix even when single predictor
  if (!is.matrix(coef_sims)) coef_sims <- matrix(coef_sims, ncol = 1)

  # Design matrix
  if (has_quad) {
    X_mat <- cbind(x_seq, x_seq^2)
  } else {
    X_mat <- matrix(x_seq, ncol = 1)
  }

  sims <- X_mat %*% t(coef_sims)   # n_seq x nsim

  # ── Assemble data frames ─────────────────────────────────────────────────────
  fit_df <- data.frame(
    x   = x_seq,
    fit = fitted_line,
    lwr = apply(sims, 1, quantile, 0.025),
    upr = apply(sims, 1, quantile, 0.975)
  )

  res_df <- data.frame(
    x   = xvar,
    res = partial_res
  )

  # ── Turning point (quadratic only) ───────────────────────────────────────────
  if (has_quad) {
    inflection <- unname(-coefs[var] / (2 * coefs[quad_name]))
    cat("Turning point (", var, "):", round(inflection, 3), "\n")
  }

  # ── x axis label ─────────────────────────────────────────────────────────────
  if (is.null(x_label)) x_label <- var

  # ── Plot ──────────────────────────────────────────────────────────────────────
  p <- ggplot() +
    geom_point(data = res_df,
               aes(x = x, y = res),
               size = point_size, alpha = point_alpha, color = "grey50") +
    geom_ribbon(data = fit_df,
                aes(x = x, ymin = lwr, ymax = upr),
                fill = line_color, alpha = 0.25) +
    geom_line(data = fit_df,
              aes(x = x, y = fit),
              color = line_color, linewidth = 1) +
    geom_hline(yintercept = 0,
               linetype = "dashed", color = "grey40", linewidth = 0.4) +
    labs(x = x_label, y = y_label) +
    theme_classic(base_size = 13) +
    theme(
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 11)
    )

  # Add turning-point marker only for quadratic (only when inside the EH range)
  if (has_quad) {
    if (inflection < max(xvar) & inflection > min(xvar)){
      p <- p +
        geom_vline(xintercept = inflection,
                   linetype = "dotted", color = "black", linewidth = 0.6) +
        annotate("text",
                 x     = -Inf,
                 y     = Inf,
                 label = paste0("atop(R^2 == ", round(R2, 2),
                                ", Turning~point == ", round(inflection, 2), ")"),
                 hjust = -0.1,    # small positive = slight right offset from left edge
                 vjust = 1.5,     # positive = slight down offset from top edge
                 size = 5, color = "black",
                 parse = TRUE)
    }else{
      p <- p +
        annotate("text",
                 x     = -Inf,
                 y     = Inf,
                 label = paste0("atop(R^2 == ", round(R2, 2),")"),
                 hjust = -0.1,
                 vjust = 1.5,
                 size = 5, color = "black",
                 parse = TRUE)
    }
  }

  return(p)
}
