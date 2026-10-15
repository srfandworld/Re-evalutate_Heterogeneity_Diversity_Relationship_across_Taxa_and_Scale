# =============================================================================
# fun_sar_partial_plot.R - partial-residual plot + turning-point bootstrap
#
# Source   : EmbeddingBio/analysis/fun_sar_partial_plot.R, lines 45-end (verbatim).
#            Lines 1-44 (old sar_modeling_* using log(richness/effort)) were
#            dropped; model fitting now lives in 05_HDR_models/01_fit_HDR_models.R.
# Requires : per-scale tables FIA_table_<s>, BBS_table_<s>, mammal_table_<s> in
#            the calling environment (looked up with get()).
# =============================================================================

library(dplyr)
library(tidyr)
library(glue)
library(ggplot2)
library(performance)

sar_partial_plot <- function(model,        # spautolm object
                             var,          # predictor name e.g. "EH_400m_t"
                             scale,
                             taxa,
                             quadratic  = FALSE,  # TRUE if var has I(var^2) term
                             nsim       = 500,
                             n_seq      = 200,
                             seed       = 123,
                             point_size = 0.4,
                             point_alpha= 0.5,
                             line_color = "steelblue",
                             x_label    = NULL,
                             y_label    = "Partial residual (log MSA)") {
  
  # scale, data
  if(taxa == "bird"){
    data <- get(glue("BBS_table_{scale}")) %>% drop_na()
    
    sar_coef_table <- summary(model)$Coef
    
    y <- model$Y
    y_hat <- model$fit$fitted.values
    R2 <- cor(y, y_hat)^2
    
  }else if(taxa == "mammal"){
    data <- get(glue("mammal_table_{scale}")) %>% 
      drop_na() %>% 
      group_by(Camera_Trap_Array) %>% 
      summarise(across(c(richness_Estimator, abundance, MSA, contains("_t")), 
                       ~ mean(.x, na.rm = TRUE)),
                .groups = "drop")
    
    sar_coef_table <- summary(model)$Coef
    
    y <- model$Y
    y_hat <- model$fit$fitted.values
    R2 <- cor(y, y_hat)^2
    
  }else{
    data <- get(glue('FIA_table_{scale}')) %>% 
      drop_na()
    
    sar_coef_table <- summary(model)$coef
    R2 <- unname(r2(model)[[2]])
  }
  
  
  # ── Extract model components ────────────────────────────────────────────────
  coefs  <- coef(model)
  resids <- residuals(model)
  
  
  # ── Predictor sequence ──────────────────────────────────────────────────────
  x_seq <- seq(min(data[[var]], na.rm = TRUE),
               max(data[[var]], na.rm = TRUE),
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
      coefs[var]       * data[[var]] +
      coefs[quad_name] * data[[var]]^2
  } else {
    partial_res <- resids + coefs[var] * data[[var]]
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
    x   = data[[var]],
    res = partial_res
  )
  
  # ── Inflection point (quadratic only) ────────────────────────────────────────
  if (has_quad) {
    inflection <- -coefs[var] / (2 * coefs[quad_name])
    cat("Inflection point (", var, "):", round(inflection, 3), "\n")
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
  
  # Add inflection marker only for quadratic
  if (has_quad) {
    if (inflection < max(data[[var]]) & inflection > min(data[[var]])){
      p <- p +
        
        geom_vline(xintercept = inflection,
                   linetype = "dotted", color = "black", linewidth = 0.6) +
        annotate("text",
                 x     = -Inf,
                 y     = Inf,
                 label = paste0("atop(R^2 == ", round(R2, 2), 
                                ", Inflection == ", round(inflection, 2), ")"),
                 hjust = -0.1,    # small positive = slight right offset from left edge
                 vjust = 1.5,     # positive = slight down offset from top edge
                 size = 5, color = "black",
                 parse = TRUE)
    }else{
        p <-p + 
          annotate("text",
                   x     = -Inf,
                   y     = Inf,
                   label = paste0("atop(R^2 == ", round(R2, 2),")"),
                   hjust = -0.1,    # small positive = slight right offset from left edge
                   vjust = 1.5,     # positive = slight down offset from top edge
                   size = 5, color = "black",
                   parse = TRUE)
      }
  }
  
  return(p)
}





#Bootstraps the entire SAR model to get empirical distribution of z:
inflection_boot <- function(model,taxa, nsim = 500, seed = 123) {
  set.seed(seed)
  
  
  
  # scale, data
  if(taxa == "bird"){
    scale <- str_extract(colnames(model$X)[2], "\\d+(km|m)")
    EH_var <- glue("EH_{scale}")
    
    data <- get(glue("BBS_table_{scale}")) %>% drop_na()
    
    #sar model
    sar_table  <- summary(model)$Coef
  }else if(taxa == "mammal"){
    scale <- str_extract(colnames(model$X)[2], "\\d+(km|m)")
    EH_var <- glue("EH_{scale}")
    
    data <- get(glue("mammal_table_{scale}")) %>% 
      drop_na() %>% 
      group_by(Camera_Trap_Array) %>% 
      summarise(across(where(is.numeric), 
                       ~ mean(.x, na.rm = TRUE)),
                .groups = "drop")
    
    #sar model
    sar_table  <- summary(model)$Coef
  }else{
    scale <- str_extract(colnames(model$model)[2], "\\d+(km|m)")
    EH_var <- glue("EH_{scale}")
    
    data <- get(glue('FIA_table_{scale}')) %>% 
      drop_na()
    
    #lm model
    sar_table  <- summary(model)$coef
    
  }
  
  
    # Simulate coefficients from their sampling distribution
  coef_names <- rownames(sar_table)
  means      <- sar_table[, "Estimate"]
  ses        <- sar_table[, "Std. Error"]
  
  # Draw nsim coefficient sets
  z_sims <- replicate(nsim, {
    b1 <- rnorm(1, means[2], ses[2])
    b2 <- rnorm(1, means[3], ses[3])
    -b1 / (2 * b2)
  })
  
  z_obs <- -means[2] / (2 * means[3])
  
  #orignial scale mean and sd
  mu  <- mean(data[[EH_var]], na.rm = TRUE)
  sdv <- sd(data[[EH_var]], na.rm = TRUE)
  
  data.frame(
    z      = mu + sdv*z_obs,
    se_z   = mu + sdv*sd(z_sims),
    lwr_95 = mu + sdv*quantile(z_sims, 0.025),
    upr_95 = mu + sdv*quantile(z_sims, 0.975)
  )
}



