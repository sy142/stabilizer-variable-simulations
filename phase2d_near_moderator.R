suppressPackageStartupMessages({
  library(parallel)
  library(boot)
})

source("phase1_simulation.R")

generate_data_near_moderator <- function(n_groups, n_per_group, mi_severity,
                                         interaction_strength = 0.05) {

  group_params <- vector("list", n_groups)
  all_data <- vector("list", n_groups)

  for (k in 1:n_groups) {
    n <- n_per_group

    beta_base <- rnorm(1, -0.20, 0.08)

    X <- rnorm(n, 0, 1)
    Z <- rnorm(n, 0, 1)

    Y_true <- beta_base * X + interaction_strength * X * Z + rnorm(n, 0, 0.50)

    mi_mod <- ifelse(mi_severity > 0.3, 0.20, mi_severity * 0.5)
    Y_measured <- Y_true + rnorm(n, 0, mi_mod)

    fit_baseline <- lm(Y_measured ~ X)
    beta_baseline <- coef(fit_baseline)[2]

    resid_Y <- residuals(lm(Y_measured ~ Z))
    resid_X <- residuals(lm(X ~ Z))

    if (sd(resid_X) < 1e-10) {
      beta_adjusted <- beta_baseline
    } else {
      beta_adjusted <- coef(lm(resid_Y ~ resid_X - 1))[1]
    }

    group_params[[k]] <- list(
      group = k,
      n = n,
      beta_baseline = as.numeric(beta_baseline),
      beta_adjusted = as.numeric(beta_adjusted)
    )

    all_data[[k]] <- data.frame(
      group = k,
      X = X,
      Z = Z,
      Y = Y_measured
    )
  }

  list(
    scenario = paste0("NearMod_", interaction_strength),
    group_params = group_params,
    data = do.call(rbind, all_data)
  )
}

run_single_near_moderator <- function(K, n, mi_severity,
                                       interaction_strength, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  sim_data <- generate_data_near_moderator(K, n, mi_severity,
                                            interaction_strength)

  metrics <- compute_metrics(sim_data$group_params)
  stats   <- bootstrap_test(sim_data$group_params, n_boot = 500)

  data.frame(
    interaction   = interaction_strength,
    K             = K,
    n             = n,
    mi_severity   = mi_severity,
    delta_log     = metrics$delta_log,
    cv_0          = metrics$cv_0,
    cv_1          = metrics$cv_1,
    var_reduction = metrics$var_reduction_pct,
    orient_share  = metrics$orientation_share,
    dist_red_pct  = metrics$dist_reduction_pct,
    p_value       = stats$p_value,
    z_stat        = stats$z_stat,
    se            = stats$se,
    ci_lower      = stats$ci_lower,
    ci_upper      = stats$ci_upper,
    bootstrap_pass      = stats$bootstrap_pass,
    binom_var_pass      = stats$binom_variance_pass,
    binom_align_pass    = stats$binom_alignment_pass,
    mechanism           = stats$mechanism_detected,
    binom_p_variance    = stats$binom_p_variance,
    binom_p_alignment   = stats$binom_p_alignment,
    stringsAsFactors    = FALSE
  )
}

cat("Phase 2D functions loaded.\n")
cat("Ready: run_single_near_moderator(K, n, mi_severity, interaction_strength, seed)\n")
