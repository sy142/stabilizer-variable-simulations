suppressPackageStartupMessages({
  library(parallel)
  library(boot)
})

set.seed(9186)

generate_data_typeA <- function(n_groups, n_per_group, mi_severity, 
                                true_beta = -0.25, baseline_beta_mean = 0.25,
                                baseline_beta_sd = 0.08, sign_inconsistency_prob = 0.20) {
  
  group_params <- vector("list", n_groups)
  all_data <- vector("list", n_groups)
  
  for (k in 1:n_groups) {
    n <- n_per_group
    
    U <- rnorm(n, 0, 1)
    X <- 0.70 * U + rnorm(n, 0, sqrt(1 - 0.70^2))
    Z <- -0.70 * U + rnorm(n, 0, sqrt(1 - 0.70^2))
    
    Y_true <- true_beta * X + 0.10 * U + rnorm(n, 0, 0.35)
    
    b_k <- rnorm(1, 0, mi_severity * 0.6)
    
    Y_measured <- Y_true + b_k * Z + rnorm(n, 0, mi_severity * 0.15)
    
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
      beta_adjusted = as.numeric(beta_adjusted),
      b_k = b_k
    )
    
    all_data[[k]] <- data.frame(
      group = k,
      X = X,
      Z = Z,
      Y = Y_measured
    )
  }
  
  list(
    scenario = "TypeA",
    group_params = group_params,
    data = do.call(rbind, all_data)
  )
}

generate_data_typeB <- function(n_groups, n_per_group, mi_severity,
                                true_beta = -0.15, z_x_correlation = 0.40) {
  
  group_params <- vector("list", n_groups)
  all_data <- vector("list", n_groups)
  
  gamma_const <- abs(0.45 * mi_severity + 0.15)
  
  for (k in 1:n_groups) {
    n <- n_per_group
    
    U <- rnorm(n, 0, 1)
    rho <- z_x_correlation
    X <- rho * U + rnorm(n, 0, sqrt(1 - rho^2))
    Z <- rho * U + rnorm(n, 0, sqrt(1 - rho^2))
    
    beta_k <- rnorm(1, mean = true_beta, sd = 0.05)
    gamma_k <- gamma_const
    
    Y_true <- beta_k * X + rnorm(n, 0, 0.35)
    eps_meas <- rnorm(n, 0, mi_severity * 0.04 + 0.085)
    Y_measured <- Y_true + gamma_k * Z + eps_meas
    
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
      beta_adjusted = as.numeric(beta_adjusted),
      gamma_k = gamma_k
    )
    
    all_data[[k]] <- data.frame(
      group = k,
      X = X,
      Z = Z,
      Y = Y_measured
    )
  }
  
  list(
    scenario = "TypeB",
    group_params = group_params,
    data = do.call(rbind, all_data)
  )
}

generate_data_typeAB <- function(n_groups, n_per_group, mi_severity,
                                 true_beta = -0.20, z_x_correlation = 0.55) {
  
  group_params <- vector("list", n_groups)
  all_data <- vector("list", n_groups)
  
  for (k in 1:n_groups) {
    n <- n_per_group
    
    U <- rnorm(n, 0, 1)
    rho <- z_x_correlation
    X <- rho * U + rnorm(n, 0, sqrt(1 - rho^2))
    Z <- rho * U + rnorm(n, 0, sqrt(1 - rho^2))
    
    b_k <- rnorm(1, 0, mi_severity * 0.75)
    
    gamma_base <- 0.55 * mi_severity + 0.30
    if (runif(1) < 0.15) {
      gamma_k <- -gamma_base + rnorm(1, 0, 0.05)
    } else {
      gamma_k <-  gamma_base + rnorm(1, 0, 0.05)
    }
    
    lambda_k <- rbeta(1, 3, 7)
    mix_coef <- lambda_k * b_k + (1 - lambda_k) * gamma_k
    
    Y_true <- true_beta * X + 0.1 * U + rnorm(n, 0, 0.35)
    eps_meas <- rnorm(n, 0, mi_severity * 0.06 + 0.10)
    Y_measured <- Y_true + mix_coef * Z + eps_meas
    
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
      group = k, n = n,
      beta_baseline = as.numeric(beta_baseline),
      beta_adjusted = as.numeric(beta_adjusted),
      mix_coef = mix_coef, lambda_k = lambda_k,
      b_k = b_k, gamma_k = gamma_k
    )
    
    all_data[[k]] <- data.frame(group = k, X = X, Z = Z, Y = Y_measured)
  }
  
  list(scenario = "TypeAB",
       group_params = group_params,
       data = do.call(rbind, all_data))
}




generate_data_null <- function(n_groups, n_per_group, mi_severity) {
  
  group_params <- vector("list", n_groups)
  all_data <- vector("list", n_groups)
  
  for (k in 1:n_groups) {
    n <- n_per_group
    
    beta_true <- rnorm(1, -0.20, 0.10)
    
    X <- rnorm(n, 0, 1)
    Z <- rnorm(n, 0, 1)
    
    Y_true <- beta_true * X + rnorm(n, 0, 0.7)
    
    mi_null <- ifelse(mi_severity > 0.3, 0.20, mi_severity * 0.5)
    Y_measured <- Y_true + rnorm(n, 0, mi_null)
    
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
    scenario = "Null",
    group_params = group_params,
    data = do.call(rbind, all_data)
  )
}

generate_data_moderator <- function(n_groups, n_per_group, mi_severity) {
  
  group_params <- vector("list", n_groups)
  all_data <- vector("list", n_groups)
  
  for (k in 1:n_groups) {
    n <- n_per_group
    
    beta_base <- rnorm(1, -0.20, 0.08)
    
    X <- rnorm(n, 0, 1)
    Z <- rnorm(n, 0, 1)
    
    Y_true <- beta_base * X + 0.25 * X * Z + rnorm(n, 0, 0.5)
    
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
    scenario = "Moderator",
    group_params = group_params,
    data = do.call(rbind, all_data)
  )
}

compute_metrics <- function(group_params) {
  
  n_groups <- length(group_params)
  beta_0 <- sapply(group_params, function(x) x$beta_baseline)
  beta_1 <- sapply(group_params, function(x) x$beta_adjusted)
  weights <- sapply(group_params, function(x) x$n)
  
  mean_0 <- weighted.mean(beta_0, weights)
  mean_1 <- weighted.mean(beta_1, weights)
  
  var_0 <- sum(weights * (beta_0 - mean_0)^2) / sum(weights)
  var_1 <- sum(weights * (beta_1 - mean_1)^2) / sum(weights)
  
  sd_0 <- sqrt(var_0)
  sd_1 <- sqrt(var_1)
  
  epsilon_mu <- 0.0001
  mean_0_protected <- sign(mean_0) * max(abs(mean_0), epsilon_mu)
  mean_1_protected <- sign(mean_1) * max(abs(mean_1), epsilon_mu)
  
  cv_0 <- (sd_0 / abs(mean_0_protected)) * 100
  cv_1 <- (sd_1 / abs(mean_1_protected)) * 100
  
  eps_log <- 1e-8
  delta_log <- log(cv_0 + eps_log) - log(cv_1 + eps_log)
  pct_reduction <- 100 * (1 - exp(-delta_log))
  
  var_reduction_pct <- ifelse(var_0 > 0, ((var_0 - var_1) / var_0) * 100, 0)
  
  delta_sigma <- sd_0 - sd_1
  delta_mu <- abs(mean_1) - abs(mean_0)
  
  eps <- 1e-8
  dlog_sigma <- log(sd_1 + eps) - log(sd_0 + eps)
  dlog_mu <- log(abs(mean_1) + eps) - log(abs(mean_0) + eps)
  orientation_share <- abs(dlog_mu) / (abs(dlog_mu) + abs(dlog_sigma) + eps)
  
  dist_0 <- abs(beta_0 - mean_0)
  dist_1 <- abs(beta_1 - mean_1)
  dist_reduction_count <- sum(dist_1 < dist_0)
  dist_reduction_pct <- dist_reduction_count / n_groups
  
  sign_changes <- sum(sign(beta_0) != sign(beta_1))
  
  list(
    mean_0 = mean_0,
    mean_1 = mean_1,
    sd_0 = sd_0,
    sd_1 = sd_1,
    cv_0 = cv_0,
    cv_1 = cv_1,
    delta_log = delta_log,
    pct_reduction = pct_reduction,
    var_reduction_pct = var_reduction_pct,
    orientation_share = orientation_share,
    dist_reduction_pct = dist_reduction_pct,
    sign_changes = sign_changes,
    n_groups = n_groups
  )
}

bootstrap_test <- function(group_params, n_boot = 200,
                           min_delta = 0.05,
                           min_dlog_mu = 0.02) {
  metrics <- compute_metrics(group_params)
  delta_log_obs <- metrics$delta_log
  
  n_groups <- length(group_params)
  beta_0_vec <- sapply(group_params, function(x) x$beta_baseline)
  beta_1_vec <- sapply(group_params, function(x) x$beta_adjusted)
  
  delta_beta   <- beta_1_vec - beta_0_vec
  delta_beta_c <- delta_beta - mean(delta_beta)
  
  boot_delta_log <- numeric(n_boot)
  boot_dlog_mu   <- numeric(n_boot)
  
  for (b in 1:n_boot) {
    sgn <- sample(c(-1, 1), n_groups, replace = TRUE)
    b0 <- beta_0_vec
    b1 <- beta_0_vec + sgn * delta_beta_c
    
    gp <- group_params
    for (i in 1:n_groups) {
      gp[[i]]$beta_baseline <- b0[i]
      gp[[i]]$beta_adjusted <- b1[i]
    }
    
    m <- compute_metrics(gp)
    boot_delta_log[b] <- m$delta_log
    boot_dlog_mu[b] <- log(abs(m$mean_1) + 1e-8) - log(abs(m$mean_0) + 1e-8)
  }
  
  p_value      <- mean(boot_delta_log >= delta_log_obs)
  se_delta_log <- sd(boot_delta_log, na.rm = TRUE)
  z_stat       <- delta_log_obs / (se_delta_log + 1e-12)
  ci           <- quantile(boot_delta_log, c(0.025, 0.975), na.rm = TRUE)
  
  n_closer <- sum(abs(beta_1_vec - metrics$mean_1) < abs(beta_0_vec - metrics$mean_0))
  binom_p_variance <- binom.test(n_closer, n_groups, p = 0.5,
                                 alternative = "greater")$p.value
  
  dlog_mu_obs <- log(abs(metrics$mean_1) + 1e-8) - log(abs(metrics$mean_0) + 1e-8)
  
  if (abs(metrics$mean_0) < 1e-6 && abs(metrics$mean_1) < 1e-6) {
    binom_p_alignment <- 1.0
  } else {
    binom_p_alignment <- mean(boot_dlog_mu >= dlog_mu_obs)
  }
  
  bootstrap_pass      <- (delta_log_obs > min_delta && p_value < 0.05)
  binom_variance_pass <- (binom_p_variance < 0.05)
  binom_alignment_pass<- (dlog_mu_obs > min_dlog_mu && binom_p_alignment < 0.05)
  
  mechanism <- "None"
  if (bootstrap_pass && binom_variance_pass && binom_alignment_pass) {
    mechanism <- "TypeAB"
  } else if (bootstrap_pass && binom_variance_pass) {
    mechanism <- "TypeA"
  } else if (bootstrap_pass && binom_alignment_pass) {
    mechanism <- "TypeB"
  } else if (bootstrap_pass) {
    mechanism <- "Marginal"
  }
  
  list(
    delta_log = delta_log_obs,
    p_value = p_value,
    z_stat = z_stat,
    se = se_delta_log,
    ci_lower = ci[1],
    ci_upper = ci[2],
    binom_p_variance = binom_p_variance,
    binom_p_alignment = binom_p_alignment,
    dlog_mu = dlog_mu_obs,
    bootstrap_pass = bootstrap_pass,
    binom_variance_pass = binom_variance_pass,
    binom_alignment_pass = binom_alignment_pass,
    mechanism_detected = mechanism
  )
}

run_single_sim <- function(scenario, K, n, mi_severity, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  sim_data <- switch(
    scenario,
    "TypeA" = generate_data_typeA(K, n, mi_severity),
    "TypeB" = generate_data_typeB(K, n, mi_severity),
    "TypeAB" = generate_data_typeAB(K, n, mi_severity),
    "Null" = generate_data_null(K, n, mi_severity),
    "Moderator" = generate_data_moderator(K, n, mi_severity)
  )
  
  metrics <- compute_metrics(sim_data$group_params)
  stats <- bootstrap_test(sim_data$group_params, n_boot = 200)
  
  c(
    scenario = scenario,
    K = K,
    n = n,
    mi_severity = mi_severity,
    as.list(metrics),
    as.list(stats)
  )
}

run_phase1_parallel <- function(n_cores = 14, n_reps = 1000) {
  
  scenarios <- c("TypeA", "TypeB", "TypeAB", "Null", "Moderator")
  K_levels <- c(5, 6, 7, 8, 9, 10, 15, 20)
  n_levels <- c(50, 100, 200, 500, 1000)
  mi_levels <- c(0.2, 0.3, 0.45, 0.65)
  
  conditions <- expand.grid(
    scenario = scenarios,
    K = K_levels,
    n = n_levels,
    mi = mi_levels,
    rep = 1:n_reps,
    stringsAsFactors = FALSE
  )
  
  cat("Total simulations:", nrow(conditions), "\n")
  cat("Estimated time: 18-22 hours with", n_cores, "cores\n\n")
  
  cl <- makeCluster(n_cores)
  
  clusterEvalQ(cl, {
    library(boot)
  })
  
  clusterExport(cl, ls(.GlobalEnv), envir = .GlobalEnv)
  
  cat("Starting parallel execution...\n")
  start_time <- Sys.time()
  
  results <- parLapply(cl, 1:nrow(conditions), function(i) {
    row <- conditions[i, ]
    tryCatch({
      run_single_sim(
        scenario = row$scenario,
        K = row$K,
        n = row$n,
        mi_severity = row$mi,
        seed = 9186 + i
      )
    }, error = function(e) {
      list(scenario = row$scenario, K = row$K, n = row$n, 
           mi_severity = row$mi, error = as.character(e))
    })
  })
  
  stopCluster(cl)
  
  end_time <- Sys.time()
  cat("\nCompleted in:", difftime(end_time, start_time, units = "hours"), "hours\n")
  
  results_df <- do.call(rbind, lapply(results, function(x) {
    if (is.null(x$error)) {
      as.data.frame(x, stringsAsFactors = FALSE)
    } else {
      data.frame(scenario = x$scenario, K = x$K, n = x$n, 
                 mi_severity = x$mi_severity, error = x$error)
    }
  }))
  
  return(results_df)
}

cat("Phase 1 simulation functions loaded successfully.\n")
cat("Ready to run: results <- run_phase1_parallel(n_cores = 14, n_reps = 1000)\n")