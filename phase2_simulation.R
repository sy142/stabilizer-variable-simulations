suppressPackageStartupMessages({
  library(parallel)
  library(boot)
})

source("phase1_simulation.R")

generate_data_typeAB_custom <- function(n_groups, n_per_group, mi_severity,
                                        true_beta = -0.20, z_x_correlation = 0.55,
                                        noise_level = NULL) {
  
  group_params <- vector("list", n_groups)
  all_data <- vector("list", n_groups)
  
  base_noise <- ifelse(is.null(noise_level), 0.35, noise_level)
  
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
    
    Y_true <- true_beta * X + 0.1 * U + rnorm(n, 0, base_noise)
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

run_single_sim_custom <- function(scenario, K, n, mi_severity, 
                                  n_boot = 200, seed = NULL,
                                  noise_level = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  if (scenario == "TypeAB" && !is.null(noise_level)) {
    sim_data <- generate_data_typeAB_custom(K, n, mi_severity, noise_level = noise_level)
  } else {
    sim_data <- switch(
      scenario,
      "TypeA" = generate_data_typeA(K, n, mi_severity),
      "TypeB" = generate_data_typeB(K, n, mi_severity),
      "TypeAB" = generate_data_typeAB(K, n, mi_severity),
      "Null" = generate_data_null(K, n, mi_severity),
      "Moderator" = generate_data_moderator(K, n, mi_severity)
    )
  }
  
  metrics <- compute_metrics(sim_data$group_params)
  stats <- bootstrap_test(sim_data$group_params, n_boot = n_boot)
  
  c(
    scenario = scenario,
    K = K,
    n = n,
    mi_severity = mi_severity,
    n_boot = n_boot,
    noise_level = ifelse(is.null(noise_level), NA, noise_level),
    as.list(metrics),
    as.list(stats)
  )
}

run_phase2A_bootstrap <- function(n_cores = 14) {
  
  cat("\n==============================================\n")
  cat("PHASE 2A: Bootstrap Convergence Analysis\n")
  cat("==============================================\n\n")
  
  set.seed(9186)
  
  scenarios <- c("TypeAB", "Null")
  bootstrap_levels <- c(500, 1000, 2000)
  
  conditions <- expand.grid(
    scenario = scenarios,
    K = 20,
    n = 200,
    mi = 0.45,
    n_boot = bootstrap_levels,
    rep = 1:500,
    stringsAsFactors = FALSE
  )
  
  cat("Total simulations:", nrow(conditions), "\n")
  cat("Running on", n_cores, "cores...\n\n")
  
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl, { library(boot) })
  clusterExport(cl, ls(.GlobalEnv), envir = .GlobalEnv)
  
  start_time <- Sys.time()
  
  results <- parLapply(cl, 1:nrow(conditions), function(i) {
    row <- conditions[i, ]
    tryCatch({
      run_single_sim_custom(
        scenario = row$scenario,
        K = row$K,
        n = row$n,
        mi_severity = row$mi,
        n_boot = row$n_boot,
        seed = 9186 + i
      )
    }, error = function(e) {
      list(scenario = row$scenario, error = as.character(e))
    })
  })
  
  stopCluster(cl)
  
  end_time <- Sys.time()
  cat("Completed in:", difftime(end_time, start_time, units = "mins"), "minutes\n")
  
  results_df <- do.call(rbind, lapply(results, function(x) {
    if (is.null(x$error)) {
      as.data.frame(x, stringsAsFactors = FALSE)
    } else {
      data.frame(scenario = x$scenario, error = x$error)
    }
  }))
  
  return(results_df)
}

run_phase2B_noise <- function(n_cores = 14) {
  
  cat("\n==============================================\n")
  cat("PHASE 2B: Noise Trajectory Analysis\n")
  cat("==============================================\n\n")
  
  set.seed(9186)
  
  noise_levels <- seq(0.20, 0.70, by = 0.05)
  
  conditions <- expand.grid(
    scenario = "TypeAB",
    K = 20,
    n = 200,
    mi = 0.45,
    noise = noise_levels,
    rep = 1:300,
    stringsAsFactors = FALSE
  )
  
  cat("Total simulations:", nrow(conditions), "\n")
  cat("Noise levels:", length(noise_levels), "(", min(noise_levels), "to", max(noise_levels), ")\n")
  cat("Running on", n_cores, "cores...\n\n")
  
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl, { library(boot) })
  clusterExport(cl, ls(.GlobalEnv), envir = .GlobalEnv)
  
  start_time <- Sys.time()
  
  results <- parLapply(cl, 1:nrow(conditions), function(i) {
    row <- conditions[i, ]
    tryCatch({
      run_single_sim_custom(
        scenario = row$scenario,
        K = row$K,
        n = row$n,
        mi_severity = row$mi,
        n_boot = 1000,
        noise_level = row$noise,
        seed = 10000 + i
      )
    }, error = function(e) {
      list(scenario = row$scenario, error = as.character(e))
    })
  })
  
  stopCluster(cl)
  
  end_time <- Sys.time()
  cat("Completed in:", difftime(end_time, start_time, units = "mins"), "minutes\n")
  
  results_df <- do.call(rbind, lapply(results, function(x) {
    if (is.null(x$error)) {
      as.data.frame(x, stringsAsFactors = FALSE)
    } else {
      data.frame(scenario = x$scenario, error = x$error)
    }
  }))
  
  return(results_df)
}

run_phase2C_mi <- function(n_cores = 14) {
  
  cat("\n==============================================\n")
  cat("PHASE 2C: MI Trajectory Analysis\n")
  cat("==============================================\n\n")
  
  set.seed(9186)
  
  mi_levels <- seq(0.15, 0.70, by = 0.05)
  
  conditions <- expand.grid(
    scenario = "TypeAB",
    K = 20,
    n = 200,
    mi = mi_levels,
    rep = 1:300,
    stringsAsFactors = FALSE
  )
  
  cat("Total simulations:", nrow(conditions), "\n")
  cat("MI levels:", length(mi_levels), "(", min(mi_levels), "to", max(mi_levels), ")\n")
  cat("Running on", n_cores, "cores...\n\n")
  
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl, { library(boot) })
  clusterExport(cl, ls(.GlobalEnv), envir = .GlobalEnv)
  
  start_time <- Sys.time()
  
  results <- parLapply(cl, 1:nrow(conditions), function(i) {
    row <- conditions[i, ]
    tryCatch({
      run_single_sim_custom(
        scenario = row$scenario,
        K = row$K,
        n = row$n,
        mi_severity = row$mi,
        n_boot = 1000,
        noise_level = 0.35,
        seed = 20000 + i
      )
    }, error = function(e) {
      list(scenario = row$scenario, error = as.character(e))
    })
  })
  
  stopCluster(cl)
  
  end_time <- Sys.time()
  cat("Completed in:", difftime(end_time, start_time, units = "mins"), "minutes\n")
  
  results_df <- do.call(rbind, lapply(results, function(x) {
    if (is.null(x$error)) {
      as.data.frame(x, stringsAsFactors = FALSE)
    } else {
      data.frame(scenario = x$scenario, error = x$error)
    }
  }))
  
  return(results_df)
}

cat("Phase 2 simulation functions loaded successfully.\n")
cat("Ready to run:\n")
cat("  results_2A <- run_phase2A_bootstrap(n_cores = 14)\n")
cat("  results_2B <- run_phase2B_noise(n_cores = 14)\n")
cat("  results_2C <- run_phase2C_mi(n_cores = 14)\n")