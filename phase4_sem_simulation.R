suppressPackageStartupMessages({
  library(parallel)
  library(boot)
  library(lavaan)
})

source("phase1_simulation.R")

generate_sem_data <- function(n_groups, n_per_group, mi_severity,
                              true_beta = -0.09, n_ind_x = 6,
                              scenario = "TypeAB") {
  
  base_lx <- seq(0.57, 0.82, length.out = n_ind_x)
  res_x <- 1 - base_lx^2
  group_data_list <- vector("list", n_groups)
  
  for (k in 1:n_groups) {
    n <- n_per_group
    U_k <- rnorm(1, 0, 1)
    dl_x <- rnorm(n_ind_x, mi_severity * 0.25 * U_k, mi_severity * 0.12)
    di_x <- rnorm(n_ind_x, mi_severity * 0.35 * U_k, mi_severity * 0.18)
    lx_k <- pmax(0.25, base_lx + dl_x)
    tx_k <- di_x
    xi <- rnorm(n, 0, 1)
    Xmat <- matrix(NA, n, n_ind_x)
    for (j in 1:n_ind_x) {
      Xmat[, j] <- tx_k[j] + lx_k[j] * xi + rnorm(n, 0, sqrt(res_x[j]))
    }
    if (scenario == "TypeAB") {
      rho_k <- -0.41 + rnorm(1, 0, 0.06)
      rho_k <- max(-0.55, min(-0.25, rho_k))
      Z <- rho_k * xi + rnorm(n, 0, sqrt(1 - rho_k^2))
      gamma_k <- rnorm(1, -0.17, 0.055)
      Y <- true_beta * xi + gamma_k * Z + rnorm(n, 0, 0.90)
    } else if (scenario == "Moderator") {
      Z <- rnorm(n, 0, 1)
      Y <- true_beta * xi + 0.15 * xi * Z + rnorm(n, 0, 0.90)
    } else {
      Z <- rnorm(n, 0, 1)
      beta_k <- rnorm(1, true_beta, 0.04)
      Y <- beta_k * xi + rnorm(n, 0, 0.90)
    }
    df <- data.frame(group = rep(k, n))
    for (j in 1:n_ind_x) df[[paste0("x", j)]] <- Xmat[, j]
    df$Y <- Y
    df$Z <- Z
    group_data_list[[k]] <- df
  }
  do.call(rbind, group_data_list)
}

run_single_sem_sim <- function(scenario, K, n, mi_severity, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  gen_scenario <- sub("^CFA_", "", scenario)
  
  full_data <- generate_sem_data(K, n, mi_severity, scenario = gen_scenario)
  
  x_names <- paste0("x", 1:6)
  mod0 <- paste0("XI =~ ", paste(x_names, collapse = " + "), "\nY ~ beta*XI")
  mod1 <- paste0("XI =~ ", paste(x_names, collapse = " + "), "\nY ~ beta*XI + gamma*Z")
  
  b0_vec <- rep(NA_real_, K)
  b1_vec <- rep(NA_real_, K)
  
  for (g in 1:K) {
    gd <- full_data[full_data$group == g, ]
    f0 <- tryCatch(sem(mod0, data = gd, std.lv = TRUE, estimator = "ML"),
                   error = function(e) NULL)
    f1 <- tryCatch(sem(mod1, data = gd, std.lv = TRUE, estimator = "ML"),
                   error = function(e) NULL)
    if (is.null(f0) || is.null(f1)) next
    if (!lavInspect(f0, "converged") || !lavInspect(f1, "converged")) next
    s0 <- standardizedSolution(f0)
    s1 <- standardizedSolution(f1)
    b0_vec[g] <- s0[s0$label == "beta", "est.std"]
    b1_vec[g] <- s1[s1$label == "beta", "est.std"]
  }
  
  ok <- !is.na(b0_vec) & !is.na(b1_vec)
  if (sum(ok) < 3) {
    return(data.frame(
      scenario = scenario, K = K, n = n, mi_severity = mi_severity,
      n_converged = sum(ok), converge_rate = sum(ok) / K,
      delta_log = NA, bootstrap_pass = NA, dual_criterion_pass = NA,
      mechanism = "FAILED", stringsAsFactors = FALSE))
  }
  
  gp <- vector("list", sum(ok))
  idx <- 1
  for (g in which(ok)) {
    gp[[idx]] <- list(group = g, n = as.integer(n),
                      beta_baseline = b0_vec[g],
                      beta_adjusted = b1_vec[g])
    idx <- idx + 1
  }
  
  metrics <- compute_metrics(gp)
  stats   <- bootstrap_test(gp, n_boot = 500)
  dual <- as.numeric(stats$bootstrap_pass &
                       (stats$binom_variance_pass | stats$binom_alignment_pass))
  
  data.frame(
    scenario = scenario, K = K, n = n, mi_severity = mi_severity,
    n_converged = sum(ok), converge_rate = round(sum(ok) / K, 3),
    mean_0 = metrics$mean_0, mean_1 = metrics$mean_1,
    sd_0 = metrics$sd_0, sd_1 = metrics$sd_1,
    cv_0 = metrics$cv_0, cv_1 = metrics$cv_1,
    delta_log = metrics$delta_log,
    pct_reduction = metrics$pct_reduction,
    var_reduction_pct = metrics$var_reduction_pct,
    orientation_share = metrics$orientation_share,
    dist_reduction_pct = metrics$dist_reduction_pct,
    p_value = stats$p_value, z_stat = stats$z_stat,
    bootstrap_pass = as.numeric(stats$bootstrap_pass),
    binom_var_pass = as.numeric(stats$binom_variance_pass),
    binom_align_pass = as.numeric(stats$binom_alignment_pass),
    dual_criterion_pass = dual,
    mechanism = stats$mechanism_detected,
    stringsAsFactors = FALSE)
}

cat("Phase 4.\n")