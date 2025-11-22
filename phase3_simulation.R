suppressPackageStartupMessages({
  library(parallel)
  library(boot)
  library(dplyr)
})

source("phase1_simulation.R")

run_phase3_boundary <- function(n_cores = 14) {
  
  cat("\n==============================================\n")
  cat("PHASE 3: Boundary Condition Analysis\n")
  cat("==============================================\n\n")
  
  set.seed(9186)
  
  configs <- data.frame(
    config_id = 1:18,
    dimension = c(
      rep("Extreme_K_min", 3),
      rep("Extreme_K_max", 3),
      rep("Extreme_MI", 4),
      rep("Minimal_n", 4),
      rep("Extreme_Combo", 4)
    ),
    K = c(3, 3, 3, 50, 50, 50, 20, 20, 10, 5, 10, 10, 20, 20, 3, 50, 3, 50),
    n = c(100, 200, 500, 100, 200, 500, 200, 200, 200, 200, 30, 30, 30, 30, 30, 30, 500, 500),
    mi_severity = c(0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.10, 0.90, 0.90, 0.90, 0.45, 0.65, 0.45, 0.65, 0.90, 0.45, 0.90, 0.10),
    description = c(
      "K=3, n=100 (minimum groups)",
      "K=3, n=200 (minimum groups)",
      "K=3, n=500 (minimum groups)",
      "K=50, n=100 (maximum groups)",
      "K=50, n=200 (maximum groups)",
      "K=50, n=500 (maximum groups)",
      "MI=0.10 (very weak violation)",
      "MI=0.90 (extreme violation)",
      "K=10, MI=0.90",
      "K=5, MI=0.90",
      "K=10, n=30 (minimal sample)",
      "K=10, n=30, MI=0.65",
      "K=20, n=30 (minimal sample)",
      "K=20, n=30, MI=0.65",
      "Worst case: K=3, n=30, MI=0.90",
      "Many groups, small n: K=50, n=30",
      "Few groups, extreme MI, large n",
      "Best case: K=50, n=500, MI=0.10"
    ),
    stringsAsFactors = FALSE
  )
  
  cat("Configuration Summary:\n")
  print(configs[, c("config_id", "dimension", "K", "n", "mi_severity")])
  cat("\n")
  
  conditions <- do.call(rbind, lapply(1:nrow(configs), function(i) {
    data.frame(
      config_id = configs$config_id[i],
      dimension = configs$dimension[i],
      K = configs$K[i],
      n = configs$n[i],
      mi_severity = configs$mi_severity[i],
      description = configs$description[i],
      rep = 1:200,
      stringsAsFactors = FALSE
    )
  }))
  
  cat("Total simulations:", nrow(conditions), "\n")
  cat("Running on", n_cores, "cores...\n\n")
  
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl, { library(boot) })
  clusterExport(cl, ls(.GlobalEnv), envir = .GlobalEnv)
  
  start_time <- Sys.time()
  
  results <- parLapply(cl, 1:nrow(conditions), function(i) {
    row <- conditions[i, ]
    tryCatch({
      sim_result <- run_single_sim(
        scenario = "TypeAB",
        K = row$K,
        n = row$n,
        mi_severity = row$mi_severity,
        seed = 30000 + i
      )
      
      c(
        config_id = row$config_id,
        dimension = row$dimension,
        description = row$description,
        sim_result
      )
    }, error = function(e) {
      list(
        config_id = row$config_id,
        dimension = row$dimension,
        K = row$K,
        n = row$n,
        mi_severity = row$mi_severity,
        error = as.character(e)
      )
    })
  })
  
  stopCluster(cl)
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "mins")
  
  cat("\nCompleted in:", round(as.numeric(elapsed), 2), "minutes\n")
  
  results_df <- do.call(rbind, lapply(results, function(x) {
    if (is.null(x$error)) {
      as.data.frame(x, stringsAsFactors = FALSE)
    } else {
      data.frame(
        config_id = x$config_id,
        dimension = x$dimension,
        K = x$K,
        n = x$n,
        mi_severity = x$mi_severity,
        error = x$error,
        stringsAsFactors = FALSE
      )
    }
  }))
  
  cat("\n==============================================\n")
  cat("RESULTS SUMMARY BY CONFIGURATION\n")
  cat("==============================================\n\n")
  
  if (!"error" %in% names(results_df) || all(is.na(results_df$error))) {
    results_df$bootstrap_pass <- as.numeric(results_df$bootstrap_pass)
    results_df$delta_log <- as.numeric(results_df$delta_log)
    results_df$K <- as.numeric(results_df$K)
    results_df$n <- as.numeric(results_df$n)
    results_df$mi_severity <- as.numeric(results_df$mi_severity)
    
    summary_by_config <- results_df %>%
      group_by(config_id, dimension, K, n, mi_severity) %>%
      summarise(
        N = n(),
        Power = mean(bootstrap_pass, na.rm = TRUE) * 100,
        Mean_DeltaLog = mean(delta_log, na.rm = TRUE),
        SD_DeltaLog = sd(delta_log, na.rm = TRUE),
        .groups = "drop"
      )
    
    print(summary_by_config, n = 100)
  }
  
  return(results_df)
}

cat("Phase 3 simulation functions loaded successfully.\n")
cat("Ready to run: results_3 <- run_phase3_boundary(n_cores = 14)\n")