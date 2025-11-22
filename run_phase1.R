setwd("C:/Users/Salim/Desktop/Yildiz Teknik YL Tez/Makale/EC Revs/Rev2/JASA/Kodlar RDS grafikler/SimulationCodes")

source("phase1_simulation.R")

suppressPackageStartupMessages({
  library(foreach)
  library(doParallel)
})

cat("\n==============================================\n")
cat("PHASE 1 FULL SIMULATION - SVT\n")
cat("==============================================\n\n")

set.seed(9186)

scenarios <- c("TypeA", "TypeB", "TypeAB", "Null", "Moderator")
K_levels <- c(5, 6, 7, 8, 9, 10, 15, 20)
n_levels <- c(50, 100, 200, 500, 1000)
mi_levels <- c(0.2, 0.3, 0.45, 0.65)

conditions <- expand.grid(
  scenario = scenarios,
  K = K_levels,
  n = n_levels,
  mi = mi_levels,
  rep = 1:1000,
  stringsAsFactors = FALSE
)

cat("Total simulations:", nrow(conditions), "\n")
cat("Estimated time: 18-22 hours with 12 cores\n\n")

n_cores <- 14
cl <- makeCluster(n_cores)
registerDoParallel(cl)

cat("Starting parallel execution (foreach)...\n")
start_time <- Sys.time()

results <- foreach(
  i = 1:nrow(conditions),
  .combine = rbind,
  .packages = c("boot"),
  .export = c("run_single_sim", "generate_data_typeA", "generate_data_typeB",
              "generate_data_typeAB", "generate_data_null", "generate_data_moderator",
              "compute_metrics", "bootstrap_test"),
  .errorhandling = "pass",
  .verbose = FALSE
) %dopar% {
  
  row <- conditions[i, ]
  
  result <- tryCatch({
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
  
  if (is.null(result$error)) {
    as.data.frame(result, stringsAsFactors = FALSE)
  } else {
    data.frame(scenario = result$scenario, K = result$K, n = result$n, 
               mi_severity = result$mi_severity, error = result$error,
               stringsAsFactors = FALSE)
  }
}

stopCluster(cl)

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "hours")
cat("\nCompleted in:", round(as.numeric(elapsed), 2), "hours\n")

output_dir <- "C:/Users/Salim/Desktop/Yildiz Teknik YL Tez/Makale/EC Revs/Rev2/JASA/Kodlar RDS grafikler/SimulationRecords"

saveRDS(results, file.path(output_dir, "phase1_results_full.rds"))
write.csv(results, file.path(output_dir, "phase1_results_full.csv"), row.names = FALSE)

cat("\nResults saved to:\n")
cat("- phase1_results_full.rds\n")
cat("- phase1_results_full.csv\n\n")

cat("Summary:\n")
cat("Total rows:", nrow(results), "\n")
if ("error" %in% names(results)) {
  cat("Errors:", sum(!is.na(results$error)), "\n")
}

cat("\n==============================================\n")
cat("PHASE 1 COMPLETE\n")
cat("==============================================\n")