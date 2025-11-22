setwd("C:/Users/Salim/Desktop/Yildiz Teknik YL Tez/Makale/EC Revs/Rev2/JASA/Kodlar RDS grafikler/SimulationCodes")

library(dplyr)
source("phase3_simulation.R")

output_dir <- "../SimulationRecords"

cat("\n==============================================\n")
cat("PHASE 3: BOUNDARY CONDITION ANALYSIS\n")
cat("==============================================\n\n")

cat("Testing SVT at parameter-space extremes:\n")
cat("- Extreme K: 3 and 50 groups\n")
cat("- Extreme MI: 0.10 and 0.90\n")
cat("- Minimal n: 30 per group\n")
cat("- Worst-case combinations\n\n")

results_3 <- run_phase3_boundary(n_cores = 14)

saveRDS(results_3, file.path(output_dir, "phase3_boundary_conditions.rds"))
write.csv(results_3, file.path(output_dir, "phase3_boundary_conditions.csv"), row.names = FALSE)

cat("\n??? Phase 3 saved\n")

if (!"error" %in% names(results_3) || all(is.na(results_3$error))) {
  
  cat("\n==============================================\n")
  cat("DIAGNOSTIC SUMMARY\n")
  cat("==============================================\n\n")
  
  results_3$bootstrap_pass <- as.numeric(results_3$bootstrap_pass)
  results_3$delta_log <- as.numeric(results_3$delta_log)
  
  dimension_summary <- results_3 %>%
    group_by(dimension) %>%
    summarise(
      Configs = length(unique(config_id)),
      Mean_Power = mean(bootstrap_pass, na.rm = TRUE) * 100,
      Mean_Effect = mean(delta_log, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("Summary by Dimension:\n")
  print(dimension_summary)
  
  problem_configs <- results_3 %>%
    group_by(config_id, K, n, mi_severity) %>%
    summarise(Power = mean(bootstrap_pass, na.rm = TRUE) * 100, .groups = "drop") %>%
    filter(Power < 50) %>%
    arrange(Power)
  
  if (nrow(problem_configs) > 0) {
    cat("\n??????  Low-power configurations (< 50%):\n")
    print(problem_configs)
  } else {
    cat("\n All configurations achieve > 50% power\n")
  }
}

cat("\n==============================================\n")
cat("PHASE 3 COMPLETE\n")
cat("==============================================\n")
cat("Total simulations:", nrow(results_3), "\n")
cat("Files saved to:", output_dir, "\n")