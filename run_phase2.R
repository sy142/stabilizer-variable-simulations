setwd("C:/Users/Salim/Desktop/Yildiz Teknik YL Tez/Makale/EC Revs/Rev2/JASA/Kodlar RDS grafikler/SimulationCodes")

source("phase2_simulation.R")

output_dir <- "../SimulationRecords"

cat("\n==============================================\n")
cat("PHASE 2: SENSITIVITY ANALYSIS\n")
cat("==============================================\n\n")

results_2A <- run_phase2A_bootstrap(n_cores = 14)
saveRDS(results_2A, file.path(output_dir, "phase2A_bootstrap_convergence.rds"))
write.csv(results_2A, file.path(output_dir, "phase2A_bootstrap_convergence.csv"), row.names = FALSE)
cat("n Phase 2A saved\n")

results_2B <- run_phase2B_noise(n_cores = 14)
saveRDS(results_2B, file.path(output_dir, "phase2B_noise_trajectory.rds"))
write.csv(results_2B, file.path(output_dir, "phase2B_noise_trajectory.csv"), row.names = FALSE)
cat("n Phase 2B saved\n")

results_2C <- run_phase2C_mi(n_cores = 14)
saveRDS(results_2C, file.path(output_dir, "phase2C_mi_trajectory.rds"))
write.csv(results_2C, file.path(output_dir, "phase2C_mi_trajectory.csv"), row.names = FALSE)
cat("n Phase 2C saved\n")

cat("\n==============================================\n")
cat("PHASE 2 COMPLETE\n")
cat("==============================================\n")
cat("Total simulations:", nrow(results_2A) + nrow(results_2B) + nrow(results_2C), "\n")