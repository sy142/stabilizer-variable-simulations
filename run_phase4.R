setwd("C:/Users/Salim/OneDrive/Desktop/SVT SEM/Rev1/Codes")
source("phase4_sem_simulation.R")

suppressPackageStartupMessages({
  library(foreach)
  library(doParallel)
})

set.seed(9186)

scenarios <- c("CFA_TypeAB", "CFA_Null", "CFA_Moderator")
K_levels  <- c(5, 10, 15, 20)
n_levels  <- c(100, 200, 300, 400, 500)
mi_levels <- c(0.15, 0.30, 0.45, 0.65)
n_reps    <- 100

conditions <- expand.grid(
  scenario = scenarios,
  K        = K_levels,
  n        = n_levels,
  mi       = mi_levels,
  rep      = 1:n_reps,
  stringsAsFactors = FALSE
)

cat("Phase 4: Full SEM-Estimated SVT Validation\n")
cat("Scenarios:", paste(scenarios, collapse = ", "), "\n")
cat("K:", paste(K_levels, collapse = ", "), "\n")
cat("n:", paste(n_levels, collapse = ", "), "\n")
cat("MI:", paste(mi_levels, collapse = ", "), "\n")
cat("Reps:", n_reps, "\n")
cat("Total simulations:", nrow(conditions), "\n")
cat("Estimated time: 4-8 hours with 14 cores\n\n")

n_cores <- 14
cl <- makeCluster(n_cores)
registerDoParallel(cl)

start_time <- Sys.time()

results <- foreach(
  i = 1:nrow(conditions),
  .combine = function(...) do.call(rbind, Filter(Negate(is.null), list(...))),
  .packages = c("boot", "lavaan"),
  .export = c("generate_sem_data", "run_single_sem_sim",
              "compute_metrics", "bootstrap_test"),
  .errorhandling = "remove"
) %dopar% {
  row <- conditions[i, ]
  run_single_sem_sim(row$scenario, row$K, row$n, row$mi, seed = 60000 + i)
}

stopCluster(cl)

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "hours")
cat("Completed in:", round(as.numeric(elapsed), 2), "hours\n\n")

output_dir <- "C:/Users/Salim/OneDrive/Desktop/SVT SEM/Phase4"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

saveRDS(results, file.path(output_dir, "phase4_full_results.rds"))
write.csv(results, file.path(output_dir, "phase4_full_results.csv"), row.names = FALSE)

cat("Rows:", nrow(results), "\n\n")

valid <- results[results$mechanism != "FAILED", ]

valid$bootstrap_pass    <- as.numeric(valid$bootstrap_pass)
valid$dual_criterion_pass <- as.numeric(valid$dual_criterion_pass)

cat("=== GLOBAL BY SCENARIO ===\n\n")
for (sc in scenarios) {
  sub <- valid[valid$scenario == sc, ]
  if (nrow(sub) == 0) next
  cat(sprintf("%s (n=%d): Boot=%.1f%% Dual=%.1f%% DLog=%.3f SD=%.3f\n",
              sc, nrow(sub),
              mean(sub$bootstrap_pass) * 100,
              mean(sub$dual_criterion_pass) * 100,
              mean(sub$delta_log, na.rm = TRUE),
              sd(sub$delta_log, na.rm = TRUE)))
}

cat("\n=== CFA_TypeAB BY CONDITION ===\n\n")
typeab <- valid[valid$scenario == "CFA_TypeAB", ]
if (nrow(typeab) > 0) {
  tab <- aggregate(
    cbind(bootstrap_pass, dual_criterion_pass, delta_log) ~ K + mi_severity,
    data = typeab,
    FUN  = function(x) round(mean(x, na.rm = TRUE), 3)
  )
  tab <- tab[order(tab$mi_severity, tab$K), ]
  cat(sprintf("%-4s %-6s %10s %10s %10s\n", "K", "MI", "Boot%", "Dual%", "DLog"))
  for (i in 1:nrow(tab)) {
    cat(sprintf("%-4d %-6.2f %9.1f%% %9.1f%% %10.3f\n",
                tab$K[i], tab$mi_severity[i],
                tab$bootstrap_pass[i] * 100,
                tab$dual_criterion_pass[i] * 100,
                tab$delta_log[i]))
  }
}

cat("\n=== MECHANISM (CFA_TypeAB) ===\n\n")
if (nrow(typeab) > 0) {
  mech <- table(typeab$mechanism)
  pct <- round(prop.table(mech) * 100, 1)
  for (m in names(pct)) cat(sprintf("  %-12s: %.1f%%\n", m, pct[m]))
}

cat("\n=== FALSE POSITIVE RATES ===\n\n")
for (sc in c("CFA_Null", "CFA_Moderator")) {
  sub <- valid[valid$scenario == sc, ]
  if (nrow(sub) == 0) next
  fpr_boot <- mean(sub$bootstrap_pass) * 100
  fpr_dual <- mean(sub$dual_criterion_pass) * 100
  cat(sprintf("%s: Boot FPR=%.1f%% Dual FPR=%.1f%%\n", sc, fpr_boot, fpr_dual))
}

cat("\n=== COMPARISON: Phase 4 (CFA) vs Phase 1 (Regression) ===\n\n")
if (nrow(typeab) > 0) {
  null_sub <- valid[valid$scenario == "CFA_Null", ]
  cat(sprintf("Phase 4 (CFA): Power=%.1f%% DLog=%.3f FPR=%.1f%%\n",
              mean(typeab$dual_criterion_pass) * 100,
              mean(typeab$delta_log, na.rm = TRUE),
              mean(null_sub$dual_criterion_pass) * 100))
  cat("Phase 1 (Reg): Power~99%  DLog~2.25  FPR~1.4%\n")
}

saveRDS(list(
  global = aggregate(cbind(bootstrap_pass, dual_criterion_pass, delta_log) ~ scenario,
                     data = valid, FUN = function(x) round(mean(x, na.rm = TRUE), 4)),
  detail = if (nrow(typeab) > 0) aggregate(
    cbind(bootstrap_pass, dual_criterion_pass, delta_log) ~ K + mi_severity,
    data = typeab, FUN = function(x) round(mean(x, na.rm = TRUE), 4)) else NULL
), file.path(output_dir, "phase4_summary.rds"))

cat("Phase 4 CFA-based SVT functions loaded.\n")