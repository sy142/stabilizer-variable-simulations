setwd("C:/Users/Salim/Desktop/Yildiz Teknik YL Tez/Makale/EC Revs/Rev2/Mathematics/Codes")
source("phase0_simulation.R")

suppressPackageStartupMessages({
  library(foreach)
  library(doParallel)
})

cat("\n==============================================\n")
cat("PHASE 0: ADAPTIVE MI WEIGHT VALIDATION\n")
cat("abs() vs max(0,...) Discriminant Power\n")
cat("==============================================\n\n")

set.seed(9186)

mi_levels <- c(0.10, 0.20, 0.30, 0.45, 0.65, 0.90)
n_levels  <- c(100, 200, 500)
n_reps    <- 200

conditions <- expand.grid(
  mi  = mi_levels,
  n   = n_levels,
  rep = 1:n_reps,
  stringsAsFactors = FALSE
)

cat("Conditions: ", length(mi_levels), " MI x ", length(n_levels), " n x ", n_reps, " reps\n", sep = "")
cat("Total simulations:", nrow(conditions), "\n")
cat("Each replication: 10 moderators x 3 CFA fits = 30 lavaan fits\n")
cat("Total CFA fits: ~", nrow(conditions) * 30, "\n")
cat("Estimated time: 4-8 hours with 14 cores\n\n")

n_cores <- 14
cl <- makeCluster(n_cores)
registerDoParallel(cl)

cat("Starting parallel execution...\n")
start_time <- Sys.time()

results <- foreach(
  i = 1:nrow(conditions),
  .combine = function(...) {
    args <- list(...)
    valid <- Filter(Negate(is.null), args)
    if (length(valid) == 0) return(NULL)
    do.call(rbind, valid)
  },
  .packages = c("lavaan"),
  .export = c("generate_cfa_groups", "run_mi_assessment",
              "compute_weights_both", "compute_auc",
              "run_phase0_single"),
  .errorhandling = "remove",
  .verbose = FALSE
) %dopar% {
  row <- conditions[i, ]
  tryCatch(
    run_phase0_single(
      mi_severity = row$mi,
      n_per_group = row$n,
      seed = 9186 + i
    ),
    error = function(e) NULL
  )
}

stopCluster(cl)

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "hours")
cat("\nCompleted in:", round(as.numeric(elapsed), 2), "hours\n")

output_dir <- "C:/Users/Salim/Desktop/Yildiz Teknik YL Tez/Makale/EC Revs/Rev2/Mathematics/Ciktilar"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

saveRDS(results, file.path(output_dir, "phase0_results.rds"))
write.csv(results, file.path(output_dir, "phase0_results.csv"), row.names = FALSE)

cat("\nResults saved to:\n")
cat("  phase0_results.rds\n")
cat("  phase0_results.csv\n")
cat("Total rows:", nrow(results), "\n\n")

cat("==============================================\n")
cat("SUMMARY: abs() vs max(0,...) COMPARISON\n")
cat("==============================================\n\n")

agg <- aggregate(
  cbind(chen_accuracy, chen_sensitivity, chen_specificity,
        auc_abs, auc_max0, auc_advantage, d_abs, d_max0, max0_fallback) ~
    mi_severity + n_per_group,
  data = results,
  FUN = function(x) round(mean(x, na.rm = TRUE), 4)
)
agg <- agg[order(agg$mi_severity, agg$n_per_group), ]

cat("--- Chen Threshold Performance ---\n")
for (i in 1:nrow(agg)) {
  cat(sprintf("MI=%.2f n=%3d | Acc=%.3f Sens=%.3f Spec=%.3f\n",
              agg$mi_severity[i], agg$n_per_group[i],
              agg$chen_accuracy[i], agg$chen_sensitivity[i], agg$chen_specificity[i]))
}

cat("\n--- AUC Comparison (higher = better) ---\n")
cat(sprintf("%-12s %-6s %8s %8s %10s\n", "MI", "n", "abs()", "max(0,)", "advantage"))
for (i in 1:nrow(agg)) {
  cat(sprintf("MI=%.2f      n=%3d %8.4f %8.4f %+10.4f\n",
              agg$mi_severity[i], agg$n_per_group[i],
              agg$auc_abs[i], agg$auc_max0[i], agg$auc_advantage[i]))
}

cat("\n--- Cohen's d Comparison (higher = better separation) ---\n")
cat(sprintf("%-12s %-6s %8s %8s\n", "MI", "n", "abs()", "max(0,)"))
for (i in 1:nrow(agg)) {
  cat(sprintf("MI=%.2f      n=%3d %8.3f %8.3f\n",
              agg$mi_severity[i], agg$n_per_group[i],
              agg$d_abs[i], agg$d_max0[i]))
}

cat("\n--- max(0,...) Fallback Rate ---\n")
fb <- aggregate(max0_fallback ~ mi_severity, data = results,
                FUN = function(x) round(mean(x, na.rm = TRUE) * 100, 1))
for (i in 1:nrow(fb)) {
  cat(sprintf("MI=%.2f: %.1f%% fallback to equal weights\n",
              fb$mi_severity[i], fb$max0_fallback[i]))
}

cat("\n--- Mean Adaptive Weights ---\n")
w_agg <- aggregate(
  cbind(w_cfi_abs, w_tli_abs, w_rmsea_abs, w_srmr_abs,
        w_cfi_max0, w_tli_max0, w_rmsea_max0, w_srmr_max0) ~ mi_severity,
  data = results,
  FUN = function(x) round(mean(x, na.rm = TRUE), 4)
)
cat("\nabs() weights:\n")
cat(sprintf("%-8s %8s %8s %8s %8s\n", "MI", "CFI", "TLI", "RMSEA", "SRMR"))
for (i in 1:nrow(w_agg)) {
  cat(sprintf("MI=%.2f  %8.4f %8.4f %8.4f %8.4f\n",
              w_agg$mi_severity[i],
              w_agg$w_cfi_abs[i], w_agg$w_tli_abs[i],
              w_agg$w_rmsea_abs[i], w_agg$w_srmr_abs[i]))
}
cat("\nmax(0,...) weights:\n")
cat(sprintf("%-8s %8s %8s %8s %8s\n", "MI", "CFI", "TLI", "RMSEA", "SRMR"))
for (i in 1:nrow(w_agg)) {
  cat(sprintf("MI=%.2f  %8.4f %8.4f %8.4f %8.4f\n",
              w_agg$mi_severity[i],
              w_agg$w_cfi_max0[i], w_agg$w_tli_max0[i],
              w_agg$w_rmsea_max0[i], w_agg$w_srmr_max0[i]))
}

cat("\n--- Mean Delta Values by True Status ---\n")
d_agg <- aggregate(
  cbind(delta_cfi_inv, delta_tli_inv, delta_rmsea_inv, delta_srmr_inv,
        delta_cfi_nv, delta_tli_nv, delta_rmsea_nv, delta_srmr_nv) ~ mi_severity,
  data = results,
  FUN = function(x) round(mean(x, na.rm = TRUE), 4)
)
cat(sprintf("%-8s | %-36s | %-36s\n", "MI", "Invariant (CFI TLI RMSEA SRMR)", "Non-invariant (CFI TLI RMSEA SRMR)"))
for (i in 1:nrow(d_agg)) {
  cat(sprintf("MI=%.2f  | %.4f %.4f %.4f %.4f | %.4f %.4f %.4f %.4f\n",
              d_agg$mi_severity[i],
              d_agg$delta_cfi_inv[i], d_agg$delta_tli_inv[i],
              d_agg$delta_rmsea_inv[i], d_agg$delta_srmr_inv[i],
              d_agg$delta_cfi_nv[i], d_agg$delta_tli_nv[i],
              d_agg$delta_rmsea_nv[i], d_agg$delta_srmr_nv[i]))
}

summary_agg <- agg
saveRDS(summary_agg, file.path(output_dir, "phase0_summary.rds"))
write.csv(summary_agg, file.path(output_dir, "phase0_summary.csv"), row.names = FALSE)

cat("\nSummary saved to phase0_summary.rds / .csv\n")

cat("\n==============================================\n")
cat("PHASE 0 COMPLETE\n")
cat("==============================================\n")

