setwd("C:/Users/Salim/Desktop/Yildiz Teknik YL Tez/Makale/EC Revs/Rev2/Mathematics/Codes")
source("phase2d_near_moderator.R")

suppressPackageStartupMessages({
  library(foreach)
  library(doParallel)
})

set.seed(9186)

interaction_levels <- c(0.00, 0.02, 0.05, 0.10, 0.15, 0.25)

conditions <- expand.grid(
  interaction = interaction_levels,
  K           = c(5, 10, 15, 20),
  n           = c(50, 100, 200),
  mi          = c(0.15, 0.30, 0.45),
  rep         = 1:500,
  stringsAsFactors = FALSE
)

cat("Phase 2D: Near-Moderator Robustness Analysis\n")
cat("Interaction levels:", paste(interaction_levels, collapse = ", "), "\n")
cat("Fixed: K=10, n=200, MI=0.45\n")
cat("Reps per level: 500\n")
cat("Total simulations:", nrow(conditions), "\n\n")

n_cores <- 14
cl <- makeCluster(n_cores)
registerDoParallel(cl)

start_time <- Sys.time()

results <- foreach(
  i = 1:nrow(conditions),
  .combine = function(...) {
    args <- list(...)
    valid <- Filter(Negate(is.null), args)
    if (length(valid) == 0) return(NULL)
    do.call(rbind, valid)
  },
  .packages = c("boot"),
  .export = c("generate_data_near_moderator", "run_single_near_moderator",
              "compute_metrics", "bootstrap_test"),
  .errorhandling = "remove",
  .verbose = FALSE
) %dopar% {
  row <- conditions[i, ]
  tryCatch(
    run_single_near_moderator(
      K                    = row$K,
      n                    = row$n,
      mi_severity          = row$mi,
      interaction_strength = row$interaction,
      seed                 = 40000 + i
    ),
    error = function(e) NULL
  )
}

stopCluster(cl)

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")
cat("Completed in:", round(as.numeric(elapsed), 2), "minutes\n\n")

output_dir <- "C:/Users/Salim/Desktop/Yildiz Teknik YL Tez/Makale/EC Revs/Rev2/Mathematics/Ciktilar"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

saveRDS(results, file.path(output_dir, "phase2d_near_moderator_results.rds"))
write.csv(results, file.path(output_dir, "phase2d_near_moderator_results.csv"), row.names = FALSE)

cat("Rows:", nrow(results), "\n\n")

results$bootstrap_pass     <- as.numeric(results$bootstrap_pass)
results$binom_var_pass     <- as.numeric(results$binom_var_pass)
results$binom_align_pass   <- as.numeric(results$binom_align_pass)

dual_pass <- as.numeric(results$bootstrap_pass == 1 &
                        (results$binom_var_pass == 1 | results$binom_align_pass == 1))
results$dual_criterion_pass <- dual_pass

agg <- aggregate(
  cbind(bootstrap_pass, dual_criterion_pass, binom_var_pass, binom_align_pass,
        delta_log, var_reduction, dist_red_pct) ~ interaction,
  data = results,
  FUN  = function(x) round(mean(x, na.rm = TRUE), 4)
)

cat("False Positive Rates by Interaction Strength\n")
cat(sprintf("%-12s %12s %12s %12s %12s %10s %10s\n",
            "Interaction", "Bootstrap", "Dual-Crit", "Var-Binom", "Align-Binom",
            "Mean DLog", "Var Red%"))

for (i in 1:nrow(agg)) {
  cat(sprintf("%-12.2f %11.1f%% %11.1f%% %11.1f%% %11.1f%% %10.4f %9.2f%%\n",
              agg$interaction[i],
              agg$bootstrap_pass[i] * 100,
              agg$dual_criterion_pass[i] * 100,
              agg$binom_var_pass[i] * 100,
              agg$binom_align_pass[i] * 100,
              agg$delta_log[i],
              agg$var_reduction[i]))
}

cat("\nMechanism Classification by Interaction Strength\n")
mech_tab <- table(results$interaction, results$mechanism)
print(mech_tab)

mech_pct <- round(prop.table(mech_tab, margin = 1) * 100, 1)
cat("\nMechanism Percentages\n")
print(mech_pct)

cat("\nInterpretation Guide\n")
cat("interaction = 0.00 : Pure null (no interaction, no stabilizer) -> baseline FPR\n")
cat("interaction = 0.02 : Negligible interaction -> should behave like null\n")
cat("interaction = 0.05 : Small interaction -> safe zone boundary\n")
cat("interaction = 0.10 : Moderate interaction -> practical threshold\n")
cat("interaction = 0.15 : Notable interaction -> transition zone\n")
cat("interaction = 0.25 : Full moderation (same as Phase 1 Moderator scenario)\n")

cat("\nSafe Zone Assessment\n")
for (i in 1:nrow(agg)) {
  fpr <- agg$dual_criterion_pass[i] * 100
  status <- ""
  if (fpr <= 5.0) status <- "SAFE (FPR within nominal alpha)"
  else if (fpr <= 10.0) status <- "CAUTION (FPR mildly elevated)"
  else status <- "UNSAFE (FPR exceeds acceptable bounds)"
  cat(sprintf("  interaction = %.2f : FPR = %.1f%% -> %s\n",
              agg$interaction[i], fpr, status))
}

saveRDS(agg, file.path(output_dir, "phase2d_summary.rds"))
write.csv(agg, file.path(output_dir, "phase2d_summary.csv"), row.names = FALSE)

cat("\nPhase 2D complete. Results saved.\n")
