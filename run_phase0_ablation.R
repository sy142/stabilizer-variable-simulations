setwd("C:/Users/Salim/Desktop/Yildiz Teknik YL Tez/Makale/EC Revs/Rev2/Mathematics/Codes")
source("phase0_ablation.R")

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

cat("Total simulations:", nrow(conditions), "\n")
cat("Same seeds as Phase 0 -- reproducing identical data, recomputing weights\n\n")

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
  .packages = c("lavaan"),
  .export = c("generate_cfa_groups", "run_mi_assessment",
              "compute_ablation_scores", "run_ablation_single"),
  .errorhandling = "remove",
  .verbose = FALSE
) %dopar% {
  row <- conditions[i, ]
  tryCatch(
    run_ablation_single(
      mi_severity = row$mi,
      n_per_group = row$n,
      seed        = 9186 + i
    ),
    error = function(e) NULL
  )
}

stopCluster(cl)

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "hours")
cat("Completed in:", round(as.numeric(elapsed), 2), "hours\n")

output_dir <- "C:/Users/Salim/Desktop/Yildiz Teknik YL Tez/Makale/EC Revs/Rev2/Mathematics/Ciktilar"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

saveRDS(results, file.path(output_dir, "phase0_ablation_results.rds"))
write.csv(results, file.path(output_dir, "phase0_ablation_results.csv"), row.names = FALSE)

cat("Rows:", nrow(results), "\n\n")

vnames <- c("full", "dp", "rp_dp", "vs_dp", "equal")
vlabels <- c("RP x VS x DP", "DP only", "RP x DP", "VS x DP", "Equal")

auc_agg <- aggregate(
  cbind(auc_full, auc_dp, auc_rp_dp, auc_vs_dp, auc_equal) ~
    mi_severity + n_per_group,
  data = results,
  FUN  = function(x) round(mean(x, na.rm = TRUE), 4)
)
auc_agg <- auc_agg[order(auc_agg$mi_severity, auc_agg$n_per_group), ]

d_agg <- aggregate(
  cbind(d_full, d_dp, d_rp_dp, d_vs_dp, d_equal) ~
    mi_severity + n_per_group,
  data = results,
  FUN  = function(x) round(mean(x, na.rm = TRUE), 3)
)
d_agg <- d_agg[order(d_agg$mi_severity, d_agg$n_per_group), ]

cat("AUC by MI Severity and Sample Size\n")
cat(sprintf("%-6s %-5s %10s %10s %10s %10s %10s\n",
            "MI", "n", vlabels[1], vlabels[2], vlabels[3], vlabels[4], vlabels[5]))
for (i in 1:nrow(auc_agg)) {
  cat(sprintf("%-6.2f %-5d %10.4f %10.4f %10.4f %10.4f %10.4f\n",
              auc_agg$mi_severity[i], auc_agg$n_per_group[i],
              auc_agg$auc_full[i], auc_agg$auc_dp[i],
              auc_agg$auc_rp_dp[i], auc_agg$auc_vs_dp[i],
              auc_agg$auc_equal[i]))
}

cat("\nCohen's d by MI Severity and Sample Size\n")
cat(sprintf("%-6s %-5s %10s %10s %10s %10s %10s\n",
            "MI", "n", vlabels[1], vlabels[2], vlabels[3], vlabels[4], vlabels[5]))
for (i in 1:nrow(d_agg)) {
  cat(sprintf("%-6.2f %-5d %10.3f %10.3f %10.3f %10.3f %10.3f\n",
              d_agg$mi_severity[i], d_agg$n_per_group[i],
              d_agg$d_full[i], d_agg$d_dp[i],
              d_agg$d_rp_dp[i], d_agg$d_vs_dp[i],
              d_agg$d_equal[i]))
}

cat("\nOverall Mean AUC Across All Conditions\n")
for (v in seq_along(vnames)) {
  col <- paste0("auc_", vnames[v])
  m <- round(mean(results[[col]], na.rm = TRUE), 4)
  cat(sprintf("  %-12s: %.4f\n", vlabels[v], m))
}

cat("\nOverall Mean Cohen's d Across All Conditions\n")
for (v in seq_along(vnames)) {
  col <- paste0("d_", vnames[v])
  m <- round(mean(results[[col]], na.rm = TRUE), 3)
  cat(sprintf("  %-12s: %.3f\n", vlabels[v], m))
}

cat("\nAUC Advantage of Full Model Over Each Reduced Variant\n")
for (v in 2:length(vnames)) {
  col <- paste0("auc_", vnames[v])
  diff <- results$auc_full - results[[col]]
  m <- round(mean(diff, na.rm = TRUE), 4)
  win  <- sum(diff > 0.001, na.rm = TRUE)
  tie  <- sum(abs(diff) <= 0.001, na.rm = TRUE)
  loss <- sum(diff < -0.001, na.rm = TRUE)
  cat(sprintf("  vs %-12s: mean diff = %+.4f  (win/tie/loss = %d/%d/%d)\n",
              vlabels[v], m, win, tie, loss))
}

saveRDS(list(auc = auc_agg, d = d_agg), file.path(output_dir, "phase0_ablation_summary.rds"))
write.csv(auc_agg, file.path(output_dir, "phase0_ablation_auc.csv"), row.names = FALSE)
write.csv(d_agg, file.path(output_dir, "phase0_ablation_d.csv"), row.names = FALSE)

cat("\nAblation complete. Results saved.\n")
