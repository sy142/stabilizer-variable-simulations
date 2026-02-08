suppressPackageStartupMessages({
  library(lavaan)
})

set.seed(9186)

generate_cfa_groups <- function(n_per_group, n_groups, base_loadings,
                                residual_vars, loading_shifts, intercept_shifts) {
  p <- length(base_loadings)
  all_data <- vector("list", n_groups)

  for (g in 1:n_groups) {
    xi <- rnorm(n_per_group)
    lambdas <- pmax(0.20, base_loadings + loading_shifts[g, ])
    taus <- intercept_shifts[g, ]

    X <- matrix(NA, n_per_group, p)
    for (j in 1:p) {
      X[, j] <- taus[j] + lambdas[j] * xi + rnorm(n_per_group, 0, sqrt(residual_vars[j]))
    }

    df <- as.data.frame(X)
    names(df) <- paste0("X", 1:p)
    df$group <- g
    all_data[[g]] <- df
  }

  do.call(rbind, all_data)
}

run_mi_assessment <- function(data, cfa_model) {
  data$group <- factor(data$group)

  fit_config <- tryCatch(
    cfa(cfa_model, data = data, group = "group", std.lv = TRUE, estimator = "MLR"),
    error = function(e) NULL
  )
  if (is.null(fit_config)) return(NULL)

  fit_metric <- tryCatch(
    cfa(cfa_model, data = data, group = "group", group.equal = "loadings",
        std.lv = TRUE, estimator = "MLR"),
    error = function(e) NULL
  )
  if (is.null(fit_metric)) return(NULL)

  fit_scalar <- tryCatch(
    cfa(cfa_model, data = data, group = "group",
        group.equal = c("loadings", "intercepts"),
        std.lv = TRUE, estimator = "MLR"),
    error = function(e) NULL
  )
  if (is.null(fit_scalar)) return(NULL)

  fm_cfg <- fitMeasures(fit_config, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr"))
  fm_met <- fitMeasures(fit_metric, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr"))
  fm_sca <- fitMeasures(fit_scalar, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr"))

  if (any(is.na(c(fm_cfg, fm_met, fm_sca)))) return(NULL)

  d_met <- c(
    CFI   = abs(as.numeric(fm_cfg["cfi.scaled"]   - fm_met["cfi.scaled"])),
    TLI   = abs(as.numeric(fm_cfg["tli.scaled"]   - fm_met["tli.scaled"])),
    RMSEA = abs(as.numeric(fm_met["rmsea.scaled"] - fm_cfg["rmsea.scaled"])),
    SRMR  = abs(as.numeric(fm_met["srmr"]         - fm_cfg["srmr"]))
  )

  d_sca <- c(
    CFI   = abs(as.numeric(fm_met["cfi.scaled"]   - fm_sca["cfi.scaled"])),
    TLI   = abs(as.numeric(fm_met["tli.scaled"]   - fm_sca["tli.scaled"])),
    RMSEA = abs(as.numeric(fm_sca["rmsea.scaled"] - fm_met["rmsea.scaled"])),
    SRMR  = abs(as.numeric(fm_sca["srmr"]         - fm_met["srmr"]))
  )

  d_max <- pmax(d_met, d_sca)

  is_inv <- (d_met["CFI"] <= 0.010 & d_met["TLI"] <= 0.010 &
             d_met["RMSEA"] <= 0.015 & d_met["SRMR"] <= 0.030 &
             d_sca["CFI"] <= 0.010 & d_sca["TLI"] <= 0.010 &
             d_sca["RMSEA"] <= 0.015 & d_sca["SRMR"] <= 0.030)

  list(
    delta_metric = d_met,
    delta_scalar = d_sca,
    delta_max = d_max,
    is_invariant_chen = as.logical(is_inv),
    fit_config = fm_cfg,
    fit_metric = fm_met,
    fit_scalar = fm_sca
  )
}

compute_weights_both <- function(delta_mat, detected_inv) {

  norm_mat <- cbind(
    CFI   = pmin(delta_mat[, "CFI"]   / 0.010, 1),
    TLI   = pmin(delta_mat[, "TLI"]   / 0.010, 1),
    RMSEA = pmin(delta_mat[, "RMSEA"] / 0.015, 1),
    SRMR  = pmin(delta_mat[, "SRMR"]  / 0.030, 1)
  )

  cv_vals <- apply(norm_mat, 2, function(x) sd(x) / (mean(x) + 0.001))
  vs <- 1 / (1 + cv_vals)

  cor_m <- tryCatch(cor(delta_mat), error = function(e) diag(4))
  mean_r <- apply(abs(cor_m), 1, function(x) mean(x[x < 1]))
  rp <- 1 / (1 + mean_r)

  indices <- c("CFI", "TLI", "RMSEA", "SRMR")
  dp_abs <- dp_max0 <- setNames(numeric(4), indices)

  n_inv <- sum(detected_inv)
  n_noninv <- sum(!detected_inv)

  for (idx in indices) {
    if (n_inv > 0 && n_noninv > 0) {
      inv_m <- mean(norm_mat[detected_inv, idx])
      noninv_m <- mean(norm_mat[!detected_inv, idx])
      dp_abs[idx] <- abs(noninv_m - inv_m)
      dp_max0[idx] <- max(0, noninv_m - inv_m)
    } else {
      dp_abs[idx] <- sd(norm_mat[, idx])
      dp_max0[idx] <- sd(norm_mat[, idx])
    }
  }

  make_weights <- function(dp) {
    raw <- dp * rp * vs
    if (sum(raw) > 0) raw / sum(raw) else rep(0.25, 4)
  }

  w_abs  <- make_weights(dp_abs)
  w_max0 <- make_weights(dp_max0)

  scores_abs  <- as.numeric(norm_mat %*% w_abs)
  scores_max0 <- as.numeric(norm_mat %*% w_max0)

  list(
    w_abs = w_abs, w_max0 = w_max0,
    dp_abs = dp_abs, dp_max0 = dp_max0,
    scores_abs = scores_abs, scores_max0 = scores_max0,
    norm_mat = norm_mat,
    max0_fallback = (sum(dp_max0) == 0)
  )
}

compute_auc <- function(scores, positive_labels) {
  n_pos <- sum(positive_labels)
  n_neg <- sum(!positive_labels)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  ranks <- rank(scores)
  (sum(ranks[positive_labels]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

run_phase0_single <- function(mi_severity, n_per_group, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n_ind <- 8
  n_mod <- 10
  n_noninv <- 5
  n_grp <- 2

  base_load <- seq(0.65, 0.85, length.out = n_ind)
  res_var <- 1 - base_load^2
  cfa_mod <- paste0("F =~ ", paste(paste0("X", 1:n_ind), collapse = " + "))

  true_noninvariant <- c(rep(TRUE, n_noninv), rep(FALSE, n_mod - n_noninv))

  delta_list <- vector("list", n_mod)
  chen_inv <- logical(n_mod)
  fit_details <- vector("list", n_mod)

  for (m in 1:n_mod) {
    l_shift <- matrix(0, n_grp, n_ind)
    i_shift <- matrix(0, n_grp, n_ind)

    if (true_noninvariant[m]) {
      for (g in 2:n_grp) {
        l_shift[g, ] <- rnorm(n_ind, 0, mi_severity * 0.25)
        i_shift[g, ] <- rnorm(n_ind, 0, mi_severity * 0.40)
      }
    }

    mod_data <- generate_cfa_groups(n_per_group, n_grp, base_load,
                                     res_var, l_shift, i_shift)

    mi_res <- run_mi_assessment(mod_data, cfa_mod)

    if (!is.null(mi_res)) {
      delta_list[[m]] <- mi_res$delta_max
      chen_inv[m] <- mi_res$is_invariant_chen
      fit_details[[m]] <- mi_res
    } else {
      delta_list[[m]] <- rep(NA_real_, 4)
      chen_inv[m] <- NA
    }
  }

  valid <- sapply(delta_list, function(x) !any(is.na(x)))
  if (sum(valid) < 4) return(NULL)

  d_mat <- do.call(rbind, delta_list[valid])
  colnames(d_mat) <- c("CFI", "TLI", "RMSEA", "SRMR")
  true_nv <- true_noninvariant[valid]
  chen_v <- chen_inv[valid]

  chen_separates <- (sum(chen_v) > 0 && sum(!chen_v) > 0)

  comp <- compute_weights_both(d_mat, chen_v)

  auc_abs  <- compute_auc(comp$scores_abs, true_nv)
  auc_max0 <- compute_auc(comp$scores_max0, true_nv)

  chen_acc <- mean(chen_v == !true_nv)

  chen_sens <- NA_real_
  chen_spec <- NA_real_
  if (sum(true_nv) > 0)  chen_sens <- mean(!chen_v[true_nv])
  if (sum(!true_nv) > 0) chen_spec <- mean(chen_v[!true_nv])

  pool_sd <- function(x, y) {
    if (length(x) < 2 || length(y) < 2) return(NA_real_)
    vx <- var(x); vy <- var(y)
    if (is.na(vx) || is.na(vy)) return(NA_real_)
    sqrt((vx + vy) / 2)
  }

  s_abs_inv  <- comp$scores_abs[!true_nv]
  s_abs_nv   <- comp$scores_abs[true_nv]
  s_max0_inv <- comp$scores_max0[!true_nv]
  s_max0_nv  <- comp$scores_max0[true_nv]

  psd_abs  <- pool_sd(s_abs_inv, s_abs_nv)
  psd_max0 <- pool_sd(s_max0_inv, s_max0_nv)

  d_abs  <- if (!is.na(psd_abs) && psd_abs > 0)   (mean(s_abs_nv) - mean(s_abs_inv)) / psd_abs     else NA_real_
  d_max0 <- if (!is.na(psd_max0) && psd_max0 > 0) (mean(s_max0_nv) - mean(s_max0_inv)) / psd_max0  else NA_real_

  mean_delta_inv <- colMeans(d_mat[!true_nv, , drop = FALSE])
  mean_delta_nv  <- colMeans(d_mat[true_nv, , drop = FALSE])

  data.frame(
    mi_severity      = mi_severity,
    n_per_group      = n_per_group,
    n_valid          = sum(valid),
    chen_separates   = chen_separates,
    chen_accuracy    = chen_acc,
    chen_sensitivity = chen_sens,
    chen_specificity = chen_spec,
    auc_abs          = auc_abs,
    auc_max0         = auc_max0,
    auc_advantage    = auc_max0 - auc_abs,
    d_abs            = d_abs,
    d_max0           = d_max0,
    mean_inv_abs     = mean(s_abs_inv),
    mean_noninv_abs  = mean(s_abs_nv),
    mean_inv_max0    = mean(s_max0_inv),
    mean_noninv_max0 = mean(s_max0_nv),
    w_cfi_abs        = comp$w_abs["CFI"],
    w_tli_abs        = comp$w_abs["TLI"],
    w_rmsea_abs      = comp$w_abs["RMSEA"],
    w_srmr_abs       = comp$w_abs["SRMR"],
    w_cfi_max0       = comp$w_max0["CFI"],
    w_tli_max0       = comp$w_max0["TLI"],
    w_rmsea_max0     = comp$w_max0["RMSEA"],
    w_srmr_max0      = comp$w_max0["SRMR"],
    dp_cfi_abs       = comp$dp_abs["CFI"],
    dp_tli_abs       = comp$dp_abs["TLI"],
    dp_rmsea_abs     = comp$dp_abs["RMSEA"],
    dp_srmr_abs      = comp$dp_abs["SRMR"],
    dp_cfi_max0      = comp$dp_max0["CFI"],
    dp_tli_max0      = comp$dp_max0["TLI"],
    dp_rmsea_max0    = comp$dp_max0["RMSEA"],
    dp_srmr_max0     = comp$dp_max0["SRMR"],
    max0_fallback    = comp$max0_fallback,
    delta_cfi_inv    = mean_delta_inv["CFI"],
    delta_tli_inv    = mean_delta_inv["TLI"],
    delta_rmsea_inv  = mean_delta_inv["RMSEA"],
    delta_srmr_inv   = mean_delta_inv["SRMR"],
    delta_cfi_nv     = mean_delta_nv["CFI"],
    delta_tli_nv     = mean_delta_nv["TLI"],
    delta_rmsea_nv   = mean_delta_nv["RMSEA"],
    delta_srmr_nv    = mean_delta_nv["SRMR"],
    stringsAsFactors = FALSE
  )
}

cat("Phase 0 simulation functions loaded.\n")
cat("Ready: run_phase0_single(mi_severity, n_per_group, seed)\n")
