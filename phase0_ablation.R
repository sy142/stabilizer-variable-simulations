suppressPackageStartupMessages({
  library(lavaan)
  library(foreach)
  library(doParallel)
})

source("phase0_simulation.R")

compute_ablation_scores <- function(norm_mat, detected_inv, delta_mat) {

  cv_vals <- apply(norm_mat, 2, function(x) sd(x) / (mean(x) + 0.001))
  vs <- 1 / (1 + cv_vals)

  cor_m <- tryCatch(cor(delta_mat), error = function(e) diag(4))
  mean_r <- apply(abs(cor_m), 1, function(x) mean(x[x < 1]))
  rp <- 1 / (1 + mean_r)

  indices <- c("CFI", "TLI", "RMSEA", "SRMR")
  dp <- setNames(numeric(4), indices)

  n_inv <- sum(detected_inv)
  n_noninv <- sum(!detected_inv)

  for (idx in indices) {
    if (n_inv > 0 && n_noninv > 0) {
      inv_m <- mean(norm_mat[detected_inv, idx])
      noninv_m <- mean(norm_mat[!detected_inv, idx])
      dp[idx] <- max(0, noninv_m - inv_m)
    } else {
      dp[idx] <- sd(norm_mat[, idx])
    }
  }

  make_w <- function(raw) {
    if (sum(raw) > 0) raw / sum(raw) else rep(0.25, 4)
  }

  w_full    <- make_w(rp * vs * dp)
  w_dp      <- make_w(dp)
  w_rp_dp   <- make_w(rp * dp)
  w_vs_dp   <- make_w(vs * dp)
  w_equal   <- rep(0.25, 4)

  list(
    scores_full    = as.numeric(norm_mat %*% w_full),
    scores_dp      = as.numeric(norm_mat %*% w_dp),
    scores_rp_dp   = as.numeric(norm_mat %*% w_rp_dp),
    scores_vs_dp   = as.numeric(norm_mat %*% w_vs_dp),
    scores_equal   = as.numeric(norm_mat %*% w_equal),
    w_full = w_full, w_dp = w_dp,
    w_rp_dp = w_rp_dp, w_vs_dp = w_vs_dp, w_equal = w_equal
  )
}

run_ablation_single <- function(mi_severity, n_per_group, seed = NULL) {

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

  norm_mat <- cbind(
    CFI   = pmin(d_mat[, "CFI"]   / 0.010, 1),
    TLI   = pmin(d_mat[, "TLI"]   / 0.010, 1),
    RMSEA = pmin(d_mat[, "RMSEA"] / 0.015, 1),
    SRMR  = pmin(d_mat[, "SRMR"]  / 0.030, 1)
  )

  abl <- compute_ablation_scores(norm_mat, chen_v, d_mat)

  auc_fn <- function(scores, labels) {
    n_pos <- sum(labels)
    n_neg <- sum(!labels)
    if (n_pos == 0 || n_neg == 0) return(NA_real_)
    ranks <- rank(scores)
    (sum(ranks[labels]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
  }

  d_fn <- function(scores, labels) {
    s_inv <- scores[!labels]
    s_nv  <- scores[labels]
    if (length(s_inv) < 2 || length(s_nv) < 2) return(NA_real_)
    vx <- var(s_inv); vy <- var(s_nv)
    if (is.na(vx) || is.na(vy)) return(NA_real_)
    psd <- sqrt((vx + vy) / 2)
    if (psd <= 0) return(NA_real_)
    (mean(s_nv) - mean(s_inv)) / psd
  }

  vnames <- c("full", "dp", "rp_dp", "vs_dp", "equal")
  slist  <- list(abl$scores_full, abl$scores_dp, abl$scores_rp_dp,
                 abl$scores_vs_dp, abl$scores_equal)

  row <- data.frame(
    mi_severity = mi_severity,
    n_per_group = n_per_group,
    n_valid     = sum(valid),
    stringsAsFactors = FALSE
  )

  for (v in seq_along(vnames)) {
    row[[paste0("auc_", vnames[v])]]  <- auc_fn(slist[[v]], true_nv)
    row[[paste0("d_", vnames[v])]]    <- d_fn(slist[[v]], true_nv)
  }

  row$chen_accuracy <- mean(chen_v == !true_nv)

  return(row)
}

cat("Ablation functions loaded.\n")
cat("Ready: run_ablation_single(mi_severity, n_per_group, seed)\n")
