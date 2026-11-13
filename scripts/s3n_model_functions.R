# ==============================================================================
# Shared S3N Backend Loading and Domain Fitting / Cross-Validation Engine
# ==============================================================================
# Unified modeling engine shared by run_s3n_determinand.R and 02_fit_s3n_models.R.
# ==============================================================================

load_s3n_backend <- function(base_dir = ".") {
  cand_libs <- c(
    Sys.getenv("R_LIBS_USER"),
    file.path(Sys.getenv("HOME"), "R/x86_64-redhat-linux-gnu-library/4.6"),
    file.path(Sys.getenv("HOME"), "R/x86_64-redhat-linux-gnu-library/4.5")
  )
  for (cl in cand_libs) {
    if (nzchar(cl) && dir.exists(cl) && !(normalizePath(cl, mustWork = FALSE) %in% normalizePath(.libPaths(), mustWork = FALSE))) {
      .libPaths(c(cl, .libPaths()))
    }
  }
  if (!requireNamespace("S3N", quietly = TRUE)) {
    stop("Required R package 'S3N' could not be found in .libPaths() or R_LIBS_USER.")
  }
  suppressPackageStartupMessages(library(S3N))

  so_path <- rel_path(base_dir, "src", "s3n_fast_neighbors.so")
  if (!is.loaded("C_compute_s3n_neighbors")) {
    if (!file.exists(so_path)) {
      stop(sprintf("Compiled C++ neighbor library not found at '%s'. Run: R CMD SHLIB -o %s %s",
                   so_path, so_path, rel_path(base_dir, "src", "s3n_fast_neighbors.cpp")))
    }
    dyn.load(so_path)
  }
  invisible(TRUE)
}

parse_n_omp_arg <- function(args = character(0), default = 1L) {
  thr_arg <- grep("^--(threads|n-omp|omp)=", args, value = TRUE)
  if (length(thr_arg) > 0) {
    val <- suppressWarnings(as.integer(sub("^--(threads|n-omp|omp)=", "", thr_arg[1])))
    if (is.finite(val) && val >= 1L) return(val)
  }
  as.integer(default)
}

fit_s3n_domain <- function(df,
                           domain_id,
                           domain_label,
                           trans_obj,
                           var_specs = NULL,
                           m_neighbors = 15L,
                           n_boot = 50L,
                           n_omp = 1L,
                           use_meas_error = FALSE,
                           net_meta = NULL,
                           obs_col = "obs_val",
                           legacy_nitrate_mode = FALSE,
                           out_root = NULL) {
  obs_idx <- which(df$has_obs == 1L)
  n_obs   <- length(obs_idx)
  if (n_obs < 15L) {
    cat(sprintf("  Skipping separate fit for %s (only %d observed reaches).\n", domain_label, n_obs))
    return(NULL)
  }

  if (isTRUE(legacy_nitrate_mode)) {
    cat(sprintf("\n----------------------------------------------------------------------\n"))
    cat(sprintf("  Domain: %s (%d reaches, %d observed)\n", domain_label, nrow(df), n_obs))
    cat(sprintf("----------------------------------------------------------------------\n"))
  } else {
    cat(sprintf("\n--- Fitting %s (%d reaches, %d observed) ---\n", domain_label, nrow(df), n_obs))
  }

  obs_idx  <- obs_idx[order(df$CatchID[obs_idx], df$obs_up_dist_km[obs_idx])]
  pred_idx <- which(df$has_obs == 0L)
  df_ord   <- df[c(obs_idx, pred_idx), ]
  q_tot    <- nrow(df_ord)
  m_use    <- min(as.integer(m_neighbors), max(2L, n_obs - 2L))

  raw_obs_vals <- df_ord[[obs_col]][1:n_obs]
  y_vec_obs <- if ("obs_y" %in% names(df_ord) && !isTRUE(legacy_nitrate_mode)) {
    as.numeric(df_ord$obs_y[1:n_obs])
  } else {
    as.numeric(trans_obj$fwd(raw_obs_vals))
  }
  y_obs <- as.matrix(y_vec_obs)

  # Heteroskedastic measurement error check (Approach A)
  v_meas <- numeric(n_obs)
  has_me <- FALSE
  if ("meas_uncertainty" %in% names(df_ord)) {
    u_raw <- df_ord$meas_uncertainty[1:n_obs]
    if (any(!is.na(u_raw) & u_raw > 0)) {
      v_meas <- calc_response_meas_variance(raw_obs_vals, u_raw, trans_obj)
      has_gam <- ("meas_method" %in% names(df_ord)) && isTRUE(any(tolower(df_ord$meas_method[1:n_obs]) == "gam", na.rm = TRUE))
      has_meta_gam <- !is.null(net_meta) && !is.null(net_meta$meas_method) && identical(tolower(net_meta$meas_method), "gam")
      has_me <- isTRUE((mean(v_meas) > 1e-8) && (isTRUE(use_meas_error) || has_gam || has_meta_gam))
    }
  }

  if (has_me) {
    v_bar   <- mean(v_meas)
    mean_se <- mean(sqrt(v_meas))
    cat("  [Approach A: Measurement Error Kriging] Active: accounting for observation uncertainty\n")
    cat(sprintf("    Mean measurement SE: %.4f (%s scale), Mean error variance: %.6f\n",
                mean_se, trans_obj$short_label, v_bar))
  }

  # Integer topology lookup (avoids as.character string hashing on 90k+ reaches)
  par_idx  <- match(df_ord$parent_dn, df_ord$RivID, nomatch = 0L)
  parent_0 <- ifelse(df_ord$parent_dn == 0 | par_idx == 0L, -1L, par_idx - 1L)
  up_dist  <- ifelse(df_ord$has_obs == 1L, df_ord$obs_up_dist_km, df_ord$up_dist_mid_km)

  X_all <- build_design_matrix(df_ord, include_peat = (domain_id == "england"), var_specs = var_specs)
  cat(sprintf("  Fixed-effect covariates (%d): %s\n", ncol(X_all), paste(colnames(X_all), collapse = ", ")))
  X_obs <- X_all[1:n_obs, , drop = FALSE]
  vy    <- max(1e-4, var(y_vec_obs))

  t0 <- proc.time()
  tt <- .Call("C_build_tree_tin_tout", as.integer(parent_0))
  build_neighbors <- function(model_int) {
    .Call("C_compute_s3n_neighbors",
          as.integer(n_obs), as.integer(m_use), as.integer(parent_0),
          as.integer(df_ord$CatchID), as.integer(tt[[1]]), as.integer(tt[[2]]),
          as.double(df_ord$Dist2Mth_km), as.double(up_dist), as.double(df_ord$log_AFV),
          as.double(df_ord$mid_X), as.double(df_ord$mid_Y), as.integer(model_int))
  }
  nb_tu <- build_neighbors(0L)
  nb_td <- build_neighbors(1L)
  t_nb  <- (proc.time() - t0)[3]

  t1 <- proc.time()
  fit_tu <- BRISC_estimation_stream(
    coords = as.matrix(1:n_obs), y = y_obs, x = X_obs,
    neighbor = nb_tu$obs_neighbors, sigma.sq = vy * 0.7, tau.sq = vy * 0.1, phi = 0.05,
    cov.model = "exponential", n_omp = n_omp, nugget_status = 1, verbose = FALSE
  )
  fit_td <- BRISC_estimation_stream(
    coords = as.matrix(1:n_obs), y = y_obs, x = X_obs,
    neighbor = nb_td$obs_neighbors, sigma.sq = vy * 0.7, tau.sq = vy * 0.1, phi = 0.05,
    cov.model = "exponential", n_omp = n_omp, nugget_status = 1, verbose = FALSE
  )
  t_est <- (proc.time() - t1)[3]

  if (has_me) {
    tau_micro_tu <- max(1e-6 * fit_tu$Theta["sigma.sq"], fit_tu$Theta["tau.sq"] - v_bar)
    theta_tu_me  <- c(fit_tu$Theta["sigma.sq"], tau_micro_tu, fit_tu$Theta["phi"])
    tau_micro_td <- max(1e-6 * fit_td$Theta["sigma.sq"], fit_td$Theta["tau.sq"] - v_bar)
    theta_td_me  <- c(fit_td$Theta["sigma.sq"], tau_micro_td, fit_td$Theta["phi"])
    cat(sprintf("    Tail-Up Decomposition: Total Nugget = %.4f -> Microscale = %.4f + Meas Error = %.4f\n",
                fit_tu$Theta["tau.sq"], tau_micro_tu, v_bar))
    cat(sprintf("    Tail-Down Decomposition: Total Nugget = %.4f -> Microscale = %.4f + Meas Error = %.4f\n",
                fit_td$Theta["tau.sq"], tau_micro_td, v_bar))
  }

  t2 <- proc.time()
  bs_tu <- BRISC_bootstrap(fit_tu, n_boot = n_boot, n_omp = n_omp, verbose = FALSE)
  bs_td <- BRISC_bootstrap(fit_td, n_boot = n_boot, n_omp = n_omp, verbose = FALSE)
  t_bs  <- (proc.time() - t2)[3]

  t3 <- proc.time()
  predict_component <- function(fit_obj, nb_obj, theta_me) {
    if (has_me) {
      .Call("C_s3n_predict_meas_error",
            X_obs, y_obs, v_meas, X_all,
            fit_obj$Beta, theta_me,
            as.integer(m_use),
            nb_obj$pred_neighbors$nnIndx, nb_obj$pred_neighbors$nnDist, nb_obj$pred_neighbors$nnWght,
            nb_obj$obsobs_dist, nb_obj$obsobs_wt)
    } else {
      BRISC_prediction_stream(
        BRISC_Out = fit_obj, coords.0 = as.matrix(1:q_tot), X.0 = X_all,
        neighbor = nb_obj$pred_neighbors, obsobs_dist = nb_obj$obsobs_dist, obsobs_wt = nb_obj$obsobs_wt,
        n_omp = n_omp, verbose = FALSE
      )
    }
  }
  pr_tu  <- predict_component(fit_tu, nb_tu, if (has_me) theta_tu_me else NULL)
  pr_td  <- predict_component(fit_td, nb_td, if (has_me) theta_td_me else NULL)
  t_pred <- (proc.time() - t3)[3]

  se_tu <- (pr_tu$prediction.ci[, 2] - pr_tu$prediction) / 1.959964
  se_td <- (pr_td$prediction.ci[, 2] - pr_td$prediction) / 1.959964
  sill_nug_tu <- sum(if (has_me) theta_tu_me[1:2] else fit_tu$Theta[1:2])
  sill_nug_td <- sum(if (has_me) theta_td_me[1:2] else fit_td$Theta[1:2])
  info_tu <- pmax(1e-4, 1 - (se_tu^2) / sill_nug_tu)
  info_td <- pmax(1e-4, 1 - (se_td^2) / sill_nug_td)

  pred_comb <- (info_tu * pr_tu$prediction + info_td * pr_td$prediction) / (info_tu + info_td)
  se_comb   <- sqrt((info_tu^2 * se_tu^2 + info_td^2 * se_td^2) / ((info_tu + info_td)^2))

  if (isTRUE(legacy_nitrate_mode)) {
    pred_df <- data.frame(
      RivID              = df_ord$RivID,
      CatchID            = df_ord$CatchID,
      parent_dn          = df_ord$parent_dn,
      Shreve             = df_ord$Shreve,
      Strahler           = df_ord$Strahler,
      Dist2Mth_km        = df_ord$Dist2Mth_km,
      US_Accum_km        = df_ord$US_Accum_km,
      start_X            = df_ord$start_X,
      start_Y            = df_ord$start_Y,
      mid_X              = df_ord$mid_X,
      mid_Y              = df_ord$mid_Y,
      end_X              = df_ord$end_X,
      end_Y              = df_ord$end_Y,
      Elevation          = df_ord$Elevation,
      Is_Aquifer         = df_ord$Is_Aquifer,
      has_obs            = df_ord$has_obs,
      obs_nitrate_mgL    = df_ord[[obs_col]],
      obs_value          = ifelse(df_ord$has_obs == 1, df_ord[[obs_col]], NA_real_),
      obs_model_sc       = c(round(y_vec_obs, 5), rep(NA_real_, q_tot - n_obs)),
      pred_tu_sc         = round(pr_tu$prediction, 5),
      se_tu_sc           = round(se_tu, 5),
      pred_td_sc         = round(pr_td$prediction, 5),
      se_td_sc           = round(se_td, 5),
      pred_model_sc      = round(pred_comb, 5),
      se_model_sc        = round(se_comb, 5),
      pred_tu            = round(trans_obj$inv(pr_tu$prediction), 4),
      pred_td            = round(trans_obj$inv(pr_td$prediction), 4),
      pred_combined      = round(trans_obj$inv(pred_comb), 4),
      obs_log_nitrate    = c(round(y_vec_obs, 5), rep(NA_real_, q_tot - n_obs)),
      pred_log_tu        = round(pr_tu$prediction, 5),
      se_log_tu          = round(se_tu, 5),
      pred_mgL_tu        = round(trans_obj$inv(pr_tu$prediction), 4),
      pred_log_td        = round(pr_td$prediction, 5),
      se_log_td          = round(se_td, 5),
      pred_mgL_td        = round(trans_obj$inv(pr_td$prediction), 4),
      pred_log_combined  = round(pred_comb, 5),
      se_log_combined    = round(se_comb, 5),
      pred_mgL_combined  = round(trans_obj$inv(pred_comb), 4)
    )
  } else {
    pred_df <- data.frame(
      RivID         = df_ord$RivID,
      CatchID       = df_ord$CatchID,
      Shreve        = df_ord$Shreve,
      Strahler      = df_ord$Strahler,
      start_X       = df_ord$start_X,
      start_Y       = df_ord$start_Y,
      mid_X         = df_ord$mid_X,
      mid_Y         = df_ord$mid_Y,
      end_X         = df_ord$end_X,
      end_Y         = df_ord$end_Y,
      Elevation     = df_ord$Elevation,
      Is_Aquifer    = df_ord$Is_Aquifer,
      has_obs       = df_ord$has_obs,
      obs_value     = df_ord[[obs_col]],
      obs_model_sc  = df_ord$obs_y,
      pred_tu_sc    = pr_tu$prediction,
      pred_td_sc    = pr_td$prediction,
      pred_model_sc = pred_comb,
      se_model_sc   = se_comb,
      pred_tu       = trans_obj$inv(pr_tu$prediction),
      pred_td       = trans_obj$inv(pr_td$prediction),
      pred_combined = trans_obj$inv(pred_comb)
    )
  }
  if ("meas_uncertainty" %in% names(df_ord)) pred_df$meas_uncertainty <- df_ord$meas_uncertainty
  if ("meas_method" %in% names(df_ord))      pred_df$meas_method      <- df_ord$meas_method
  if (isTRUE(legacy_nitrate_mode) && !is.null(out_root) && nzchar(out_root)) {
    pred_df <- add_normalised_columns(pred_df)
    write.csv(pred_df, file.path(out_root, sprintf("predictions/%s_nitrate_predictions.csv", domain_id)), row.names = FALSE)
  }

  # 5-Fold Cross-Validation (using topology vectors directly without copying df_ord)
  set.seed(2026)
  folds   <- sample(rep(1:5, length.out = n_obs))
  cv_ols  <- numeric(n_obs)
  cv_tu   <- numeric(n_obs)
  cv_td   <- numeric(n_obs)
  cv_comb <- numeric(n_obs)
  y_vec   <- y_vec_obs

  riv_vec   <- df_ord$RivID
  par_vec   <- df_ord$parent_dn
  cid_vec   <- df_ord$CatchID
  d2m_vec   <- df_ord$Dist2Mth_km
  afv_vec   <- df_ord$log_AFV
  mx_vec    <- df_ord$mid_X
  my_vec    <- df_ord$mid_Y
  unobs_tail <- if (q_tot > n_obs) (n_obs + 1L):q_tot else integer(0)

  for (k in 1:5) {
    te <- which(folds == k)
    tr <- which(folds != k)
    if (length(te) == 0L || length(tr) < ncol(X_obs)) next

    beta_ols <- tryCatch(
      qr.solve(X_obs[tr, , drop = FALSE], y_vec[tr]),
      error = function(e) {
        b <- suppressWarnings(lm.fit(X_obs[tr, , drop = FALSE], y_vec[tr])$coefficients)
        b[!is.finite(b)] <- 0.0
        b
      }
    )
    cv_ols[te] <- as.numeric(X_obs[te, , drop = FALSE] %*% beta_ols)

    ord_cv     <- c(tr, te, unobs_tail)
    riv_cv     <- riv_vec[ord_cv]
    par_raw_cv <- par_vec[ord_cv]
    par_idx_cv <- match(par_raw_cv, riv_cv, nomatch = 0L)
    par_cv     <- ifelse(par_raw_cv == 0 | par_idx_cv == 0L, -1L, par_idx_cv - 1L)
    n_tr       <- length(tr)
    n_eval     <- n_tr + length(te)
    m_cv       <- min(m_use, max(2L, n_tr - 2L))
    up_cv      <- up_dist[ord_cv]
    cid_cv     <- cid_vec[ord_cv]
    d2m_cv     <- d2m_vec[ord_cv]
    afv_cv     <- afv_vec[ord_cv]
    mx_cv      <- mx_vec[ord_cv]
    my_cv      <- my_vec[ord_cv]
    tt_cv      <- .Call("C_build_tree_tin_tout", as.integer(par_cv))

    X_eval <- X_obs[c(tr, te), , drop = FALSE]
    X_tr   <- X_obs[tr, , drop = FALSE]
    y_tr   <- as.matrix(y_vec[tr])

    run_cv_component <- function(model_int) {
      nb_cv <- .Call("C_compute_s3n_neighbors",
                     as.integer(n_tr), as.integer(m_cv), as.integer(par_cv),
                     as.integer(cid_cv), as.integer(tt_cv[[1]]), as.integer(tt_cv[[2]]),
                     as.double(d2m_cv), as.double(up_cv), as.double(afv_cv),
                     as.double(mx_cv), as.double(my_cv), as.integer(model_int))
      p_cv <- list(
        nnIndx = nb_cv$pred_neighbors$nnIndx[1:(n_eval * m_cv)],
        nnDist = nb_cv$pred_neighbors$nnDist[1:(n_eval * m_cv)],
        nnWght = nb_cv$pred_neighbors$nnWght[1:(n_eval * m_cv)]
      )
      f_cv <- BRISC_estimation_stream(
        coords = as.matrix(1:n_tr), y = y_tr, x = X_tr,
        neighbor = nb_cv$obs_neighbors, sigma.sq = vy * 0.7, tau.sq = vy * 0.1, phi = 0.05,
        cov.model = "exponential", n_omp = n_omp, nugget_status = 1, verbose = FALSE
      )
      if (has_me) {
        v_tr       <- v_meas[tr]
        v_bar_tr   <- mean(v_tr)
        tau_mic_cv <- max(1e-6 * f_cv$Theta["sigma.sq"], f_cv$Theta["tau.sq"] - v_bar_tr)
        theta_me   <- c(f_cv$Theta["sigma.sq"], tau_mic_cv, f_cv$Theta["phi"])
        pr_raw     <- .Call("C_s3n_predict_meas_error",
                            X_tr, y_tr, v_tr, X_eval,
                            f_cv$Beta, theta_me, as.integer(m_cv),
                            p_cv$nnIndx, p_cv$nnDist, p_cv$nnWght,
                            nb_cv$obsobs_dist, nb_cv$obsobs_wt)
        denom <- sum(theta_me[1:2])
      } else {
        pr_raw <- BRISC_prediction_stream(
          BRISC_Out = f_cv, coords.0 = as.matrix(1:n_eval), X.0 = X_eval,
          neighbor = p_cv, obsobs_dist = nb_cv$obsobs_dist, obsobs_wt = nb_cv$obsobs_wt,
          n_omp = n_omp, verbose = FALSE
        )
        denom <- sum(f_cv$Theta[1:2])
      }
      pred_te <- pr_raw$prediction[(n_tr + 1L):n_eval]
      se_te   <- (pr_raw$prediction.ci[(n_tr + 1L):n_eval, 2] - pred_te) / 1.959964
      inf_te  <- pmax(1e-4, 1 - (se_te^2) / denom)
      list(pred = pred_te, inf = inf_te)
    }

    res_tu_cv  <- run_cv_component(0L)
    res_td_cv  <- run_cv_component(1L)
    cv_tu[te]  <- res_tu_cv$pred
    cv_td[te]  <- res_td_cv$pred
    cv_comb[te] <- (res_tu_cv$inf * cv_tu[te] + res_td_cv$inf * cv_td[te]) / (res_tu_cv$inf + res_td_cv$inf)
  }

  rmse <- function(a, b) sqrt(mean((a - b)^2))
  mae  <- function(a, b) mean(abs(a - b))
  r2   <- function(a, b) {
    ss_tot <- sum((a - mean(a))^2)
    if (!is.finite(ss_tot) || ss_tot < 1e-12) return(0.0)
    1 - sum((a - b)^2) / ss_tot
  }

  cv_summary <- data.frame(
    Domain         = domain_label,
    n_obs          = n_obs,
    n_reaches      = q_tot,
    Model          = c("OLS (Non-Spatial)", "S3N Tail-Down (Stream Distance)", "S3N Tail-Up (Flow-Connected + AFV)", "S3N Combined (Tail-Up + Tail-Down)"),
    Spatial_Method = if (has_me) "Measurement Error Kriging" else "Standard S3N"
  )
  rmse_vals <- round(c(rmse(y_vec, cv_ols), rmse(y_vec, cv_td), rmse(y_vec, cv_tu), rmse(y_vec, cv_comb)), 4)
  mae_vals  <- round(c(mae(y_vec, cv_ols),  mae(y_vec, cv_td),  mae(y_vec, cv_tu),  mae(y_vec, cv_comb)), 4)
  r2_vals   <- round(c(r2(y_vec, cv_ols),   r2(y_vec, cv_td),   r2(y_vec, cv_tu),   r2(y_vec, cv_comb)), 4)
  if (isTRUE(legacy_nitrate_mode)) {
    cv_summary$RMSE_log <- rmse_vals
    cv_summary$MAE_log  <- mae_vals
    cv_summary$R2_log   <- r2_vals
  } else {
    cv_summary$RMSE <- rmse_vals
    cv_summary$MAE  <- mae_vals
    cv_summary$R2   <- r2_vals
  }
  print(cv_summary, row.names = FALSE)

  param_names <- c(colnames(X_obs), "sigma.sq", "tau.sq", "phi")
  param_df <- rbind(
    data.frame(Domain = domain_label, Model = "S3N Tail-Up (FC + AFV)", Parameter = param_names,
               Estimate = round(c(fit_tu$Beta, fit_tu$Theta), 4),
               CI_Lower = round(bs_tu$confidence.interval[1, ], 4),
               CI_Upper = round(bs_tu$confidence.interval[2, ], 4)),
    data.frame(Domain = domain_label, Model = "S3N Tail-Down (Stream)", Parameter = param_names,
               Estimate = round(c(fit_td$Beta, fit_td$Theta), 4),
               CI_Lower = round(bs_td$confidence.interval[1, ], 4),
               CI_Upper = round(bs_td$confidence.interval[2, ], 4))
  )
  if (has_me) {
    param_df <- rbind(
      param_df,
      data.frame(Domain = domain_label, Model = "S3N Tail-Up (FC + AFV)", Parameter = "tau.sq (microscale)",
                 Estimate = round(tau_micro_tu, 4), CI_Lower = NA, CI_Upper = NA),
      data.frame(Domain = domain_label, Model = "S3N Tail-Down (Stream)", Parameter = "tau.sq (microscale)",
                 Estimate = round(tau_micro_td, 4), CI_Lower = NA, CI_Upper = NA),
      data.frame(Domain = domain_label, Model = "Measurement Error", Parameter = "mean_meas_variance",
                 Estimate = round(v_bar, 6), CI_Lower = NA, CI_Upper = NA)
    )
  }

  cv_pts <- data.frame(
    RivID   = df_ord$RivID[1:n_obs],
    CatchID = df_ord$CatchID[1:n_obs]
  )
  if (isTRUE(legacy_nitrate_mode)) cv_pts$obs_log <- y_vec else cv_pts$obs_y <- y_vec
  cv_pts$cv_ols  <- cv_ols
  cv_pts$cv_tu   <- cv_tu
  cv_pts$cv_td   <- cv_td
  cv_pts$cv_comb <- cv_comb

  list(
    pred_df    = pred_df,
    cv_summary = cv_summary,
    param_df   = param_df,
    cv_pts     = cv_pts,
    runtimes   = c(neighbor_sec = t_nb, est_sec = t_est, boot_sec = t_bs, pred_sec = t_pred)
  )
}
