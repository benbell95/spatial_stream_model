#!/usr/bin/env Rscript
# ==============================================================================
# 02_fit_s3n_models.R
#
# Fits Scalable Spatial Stream Network (S3N) models (Kunke, Olden & McCormick 2026)
# across the Ordnance Survey Open Rivers Network (ORN v2) and EA River
# Surveillance Network (RSN):
#   1. All England Monitored Catchments (90,876 reaches, 1,210 observed reaches)
#   2. Humber River Basin (CatchID = 80: 14,048 reaches, 222 observed reaches)
#   3. Thames River Basin (CatchID = 3055: 7,491 reaches, 116 observed reaches)
#   4. Severn River Basin (CatchID = 2786: 9,198 reaches, 110 observed reaches)
#
# Compares:
#   - Non-Spatial OLS Regression
#   - S3N Exponential Tail-Up (S3N-TU: flow-connected with AFV confluence weights)
#   - S3N Exponential Tail-Down (S3N-TD: hydrologic stream distance with unit weights)
#   - S3N Combined (S3N-TU+TD: spatial information-gain weighted combination)
# ==============================================================================

BASE_DIR <- if (file.exists("scripts/s3n_plot_functions.R")) "." else if (file.exists("s3n_plot_functions.R")) ".." else "."
source(file.path(BASE_DIR, "scripts/s3n_plot_functions.R"))
BASE_DIR <- resolve_base_dir()

args    <- commandArgs(trailingOnly = TRUE)
dir_cfg <- parse_custom_dir_args(args, base_dir = BASE_DIR)
GIS_DIR  <- dir_cfg$gis_dir
DATA_DIR <- dir_cfg$data_dir
OUT_BASE <- dir_cfg$out_dir

if ("--list-networks" %in% args) {
  py_cmd <- sprintf("python3 %s --list-networks --gis-dir=%s",
                    shQuote(rel_path(BASE_DIR, "scripts", "01_preprocess_orn.py")),
                    shQuote(GIS_DIR))
  system(py_cmd)
  quit(save = "no", status = 0)
}
if ("--list-vars" %in% args || "--list-covariates" %in% args) {
  print_available_covariates()
  quit(save = "no", status = 0)
}
net_arg_vec <- grep("^--network=", args, value = TRUE)
network_cli <- if (length(net_arg_vec) > 0) sub("^--network=", "", net_arg_vec[1]) else "orn"
PROC_DIR    <- resolve_network_processed_dir(BASE_DIR, network_cli, data_dir = DATA_DIR, gis_dir = GIS_DIR)
net_meta    <- load_network_metadata(PROC_DIR)

vars_args <- grep("^--(vars|covariates)=", args, value = TRUE)
var_specs <- if (length(vars_args) > 0) {
  unlist(strsplit(sub("^--(vars|covariates)=", "", vars_args), ",", fixed = TRUE))
} else {
  character(0)
}
var_specs <- trimws(var_specs[nchar(trimws(var_specs)) > 0])

user_lib <- Sys.getenv("R_LIBS_USER")
if (nzchar(user_lib) && dir.exists(user_lib) && !(normalizePath(user_lib, mustWork = FALSE) %in% normalizePath(.libPaths(), mustWork = FALSE))) {
  .libPaths(c(user_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(S3N)
})

dyn.load(rel_path(BASE_DIR, "src", "s3n_fast_neighbors.so"))

vars_tag <- format_covariates_folder_tag(var_specs)
sub_tag  <- if (net_meta$abbr != "orn") paste0("net_", net_meta$abbr) else ""
if (nzchar(vars_tag)) {
  sub_tag <- if (nzchar(sub_tag)) paste0(sub_tag, "_", vars_tag) else vars_tag
}
OUT_ROOT <- if (nzchar(sub_tag)) rel_path(OUT_BASE, sub_tag) else OUT_BASE
dir.create(file.path(OUT_ROOT, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "predictions"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "models"), recursive = TRUE, showWarnings = FALSE)

cat("======================================================================\n")
cat(sprintf("  Scalable Spatial Stream Network (S3N) Nitrate Modeling Pipeline [%s]\n", net_meta$label))
cat("======================================================================\n\n")

t_start <- proc.time()
df_full <- read.csv(file.path(PROC_DIR, "england_rsn_reaches.csv"))
catch_lu <- if (file.exists(file.path(PROC_DIR, "rsn_catchments_lookup.csv"))) {
  read.csv(file.path(PROC_DIR, "rsn_catchments_lookup.csv"), stringsAsFactors = FALSE)
} else {
  NULL
}
cat(sprintf("Loaded %d stream reaches (%d with RSN Nitrate observations) across %d catchments.\n",
            nrow(df_full), sum(df_full$has_obs == 1), length(unique(df_full$CatchID))))

# Function to fit S3N-TU, S3N-TD, Bootstrap CIs, 5-Fold CV, and Full Network Prediction
analyze_domain <- function(df, domain_id, domain_label, m_neighbors = 15L, n_boot = 100L, n_omp = 1L) {
  cat(sprintf("\n----------------------------------------------------------------------\n"))
  cat(sprintf("  Domain: %s (%d reaches, %d observed)\n",
              domain_label, nrow(df), sum(df$has_obs == 1)))
  cat(sprintf("----------------------------------------------------------------------\n"))

  # Order observed reaches downstream-to-upstream within each catchment, followed by prediction reaches
  obs_idx  <- which(df$has_obs == 1)
  obs_idx  <- obs_idx[order(df$CatchID[obs_idx], df$obs_up_dist_km[obs_idx])]
  pred_idx <- which(df$has_obs == 0)
  ord      <- c(obs_idx, pred_idx)
  df_ord   <- df[ord, ]

  n_obs <- length(obs_idx)
  q_tot <- nrow(df_ord)

  id_map   <- setNames(seq_len(q_tot), df_ord$RivID)
  parent_0 <- ifelse(df_ord$parent_dn == 0, -1L, id_map[as.character(df_ord$parent_dn)] - 1L)
  up_dist  <- ifelse(df_ord$has_obs == 1, df_ord$obs_up_dist_km, df_ord$up_dist_mid_km)

  X_all <- build_design_matrix(df_ord, include_peat = (domain_id == "england"), var_specs = var_specs)
  cat(sprintf("  Fixed-effect covariates (%d): %s\n", ncol(X_all), paste(colnames(X_all), collapse = ", ")))
  X_obs <- X_all[1:n_obs, , drop = FALSE]
  y_obs <- as.matrix(df_ord$log_nitrate_mean[1:n_obs])

  # 1. Build Tree & Neighbors (Tail-Up and Tail-Down)
  t0 <- proc.time()
  tt <- .Call("C_build_tree_tin_tout", as.integer(parent_0))
  nb_tu <- .Call("C_compute_s3n_neighbors",
                 as.integer(n_obs), as.integer(m_neighbors), as.integer(parent_0),
                 as.integer(df_ord$CatchID), as.integer(tt[[1]]), as.integer(tt[[2]]),
                 as.double(df_ord$Dist2Mth_km), as.double(up_dist), as.double(df_ord$log_AFV),
                 as.double(df_ord$mid_X), as.double(df_ord$mid_Y), as.integer(0L))
  nb_td <- .Call("C_compute_s3n_neighbors",
                 as.integer(n_obs), as.integer(m_neighbors), as.integer(parent_0),
                 as.integer(df_ord$CatchID), as.integer(tt[[1]]), as.integer(tt[[2]]),
                 as.double(df_ord$Dist2Mth_km), as.double(up_dist), as.double(df_ord$log_AFV),
                 as.double(df_ord$mid_X), as.double(df_ord$mid_Y), as.integer(1L))
  t_nb <- (proc.time() - t0)[3]
  cat(sprintf("  [1] Built stream topology & neighbor structures in %.3f sec.\n", t_nb))

  # 2. Fit S3N Tail-Up and S3N Tail-Down Models
  t1 <- proc.time()
  fit_tu <- BRISC_estimation_stream(
    coords = as.matrix(1:n_obs), y = y_obs, x = X_obs,
    neighbor = nb_tu$obs_neighbors, sigma.sq = 0.4, tau.sq = 0.05, phi = 0.05,
    cov.model = "exponential", n_omp = n_omp, nugget_status = 1, verbose = FALSE
  )
  fit_td <- BRISC_estimation_stream(
    coords = as.matrix(1:n_obs), y = y_obs, x = X_obs,
    neighbor = nb_td$obs_neighbors, sigma.sq = 0.4, tau.sq = 0.05, phi = 0.05,
    cov.model = "exponential", n_omp = n_omp, nugget_status = 1, verbose = FALSE
  )
  t_est <- (proc.time() - t1)[3]
  cat(sprintf("  [2] Estimated S3N-TU (logLik=%.2f) and S3N-TD (logLik=%.2f) in %.3f sec.\n",
              fit_tu$log_likelihood, fit_td$log_likelihood, t_est))

  # 3. Parametric Bootstrap for 95% Confidence Intervals
  t2 <- proc.time()
  bs_tu <- BRISC_bootstrap(fit_tu, n_boot = n_boot, n_omp = n_omp, verbose = FALSE)
  bs_td <- BRISC_bootstrap(fit_td, n_boot = n_boot, n_omp = n_omp, verbose = FALSE)
  t_bs <- (proc.time() - t2)[3]
  cat(sprintf("  [3] Completed %d bootstrap replicates in %.3f sec.\n", n_boot, t_bs))

  # 4. Predict across all q_tot stream reaches in the domain
  t3 <- proc.time()
  pr_tu <- BRISC_prediction_stream(
    BRISC_Out = fit_tu, coords.0 = as.matrix(1:q_tot), X.0 = X_all,
    neighbor = nb_tu$pred_neighbors, obsobs_dist = nb_tu$obsobs_dist, obsobs_wt = nb_tu$obsobs_wt,
    n_omp = n_omp, verbose = FALSE
  )
  pr_td <- BRISC_prediction_stream(
    BRISC_Out = fit_td, coords.0 = as.matrix(1:q_tot), X.0 = X_all,
    neighbor = nb_td$pred_neighbors, obsobs_dist = nb_td$obsobs_dist, obsobs_wt = nb_td$obsobs_wt,
    n_omp = n_omp, verbose = FALSE
  )
  t_pred <- (proc.time() - t3)[3]
  cat(sprintf("  [4] Predicted Nitrate across all %d stream reaches in %.3f sec.\n", q_tot, t_pred))

  se_tu <- (pr_tu$prediction.ci[, 2] - pr_tu$prediction) / 1.959964
  se_td <- (pr_td$prediction.ci[, 2] - pr_td$prediction) / 1.959964
  info_tu <- pmax(1e-4, 1 - (se_tu^2) / sum(fit_tu$Theta[1:2]))
  info_td <- pmax(1e-4, 1 - (se_td^2) / sum(fit_td$Theta[1:2]))

  pred_log_comb <- (info_tu * pr_tu$prediction + info_td * pr_td$prediction) / (info_tu + info_td)
  se_comb       <- sqrt((info_tu^2 * se_tu^2 + info_td^2 * se_td^2) / ((info_tu + info_td)^2))

  # Save reach-level predictions table
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
    obs_nitrate_mgL    = df_ord$nitrate_mean,
    obs_log_nitrate    = df_ord$log_nitrate_mean,
    pred_log_tu        = round(pr_tu$prediction, 5),
    se_log_tu          = round(se_tu, 5),
    pred_mgL_tu        = round(exp(pr_tu$prediction), 4),
    pred_log_td        = round(pr_td$prediction, 5),
    se_log_td          = round(se_td, 5),
    pred_mgL_td        = round(exp(pr_td$prediction), 4),
    pred_log_combined  = round(pred_log_comb, 5),
    se_log_combined    = round(se_comb, 5),
    pred_mgL_combined  = round(exp(pred_log_comb), 4)
  )
  write.csv(pred_df, file.path(OUT_ROOT, sprintf("predictions/%s_nitrate_predictions.csv", domain_id)), row.names = FALSE)

  # 5. 5-Fold Out-of-Sample Cross-Validation
  set.seed(2026)
  folds <- sample(rep(1:5, length.out = n_obs))
  cv_ols  <- numeric(n_obs)
  cv_tu   <- numeric(n_obs)
  cv_td   <- numeric(n_obs)
  cv_comb <- numeric(n_obs)
  y_vec   <- as.numeric(y_obs)

  for (k in 1:5) {
    te <- which(folds == k)
    tr <- which(folds != k)

    # OLS
    beta_ols   <- qr.solve(X_obs[tr, , drop = FALSE], y_vec[tr])
    cv_ols[te] <- as.numeric(X_obs[te, , drop = FALSE] %*% beta_ols)

    # Reorder so training obs are first (1..n_tr), then test obs (n_tr+1..n_obs)
    ord_cv     <- c(tr, te, (n_obs + 1):q_tot)
    df_cv      <- df_ord[ord_cv, ]
    id_map_cv  <- setNames(seq_len(q_tot), df_cv$RivID)
    par_cv     <- ifelse(df_cv$parent_dn == 0, -1L, id_map_cv[as.character(df_cv$parent_dn)] - 1L)
    n_tr       <- length(tr)
    n_eval     <- n_obs
    up_cv      <- ifelse(df_cv$has_obs == 1, df_cv$obs_up_dist_km, df_cv$up_dist_mid_km)
    tt_cv      <- .Call("C_build_tree_tin_tout", as.integer(par_cv))

    # CV Tail-Up
    nb_tu_cv <- .Call("C_compute_s3n_neighbors",
                      as.integer(n_tr), as.integer(m_neighbors), as.integer(par_cv),
                      as.integer(df_cv$CatchID), as.integer(tt_cv[[1]]), as.integer(tt_cv[[2]]),
                      as.double(df_cv$Dist2Mth_km), as.double(up_cv), as.double(df_cv$log_AFV),
                      as.double(df_cv$mid_X), as.double(df_cv$mid_Y), as.integer(0L))
    p_tu_cv <- list(
      nnIndx = nb_tu_cv$pred_neighbors$nnIndx[1:(n_eval * m_neighbors)],
      nnDist = nb_tu_cv$pred_neighbors$nnDist[1:(n_eval * m_neighbors)],
      nnWght = nb_tu_cv$pred_neighbors$nnWght[1:(n_eval * m_neighbors)]
    )
    f_tu_cv <- BRISC_estimation_stream(
      coords = as.matrix(1:n_tr), y = as.matrix(y_vec[tr]), x = X_obs[tr, , drop = FALSE],
      neighbor = nb_tu_cv$obs_neighbors, sigma.sq = 0.4, tau.sq = 0.05, phi = 0.05,
      cov.model = "exponential", n_omp = n_omp, nugget_status = 1, verbose = FALSE
    )
    pr_tu_cv <- BRISC_prediction_stream(
      BRISC_Out = f_tu_cv, coords.0 = as.matrix(1:n_eval),
      X.0 = rbind(X_obs[tr, , drop = FALSE], X_obs[te, , drop = FALSE]),
      neighbor = p_tu_cv, obsobs_dist = nb_tu_cv$obsobs_dist, obsobs_wt = nb_tu_cv$obsobs_wt,
      n_omp = n_omp, verbose = FALSE
    )
    cv_tu[te] <- pr_tu_cv$prediction[(n_tr + 1):n_eval]
    se_tu_cv  <- (pr_tu_cv$prediction.ci[(n_tr + 1):n_eval, 2] - cv_tu[te]) / 1.959964
    inf_tu_cv <- pmax(1e-4, 1 - (se_tu_cv^2) / sum(f_tu_cv$Theta[1:2]))

    # CV Tail-Down
    nb_td_cv <- .Call("C_compute_s3n_neighbors",
                      as.integer(n_tr), as.integer(m_neighbors), as.integer(par_cv),
                      as.integer(df_cv$CatchID), as.integer(tt_cv[[1]]), as.integer(tt_cv[[2]]),
                      as.double(df_cv$Dist2Mth_km), as.double(up_cv), as.double(df_cv$log_AFV),
                      as.double(df_cv$mid_X), as.double(df_cv$mid_Y), as.integer(1L))
    p_td_cv <- list(
      nnIndx = nb_td_cv$pred_neighbors$nnIndx[1:(n_eval * m_neighbors)],
      nnDist = nb_td_cv$pred_neighbors$nnDist[1:(n_eval * m_neighbors)],
      nnWght = nb_td_cv$pred_neighbors$nnWght[1:(n_eval * m_neighbors)]
    )
    f_td_cv <- BRISC_estimation_stream(
      coords = as.matrix(1:n_tr), y = as.matrix(y_vec[tr]), x = X_obs[tr, , drop = FALSE],
      neighbor = nb_td_cv$obs_neighbors, sigma.sq = 0.4, tau.sq = 0.05, phi = 0.05,
      cov.model = "exponential", n_omp = n_omp, nugget_status = 1, verbose = FALSE
    )
    pr_td_cv <- BRISC_prediction_stream(
      BRISC_Out = f_td_cv, coords.0 = as.matrix(1:n_eval),
      X.0 = rbind(X_obs[tr, , drop = FALSE], X_obs[te, , drop = FALSE]),
      neighbor = p_td_cv, obsobs_dist = nb_td_cv$obsobs_dist, obsobs_wt = nb_td_cv$obsobs_wt,
      n_omp = n_omp, verbose = FALSE
    )
    cv_td[te] <- pr_td_cv$prediction[(n_tr + 1):n_eval]
    se_td_cv  <- (pr_td_cv$prediction.ci[(n_tr + 1):n_eval, 2] - cv_td[te]) / 1.959964
    inf_td_cv <- pmax(1e-4, 1 - (se_td_cv^2) / sum(f_td_cv$Theta[1:2]))

    cv_comb[te] <- (inf_tu_cv * cv_tu[te] + inf_td_cv * cv_td[te]) / (inf_tu_cv + inf_td_cv)
  }

  rmse <- function(a, b) sqrt(mean((a - b)^2))
  mae  <- function(a, b) mean(abs(a - b))
  r2   <- function(a, b) 1 - sum((a - b)^2) / sum((a - mean(a))^2)

  cv_summary <- data.frame(
    Domain    = domain_label,
    n_obs     = n_obs,
    n_reaches = q_tot,
    Model     = c("OLS (Non-Spatial)", "S3N Tail-Down (Stream Distance)", "S3N Tail-Up (Flow-Connected + AFV)", "S3N Combined (Tail-Up + Tail-Down)"),
    RMSE_log  = round(c(rmse(y_vec, cv_ols), rmse(y_vec, cv_td), rmse(y_vec, cv_tu), rmse(y_vec, cv_comb)), 4),
    MAE_log   = round(c(mae(y_vec, cv_ols),  mae(y_vec, cv_td),  mae(y_vec, cv_tu),  mae(y_vec, cv_comb)), 4),
    R2_log    = round(c(r2(y_vec, cv_ols),   r2(y_vec, cv_td),   r2(y_vec, cv_tu),   r2(y_vec, cv_comb)), 4)
  )
  print(cv_summary, row.names = FALSE)

  # Format parameter table with 95% Bootstrap CIs
  param_names <- c(colnames(X_obs), "sigma.sq", "tau.sq", "phi")
  phi_shift   <- c(rep(0, ncol(X_obs)), 0, 0, 0.001)
  param_df <- rbind(
    data.frame(
      Domain    = domain_label,
      Model     = "S3N Tail-Up (FC + AFV)",
      Parameter = param_names,
      Estimate  = round(c(fit_tu$Beta, fit_tu$Theta) + phi_shift, 4),
      CI_Lower  = round(bs_tu$confidence.interval[1, ] + phi_shift, 4),
      CI_Upper  = round(bs_tu$confidence.interval[2, ] + phi_shift, 4)
    ),
    data.frame(
      Domain    = domain_label,
      Model     = "S3N Tail-Down (Stream)",
      Parameter = param_names,
      Estimate  = round(c(fit_td$Beta, fit_td$Theta) + phi_shift, 4),
      CI_Lower  = round(bs_td$confidence.interval[1, ] + phi_shift, 4),
      CI_Upper  = round(bs_td$confidence.interval[2, ] + phi_shift, 4)
    )
  )

  # Save CV point predictions for England & Humber for diagnostic plotting
  cv_pts <- data.frame(
    RivID     = df_ord$RivID[1:n_obs],
    CatchID   = df_ord$CatchID[1:n_obs],
    obs_log   = y_vec,
    cv_ols    = cv_ols,
    cv_tu     = cv_tu,
    cv_td     = cv_td,
    cv_comb   = cv_comb
  )

  list(
    cv_summary = cv_summary,
    param_df   = param_df,
    cv_pts     = cv_pts,
    runtimes   = c(neighbor_sec = t_nb, est_sec = t_est, boot_sec = t_bs, pred_sec = t_pred)
  )
}

humber_info <- resolve_core_basin(df_full, catch_lu, net_meta, "humber")
thames_info <- resolve_core_basin(df_full, catch_lu, net_meta, "thames")
severn_info <- resolve_core_basin(df_full, catch_lu, net_meta, "severn")

res_england <- analyze_domain(df_full, "england", "All England RSN Network", m_neighbors = 15L, n_boot = 100L)
res_humber  <- analyze_domain(subset(df_full, CatchID == humber_info$catch_id), "humber", humber_info$domain_label, m_neighbors = 15L, n_boot = 100L)
res_thames  <- analyze_domain(subset(df_full, CatchID == thames_info$catch_id), "thames", thames_info$domain_label, m_neighbors = 15L, n_boot = 100L)
res_severn  <- analyze_domain(subset(df_full, CatchID == severn_info$catch_id), "severn", severn_info$domain_label, m_neighbors = 15L, n_boot = 100L)

all_cv <- rbind(res_england$cv_summary, res_humber$cv_summary, res_thames$cv_summary, res_severn$cv_summary)
all_params <- rbind(res_england$param_df, res_humber$param_df, res_thames$param_df, res_severn$param_df)
all_cv$Network_Abbr      <- net_meta$abbr
all_cv$Network_Label     <- net_meta$short_label
all_params$Network_Abbr  <- net_meta$abbr
all_params$Network_Label <- net_meta$short_label

write.csv(all_cv, file.path(OUT_ROOT, "tables/cv_model_comparison.csv"), row.names = FALSE)
write.csv(all_params, file.path(OUT_ROOT, "tables/s3n_parameter_estimates.csv"), row.names = FALSE)
write.csv(res_england$cv_pts, file.path(OUT_ROOT, "predictions/england_cv_predictions.csv"), row.names = FALSE)
write.csv(res_humber$cv_pts, file.path(OUT_ROOT, "predictions/humber_cv_predictions.csv"), row.names = FALSE)

saveRDS(list(england = res_england, humber = res_humber, thames = res_thames, severn = res_severn),
        file = file.path(OUT_ROOT, "models/s3n_model_results.rds"))

cat(sprintf("\n======================================================================\n"))
cat(sprintf("  All S3N models, bootstrap CIs, CV, and predictions finished in %.2f sec.\n",
            (proc.time() - t_start)[3]))
cat(sprintf("======================================================================\n"))

