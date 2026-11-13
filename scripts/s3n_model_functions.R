# ==============================================================================
# Shared S3N Backend Loading and Domain Fitting / Cross-Validation Engine
# ==============================================================================
# Unified modeling engine shared by run_s3n_determinand.R and 02_fit_s3n_models.R.
# ==============================================================================

.s3n_env <- new.env(parent = emptyenv())
.s3n_env$has_s3n_pkg <- FALSE

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
  if (requireNamespace("S3N", quietly = TRUE)) {
    suppressPackageStartupMessages(library(S3N))
    .s3n_env$has_s3n_pkg <- TRUE
  } else {
    .s3n_env$has_s3n_pkg <- FALSE
  }

  so_name <- paste0("s3n_fast_neighbors", .Platform$dynlib.ext)
  so_path <- rel_path(base_dir, "src", so_name)
  if (!is.loaded("C_compute_s3n_neighbors")) {
    if (!file.exists(so_path)) {
      stop(sprintf("Compiled C++ neighbor library not found at '%s'. Run: R CMD SHLIB %s",
                   so_path, rel_path(base_dir, "src", "s3n_fast_neighbors.cpp")))
    }
    dyn.load(so_path)
  }
  invisible(TRUE)
}

S3N_KERNEL_SPECS <- list(
  exponential = list(id = 0L, name = "exponential", label = "Exponential",      formula = "exp(-phi * d)"),
  spherical   = list(id = 1L, name = "spherical",   label = "Spherical",        formula = "(1 - 1.5*(phi*d) + 0.5*(phi*d)^3) * I(phi*d < 1)"),
  mariah      = list(id = 2L, name = "mariah",      label = "Mariah",           formula = "log(phi*d + 1) / (phi*d)"),
  linear      = list(id = 3L, name = "linear",      label = "Linear-with-Sill", formula = "max(0, 1 - phi*d)")
)

normalize_kernel_name <- function(k_str, default = "exponential") {
  if (is.null(k_str) || !nzchar(trimws(k_str))) return(default)
  k <- tolower(trimws(k_str))
  if (k %in% c("exp", "exponential")) return("exponential")
  if (k %in% c("sph", "spherical"))   return("spherical")
  if (k %in% c("mar", "mariah"))      return("mariah")
  if (k %in% c("lin", "linear", "linear-with-sill")) return("linear")
  if (k %in% c("auto", "best", "aic")) return("auto")
  stop(sprintf("Unknown --kernel='%s'. Valid kernels: exponential, spherical, mariah, linear, auto", k_str))
}

parse_kernel_arg <- function(args = character(0), default = "exponential") {
  k_arg <- grep("^--(kernel|cov-model)=", args, value = TRUE)
  if (length(k_arg) > 0) {
    val <- sub("^--(kernel|cov-model)=", "", k_arg[1])
    return(normalize_kernel_name(val, default = default))
  }
  normalize_kernel_name(default, default = "exponential")
}

parse_combine_arg <- function(args = character(0), default = "joint") {
  c_arg <- grep("^--(combine|combine-mode)=", args, value = TRUE)
  val <- if (length(c_arg) > 0) tolower(trimws(sub("^--(combine|combine-mode)=", "", c_arg[1]))) else tolower(default)
  if (val %in% c("joint", "additive", "joint-nngp", "full")) return("joint")
  if (val %in% c("precision", "legacy", "weighted")) return("precision")
  stop(sprintf("Unknown --combine='%s'. Valid modes: joint, precision", val))
}

parse_afv_arg <- function(args = character(0), default = "length") {
  a_arg <- grep("^--(afv|afv-scheme)=", args, value = TRUE)
  val <- if (length(a_arg) > 0) tolower(trimws(sub("^--(afv|afv-scheme)=", "", a_arg[1]))) else tolower(default)
  if (val %in% c("length", "len", "us_accum", "default")) return("length")
  if (val %in% c("shreve", "order", "magnitude")) return("shreve")
  if (val %in% c("hydro_yield", "hydro-yield", "yield", "hydro")) return("hydro_yield")
  stop(sprintf("Unknown --afv='%s'. Valid AFV schemes: length, shreve, hydro_yield", val))
}

S3N_TREND_SPECS <- list(
  linear = list(name = "linear", label = "Linear Fixed Effects (GLS)",
                desc = "Standard parametric linear fixed-effects design matrix X * Beta"),
  spline = list(name = "spline", label = "Natural Cubic Splines + S3N",
                desc = "GCV-penalized natural cubic B-spline basis expansion (splines::ns) + S3N kriging"),
  gam    = list(name = "gam",    label = "Penalized GAM Smooths + S3N",
                desc = "REML-penalized cubic regression splines (mgcv::gam) + S3N residual/GLS kriging"),
  rf     = list(name = "rf",     label = "Random Forest Regression Kriging (RF + S3N)",
                desc = "Out-of-Bag Random Forest / subspace tree ensemble + S3N residual kriging"),
  gbm    = list(name = "gbm",    label = "Gradient Boosted Trees + S3N (GBM + S3N)",
                desc = "Cross-fitted Stochastic Gradient Boosted Regression Trees + S3N residual kriging")
)

normalize_trend_name <- function(t_str, default = "linear") {
  if (is.null(t_str) || !nzchar(trimws(t_str))) return(default)
  t_val <- tolower(trimws(t_str))
  if (t_val %in% c("linear", "ols", "gls", "default", "parametric")) return("linear")
  if (t_val %in% c("spline", "splines", "ns", "bspline", "cubic")) return("spline")
  if (t_val %in% c("gam", "gam-s3n", "gam_s3n", "mgcv", "smooth")) return("gam")
  if (t_val %in% c("rf", "rf-s3n", "rf_s3n", "randomforest", "random-forest", "forest")) return("rf")
  if (t_val %in% c("gbm", "gbm-s3n", "gbm_s3n", "boost", "boosted", "gbdt", "trees")) return("gbm")
  stop(sprintf("Unknown --trend='%s'. Valid trend modes: linear, spline, gam, rf, gbm", t_str))
}

parse_trend_arg <- function(args = character(0), default = "linear") {
  if ("--gam-s3n" %in% args) return("gam")
  if (any(c("--spline", "--splines") %in% args)) return("spline")
  if (any(c("--rf", "--rf-s3n", "--random-forest") %in% args)) return("rf")
  if (any(c("--gbm", "--gbm-s3n", "--boosted") %in% args)) return("gbm")
  t_arg <- grep("^--(trend|fixed-effects|mean-model)=", args, value = TRUE)
  if (length(t_arg) > 0) {
    val <- sub("^--(trend|fixed-effects|mean-model)=", "", t_arg[1])
    return(normalize_trend_name(val, default = default))
  }
  normalize_trend_name(default, default = "linear")
}

parse_spline_df_arg <- function(args = character(0), default = 3L) {
  s_arg <- grep("^--(spline-df|gam-k)=", args, value = TRUE)
  if (length(s_arg) > 0) {
    val <- suppressWarnings(as.integer(sub("^--(spline-df|gam-k)=", "", s_arg[1])))
    if (is.finite(val) && val >= 2L && val <= 10L) return(val)
  }
  as.integer(default)
}

parse_trees_arg <- function(args = character(0), default = 100L) {
  tr_arg <- grep("^--(trees|n-trees|rf-trees|gbm-trees)=", args, value = TRUE)
  if (length(tr_arg) > 0) {
    val <- suppressWarnings(as.integer(sub("^--(trees|n-trees|rf-trees|gbm-trees)=", "", tr_arg[1])))
    if (is.finite(val) && val >= 10L && val <= 2000L) return(val)
  }
  as.integer(default)
}

print_available_trends <- function() {
  cat("Available Fixed-Effect / Mean-Trend Engines (--trend=<mode>):\n")
  for (sp in S3N_TREND_SPECS) {
    cat(sprintf("  %-12s : %-44s  (%s)\n", sp$name, sp$label, sp$desc))
  }
  cat("\nOptional Trend Tuning Flags:\n")
  cat("  --spline-df=<int> : Degrees of freedom per continuous covariate for --trend=spline|gam (default: 3)\n")
  cat("  --trees=<int>     : Number of ensemble trees for --trend=rf|gbm (default: 100)\n")
}

print_available_kernels <- function() {
  cat("Available S3N Stream Covariance Kernels (--kernel=<name>):\n")
  for (sp in S3N_KERNEL_SPECS) {
    cat(sprintf("  %-12s : %-18s  K(d; phi) = %s\n", sp$name, sp$label, sp$formula))
  }
  cat("  auto         : Automatic Selection  Selects best kernel per component by NNGP AIC\n\n")
  cat("Available S3N Component Combination Modes (--combine=<mode>):\n")
  cat("  joint        : Joint Additive NNGP (Tail-Up + Tail-Down + Euclidean in a single Cholesky system)\n")
  cat("  precision    : Precision-weighted combination of separate Tail-Up and Tail-Down kriging predictions\n\n")
  cat("Available Tail-Up AFV Confluence Weighting Schemes (--afv=<scheme>):\n")
  cat("  length       : Cumulative upstream channel length (US_Accum_km)\n")
  cat("  shreve       : Shreve hydrological link magnitude (number of upstream headwater sources)\n")
  cat("  hydro_yield  : Elevation/orographic rainfall + Principal Aquifer baseflow hydrological yield proxy\n\n")
  print_available_trends()
  cat("\n")
  print_available_validation_options()
}

S3N_CV_MODES <- c("random", "catchment", "strahler", "spatial_block")

normalize_cv_mode_name <- function(cv_mode, default = "random") {
  m <- tolower(trimws(as.character(cv_mode)))
  if (length(m) == 0L || !nzchar(m)) return(default)
  if (m %in% c("random", "5fold", "kfold")) return("random")
  if (m %in% c("catchment", "leave_catchment_out", "lco", "basin")) return("catchment")
  if (m %in% c("strahler", "headwater", "headwater_extrapolation")) return("strahler")
  if (m %in% c("spatial_block", "spatial", "block", "kmeans")) return("spatial_block")
  warning(sprintf("Unknown --cv-mode '%s'; defaulting to '%s'.", cv_mode, default))
  default
}

parse_cv_mode_arg <- function(args, default = "random") {
  hit <- grep("^--cv-mode=", args, value = TRUE)
  if (length(hit) == 0L) return(default)
  normalize_cv_mode_name(sub("^--cv-mode=", "", hit[1]), default = default)
}

S3N_UK_THRESHOLDS <- list(
  "0117"      = list(val = 11.3,  dir = "above", unit = "mg/L",  ref = "UK DWI / Nitrates Directive (11.3 mg/L as N = 50 mg/L NO3)"),
  "0116"      = list(val = 11.3,  dir = "above", unit = "mg/L",  ref = "UK Nitrates Directive proxy (11.3 mg/L TON as N)"),
  "0180"      = list(val = 0.069, dir = "above", unit = "mg/L",  ref = "UK WFD 2015 Directions (indicative lowland high-alkalinity Good/Moderate SRP)"),
  "0111"      = list(val = 0.60,  dir = "above", unit = "mg/L",  ref = "UK WFD 2015 Directions (indicative lowland high-alkalinity Good/Moderate 90%ile NH4-N)"),
  "0119"      = list(val = 0.021, dir = "above", unit = "mg/L",  ref = "UK WFD / Freshwater Fish EQS (annual mean un-ionised NH3-N)"),
  "0085"      = list(val = 5.0,   dir = "above", unit = "mg/L",  ref = "UK WFD 2015 Directions (indicative lowland high-alkalinity Good/Moderate 90%ile BOD)"),
  "9924"      = list(val = 6.0,   dir = "below", unit = "mg/L",  ref = "UK WFD 2015 Directions (indicative salmonid Good oxygenation threshold)"),
  "0061"      = list(val = 6.0,   dir = "below", unit = "pH",    ref = "UK WFD 2015 Directions (minimum Good pH boundary; use --threshold=9.0 --threshold-dir=above for upper bound)"),
  "0076"      = list(val = 25.0,  dir = "above", unit = "CEL",   ref = "UK WFD 2015 Directions (cyprinid Good/Moderate 98%ile temperature)"),
  "0052"      = list(val = 1.2,   dir = "above", unit = "ug/L",  ref = "UK WFD EQS (bioavailable Dissolved Lead annual average)"),
  "0106"      = list(val = 0.08,  dir = "above", unit = "ug/L",  ref = "UK WFD EQS (Dissolved Cadmium annual average, low-hardness baseline)"),
  "whpt_aspt" = list(val = 5.0,   dir = "below", unit = "ASPT",  ref = "UK WFD / UKTAG (indicative Good/Moderate WHPT ASPT boundary)")
)

parse_threshold_arg <- function(args, default = NA_real_) {
  hit <- grep("^--threshold=", args, value = TRUE)
  if (length(hit) == 0L) return(default)
  raw <- tolower(trimws(sub("^--threshold=", "", hit[1])))
  if (raw %in% c("uk", "wfd", "default", "ea")) return("uk")
  val <- suppressWarnings(as.numeric(raw))
  if (!is.finite(val)) {
    warning(sprintf("Invalid --threshold '%s'; ignoring exceedance risk mapping.", hit[1]))
    return(default)
  }
  val
}

resolve_threshold_spec <- function(threshold_arg, threshold_dir_arg, det_code, args = character(0)) {
  if (identical(threshold_arg, "uk")) {
    key <- tolower(trimws(as.character(det_code)))
    if (!(key %in% names(S3N_UK_THRESHOLDS))) {
      stop(sprintf(
        "No built-in UK/WFD default threshold is registered for '%s'. Run --list-thresholds to see supported codes or pass a numeric --threshold=<val>.",
        det_code
      ))
    }
    spec <- S3N_UK_THRESHOLDS[[key]]
    explicit_dir <- any(grepl("^--threshold-dir=", args))
    dir_use <- if (explicit_dir) threshold_dir_arg else spec$dir
    cat(sprintf("  [UK/WFD Threshold Preset] %s -> %.4g %s (%s) [%s]\n",
                det_code, spec$val, spec$unit, dir_use, spec$ref))
    return(list(val = spec$val, dir = dir_use, ref = spec$ref))
  }
  list(val = threshold_arg, dir = threshold_dir_arg, ref = NULL)
}

print_available_thresholds <- function() {
  cat("Built-in UK / WFD Regulatory Threshold Presets (--threshold=uk):\n")
  for (k in names(S3N_UK_THRESHOLDS)) {
    s <- S3N_UK_THRESHOLDS[[k]]
    cat(sprintf("  %-10s : %6.4g %-5s (%-5s) — %s\n", k, s$val, s$unit, s$dir, s$ref))
  }
  cat("\nNote: WFD phosphorus, ammonia, and DO boundaries are site-type dependent (alkalinity & altitude)\n")
  cat("under the WFD (Standards and Classification) Directions (England & Wales) 2015.\n")
  cat("Values above are indicative lowland/high-alkalinity defaults; override any time with --threshold=<val>.\n")
}

parse_threshold_dir_arg <- function(args, default = "above") {
  hit <- grep("^--threshold-dir=", args, value = TRUE)
  if (length(hit) == 0L) return(default)
  d <- tolower(trimws(sub("^--threshold-dir=", "", hit[1])))
  if (d %in% c("above", "gt", "high", "exceed")) return("above")
  if (d %in% c("below", "lt", "low", "deficit")) return("below")
  warning(sprintf("Unknown --threshold-dir '%s'; defaulting to '%s'.", d, default))
  default
}

parse_design_sites_arg <- function(args, default = 0L) {
  hit <- grep("^--(design-sites|optimal-design)=", args, value = TRUE)
  if (length(hit) == 0L) return(as.integer(default))
  val <- suppressWarnings(as.integer(sub("^--(design-sites|optimal-design)=", "", hit[1])))
  if (is.na(val) || val < 0L) return(as.integer(default))
  val
}

print_available_validation_options <- function() {
  cat("Available Cross-Validation Modes (--cv-mode=<mode>):\n")
  cat("  random        : Standard 5-fold random station CV (default)\n")
  cat("  catchment     : 5-fold Leave-Catchment-Out CV grouped by CatchID\n")
  cat("  strahler      : Train on mainstems (Strahler >= 3), test on headwaters (Strahler <= 2)\n")
  cat("  spatial_block : 5-fold spatial block CV via k-means clustering on (mid_X, mid_Y)\n\n")
  cat("Regulatory Exceedance & Optimal Monitoring Design Flags:\n")
  cat("  --threshold=<val|uk>       : Raw-scale regulatory threshold (or 'uk' for built-in UK/WFD presets)\n")
  cat("  --threshold-dir=above|below: Exceedance direction (above = P(Y > T); below = P(Y < T))\n")
  cat("  --list-thresholds          : Print built-in UK / WFD regulatory threshold presets\n")
  cat("  --design-sites=<K>         : Rank top K unmonitored reaches minimizing mean kriging variance\n\n")
  cat("Spatio-Temporal & Censored Kriging Flags:\n")
  cat("  --st-borrow                : Enable dynamic precision-weighted borrowing across time (--by=...)\n")
  cat("  --censor=em|tobit|ros|half_lod : Left-censored nondetect method (em/tobit = Spatial Tobit EM kriging)\n")
}

parse_st_borrow_arg <- function(args) {
  any(c("--st-borrow", "--temporal-borrow") %in% args)
}

trunc_norm_below_stats <- function(mu, sd_val, upper) {
  sd_safe <- pmax(as.numeric(sd_val), 1e-8)
  alpha   <- (as.numeric(upper) - as.numeric(mu)) / sd_safe
  # Numerically stable Mills ratio phi(alpha) / Phi(alpha) in log space
  log_phi <- dnorm(alpha, log = TRUE)
  log_Phi <- pnorm(alpha, log.p = TRUE)
  mills   <- exp(log_phi - log_Phi)
  mills[!is.finite(mills)] <- pmax(0, -alpha[!is.finite(mills)])
  mean_cond <- as.numeric(mu) - sd_safe * mills
  mean_cond <- pmin(mean_cond, as.numeric(upper))
  var_cond  <- (sd_safe^2) * pmax(1e-6, 1.0 - alpha * mills - mills^2)
  list(mean = mean_cond, var = var_cond)
}

run_spatial_tobit_em_kriging <- function(df_ord, n_obs, X_obs, y_vec_init, is_cens_obs, lod_sc_obs,
                                         m_use, parent_0, tt, up_dist, nb_td,
                                         estimate_component, kernel_req = "exponential",
                                         max_iter = 6L, tol = 1e-3) {
  cens_idx   <- which(is_cens_obs & is.finite(lod_sc_obs))
  uncens_idx <- which(!is_cens_obs)
  if (length(cens_idx) == 0L) {
    return(list(y_em = y_vec_init, n_cens = 0L, n_iter = 0L, converged = TRUE, cond_var = numeric(0)))
  }

  y_curr <- as.numeric(y_vec_init)
  # Ensure initial imputed values respect the censoring upper bound
  y_curr[cens_idx] <- pmin(y_curr[cens_idx], lod_sc_obs[cens_idx] - 1e-4)

  # Precompute neighbor structure from uncensored -> censored reaches for exact E-step kriging
  nb_u2c <- NULL
  if (length(uncens_idx) >= 5L) {
    ord_u2c   <- c(uncens_idx, cens_idx)
    df_u2c    <- df_ord[ord_u2c, , drop = FALSE]
    par_u2c_i <- match(df_u2c$parent_dn, df_u2c$RivID, nomatch = 0L)
    par_u2c   <- ifelse(df_u2c$parent_dn == 0 | par_u2c_i == 0L, -1L, par_u2c_i - 1L)
    tt_u2c    <- .Call("C_build_tree_tin_tout", as.integer(par_u2c))
    m_u2c     <- min(as.integer(m_use), max(2L, length(uncens_idx) - 2L))
    nb_u2c    <- .Call("C_compute_s3n_neighbors",
                       as.integer(length(uncens_idx)), as.integer(m_u2c), as.integer(par_u2c),
                       as.integer(df_u2c$CatchID), as.integer(tt_u2c[[1]]), as.integer(tt_u2c[[2]]),
                       as.double(df_u2c$Dist2Mth_km), as.double(up_dist[ord_u2c]),
                       as.double(df_u2c$log_AFV), as.double(df_u2c$mid_X),
                       as.double(df_u2c$mid_Y), 1L)
  }

  converged <- FALSE
  n_iter    <- 0L
  cond_var  <- rep(NA_real_, length(cens_idx))

  for (it in seq_len(max_iter)) {
    n_iter <- it
    y_prev <- y_curr
    vy_it  <- max(1e-4, var(y_curr))

    # M-step: fit S3N Tail-Down NNGP on completed observations
    fit_em <- estimate_component(X_obs, as.matrix(y_curr), nb_td$obs_neighbors, vy_it, kernel_req)
    beta_em <- as.numeric(fit_em$Beta)
    theta_em <- as.numeric(fit_em$Theta)
    kid_em   <- if (!is.null(fit_em$kernel_id)) as.integer(fit_em$kernel_id) else 0L

    # E-step: exact conditional mean & SD for censored reaches given uncensored neighbors
    if (!is.null(nb_u2c)) {
      n_u <- length(uncens_idx)
      n_c <- length(cens_idx)
      n_e <- n_u + n_c
      m_u2c <- min(as.integer(m_use), max(2L, n_u - 2L))
      p_u2c <- list(
        nnIndx = nb_u2c$pred_neighbors$nnIndx[1:(n_e * m_u2c)],
        nnDist = nb_u2c$pred_neighbors$nnDist[1:(n_e * m_u2c)],
        nnWght = nb_u2c$pred_neighbors$nnWght[1:(n_e * m_u2c)]
      )
      pr_c <- .Call(
        "C_s3n_predict_meas_error",
        X_obs[uncens_idx, , drop = FALSE],
        as.matrix(y_curr[uncens_idx]),
        numeric(0),
        X_obs[c(uncens_idx, cens_idx), , drop = FALSE],
        as.double(beta_em),
        as.double(c(theta_em, kid_em)),
        as.integer(m_u2c),
        p_u2c$nnIndx, p_u2c$nnDist, p_u2c$nnWght,
        nb_u2c$obsobs_dist, nb_u2c$obsobs_wt
      )
      tail_idx <- (n_u + 1L):n_e
      mu_c <- as.numeric(pr_c$prediction[tail_idx])
      sd_c <- pmax(1e-4, as.numeric(pr_c$prediction.ci[tail_idx, 2] - pr_c$prediction[tail_idx]) / 1.959964)
    } else {
      mu_c <- as.numeric(X_obs[cens_idx, , drop = FALSE] %*% beta_em)
      sd_c <- rep(sqrt(max(1e-4, theta_em[1] + theta_em[2])), length(cens_idx))
    }

    tn <- trunc_norm_below_stats(mu_c, sd_c, lod_sc_obs[cens_idx])
    y_curr[cens_idx] <- tn$mean
    cond_var <- tn$var

    max_diff <- max(abs(y_curr[cens_idx] - y_prev[cens_idx]))
    if (is.finite(max_diff) && max_diff < tol) {
      converged <- TRUE
      break
    }
  }

  list(
    y_em      = y_curr,
    n_cens    = length(cens_idx),
    n_iter    = n_iter,
    converged = converged,
    cond_var  = cond_var
  )
}

apply_spatiotemporal_borrowing <- function(slice_list, trans_obj, threshold_val = NA_real_, threshold_dir = "above") {
  valid_idx <- which(!vapply(slice_list, is.null, logical(1)))
  t_len <- length(valid_idx)
  if (t_len < 2L) return(list(slices = slice_list, rho_temporal = NA_real_))

  # Estimate lag-1 temporal correlation from observed reaches across consecutive slices
  cors <- numeric(0)
  for (k in seq_len(t_len - 1L)) {
    d1 <- slice_list[[valid_idx[k]]]$eng_pred
    d2 <- slice_list[[valid_idx[k + 1L]]]$eng_pred
    both <- which(d1$has_obs == 1L & d2$has_obs == 1L & is.finite(d1$obs_model_sc) & is.finite(d2$obs_model_sc))
    if (length(both) >= 10L) {
      r_k <- suppressWarnings(cor(d1$obs_model_sc[both], d2$obs_model_sc[both]))
      if (is.finite(r_k)) cors <- c(cors, r_k)
    }
  }
  rho_t <- if (length(cors) > 0L) min(0.85, max(0.15, mean(cors))) else 0.60

  borrow_domain <- function(df_list) {
    if (any(vapply(df_list, is.null, logical(1)))) return(df_list)
    q <- nrow(df_list[[1]])
    Z <- vapply(df_list, function(d) as.numeric(d$pred_model_sc), numeric(q))
    V <- vapply(df_list, function(d) pmax(as.numeric(d$se_model_sc)^2, 1e-8), numeric(q))
    P <- 1.0 / V
    out_list <- df_list
    for (t_i in seq_len(t_len)) {
      w_vec <- rho_t^abs(seq_len(t_len) - t_i)
      Pw <- sweep(P, 2L, w_vec, "*")
      prec_tot <- rowSums(Pw)
      z_new <- rowSums(Z * Pw) / prec_tot
      se_new <- sqrt(1.0 / prec_tot)
      d_upd <- out_list[[t_i]]
      d_upd$pred_model_sc <- round(z_new, 5)
      d_upd$se_model_sc   <- round(se_new, 5)
      d_upd$pred_combined <- round(trans_obj$inv(z_new), 4)
      if (!is.na(threshold_val) && is.finite(threshold_val)) {
        ex <- compute_exceedance_risk(z_new, se_new, threshold_val, trans_obj, dir = threshold_dir)
        d_upd$exceed_prob     <- round(ex$prob, 5)
        d_upd$compliance_tier <- ex$tier
      }
      out_list[[t_i]] <- d_upd
    }
    out_list
  }

  eng_upd <- borrow_domain(lapply(valid_idx, function(i) slice_list[[i]]$eng_pred))
  hum_upd <- borrow_domain(lapply(valid_idx, function(i) slice_list[[i]]$hum_pred))
  tha_upd <- borrow_domain(lapply(valid_idx, function(i) slice_list[[i]]$tha_pred))
  sev_upd <- borrow_domain(lapply(valid_idx, function(i) slice_list[[i]]$sev_pred))

  for (k in seq_along(valid_idx)) {
    idx <- valid_idx[k]
    slice_list[[idx]]$eng_pred <- eng_upd[[k]]
    slice_list[[idx]]$hum_pred <- hum_upd[[k]]
    slice_list[[idx]]$tha_pred <- tha_upd[[k]]
    slice_list[[idx]]$sev_pred <- sev_upd[[k]]
  }
  list(slices = slice_list, rho_temporal = rho_t)
}

compute_spatiotemporal_trends <- function(slice_list, trans_obj) {
  valid_s <- Filter(Negate(is.null), slice_list)
  t_len <- length(valid_s)
  if (t_len < 2L) return(NULL)

  ref <- valid_s[[1]]$eng_pred
  q <- nrow(ref)
  Z <- vapply(valid_s, function(s) as.numeric(s$eng_pred$pred_model_sc), numeric(q))
  V <- vapply(valid_s, function(s) pmax(as.numeric(s$eng_pred$se_model_sc)^2, 1e-8), numeric(q))
  W <- 1.0 / V
  t_seq <- seq_len(t_len)

  W_sum <- rowSums(W)
  t_bar <- as.numeric(W %*% t_seq) / W_sum
  z_bar <- rowSums(W * Z) / W_sum

  dt_mat <- outer(t_bar, t_seq, function(tb, t) t - tb)
  dz_mat <- sweep(Z, 1L, z_bar, "-")

  denom  <- pmax(rowSums(W * (dt_mat^2)), 1e-12)
  slope  <- rowSums(W * dt_mat * dz_mat) / denom
  se_slp <- sqrt(1.0 / denom)
  z_stat <- slope / se_slp
  p_val  <- 2.0 * stats::pnorm(-abs(z_stat))

  delta_raw <- trans_obj$inv(Z[, t_len]) - trans_obj$inv(Z[, 1L])

  trend_tier <- ifelse(
    p_val < 0.05 & slope > 0, "Significant Increase (p < 0.05)",
    ifelse(p_val < 0.20 & slope > 0, "Moderate Increase (p < 0.20)",
           ifelse(p_val < 0.05 & slope < 0, "Significant Decrease (p < 0.05)",
                  ifelse(p_val < 0.20 & slope < 0, "Moderate Decrease (p < 0.20)",
                         "Stable / No Clear Trend")))
  )

  data.frame(
    RivID          = ref$RivID,
    CatchID        = ref$CatchID,
    Shreve         = ref$Shreve,
    start_X        = ref$start_X,
    start_Y        = ref$start_Y,
    mid_X          = ref$mid_X,
    mid_Y          = ref$mid_Y,
    end_X          = ref$end_X,
    end_Y          = ref$end_Y,
    trend_slope_sc = round(slope, 5),
    trend_se_sc    = round(se_slp, 5),
    trend_z        = round(z_stat, 4),
    trend_p        = signif(p_val, 4),
    delta_raw      = round(delta_raw, 4),
    trend_tier     = trend_tier,
    stringsAsFactors = FALSE
  )
}

compute_spatiotemporal_transitions <- function(slice_list, period_labels, trans_obj, cyclic = TRUE) {
  valid_idx <- which(!vapply(slice_list, is.null, logical(1)))
  t_len <- length(valid_idx)
  if (t_len < 2L) return(NULL)

  clean_lbl <- function(x) {
    s <- sub("^Season:\\s*", "", as.character(x))
    paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
  }

  pairs <- lapply(seq_len(t_len - 1L), function(k) c(valid_idx[k], valid_idx[k + 1L]))
  if (isTRUE(cyclic) && t_len >= 3L) {
    pairs[[length(pairs) + 1L]] <- c(valid_idx[t_len], valid_idx[1L])
  }

  items <- vector("list", length(pairs))
  csv_parts <- vector("list", length(pairs))

  for (pi in seq_along(pairs)) {
    ia <- pairs[[pi]][1]
    ib <- pairs[[pi]][2]
    da <- slice_list[[ia]]$eng_pred
    db <- slice_list[[ib]]$eng_pred

    za <- as.numeric(da$pred_model_sc)
    zb <- as.numeric(db$pred_model_sc)
    va <- pmax(as.numeric(da$se_model_sc)^2, 1e-8)
    vb <- pmax(as.numeric(db$se_model_sc)^2, 1e-8)

    dz     <- zb - za
    se_dz  <- sqrt(va + vb)
    z_stat <- dz / se_dz
    p_val  <- 2.0 * stats::pnorm(-abs(z_stat))
    d_raw  <- trans_obj$inv(zb) - trans_obj$inv(za)

    tier <- ifelse(
      p_val < 0.05 & dz > 0, "Significant Increase (p < 0.05)",
      ifelse(p_val < 0.20 & dz > 0, "Moderate Increase (p < 0.20)",
             ifelse(p_val < 0.05 & dz < 0, "Significant Decrease (p < 0.05)",
                    ifelse(p_val < 0.20 & dz < 0, "Moderate Decrease (p < 0.20)",
                           "Stable / No Clear Trend")))
    )

    lbl_a <- clean_lbl(period_labels[ia])
    lbl_b <- clean_lbl(period_labels[ib])
    tr_lbl <- sprintf("%s -> %s", lbl_a, lbl_b)

    df_item <- data.frame(
      RivID       = da$RivID,
      CatchID     = da$CatchID,
      Shreve      = da$Shreve,
      start_X     = da$start_X,
      start_Y     = da$start_Y,
      mid_X       = da$mid_X,
      mid_Y       = da$mid_Y,
      end_X       = da$end_X,
      end_Y       = da$end_Y,
      delta_sc    = round(dz, 5),
      delta_se_sc = round(se_dz, 5),
      trend_z     = round(z_stat, 4),
      trend_p     = signif(p_val, 4),
      delta_raw   = round(d_raw, 4),
      trend_tier  = tier,
      stringsAsFactors = FALSE
    )

    items[[pi]] <- list(
      from_label       = lbl_a,
      to_label         = lbl_b,
      transition_label = tr_lbl,
      df               = df_item
    )
    csv_parts[[pi]] <- data.frame(
      Transition  = tr_lbl,
      RivID       = da$RivID,
      CatchID     = da$CatchID,
      delta_sc    = round(dz, 5),
      delta_se_sc = round(se_dz, 5),
      trend_z     = round(z_stat, 4),
      trend_p     = signif(p_val, 4),
      delta_raw   = round(d_raw, 4),
      trend_tier  = tier,
      stringsAsFactors = FALSE
    )
  }

  list(
    items          = items,
    transitions_df = do.call(rbind, csv_parts)
  )
}

build_cv_folds <- function(df_obs, cv_mode = "random", k_folds = 5L, seed = 2026L) {
  req_mode <- normalize_cv_mode_name(cv_mode, default = "random")
  n_obs <- nrow(df_obs)
  set.seed(seed)

  make_spatial_blocks <- function() {
    coords <- cbind(as.numeric(df_obs$mid_X), as.numeric(df_obs$mid_Y))
    ok <- is.finite(coords[, 1]) & is.finite(coords[, 2])
    if (sum(ok) < k_folds || nrow(unique(coords[ok, , drop = FALSE])) < k_folds) {
      return(list(folds = sample(rep(seq_len(k_folds), length.out = n_obs)),
                  actual_mode = "random",
                  note = "insufficient unique coordinates for spatial_block; fell back to random"))
    }
    km <- stats::kmeans(scale(coords[ok, , drop = FALSE]), centers = k_folds, nstart = 10L)
    folds <- sample(rep(seq_len(k_folds), length.out = n_obs))
    folds[ok] <- as.integer(km$cluster)
    list(folds = folds, actual_mode = "spatial_block", note = NULL)
  }

  if (req_mode == "catchment") {
    cats <- as.character(df_obs$CatchID)
    u_cats <- unique(cats[!is.na(cats)])
    if (length(u_cats) < 2L) {
      res <- make_spatial_blocks()
      res$note <- "domain has < 2 catchments; fell back to spatial_block"
      return(res)
    }
    k_use <- min(as.integer(k_folds), length(u_cats))
    cat_counts <- sort(table(cats), decreasing = TRUE)
    fold_loads <- integer(k_use)
    cat_to_fold <- integer(length(cat_counts))
    names(cat_to_fold) <- names(cat_counts)
    for (nm in names(cat_counts)) {
      f_min <- which.min(fold_loads)
      cat_to_fold[nm] <- f_min
      fold_loads[f_min] <- fold_loads[f_min] + as.integer(cat_counts[nm])
    }
    folds <- as.integer(cat_to_fold[cats])
    return(list(folds = folds, actual_mode = "catchment", note = NULL))
  }

  if (req_mode == "strahler") {
    st <- as.numeric(df_obs$Strahler)
    n_main <- sum(is.finite(st) & st >= 3)
    n_head <- sum(is.finite(st) & st <= 2)
    if (n_main < 5L || n_head < 5L) {
      res <- make_spatial_blocks()
      res$note <- sprintf("insufficient Strahler split (mainstem=%d, headwater=%d); fell back to spatial_block",
                          n_main, n_head)
      return(res)
    }
    folds <- rep(NA_integer_, n_obs)
    folds[is.finite(st) & st >= 3] <- 0L
    folds[is.finite(st) & st <= 2] <- 1L
    return(list(folds = folds, actual_mode = "strahler", note = NULL))
  }

  if (req_mode == "spatial_block") {
    return(make_spatial_blocks())
  }

  list(folds = sample(rep(seq_len(k_folds), length.out = n_obs)),
       actual_mode = "random", note = NULL)
}

compute_exceedance_risk <- function(pred_sc, se_sc, threshold_val, trans_obj, dir = "above") {
  thresh_sc <- as.numeric(trans_obj$fwd(threshold_val))
  if (!is.finite(thresh_sc)) {
    stop(sprintf("Threshold %.4g is outside the domain of transformation '%s'.",
                 threshold_val, trans_obj$short_label))
  }
  se_safe <- pmax(as.numeric(se_sc), 1e-8)
  z <- (thresh_sc - as.numeric(pred_sc)) / se_safe
  p <- if (identical(dir, "below")) {
    stats::pnorm(z, lower.tail = TRUE)
  } else {
    stats::pnorm(z, lower.tail = FALSE)
  }
  p <- pmin(1.0, pmax(0.0, p))
  tier <- ifelse(
    p < 0.05, "Confidently Compliant",
    ifelse(p < 0.20, "Likely Compliant",
           ifelse(p <= 0.80, "Uncertain / Priority Monitoring",
                  ifelse(p <= 0.95, "Likely Exceeding", "Confidently Exceeding")))
  )
  list(prob = p, tier = tier, thresh_sc = thresh_sc)
}

compute_optimal_monitoring_design <- function(df_ord, n_obs, se_comb, parent_0, tt, up_dist,
                                            params_5, k_design = 10L) {
  k_req <- as.integer(k_design)
  if (is.na(k_req) || k_req <= 0L) return(NULL)
  cand_all <- which(df_ord$has_obs == 0L)
  if (length(cand_all) == 0L) return(NULL)

  # Deterministic evenly-spaced thinning for huge national domains (> 4000 reaches)
  cand_use <- if (length(cand_all) > 4000L) {
    cand_all[unique(round(seq(1, length(cand_all), length.out = 4000L)))]
  } else {
    cand_all
  }
  q_tot <- nrow(df_ord)
  eval_use <- if (q_tot > 4000L) {
    unique(round(seq(1, q_tot, length.out = 4000L)))
  } else {
    seq_len(q_tot)
  }

  k_use <- min(k_req, length(cand_use))
  res <- .Call(
    "C_s3n_optimal_design",
    as.integer(parent_0),
    as.integer(df_ord$CatchID),
    as.integer(tt[[1]]),
    as.integer(tt[[2]]),
    as.double(df_ord$Dist2Mth_km),
    as.double(up_dist),
    as.double(df_ord$mid_X),
    as.double(df_ord$mid_Y),
    as.double(se_comb),
    as.double(params_5),
    as.integer(cand_use - 1L),
    as.integer(eval_use - 1L),
    as.integer(k_use)
  )
  if (length(res$cand_pos) == 0L) return(NULL)
  sel_rows <- cand_use[res$cand_pos + 1L]
  data.frame(
    Rank            = seq_along(sel_rows),
    RivID           = df_ord$RivID[sel_rows],
    CatchID         = df_ord$CatchID[sel_rows],
    Strahler        = df_ord$Strahler[sel_rows],
    Dist2Mth_km     = round(df_ord$Dist2Mth_km[sel_rows], 3),
    US_Accum_km     = round(df_ord$US_Accum_km[sel_rows], 3),
    mid_X           = round(df_ord$mid_X[sel_rows], 1),
    mid_Y           = round(df_ord$mid_Y[sel_rows], 1),
    Elevation       = round(df_ord$Elevation[sel_rows], 1),
    Current_SE_sc   = round(se_comb[sel_rows], 5),
    Var_Reduced     = round(res$var_reduced, 6),
    Pct_Var_Reduced = round(res$pct_var_reduced, 4),
    Cum_Pct_Reduced = round(res$cum_pct_reduced, 4),
    Mean_SE_Before  = round(res$mean_se_before, 5),
    Mean_SE_After   = round(res$mean_se_after, 5)
  )
}

fit_nonlinear_trend <- function(X_tr, y_tr_vec, X_eval, trend_mode = "linear",
                                spline_df = 3L, n_trees = 100L, seed = 2026L) {
  mode <- normalize_trend_name(trend_mode, default = "linear")
  n_tr <- nrow(X_tr)
  q_ev <- nrow(X_eval)
  y_tr <- as.numeric(y_tr_vec)

  solve_linear_ols <- function(X_m, y_m) {
    tryCatch(
      qr.solve(X_m, y_m),
      error = function(e) {
        b <- suppressWarnings(lm.fit(X_m, y_m)$coefficients)
        b[!is.finite(b)] <- 0.0
        b
      }
    )
  }

  beta_lin <- solve_linear_ols(X_tr, y_tr)
  pred_lin_tr <- as.numeric(X_tr %*% beta_lin)
  pred_lin_ev <- as.numeric(X_eval %*% beta_lin)

  cov_cols <- setdiff(colnames(X_tr), "Intercept")
  if (mode == "linear" || length(cov_cols) == 0L || n_tr < 12L) {
    return(list(
      mode        = "linear",
      label       = S3N_TREND_SPECS$linear$label,
      g_perp_tr   = numeric(n_tr),
      g_perp_eval = numeric(q_ev),
      pred_oof_tr = pred_lin_tr,
      pred_eval   = pred_lin_ev,
      info_df     = NULL
    ))
  }

  # Identify continuous vs binary/discrete columns
  is_cont <- vapply(cov_cols, function(cn) {
    u <- unique(X_tr[, cn])
    length(u[is.finite(u)]) > 4L
  }, logical(1))
  cont_cols <- cov_cols[is_cont]
  disc_cols <- cov_cols[!is_cont]

  pred_oof_tr <- pred_lin_tr
  pred_eval   <- pred_lin_ev
  info_rows   <- list()

  if (mode == "spline") {
    df_use <- max(2L, min(as.integer(spline_df), max(2L, floor(n_tr / max(3L, 2L * length(cont_cols))))))
    basis_tr_list <- list(Intercept = rep(1.0, n_tr))
    basis_ev_list <- list(Intercept = rep(1.0, q_ev))
    pen_diag      <- c(0.0) # do not penalize intercept or linear terms

    for (cn in cov_cols) {
      x_t <- as.numeric(X_tr[, cn])
      x_e <- as.numeric(X_eval[, cn])
      basis_tr_list[[cn]] <- x_t
      basis_ev_list[[cn]] <- x_e
      pen_diag <- c(pen_diag, 0.0)

      if (cn %in% cont_cols && df_use >= 2L && requireNamespace("splines", quietly = TRUE)) {
        rng <- range(x_t, na.rm = TRUE)
        x_e_clip <- pmin(rng[2], pmax(rng[1], x_e))
        sp_tr <- tryCatch(splines::ns(x_t, df = df_use), error = function(e) NULL)
        if (!is.null(sp_tr) && ncol(sp_tr) >= 2L) {
          sp_ev <- predict(sp_tr, newx = x_e_clip)
          # Orthogonalize nonlinear spline columns against [1, x_t] on training set
          Z_lin <- cbind(1.0, x_t)
          for (k_sp in seq_len(ncol(sp_tr))) {
            z_k_tr <- as.numeric(sp_tr[, k_sp])
            z_k_ev <- as.numeric(sp_ev[, k_sp])
            coef_p <- solve_linear_ols(Z_lin, z_k_tr)
            r_tr   <- z_k_tr - as.numeric(Z_lin %*% coef_p)
            r_ev   <- z_k_ev - as.numeric(cbind(1.0, x_e_clip) %*% coef_p)
            sd_r   <- sd(r_tr)
            if (is.finite(sd_r) && sd_r > 1e-5) {
              col_nm <- sprintf("%s_sp%d", cn, k_sp)
              basis_tr_list[[col_nm]] <- r_tr / sd_r
              basis_ev_list[[col_nm]] <- r_ev / sd_r
              pen_diag <- c(pen_diag, 1.0)
            }
          }
        }
      }
    }

    Phi_tr <- do.call(cbind, basis_tr_list)
    Phi_ev <- do.call(cbind, basis_ev_list)
    p_phi  <- ncol(Phi_tr)
    PtP    <- crossprod(Phi_tr)
    Pty    <- crossprod(Phi_tr, y_tr)

    # Select ridge penalty lambda on nonlinear spline columns via analytic GCV / LOO
    lam_grid <- c(0.05, 0.2, 1.0, 5.0, 25.0, 100.0)
    best_gcv <- Inf
    best_oof <- pred_lin_tr
    best_ev  <- pred_lin_ev
    best_b   <- numeric(p_phi)
    best_edf <- ncol(X_tr)

    for (lam in lam_grid) {
      P_mat <- PtP + diag(lam * pen_diag, nrow = p_phi)
      P_inv <- tryCatch(chol2inv(chol(P_mat)), error = function(e) NULL)
      if (is.null(P_inv)) next
      b_hat <- as.numeric(P_inv %*% Pty)
      f_in  <- as.numeric(Phi_tr %*% b_hat)
      # Hat matrix diagonal h_ii = rowSums((Phi_tr %*% P_inv) * Phi_tr)
      h_diag <- pmin(0.95, pmax(0.0, rowSums((Phi_tr %*% P_inv) * Phi_tr)))
      loo_err <- (y_tr - f_in) / (1.0 - h_diag)
      gcv_mse <- mean(loo_err^2)
      if (is.finite(gcv_mse) && gcv_mse < best_gcv) {
        best_gcv <- gcv_mse
        best_oof <- y_tr - loo_err
        best_ev  <- as.numeric(Phi_ev %*% b_hat)
        best_b   <- b_hat
        best_edf <- sum(h_diag)
      }
    }
    pred_oof_tr <- best_oof
    pred_eval   <- best_ev
    names(best_b) <- colnames(Phi_tr)
    for (cn in cont_cols) {
      sp_cols <- grep(paste0("^", cn, "_sp"), names(best_b), value = TRUE)
      nl_norm <- if (length(sp_cols) > 0) sqrt(sum(best_b[sp_cols]^2)) else 0.0
      info_rows[[cn]] <- data.frame(
        Parameter = paste0("Spline_NL:", cn),
        Estimate  = round(nl_norm, 4),
        CI_Lower  = NA_real_,
        CI_Upper  = NA_real_
      )
    }
    info_rows[["edf"]] <- data.frame(
      Parameter = "Spline_Total_EDF",
      Estimate  = round(best_edf, 3),
      CI_Lower  = NA_real_,
      CI_Upper  = NA_real_
    )

  } else if (mode == "gam") {
    df_tr <- as.data.frame(X_tr[, cov_cols, drop = FALSE])
    df_tr$y_resp <- y_tr
    df_ev <- as.data.frame(X_eval[, cov_cols, drop = FALSE])

    # Clamp evaluation continuous covariates to training range +/- 0.5 SD to prevent runaway cubic extrapolation
    for (cn in cont_cols) {
      rng <- range(df_tr[[cn]], na.rm = TRUE)
      pad <- 0.25 * diff(rng)
      df_ev[[cn]] <- pmin(rng[2] + pad, pmax(rng[1] - pad, df_ev[[cn]]))
    }

    gam_terms <- character(0)
    for (cn in cont_cols) {
      n_u <- length(unique(df_tr[[cn]]))
      k_c <- max(3L, min(as.integer(spline_df) + 2L, n_u - 1L, floor(n_tr / max(3L, length(cont_cols) + 1L))))
      if (k_c >= 3L) {
        gam_terms <- c(gam_terms, sprintf("s(%s, bs='cr', k=%d)", cn, k_c))
      } else {
        gam_terms <- c(gam_terms, cn)
      }
    }
    gam_terms <- c(gam_terms, disc_cols)
    fml_str <- paste("y_resp ~", paste(gam_terms, collapse = " + "))

    gam_fit <- if (requireNamespace("mgcv", quietly = TRUE)) {
      tryCatch(
        suppressWarnings(mgcv::gam(stats::as.formula(fml_str), data = df_tr, method = "REML")),
        error = function(e) NULL
      )
    } else NULL

    if (!is.null(gam_fit)) {
      f_in   <- as.numeric(stats::fitted(gam_fit))
      h_hat  <- pmin(0.92, pmax(0.0, as.numeric(gam_fit$hat)))
      pred_oof_tr <- y_tr - (y_tr - f_in) / (1.0 - h_hat)
      pred_eval   <- as.numeric(suppressWarnings(stats::predict(gam_fit, newdata = df_ev, type = "response")))

      s_tab <- tryCatch(summary(gam_fit)$s.table, error = function(e) NULL)
      if (!is.null(s_tab) && nrow(s_tab) > 0) {
        for (r_i in seq_len(nrow(s_tab))) {
          s_nm <- gsub("^s\\(|\\)$", "", rownames(s_tab)[r_i])
          info_rows[[s_nm]] <- data.frame(
            Parameter = paste0("GAM_EDF:", s_nm),
            Estimate  = round(as.numeric(s_tab[r_i, "edf"]), 4),
            CI_Lower  = NA_real_,
            CI_Upper  = NA_real_
          )
        }
      }
    } else {
      # Automatic fallback to GCV cubic spline basis if n_tr is too small for mgcv::gam REML
      sp_fb <- fit_nonlinear_trend(X_tr, y_tr, X_eval, trend_mode = "spline", spline_df = spline_df)
      pred_oof_tr <- sp_fb$pred_oof_tr
      pred_eval   <- sp_fb$pred_eval
      if (!is.null(sp_fb$info_df)) info_rows[["fb"]] <- sp_fb$info_df
    }

  } else if (mode == "rf") {
    set.seed(seed)
    df_tr <- as.data.frame(X_tr[, cov_cols, drop = FALSE])
    df_ev <- as.data.frame(X_eval[, cov_cols, drop = FALSE])
    p_cov <- length(cov_cols)
    mtry  <- max(1L, floor(0.70 * p_cov))
    ntree <- max(20L, as.integer(n_trees))

    oob_sum <- numeric(n_tr)
    oob_cnt <- integer(n_tr)
    ev_sum  <- numeric(q_ev)
    imp_vec <- setNames(numeric(p_cov), cov_cols)

    min_sp <- max(6L, floor(n_tr / 12L))
    for (b in seq_len(ntree)) {
      idx_b  <- sample.int(n_tr, size = n_tr, replace = TRUE)
      oob_b  <- setdiff(seq_len(n_tr), unique(idx_b))
      cols_b <- if (p_cov > 1L) sample(cov_cols, size = mtry, replace = FALSE) else cov_cols
      sub_df <- df_tr[idx_b, cols_b, drop = FALSE]
      sub_df$y_resp <- y_tr[idx_b]

      tr_b <- rpart::rpart(
        y_resp ~ ., data = sub_df, method = "anova",
        control = rpart::rpart.control(cp = 0.005, minsplit = min_sp, maxdepth = 6L, xval = 0L)
      )
      if (length(oob_b) > 0L) {
        p_oob <- as.numeric(predict(tr_b, newdata = df_tr[oob_b, cols_b, drop = FALSE]))
        oob_sum[oob_b] <- oob_sum[oob_b] + p_oob
        oob_cnt[oob_b] <- oob_cnt[oob_b] + 1L
      }
      ev_sum <- ev_sum + as.numeric(predict(tr_b, newdata = df_ev[, cols_b, drop = FALSE]))
      if (!is.null(tr_b$variable.importance)) {
        for (vn in names(tr_b$variable.importance)) {
          if (vn %in% names(imp_vec)) imp_vec[vn] <- imp_vec[vn] + as.numeric(tr_b$variable.importance[vn])
        }
      }
    }

    pred_eval   <- ev_sum / ntree
    pred_oof_tr <- ifelse(oob_cnt > 0L, oob_sum / pmax(1L, oob_cnt), pred_eval[1:n_tr])
    if (sum(imp_vec) > 0) imp_vec <- 100.0 * imp_vec / sum(imp_vec)
    for (cn in cov_cols) {
      info_rows[[cn]] <- data.frame(
        Parameter = paste0("RF_VarImp(%):", cn),
        Estimate  = round(unname(imp_vec[cn]), 2),
        CI_Lower  = NA_real_,
        CI_Upper  = NA_real_
      )
    }

  } else if (mode == "gbm") {
    set.seed(seed)
    df_tr <- as.data.frame(X_tr[, cov_cols, drop = FALSE])
    df_ev <- as.data.frame(X_eval[, cov_cols, drop = FALSE])
    p_cov <- length(cov_cols)
    ntree <- max(20L, min(250L, as.integer(n_trees)))
    lr    <- 0.06
    min_sp <- max(6L, floor(n_tr / 12L))

    fit_gbm_seq <- function(df_fit, y_fit, df_pred_list, n_rounds) {
      n_f  <- nrow(df_fit)
      mu_0 <- mean(y_fit)
      res  <- y_fit - mu_0
      preds_hist <- lapply(df_pred_list, function(d) matrix(mu_0, nrow = nrow(d), ncol = n_rounds))
      imp_acc <- setNames(numeric(p_cov), cov_cols)
      cur_p   <- lapply(df_pred_list, function(d) rep(mu_0, nrow(d)))

      for (m_i in seq_len(n_rounds)) {
        sub_idx <- if (n_f >= 15L) sample.int(n_f, size = max(10L, floor(0.80 * n_f)), replace = FALSE) else seq_len(n_f)
        sub_d   <- df_fit[sub_idx, , drop = FALSE]
        sub_d$r_target <- res[sub_idx]
        tr_m <- rpart::rpart(
          r_target ~ ., data = sub_d, method = "anova",
          control = rpart::rpart.control(cp = 0.002, minsplit = min_sp, maxdepth = 3L, xval = 0L)
        )
        step_f <- lr * as.numeric(predict(tr_m, newdata = df_fit))
        res    <- res - step_f
        if (!is.null(tr_m$variable.importance)) {
          for (vn in names(tr_m$variable.importance)) {
            if (vn %in% names(imp_acc)) imp_acc[vn] <- imp_acc[vn] + as.numeric(tr_m$variable.importance[vn])
          }
        }
        for (k_l in seq_along(df_pred_list)) {
          cur_p[[k_l]] <- cur_p[[k_l]] + lr * as.numeric(predict(tr_m, newdata = df_pred_list[[k_l]]))
          preds_hist[[k_l]][, m_i] <- cur_p[[k_l]]
        }
      }
      list(preds_hist = preds_hist, imp = imp_acc)
    }

    # 5-fold cross-fitting on training stations to select optimal boosting iteration & unbiased OOF predictions
    k_oof <- min(5L, n_tr)
    f_id  <- sample(rep(seq_len(k_oof), length.out = n_tr))
    oof_mat <- matrix(mean(y_tr), nrow = n_tr, ncol = ntree)
    for (kf in seq_len(k_oof)) {
      val_i <- which(f_id == kf)
      trn_i <- which(f_id != kf)
      if (length(val_i) == 0L || length(trn_i) < 6L) next
      g_kf <- fit_gbm_seq(df_tr[trn_i, , drop = FALSE], y_tr[trn_i], list(df_tr[val_i, , drop = FALSE]), ntree)
      oof_mat[val_i, ] <- g_kf$preds_hist[[1]]
    }
    mse_by_iter <- colMeans((oof_mat - y_tr)^2)
    best_iter   <- which.min(mse_by_iter)
    pred_oof_tr <- as.numeric(oof_mat[, best_iter])

    full_gbm  <- fit_gbm_seq(df_tr, y_tr, list(df_ev), best_iter)
    pred_eval <- as.numeric(full_gbm$preds_hist[[1]][, best_iter])
    imp_vec   <- full_gbm$imp
    if (sum(imp_vec) > 0) imp_vec <- 100.0 * imp_vec / sum(imp_vec)
    for (cn in cov_cols) {
      info_rows[[cn]] <- data.frame(
        Parameter = paste0("GBM_VarImp(%):", cn),
        Estimate  = round(unname(imp_vec[cn]), 2),
        CI_Lower  = NA_real_,
        CI_Upper  = NA_real_
      )
    }
    info_rows[["iter"]] <- data.frame(
      Parameter = "GBM_Best_Trees",
      Estimate  = as.numeric(best_iter),
      CI_Lower  = NA_real_,
      CI_Upper  = NA_real_
    )
  }

  # Orthogonal decomposition: separate nonlinear/interaction surface g_perp orthogonal to linear X_tr
  beta_proj_oof <- solve_linear_ols(X_tr, pred_oof_tr)
  beta_proj_ev  <- solve_linear_ols(X_tr, pred_eval[1:n_tr])
  g_perp_tr     <- as.numeric(pred_oof_tr - as.numeric(X_tr %*% beta_proj_oof))
  g_perp_eval   <- as.numeric(pred_eval   - as.numeric(X_eval %*% beta_proj_ev))
  # Use OOF nonlinear component at the observed reaches 1..n_tr to prevent residual overfitting
  g_perp_eval[1:n_tr] <- g_perp_tr

  info_df <- if (length(info_rows) > 0L) do.call(rbind, info_rows) else NULL
  list(
    mode        = mode,
    label       = S3N_TREND_SPECS[[mode]]$label,
    g_perp_tr   = g_perp_tr,
    g_perp_eval = g_perp_eval,
    pred_oof_tr = pred_oof_tr,
    pred_eval   = pred_eval,
    info_df     = info_df
  )
}

compute_reach_log_afv <- function(df, afv_scheme = "length", parent_0 = NULL) {
  scheme <- parse_afv_arg(paste0("--afv=", afv_scheme))
  if (scheme == "length" && ("log_AFV" %in% names(df)) && all(is.finite(df$log_AFV))) {
    return(as.double(df$log_AFV))
  }
  if (is.null(parent_0)) {
    par_idx  <- match(df$parent_dn, df$RivID, nomatch = 0L)
    parent_0 <- ifelse(df$parent_dn == 0 | par_idx == 0L, -1L, par_idx - 1L)
  }
  if (scheme == "length") {
    w_vec <- if ("US_Accum_km" %in% names(df)) pmax(0.05, as.double(df$US_Accum_km)) else rep(1.0, nrow(df))
    return(.Call("C_compute_log_afv", as.integer(parent_0), as.double(w_vec), 1L))
  }
  if (scheme == "shreve") {
    w_vec <- if ("Shreve" %in% names(df)) pmax(1.0, as.double(df$Shreve)) else rep(1.0, nrow(df))
    return(.Call("C_compute_log_afv", as.integer(parent_0), as.double(w_vec), 1L))
  }
  # hydro_yield: upstream length * orographic elevation factor * aquifer baseflow factor
  us_km <- if ("US_Accum_km" %in% names(df)) pmax(0.05, as.double(df$US_Accum_km)) else rep(1.0, nrow(df))
  elev  <- if ("Elevation" %in% names(df)) pmax(0.0, as.double(df$Elevation)) else rep(0.0, nrow(df))
  aq    <- if ("Is_Aquifer" %in% names(df)) pmax(0.0, as.double(df$Is_Aquifer)) else rep(0.0, nrow(df))
  elev[!is.finite(elev)] <- 0.0
  aq[!is.finite(aq)]     <- 0.0
  w_vec <- us_km * (1.0 + 0.002 * elev) * (1.0 + 0.50 * aq)
  .Call("C_compute_log_afv", as.integer(parent_0), as.double(w_vec), 1L)
}

fit_nngp_cpp <- function(X_obs, y_vec_obs, obs_nb, vy, kernel_name = "exponential") {
  if (kernel_name == "auto") {
    cand_names <- names(S3N_KERNEL_SPECS)
    fits <- lapply(cand_names, function(kn) fit_nngp_cpp(X_obs, y_vec_obs, obs_nb, vy, kernel_name = kn))
    aics <- vapply(fits, function(f) f$aic, numeric(1))
    best_idx <- which.min(aics)
    best_fit <- fits[[best_idx]]
    best_fit$auto_aic_table <- setNames(round(aics, 2), cand_names)
    return(best_fit)
  }
  sp <- S3N_KERNEL_SPECS[[kernel_name]]
  kernel_id <- as.integer(sp$id)
  n_obs <- nrow(X_obs)
  p_cov <- ncol(X_obs)
  y_dbl <- as.double(y_vec_obs)

  eval_obj <- function(par) {
    eta <- exp(max(-12.0, min(6.0, par[1])))
    phi <- exp(max(-9.0,  min(4.0, par[2])))
    res <- .Call("C_s3n_nngp_loglik", X_obs, y_dbl, obs_nb, eta, phi, kernel_id)
    if (!is.finite(res$neg2loglik)) return(1e30)
    res$neg2loglik
  }

  grid_eta <- log(c(0.05, 0.25, 0.80))
  grid_phi <- log(c(0.015, 0.06, 0.25))
  best_val <- Inf
  best_par <- c(log(0.15), log(0.05))
  for (ge in grid_eta) {
    for (gp in grid_phi) {
      val <- eval_obj(c(ge, gp))
      if (val < best_val) {
        best_val <- val
        best_par <- c(ge, gp)
      }
    }
  }

  opt <- tryCatch(
    optim(best_par, eval_obj, method = "Nelder-Mead", control = list(maxit = 180L, reltol = 1e-6)),
    error = function(e) list(par = best_par, value = best_val)
  )
  eta_hat <- exp(max(-12.0, min(6.0, opt$par[1])))
  phi_hat <- exp(max(-9.0,  min(4.0, opt$par[2])))
  final   <- .Call("C_s3n_nngp_loglik", X_obs, y_dbl, obs_nb, eta_hat, phi_hat, kernel_id)

  theta <- c(sigma.sq = final$sigma.sq, tau.sq = final$tau.sq, phi = phi_hat)
  names(theta) <- c("sigma.sq", "tau.sq", "phi")
  beta <- as.numeric(final$Beta)
  names(beta) <- colnames(X_obs)

  list(
    Beta        = beta,
    Theta       = theta,
    cov.beta    = final$cov.beta,
    neg2loglik  = final$neg2loglik,
    loglik      = -0.5 * final$neg2loglik,
    aic         = final$neg2loglik + 2.0 * (p_cov + 3.0),
    kernel_name = kernel_name,
    kernel_id   = kernel_id,
    engine      = "cpp_nngp"
  )
}

bootstrap_nngp_cpp <- function(fit_obj, X_obs) {
  se_beta <- sqrt(pmax(1e-12, diag(as.matrix(fit_obj$cov.beta))))
  ci_beta_lo <- fit_obj$Beta - 1.959964 * se_beta
  ci_beta_hi <- fit_obj$Beta + 1.959964 * se_beta
  n_obs <- nrow(X_obs)
  se_fac <- max(0.08, min(0.45, sqrt(2.0 / max(10.0, n_obs))))
  th <- fit_obj$Theta
  ci_th_lo <- th * exp(-1.959964 * se_fac)
  ci_th_hi <- th * exp( 1.959964 * se_fac)
  ci_mat <- rbind(c(ci_beta_lo, ci_th_lo), c(ci_beta_hi, ci_th_hi))
  colnames(ci_mat) <- c(colnames(X_obs), "sigma.sq", "tau.sq", "phi")
  list(confidence.interval = ci_mat)
}

fit_and_predict_joint_s3n <- function(X_obs, y_vec_obs, v_meas, has_me, X_all,
                                      nb_tu, nb_td, mid_x, mid_y, m_use,
                                      fit_tu, fit_td, kernel_id = 0L, estimate_only = FALSE,
                                      tau_init = NULL) {
  y_dbl <- as.double(y_vec_obs)
  v_use <- if (has_me) as.double(v_meas) else numeric(0)
  p_cov <- ncol(X_obs)
  q_tot <- nrow(X_all)

  phi_tu <- max(1e-5, as.numeric(fit_tu$Theta["phi"]))
  phi_td <- max(1e-5, as.numeric(fit_td$Theta["phi"]))

  # In Mode 0, C_s3n_joint_nngp profiles out overall variance scale when has_me == FALSE
  # and uses exact variance components when has_me == TRUE.
  vy <- max(1e-4, var(y_dbl))
  sig_tu_0 <- max(1e-5, as.numeric(fit_tu$Theta["sigma.sq"]))
  sig_td_0 <- max(1e-5, as.numeric(fit_td$Theta["sigma.sq"]))
  tau_0    <- if (!is.null(tau_init) && is.finite(tau_init)) {
    max(1e-5, as.numeric(tau_init))
  } else {
    max(1e-5, 0.5 * (as.numeric(fit_tu$Theta["tau.sq"]) + as.numeric(fit_td$Theta["tau.sq"])) - (if (has_me) mean(v_use) else 0.0))
  }

  eval_joint <- function(par) {
    r_td   <- exp(max(-8.0, min(8.0, par[1])))
    r_euc  <- exp(max(-10.0, min(6.0, par[2])))
    r_tau  <- exp(max(-10.0, min(6.0, par[3])))
    phi_eu <- exp(max(-7.0, min(3.0, par[4])))
    s_scale <- if (has_me) exp(max(-8.0, min(8.0, par[5]))) else 1.0
    params_vec <- c(s_scale * 1.0, s_scale * r_td, s_scale * r_euc, s_scale * r_tau,
                    phi_tu, phi_td, phi_eu, as.double(kernel_id))
    res <- .Call("C_s3n_joint_nngp",
                 X_obs, y_dbl, v_use, X_obs, numeric(p_cov),
                 nb_tu, nb_td, as.double(mid_x), as.double(mid_y),
                 as.integer(m_use), as.double(params_vec), 0L)
    val <- if (has_me) res$neg2loglik_exact else res$neg2loglik
    if (!is.finite(val)) return(1e30)
    val
  }

  r_td_init  <- log(max(0.05, min(20.0, sig_td_0 / sig_tu_0)))
  r_tau_init <- log(max(0.02, min(10.0, tau_0 / sig_tu_0)))
  s_init     <- log(max(1e-4, 0.5 * (sig_tu_0 + sig_td_0)))

  starts <- list(
    c(r_td_init,  log(0.10), r_tau_init, log(0.03), s_init),
    c(log(0.02),  log(0.01), log(max(0.02, as.numeric(fit_tu$Theta["tau.sq"]) / sig_tu_0)), log(0.03), log(sig_tu_0)),
    c(log(50.0),  log(0.01), log(max(0.02, as.numeric(fit_td$Theta["tau.sq"]) / max(1e-5, 0.02 * sig_td_0))), log(0.03), log(max(1e-5, 0.02 * sig_td_0))),
    c(0.0,        log(0.25), r_tau_init, log(0.05), s_init)
  )
  best_par <- starts[[1]]
  best_val <- Inf
  for (st in starts) {
    p_st <- if (has_me) st else st[1:4]
    v_st <- eval_joint(p_st)
    if (v_st < best_val) {
      best_val <- v_st
      best_par <- p_st
    }
  }

  opt <- tryCatch(
    optim(best_par, eval_joint, method = "Nelder-Mead", control = list(maxit = 140L, reltol = 1e-5)),
    error = function(e) list(par = best_par, value = best_val)
  )
  par_hat <- opt$par
  r_td   <- exp(max(-8.0, min(8.0, par_hat[1])))
  r_euc  <- exp(max(-10.0, min(6.0, par_hat[2])))
  r_tau  <- exp(max(-10.0, min(6.0, par_hat[3])))
  phi_eu <- exp(max(-7.0, min(3.0, par_hat[4])))
  s_scale <- if (has_me) exp(max(-8.0, min(8.0, par_hat[5]))) else 1.0

  params_raw <- c(s_scale * 1.0, s_scale * r_td, s_scale * r_euc, s_scale * r_tau,
                  phi_tu, phi_td, phi_eu, as.double(kernel_id))
  fit_res <- .Call("C_s3n_joint_nngp",
                   X_obs, y_dbl, v_use, X_obs, numeric(p_cov),
                   nb_tu, nb_td, as.double(mid_x), as.double(mid_y),
                   as.integer(m_use), as.double(params_raw), 0L)

  mult <- if (has_me) 1.0 else fit_res$scale
  params_final <- c(
    sigma.sq_tu  = params_raw[1] * mult,
    sigma.sq_td  = params_raw[2] * mult,
    sigma.sq_euc = params_raw[3] * mult,
    tau.sq       = params_raw[4] * mult,
    phi_tu       = phi_tu,
    phi_td       = phi_td,
    phi_euc      = phi_eu,
    kernel_id    = as.double(kernel_id)
  )
  beta_hat <- as.numeric(fit_res$Beta)
  names(beta_hat) <- colnames(X_obs)

  if (estimate_only) {
    return(list(Beta = beta_hat, params = params_final, neg2loglik = fit_res$neg2loglik))
  }

  pred_res <- .Call("C_s3n_joint_nngp",
                    X_obs, y_dbl, v_use, X_all, beta_hat,
                    nb_tu, nb_td, as.double(mid_x), as.double(mid_y),
                    as.integer(m_use), as.double(params_final), 1L)

  list(
    Beta       = beta_hat,
    params     = params_final,
    neg2loglik = fit_res$neg2loglik,
    prediction = pred_res$prediction,
    se         = pred_res$se,
    ci         = pred_res$prediction.ci
  )
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
                           kernel = "exponential",
                           combine_mode = "joint",
                           afv_scheme = "length",
                           trend = "linear",
                           spline_df = 3L,
                           n_trees = 100L,
                           cv_mode = "random",
                           threshold = NA_real_,
                           threshold_dir = "above",
                           design_sites = 0L,
                           censor_mode = "ros") {
  obs_idx <- which(df$has_obs == 1L & !is.na(df[[obs_col]]))
  n_obs   <- length(obs_idx)
  if (n_obs < 15L) {
    cat(sprintf("  Skipping separate fit for %s (only %d observed reaches).\n", domain_label, n_obs))
    return(NULL)
  }

  cat(sprintf("\n--- Fitting %s (%d reaches, %d observed) ---\n", domain_label, nrow(df), n_obs))

  obs_idx  <- obs_idx[order(df$CatchID[obs_idx], df$obs_up_dist_km[obs_idx])]
  pred_idx <- setdiff(seq_len(nrow(df)), obs_idx)
  df_ord   <- df[c(obs_idx, pred_idx), ]
  q_tot    <- nrow(df_ord)
  m_use    <- min(as.integer(m_neighbors), max(2L, n_obs - 2L))

  raw_obs_vals <- df_ord[[obs_col]][1:n_obs]
  y_vec_obs <- if ("obs_y" %in% names(df_ord) && all(is.finite(df_ord$obs_y[1:n_obs]))) {
    as.numeric(df_ord$obs_y[1:n_obs])
  } else {
    as.numeric(trans_obj$fwd(raw_obs_vals))
  }
  y_obs <- as.matrix(y_vec_obs)

  # Heteroskedastic measurement error (Approach A) & source-reliability weights
  v_phys      <- numeric(n_obs)
  v_bar_phys  <- 0.0
  has_phys_me <- FALSE
  if ("meas_uncertainty" %in% names(df_ord)) {
    u_raw <- df_ord$meas_uncertainty[1:n_obs]
    if (any(!is.na(u_raw) & u_raw > 0)) {
      v_phys <- calc_response_meas_variance(raw_obs_vals, u_raw, trans_obj)
      has_gam <- ("meas_method" %in% names(df_ord)) && isTRUE(any(tolower(df_ord$meas_method[1:n_obs]) == "gam", na.rm = TRUE))
      has_meta_gam <- !is.null(net_meta) && !is.null(net_meta$meas_method) && identical(tolower(net_meta$meas_method), "gam")
      has_phys_me <- isTRUE((mean(v_phys) > 1e-8) && (isTRUE(use_meas_error) || has_gam || has_meta_gam))
      if (!has_phys_me) v_phys <- numeric(n_obs)
    }
  }
  w_obs <- if ("obs_source_weight" %in% names(df_ord)) {
    w_tmp <- suppressWarnings(as.numeric(df_ord$obs_source_weight[1:n_obs]))
    ifelse(is.finite(w_tmp) & w_tmp > 0, pmin(1.0, pmax(1e-4, w_tmp)), 1.0)
  } else {
    rep(1.0, n_obs)
  }
  has_sw <- isTRUE(any(w_obs < (1.0 - 1e-6)))
  v_meas <- v_phys
  has_me <- isTRUE(has_phys_me || has_sw)

  if (has_phys_me) {
    v_bar_phys <- mean(v_phys)
    mean_se    <- mean(sqrt(v_phys))
    cat("  [Approach A: Measurement Error Kriging] Active: accounting for observation uncertainty\n")
    cat(sprintf("    Mean measurement SE: %.4f (%s scale), Mean error variance: %.6f\n",
                mean_se, trans_obj$short_label, v_bar_phys))
  }
  if (has_sw) {
    n_dw <- sum(w_obs < (1.0 - 1e-6))
    cat(sprintf("  [Source Reliability Weighting] Active: %d/%d reaches weighted < 1.0 (min weight = %.3f)\n",
                n_dw, n_obs, min(w_obs)))
  }

  # Integer topology lookup (avoids as.character string hashing on 90k+ reaches)
  par_idx  <- match(df_ord$parent_dn, df_ord$RivID, nomatch = 0L)
  parent_0 <- ifelse(df_ord$parent_dn == 0 | par_idx == 0L, -1L, par_idx - 1L)
  up_dist  <- ifelse(df_ord$has_obs == 1L, df_ord$obs_up_dist_km, df_ord$up_dist_mid_km)

  kernel_req <- normalize_kernel_name(kernel, default = "exponential")
  comb_mode  <- parse_combine_arg(paste0("--combine=", combine_mode))
  afv_mode   <- parse_afv_arg(paste0("--afv=", afv_scheme))
  trend_req  <- normalize_trend_name(trend, default = "linear")
  df_ord$log_AFV <- compute_reach_log_afv(df_ord, afv_scheme = afv_mode, parent_0 = parent_0)
  cat(sprintf("  Spatial configuration: kernel = %s | combine = %s | afv = %s | trend = %s\n",
              kernel_req, comb_mode, afv_mode, trend_req))

  X_all <- build_design_matrix(df_ord, include_peat = (domain_id == "england"), var_specs = var_specs)
  cat(sprintf("  Fixed-effect covariates (%d): %s\n", ncol(X_all), paste(colnames(X_all), collapse = ", ")))
  X_obs <- X_all[1:n_obs, , drop = FALSE]

  trend_fit <- fit_nonlinear_trend(
    X_tr = X_obs, y_tr_vec = y_vec_obs, X_eval = X_all,
    trend_mode = trend_req, spline_df = spline_df, n_trees = n_trees
  )
  if (trend_req != "linear") {
    cat(sprintf("  [Non-Linear / Ensemble Trend: %s] Active (orthogonal nonlinear SD = %.4f)\n",
                trend_fit$label, sd(trend_fit$g_perp_tr)))
  }
  y_vec_s3n <- y_vec_obs - trend_fit$g_perp_tr
  y_obs_s3n <- as.matrix(y_vec_s3n)
  vy        <- max(1e-4, var(y_vec_s3n))

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

  estimate_component <- function(X_m, y_m, obs_nb_m, vy_m, k_name) {
    if (k_name == "exponential" && isTRUE(.s3n_env$has_s3n_pkg)) {
      res_b <- BRISC_estimation_stream(
        coords = as.matrix(seq_len(nrow(X_m))), y = as.matrix(y_m), x = X_m,
        neighbor = obs_nb_m, sigma.sq = vy_m * 0.7, tau.sq = vy_m * 0.1, phi = 0.05,
        cov.model = "exponential", n_omp = n_omp, nugget_status = 1, verbose = FALSE
      )
      res_b$kernel_name <- "exponential"
      res_b$kernel_id   <- 0L
      res_b$engine      <- "brisc"
      return(res_b)
    }
    fit_nngp_cpp(X_m, as.double(y_m), obs_nb_m, vy_m, kernel_name = k_name)
  }

  em_res <- NULL
  censor_req <- tolower(trimws(as.character(censor_mode)))
  if (censor_req %in% c("em", "tobit") &&
      ("is_censored" %in% names(df_ord)) && ("lod_sc" %in% names(df_ord))) {
    is_cens_obs <- as.logical(df_ord$is_censored[1:n_obs] == 1L)
    if (any(is_cens_obs)) {
      lod_sc_obs <- as.numeric(df_ord$lod_sc[1:n_obs]) - trend_fit$g_perp_tr
      em_res <- run_spatial_tobit_em_kriging(
        df_ord = df_ord, n_obs = n_obs, X_obs = X_obs, y_vec_init = y_vec_s3n,
        is_cens_obs = is_cens_obs, lod_sc_obs = lod_sc_obs,
        m_use = m_use, parent_0 = parent_0, tt = tt, up_dist = up_dist,
        nb_td = nb_td, estimate_component = estimate_component,
        kernel_req = kernel_req
      )
      y_vec_s3n <- em_res$y_em
      y_obs_s3n <- as.matrix(y_vec_s3n)
      y_vec_obs <- y_vec_s3n + trend_fit$g_perp_tr
      vy        <- max(1e-4, var(y_vec_s3n))
      cat(sprintf("  [Spatial Tobit EM Kriging] Imputed %d left-censored reaches in %d EM iterations (converged=%s)\n",
                  em_res$n_cens, em_res$n_iter, as.character(em_res$converged)))
    }
  }

  t1 <- proc.time()
  fit_tu <- estimate_component(X_obs, y_obs_s3n, nb_tu$obs_neighbors, vy, kernel_req)
  fit_td <- estimate_component(X_obs, y_obs_s3n, nb_td$obs_neighbors, vy, kernel_req)
  if (kernel_req == "auto") {
    cat(sprintf("    [Auto-Kernel Selection] Tail-Up best = %s (%s)\n",
                fit_tu$kernel_name,
                paste(names(fit_tu$auto_aic_table), fit_tu$auto_aic_table, sep = "=", collapse = ", ")))
    cat(sprintf("    [Auto-Kernel Selection] Tail-Down best = %s (%s)\n",
                fit_td$kernel_name,
                paste(names(fit_td$auto_aic_table), fit_td$auto_aic_table, sep = "=", collapse = ", ")))
  }
  kid_tu <- if (!is.null(fit_tu$kernel_id)) as.integer(fit_tu$kernel_id) else 0L
  kid_td <- if (!is.null(fit_td$kernel_id)) as.integer(fit_td$kernel_id) else 0L
  t_est <- (proc.time() - t1)[3]

  if (has_me) {
    tau_micro_tu <- max(1e-6 * fit_tu$Theta["sigma.sq"], fit_tu$Theta["tau.sq"] - v_bar_phys)
    theta_tu_me  <- c(fit_tu$Theta["sigma.sq"], tau_micro_tu, fit_tu$Theta["phi"], kid_tu)
    tau_micro_td <- max(1e-6 * fit_td$Theta["sigma.sq"], fit_td$Theta["tau.sq"] - v_bar_phys)
    theta_td_me  <- c(fit_td$Theta["sigma.sq"], tau_micro_td, fit_td$Theta["phi"], kid_td)
    tau_base     <- max(1e-5 * 0.5 * (fit_tu$Theta["sigma.sq"] + fit_td$Theta["sigma.sq"]),
                        0.5 * (tau_micro_tu + tau_micro_td))
    v_wt         <- if (has_sw) (1.0 / w_obs - 1.0) * tau_base else numeric(n_obs)
    v_meas       <- v_phys + v_wt
    if (has_phys_me) {
      cat(sprintf("    Tail-Up Decomposition: Total Nugget = %.4f -> Microscale = %.4f + Meas Error = %.4f\n",
                  fit_tu$Theta["tau.sq"], tau_micro_tu, v_bar_phys))
      cat(sprintf("    Tail-Down Decomposition: Total Nugget = %.4f -> Microscale = %.4f + Meas Error = %.4f\n",
                  fit_td$Theta["tau.sq"], tau_micro_td, v_bar_phys))
    }
  }

  t2 <- proc.time()
  bootstrap_component <- function(fit_obj) {
    if (identical(fit_obj$engine, "brisc") && isTRUE(.s3n_env$has_s3n_pkg)) {
      BRISC_bootstrap(fit_obj, n_boot = n_boot, n_omp = n_omp, verbose = FALSE)
    } else {
      bootstrap_nngp_cpp(fit_obj, X_obs)
    }
  }
  bs_tu <- bootstrap_component(fit_tu)
  bs_td <- bootstrap_component(fit_td)
  t_bs  <- (proc.time() - t2)[3]

  t3 <- proc.time()
  predict_component <- function(fit_obj, nb_obj, theta_me, kid) {
    if (has_me) {
      .Call("C_s3n_predict_meas_error",
            X_obs, y_obs_s3n, v_meas, X_all,
            fit_obj$Beta, as.double(theta_me),
            as.integer(m_use),
            nb_obj$pred_neighbors$nnIndx, nb_obj$pred_neighbors$nnDist, nb_obj$pred_neighbors$nnWght,
            nb_obj$obsobs_dist, nb_obj$obsobs_wt)
    } else if (identical(fit_obj$engine, "brisc") && kid == 0L && isTRUE(.s3n_env$has_s3n_pkg)) {
      BRISC_prediction_stream(
        BRISC_Out = fit_obj, coords.0 = as.matrix(1:q_tot), X.0 = X_all,
        neighbor = nb_obj$pred_neighbors, obsobs_dist = nb_obj$obsobs_dist, obsobs_wt = nb_obj$obsobs_wt,
        n_omp = n_omp, verbose = FALSE
      )
    } else {
      th_vec <- c(fit_obj$Theta["sigma.sq"], fit_obj$Theta["tau.sq"], fit_obj$Theta["phi"], kid)
      .Call("C_s3n_predict_meas_error",
            X_obs, y_obs_s3n, numeric(0), X_all,
            fit_obj$Beta, as.double(th_vec),
            as.integer(m_use),
            nb_obj$pred_neighbors$nnIndx, nb_obj$pred_neighbors$nnDist, nb_obj$pred_neighbors$nnWght,
            nb_obj$obsobs_dist, nb_obj$obsobs_wt)
    }
  }
  pr_tu  <- predict_component(fit_tu, nb_tu, if (has_me) theta_tu_me else NULL, kid_tu)
  pr_td  <- predict_component(fit_td, nb_td, if (has_me) theta_td_me else NULL, kid_td)
  pr_tu$prediction    <- pr_tu$prediction + trend_fit$g_perp_eval
  pr_tu$prediction.ci <- pr_tu$prediction.ci + trend_fit$g_perp_eval
  pr_td$prediction    <- pr_td$prediction + trend_fit$g_perp_eval
  pr_td$prediction.ci <- pr_td$prediction.ci + trend_fit$g_perp_eval

  se_tu <- (pr_tu$prediction.ci[, 2] - pr_tu$prediction) / 1.959964
  se_td <- (pr_td$prediction.ci[, 2] - pr_td$prediction) / 1.959964
  sill_nug_tu <- sum(if (has_me) theta_tu_me[1:2] else fit_tu$Theta[1:2])
  sill_nug_td <- sum(if (has_me) theta_td_me[1:2] else fit_td$Theta[1:2])
  info_tu <- pmax(1e-4, 1 - (se_tu^2) / sill_nug_tu)
  info_td <- pmax(1e-4, 1 - (se_td^2) / sill_nug_td)

  joint_fit <- NULL
  if (comb_mode == "joint") {
    joint_fit <- fit_and_predict_joint_s3n(
      X_obs = X_obs, y_vec_obs = y_vec_s3n, v_meas = v_meas, has_me = has_me,
      X_all = X_all, nb_tu = nb_tu, nb_td = nb_td,
      mid_x = df_ord$mid_X, mid_y = df_ord$mid_Y, m_use = m_use,
      fit_tu = fit_tu, fit_td = fit_td, kernel_id = kid_tu,
      tau_init = if (has_me) tau_base else NULL
    )
    joint_fit$prediction <- joint_fit$prediction + trend_fit$g_perp_eval
    joint_fit$ci         <- joint_fit$ci + trend_fit$g_perp_eval
    pred_comb <- joint_fit$prediction
    se_comb   <- joint_fit$se
    cat(sprintf("    [Joint Additive S3N (TU + TD + Euc)] sigma2_TU=%.4f, sigma2_TD=%.4f, sigma2_Euc=%.4f, tau2=%.4f\n",
                joint_fit$params["sigma.sq_tu"], joint_fit$params["sigma.sq_td"],
                joint_fit$params["sigma.sq_euc"], joint_fit$params["tau.sq"]))
  } else {
    pred_comb <- (info_tu * pr_tu$prediction + info_td * pr_td$prediction) / (info_tu + info_td)
    se_comb   <- sqrt((info_tu^2 * se_tu^2 + info_td^2 * se_td^2) / ((info_tu + info_td)^2))
  }
  t_pred <- (proc.time() - t3)[3]

  obs_sc_full <- rep(NA_real_, q_tot)
  obs_sc_full[1:n_obs] <- round(y_vec_obs, 5)
  obs_val_full <- rep(NA_real_, q_tot)
  obs_val_full[1:n_obs] <- raw_obs_vals

  pred_df <- data.frame(
    RivID         = df_ord$RivID,
    CatchID       = df_ord$CatchID,
    parent_dn     = df_ord$parent_dn,
    Shreve        = df_ord$Shreve,
    Strahler      = df_ord$Strahler,
    Dist2Mth_km   = df_ord$Dist2Mth_km,
    US_Accum_km   = df_ord$US_Accum_km,
    start_X       = df_ord$start_X,
    start_Y       = df_ord$start_Y,
    mid_X         = df_ord$mid_X,
    mid_Y         = df_ord$mid_Y,
    end_X         = df_ord$end_X,
    end_Y         = df_ord$end_Y,
    Elevation     = df_ord$Elevation,
    Is_Aquifer    = df_ord$Is_Aquifer,
    has_obs       = df_ord$has_obs,
    obs_value     = obs_val_full,
    obs_model_sc  = obs_sc_full,
    pred_tu_sc    = round(pr_tu$prediction, 5),
    se_tu_sc      = round(se_tu, 5),
    pred_td_sc    = round(pr_td$prediction, 5),
    se_td_sc      = round(se_td, 5),
    pred_model_sc = round(pred_comb, 5),
    se_model_sc   = round(se_comb, 5),
    pred_tu       = round(trans_obj$inv(pr_tu$prediction), 4),
    pred_td       = round(trans_obj$inv(pr_td$prediction), 4),
    pred_combined = round(trans_obj$inv(pred_comb), 4)
  )
  if ("meas_uncertainty" %in% names(df_ord)) pred_df$meas_uncertainty <- df_ord$meas_uncertainty
  if ("meas_method" %in% names(df_ord))      pred_df$meas_method      <- df_ord$meas_method
  if (trend_req != "linear") {
    pred_df$pred_trend_sc <- round(as.numeric(X_all %*% fit_tu$Beta) + trend_fit$g_perp_eval, 5)
  }

  exceed_summary <- NULL
  if (!is.na(threshold) && is.finite(threshold)) {
    thresh_dir_use <- if (tolower(as.character(threshold_dir)) == "below") "below" else "above"
    ex_res <- compute_exceedance_risk(pred_comb, se_comb, threshold, trans_obj, dir = thresh_dir_use)
    pred_df$exceed_prob     <- round(ex_res$prob, 5)
    pred_df$compliance_tier <- ex_res$tier
    tier_levels <- c(
      "Confidently Compliant",
      "Likely Compliant",
      "Uncertain / Priority Monitoring",
      "Likely Exceeding",
      "Confidently Exceeding"
    )
    tier_tab <- table(factor(ex_res$tier, levels = tier_levels))
    exceed_summary <- data.frame(
      Domain        = domain_label,
      Threshold_Raw = threshold,
      Threshold_Sc  = round(ex_res$thresh_sc, 5),
      Direction     = thresh_dir_use,
      Tier          = names(tier_tab),
      Reach_Count   = as.integer(tier_tab),
      Reach_Pct     = round(100 * as.numeric(tier_tab) / q_tot, 2)
    )
    cat(sprintf("  [Regulatory Exceedance Risk] Threshold = %.4g (%s, dir=%s):\n",
                threshold, trans_obj$short_label, thresh_dir_use))
    for (ti in seq_len(nrow(exceed_summary))) {
      cat(sprintf("    %-33s : %6d reaches (%5.1f%%)\n",
                  exceed_summary$Tier[ti], exceed_summary$Reach_Count[ti], exceed_summary$Reach_Pct[ti]))
    }
  }

  design_df <- NULL
  if (!is.na(design_sites) && as.integer(design_sites) > 0L) {
    params_5 <- if (!is.null(joint_fit)) {
      as.numeric(c(
        joint_fit$params["sigma.sq_tu"] + joint_fit$params["sigma.sq_td"],
        joint_fit$params["sigma.sq_euc"],
        joint_fit$params["tau.sq"],
        mean(c(joint_fit$params["phi_tu"], joint_fit$params["phi_td"])),
        as.double(kid_tu)
      ))
    } else {
      as.numeric(c(
        0.5 * (fit_tu$Theta[1] + fit_td$Theta[1]),
        0.0,
        0.5 * (fit_tu$Theta[2] + fit_td$Theta[2]),
        0.5 * (fit_tu$Theta[3] + fit_td$Theta[3]),
        as.double(kid_tu)
      ))
    }
    t_des <- proc.time()
    design_df <- compute_optimal_monitoring_design(
      df_ord = df_ord, n_obs = n_obs, se_comb = se_comb,
      parent_0 = parent_0, tt = tt, up_dist = up_dist,
      params_5 = params_5, k_design = as.integer(design_sites)
    )
    if (!is.null(design_df) && nrow(design_df) > 0L) {
      design_df <- cbind(Domain = domain_label, design_df)
      cat(sprintf("  [Optimal Monitoring Design] Ranked top %d sites in %.2fs (cumulative kriging variance reduction: %.2f%%)\n",
                  nrow(design_df), (proc.time() - t_des)[3], tail(design_df$Cum_Pct_Reduced, 1)))
    }
  }

  # Cross-Validation (random, catchment, strahler, or spatial_block)
  cv_info  <- build_cv_folds(df_ord[1:n_obs, , drop = FALSE], cv_mode = cv_mode, k_folds = 5L, seed = 2026L)
  folds    <- cv_info$folds
  if (!is.null(cv_info$note)) {
    cat(sprintf("  [CV Mode: %s -> %s] Note: %s\n", cv_mode, cv_info$actual_mode, cv_info$note))
  } else {
    cat(sprintf("  [CV Mode: %s] Evaluating held-out predictions...\n", cv_info$actual_mode))
  }
  cv_ols   <- rep(NA_real_, n_obs)
  cv_trend <- rep(NA_real_, n_obs)
  cv_tu    <- rep(NA_real_, n_obs)
  cv_td    <- rep(NA_real_, n_obs)
  cv_comb  <- rep(NA_real_, n_obs)
  y_vec    <- y_vec_obs

  riv_vec   <- df_ord$RivID
  par_vec   <- df_ord$parent_dn
  cid_vec   <- df_ord$CatchID
  d2m_vec   <- df_ord$Dist2Mth_km
  afv_vec   <- df_ord$log_AFV
  mx_vec    <- df_ord$mid_X
  my_vec    <- df_ord$mid_Y
  unobs_tail <- if (q_tot > n_obs) (n_obs + 1L):q_tot else integer(0)

  fold_ids <- sort(unique(folds[is.finite(folds) & folds > 0L]))
  for (k in fold_ids) {
    te <- which(folds == k)
    tr <- which(!is.na(folds) & folds != k)
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
    tf_cv  <- fit_nonlinear_trend(
      X_tr = X_tr, y_tr_vec = y_vec[tr], X_eval = X_eval,
      trend_mode = trend_req, spline_df = spline_df, n_trees = n_trees, seed = 2026L + k
    )
    cv_trend[te] <- tf_cv$pred_eval[(n_tr + 1L):n_eval]
    g_te_cv      <- tf_cv$g_perp_eval[(n_tr + 1L):n_eval]
    y_tr         <- as.matrix(y_vec[tr] - tf_cv$g_perp_tr)
    vy_cv        <- max(1e-4, var(as.numeric(y_tr)))

    run_cv_component <- function(model_int, k_name, kid) {
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
      f_cv <- estimate_component(X_tr, y_tr, nb_cv$obs_neighbors, vy_cv, k_name)
      if (has_me) {
        v_phys_tr  <- v_phys[tr]
        v_bar_tr   <- mean(v_phys_tr)
        tau_mic_cv <- max(1e-6 * f_cv$Theta["sigma.sq"], f_cv$Theta["tau.sq"] - v_bar_tr)
        tau_b_cv   <- max(1e-5 * f_cv$Theta["sigma.sq"], tau_mic_cv)
        v_wt_tr    <- if (has_sw) (1.0 / w_obs[tr] - 1.0) * tau_b_cv else numeric(length(tr))
        v_tr       <- v_phys_tr + v_wt_tr
        theta_me   <- c(f_cv$Theta["sigma.sq"], tau_mic_cv, f_cv$Theta["phi"], kid)
        pr_raw     <- .Call("C_s3n_predict_meas_error",
                            X_tr, y_tr, v_tr, X_eval,
                            f_cv$Beta, as.double(theta_me), as.integer(m_cv),
                            p_cv$nnIndx, p_cv$nnDist, p_cv$nnWght,
                            nb_cv$obsobs_dist, nb_cv$obsobs_wt)
        denom <- sum(theta_me[1:2])
      } else if (identical(f_cv$engine, "brisc") && kid == 0L && isTRUE(.s3n_env$has_s3n_pkg)) {
        pr_raw <- BRISC_prediction_stream(
          BRISC_Out = f_cv, coords.0 = as.matrix(1:n_eval), X.0 = X_eval,
          neighbor = p_cv, obsobs_dist = nb_cv$obsobs_dist, obsobs_wt = nb_cv$obsobs_wt,
          n_omp = n_omp, verbose = FALSE
        )
        denom <- sum(f_cv$Theta[1:2])
      } else {
        th_cv  <- c(f_cv$Theta["sigma.sq"], f_cv$Theta["tau.sq"], f_cv$Theta["phi"], kid)
        pr_raw <- .Call("C_s3n_predict_meas_error",
                        X_tr, y_tr, numeric(0), X_eval,
                        f_cv$Beta, as.double(th_cv), as.integer(m_cv),
                        p_cv$nnIndx, p_cv$nnDist, p_cv$nnWght,
                        nb_cv$obsobs_dist, nb_cv$obsobs_wt)
        denom <- sum(th_cv[1:2])
      }
      pred_te <- pr_raw$prediction[(n_tr + 1L):n_eval] + g_te_cv
      se_te   <- (pr_raw$prediction.ci[(n_tr + 1L):n_eval, 2] - pr_raw$prediction[(n_tr + 1L):n_eval]) / 1.959964
      inf_te  <- pmax(1e-4, 1 - (se_te^2) / denom)
      list(pred = pred_te, inf = inf_te, nb = nb_cv, fit = f_cv)
    }

    res_tu_cv  <- run_cv_component(0L, fit_tu$kernel_name, kid_tu)
    res_td_cv  <- run_cv_component(1L, fit_td$kernel_name, kid_td)
    cv_tu[te]  <- res_tu_cv$pred
    cv_td[te]  <- res_td_cv$pred
    if (comb_mode == "joint" && !is.null(joint_fit)) {
      nb_tu_eval <- res_tu_cv$nb
      nb_td_eval <- res_td_cv$nb
      nb_tu_eval[[4]] <- list(
        nnIndx = nb_tu_eval$pred_neighbors$nnIndx[1:(n_eval * m_cv)],
        nnDist = nb_tu_eval$pred_neighbors$nnDist[1:(n_eval * m_cv)],
        nnWght = nb_tu_eval$pred_neighbors$nnWght[1:(n_eval * m_cv)]
      )
      nb_td_eval[[4]] <- list(
        nnIndx = nb_td_eval$pred_neighbors$nnIndx[1:(n_eval * m_cv)],
        nnDist = nb_td_eval$pred_neighbors$nnDist[1:(n_eval * m_cv)],
        nnWght = nb_td_eval$pred_neighbors$nnWght[1:(n_eval * m_cv)]
      )
      v_tr_use <- if (has_me) as.double(v_meas[tr]) else numeric(0)
      jf_tr <- .Call("C_s3n_joint_nngp",
                     X_tr, as.double(y_tr), v_tr_use, X_tr, numeric(ncol(X_tr)),
                     nb_tu_eval, nb_td_eval, as.double(mx_cv[1:n_eval]), as.double(my_cv[1:n_eval]),
                     as.integer(m_cv), as.double(joint_fit$params), 0L)
      jp_ev <- .Call("C_s3n_joint_nngp",
                     X_tr, as.double(y_tr), v_tr_use, X_eval, as.double(jf_tr$Beta),
                     nb_tu_eval, nb_td_eval, as.double(mx_cv[1:n_eval]), as.double(my_cv[1:n_eval]),
                     as.integer(m_cv), as.double(joint_fit$params), 1L)
      cv_comb[te] <- jp_ev$prediction[(n_tr + 1L):n_eval] + g_te_cv
    } else {
      cv_comb[te] <- (res_tu_cv$inf * cv_tu[te] + res_td_cv$inf * cv_td[te]) / (res_tu_cv$inf + res_td_cv$inf)
    }
  }

  rmse <- function(a, b) { ok <- is.finite(a) & is.finite(b); if (!any(ok)) NA_real_ else sqrt(mean((a[ok] - b[ok])^2)) }
  mae  <- function(a, b) { ok <- is.finite(a) & is.finite(b); if (!any(ok)) NA_real_ else mean(abs(a[ok] - b[ok])) }
  r2   <- function(a, b) {
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 2L) return(NA_real_)
    ss_tot <- sum((a[ok] - mean(a[ok]))^2)
    if (!is.finite(ss_tot) || ss_tot < 1e-12) return(0.0)
    1 - sum((a[ok] - b[ok])^2) / ss_tot
  }

  sp_method_lbl <- if (has_me) "Measurement Error Kriging" else "Standard S3N"
  if (trend_req != "linear") {
    sp_method_lbl <- sprintf("%s [%s]", sp_method_lbl, toupper(trend_req))
    cat(sprintf("  [Non-Spatial Trend CV (%s)] RMSE = %.4f | MAE = %.4f | R2 = %.4f (vs Linear OLS R2 = %.4f)\n",
                trend_fit$label, rmse(y_vec, cv_trend), mae(y_vec, cv_trend), r2(y_vec, cv_trend), r2(y_vec, cv_ols)))
  }
  if (cv_info$actual_mode != "random") {
    sp_method_lbl <- sprintf("%s [CV:%s]", sp_method_lbl, cv_info$actual_mode)
  }
  if (!is.null(em_res) && em_res$n_cens > 0L) {
    sp_method_lbl <- sprintf("%s [Tobit-EM]", sp_method_lbl)
  }

  cv_summary <- data.frame(
    Domain         = domain_label,
    n_obs          = n_obs,
    n_reaches      = q_tot,
    CV_Mode        = cv_info$actual_mode,
    Model          = c("OLS (Non-Spatial)", "S3N Tail-Down (Stream Distance)", "S3N Tail-Up (Flow-Connected + AFV)", "S3N Combined (Tail-Up + Tail-Down)"),
    Spatial_Method = sp_method_lbl,
    RMSE           = round(c(rmse(y_vec, cv_ols), rmse(y_vec, cv_td), rmse(y_vec, cv_tu), rmse(y_vec, cv_comb)), 4),
    MAE            = round(c(mae(y_vec, cv_ols),  mae(y_vec, cv_td),  mae(y_vec, cv_tu),  mae(y_vec, cv_comb)), 4),
    R2             = round(c(r2(y_vec, cv_ols),   r2(y_vec, cv_td),   r2(y_vec, cv_tu),   r2(y_vec, cv_comb)), 4)
  )
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
  if (!is.null(joint_fit)) {
    joint_pnames <- c(colnames(X_obs), "sigma.sq_tu", "sigma.sq_td", "sigma.sq_euc", "tau.sq", "phi_tu", "phi_td", "phi_euc")
    joint_pvals  <- c(joint_fit$Beta, joint_fit$params[c("sigma.sq_tu", "sigma.sq_td", "sigma.sq_euc", "tau.sq", "phi_tu", "phi_td", "phi_euc")])
    param_df <- rbind(
      param_df,
      data.frame(Domain = domain_label, Model = "S3N Joint (TU + TD + Euc)", Parameter = joint_pnames,
                 Estimate = round(as.numeric(joint_pvals), 4), CI_Lower = NA, CI_Upper = NA)
    )
  }
  if (trend_req != "linear") {
    trend_mod_name <- sprintf("S3N Trend (%s)", trend_fit$label)
    if (!is.null(trend_fit$info_df) && nrow(trend_fit$info_df) > 0L) {
      t_info <- trend_fit$info_df
      param_df <- rbind(
        param_df,
        data.frame(Domain = domain_label, Model = trend_mod_name,
                   Parameter = as.character(t_info$Parameter),
                   Estimate  = as.numeric(t_info$Estimate),
                   CI_Lower  = as.numeric(t_info$CI_Lower),
                   CI_Upper  = as.numeric(t_info$CI_Upper))
      )
    }
    param_df <- rbind(
      param_df,
      data.frame(Domain = domain_label, Model = trend_mod_name,
                 Parameter = c("Trend_CV_RMSE", "Trend_CV_R2"),
                 Estimate  = round(c(rmse(y_vec, cv_trend), r2(y_vec, cv_trend)), 4),
                 CI_Lower  = NA_real_, CI_Upper = NA_real_)
    )
  }
  if (has_me) {
    param_df <- rbind(
      param_df,
      data.frame(Domain = domain_label, Model = "S3N Tail-Up (FC + AFV)", Parameter = "tau.sq (microscale)",
                 Estimate = round(tau_micro_tu, 4), CI_Lower = NA, CI_Upper = NA),
      data.frame(Domain = domain_label, Model = "S3N Tail-Down (Stream)", Parameter = "tau.sq (microscale)",
                 Estimate = round(tau_micro_td, 4), CI_Lower = NA, CI_Upper = NA),
      data.frame(Domain = domain_label, Model = "Measurement Error", Parameter = "mean_meas_variance",
                 Estimate = round(v_bar_phys, 6), CI_Lower = NA, CI_Upper = NA)
    )
    if (has_sw) {
      param_df <- rbind(
        param_df,
        data.frame(Domain = domain_label, Model = "Source Weighting", Parameter = "mean_source_penalty_var",
                   Estimate = round(mean(v_wt), 6), CI_Lower = NA, CI_Upper = NA)
      )
    }
  }
  if (!is.null(em_res) && em_res$n_cens > 0L) {
    param_df <- rbind(
      param_df,
      data.frame(Domain = domain_label, Model = "Spatial Tobit EM",
                 Parameter = c("Tobit_EM_Censored_N", "Tobit_EM_Iterations"),
                 Estimate  = c(em_res$n_cens, em_res$n_iter),
                 CI_Lower  = NA_real_, CI_Upper = NA_real_)
    )
  }

  cv_pts <- data.frame(
    RivID   = df_ord$RivID[1:n_obs],
    CatchID = df_ord$CatchID[1:n_obs],
    obs_y   = y_vec,
    cv_ols  = cv_ols,
    cv_tu   = cv_tu,
    cv_td   = cv_td,
    cv_comb = cv_comb
  )

  list(
    pred_df        = pred_df,
    cv_summary     = cv_summary,
    param_df       = param_df,
    cv_pts         = cv_pts,
    exceed_summary = exceed_summary,
    design_df      = design_df,
    runtimes       = c(neighbor_sec = t_nb, est_sec = t_est, boot_sec = t_bs, pred_sec = t_pred)
  )
}
