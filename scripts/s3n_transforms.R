# ==============================================================================
# s3n_transforms.R
#
# Shared data transformations (--transform=...), Box-Cox & Yeo-Johnson MLE,
# delta-method response measurement-error variance propagation, normalisation,
# and generic/legacy schema alias resolution.
# Automatically sourced by scripts/s3n_plot_functions.R.
# ==============================================================================

estimate_boxcox_lambda <- function(y, lower = -2.0, upper = 2.0) {
  y_pos <- y[!is.na(y) & y > 0]
  n <- length(y_pos)
  if (n < 10) return(0.5)
  log_y <- log(y_pos)
  sum_log_y <- sum(log_y)

  ll_fn <- function(lam) {
    z <- if (abs(lam) < 1e-4) log_y else ((y_pos^lam) - 1.0) / lam
    z_mean <- mean(z)
    var_z  <- mean((z - z_mean)^2)
    if (var_z < 1e-12) return(-1e12)
    -0.5 * n * log(var_z) + (lam - 1.0) * sum_log_y
  }

  opt <- tryCatch(
    optimize(ll_fn, interval = c(lower, upper), maximum = TRUE),
    error = function(e) list(maximum = 0.5)
  )
  round(opt$maximum, 4)
}

estimate_yeojohnson_lambda <- function(y, lower = -2.0, upper = 2.0) {
  y_val <- y[!is.na(y) & is.finite(y)]
  n <- length(y_val)
  if (n < 10) return(0.5)
  pos <- (y_val >= 0)
  neg <- (!pos)

  ll_fn <- function(lam) {
    z <- numeric(n)
    if (any(pos)) {
      yp <- y_val[pos]
      z[pos] <- if (abs(lam) < 1e-4) log(yp + 1.0) else (((yp + 1.0)^lam) - 1.0) / lam
    }
    if (any(neg)) {
      yn <- y_val[neg]
      z[neg] <- if (abs(2.0 - lam) < 1e-4) -log(-yn + 1.0) else -(((-yn + 1.0)^(2.0 - lam)) - 1.0) / (2.0 - lam)
    }
    z_mean <- mean(z)
    var_z  <- mean((z - z_mean)^2)
    if (var_z < 1e-12) return(-1e12)
    jac <- sum(log(abs(y_val) + 1.0) * ifelse(pos, lam - 1.0, 1.0 - lam))
    -0.5 * n * log(var_z) + jac
  }

  opt <- tryCatch(
    optimize(ll_fn, interval = c(lower, upper), maximum = TRUE),
    error = function(e) list(maximum = 0.5)
  )
  round(opt$maximum, 4)
}

build_transformer <- function(transform_type = "log", y_obs = NULL, lambda = NULL, scale = NULL) {
  t_raw <- if (is.null(transform_type) || length(transform_type) == 0L || is.na(transform_type[1])) "" else tolower(trimws(as.character(transform_type[1])))
  if (!nzchar(t_raw)) t_raw <- "log"

  alias_map <- c(
    ln = "log", log = "log",
    log10 = "log10", log_10 = "log10",
    log1p = "log1p", log_1p = "log1p",
    sqrt = "sqrt", square_root = "sqrt", squareroot = "sqrt",
    cbrt = "cbrt", cube_root = "cbrt", cuberoot = "cbrt",
    fourthrt = "fourthrt", fourth_root = "fourthrt", qtrt = "fourthrt",
    reciprocal = "reciprocal", inv = "reciprocal", inverse = "reciprocal",
    boxcox = "boxcox", box_cox = "boxcox", bc = "boxcox",
    yeojohnson = "yeojohnson", yeo_johnson = "yeojohnson", yj = "yeojohnson",
    asin = "asin", arcsin = "asin", arcsine = "asin", asin_sqrt = "asin",
    logit = "logit",
    none = "none", raw = "none", identity = "none"
  )
  if (!(t_raw %in% names(alias_map))) {
    stop(sprintf("Unknown transformation '%s'. Run --list-transforms to see supported options.", transform_type))
  }
  ttype <- unname(alias_map[t_raw])

  if (ttype == "log") {
    return(list(
      type           = "log",
      label          = "Natural Logarithm (ln)",
      short_label    = "log",
      fn_prefix      = "log",
      folder_tag     = "log",
      is_transformed = TRUE,
      is_log         = TRUE,
      lambda         = NA_real_,
      scale          = 1.0,
      valid_fn       = function(y) !is.na(y) & y > 0,
      fwd            = function(y) log(pmax(1e-9, y)),
      inv            = function(z) exp(z)
    ))
  }

  if (ttype == "log10") {
    return(list(
      type           = "log10",
      label          = "Base-10 Logarithm (log10)",
      short_label    = "log10",
      fn_prefix      = "log10",
      folder_tag     = "log10",
      is_transformed = TRUE,
      is_log         = FALSE,
      lambda         = NA_real_,
      scale          = 1.0,
      valid_fn       = function(y) !is.na(y) & y > 0,
      fwd            = function(y) log10(pmax(1e-9, y)),
      inv            = function(z) 10^z
    ))
  }

  if (ttype == "log1p") {
    return(list(
      type           = "log1p",
      label          = "Shifted Logarithm ln(1+y)",
      short_label    = "log1p",
      fn_prefix      = "log1p",
      folder_tag     = "log1p",
      is_transformed = TRUE,
      is_log         = FALSE,
      lambda         = NA_real_,
      scale          = 1.0,
      valid_fn       = function(y) !is.na(y) & y >= 0,
      fwd            = function(y) log1p(pmax(0.0, y)),
      inv            = function(z) pmax(0.0, expm1(z))
    ))
  }

  if (ttype == "sqrt") {
    return(list(
      type           = "sqrt",
      label          = "Square Root (sqrt)",
      short_label    = "sqrt",
      fn_prefix      = "sqrt",
      folder_tag     = "sqrt",
      is_transformed = TRUE,
      is_log         = FALSE,
      lambda         = 0.5,
      scale          = 1.0,
      valid_fn       = function(y) !is.na(y) & y >= 0,
      fwd            = function(y) sqrt(pmax(0.0, y)),
      inv            = function(z) (pmax(0.0, z))^2
    ))
  }

  if (ttype == "cbrt") {
    return(list(
      type           = "cbrt",
      label          = "Cube Root (cbrt)",
      short_label    = "cbrt",
      fn_prefix      = "cbrt",
      folder_tag     = "cbrt",
      is_transformed = TRUE,
      is_log         = FALSE,
      lambda         = 1.0 / 3.0,
      scale          = 1.0,
      valid_fn       = function(y) !is.na(y),
      fwd            = function(y) sign(y) * (abs(y))^(1.0 / 3.0),
      inv            = function(z) z^3
    ))
  }

  if (ttype == "fourthrt") {
    return(list(
      type           = "fourthrt",
      label          = "Fourth Root (fourthrt)",
      short_label    = "fourthrt",
      fn_prefix      = "fourthrt",
      folder_tag     = "fourthrt",
      is_transformed = TRUE,
      is_log         = FALSE,
      lambda         = 0.25,
      scale          = 1.0,
      valid_fn       = function(y) !is.na(y) & y >= 0,
      fwd            = function(y) (pmax(0.0, y))^0.25,
      inv            = function(z) (pmax(0.0, z))^4
    ))
  }

  if (ttype == "reciprocal") {
    return(list(
      type           = "reciprocal",
      label          = "Negative Reciprocal (-1/y)",
      short_label    = "reciprocal",
      fn_prefix      = "recip",
      folder_tag     = "reciprocal",
      is_transformed = TRUE,
      is_log         = FALSE,
      lambda         = -1.0,
      scale          = 1.0,
      valid_fn       = function(y) !is.na(y) & y > 0,
      fwd            = function(y) -1.0 / pmax(1e-6, y),
      inv            = function(z) -1.0 / pmin(-1e-6, z)
    ))
  }

  if (ttype == "boxcox") {
    lam_val <- if (!is.null(lambda) && !is.na(as.numeric(lambda))) {
      as.numeric(lambda)
    } else if (!is.null(y_obs) && length(y_obs) > 0) {
      estimate_boxcox_lambda(y_obs)
    } else {
      0.5
    }
    lbl <- sprintf("Box-Cox (lambda = %.2f)", lam_val)
    f_tag <- sprintf("boxcox_lam%.2f", lam_val)
    return(list(
      type           = "boxcox",
      label          = lbl,
      short_label    = sprintf("boxcox[%.2f]", lam_val),
      fn_prefix      = sprintf("boxcox[lambda=%.2f]", lam_val),
      folder_tag     = f_tag,
      is_transformed = TRUE,
      is_log         = FALSE,
      lambda         = lam_val,
      scale          = 1.0,
      valid_fn       = function(y) !is.na(y) & y > 0,
      fwd            = function(y) {
        if (abs(lam_val) < 1e-4) log(pmax(1e-9, y)) else ((pmax(1e-9, y)^lam_val) - 1.0) / lam_val
      },
      inv            = function(z) {
        if (abs(lam_val) < 1e-4) {
          exp(z)
        } else {
          inner <- pmax(1e-6, lam_val * z + 1.0)
          out <- pmax(0.0, inner^(1.0 / lam_val))
          if (lam_val < 0) pmin(1e6, out) else out
        }
      }
    ))
  }

  if (ttype == "yeojohnson") {
    lam_val <- if (!is.null(lambda) && !is.na(as.numeric(lambda))) {
      as.numeric(lambda)
    } else if (!is.null(y_obs) && length(y_obs) > 0) {
      estimate_yeojohnson_lambda(y_obs)
    } else {
      0.5
    }
    lbl <- sprintf("Yeo-Johnson (lambda = %.2f)", lam_val)
    f_tag <- sprintf("yeojohnson_lam%.2f", lam_val)
    return(list(
      type           = "yeojohnson",
      label          = lbl,
      short_label    = sprintf("yeojohnson[%.2f]", lam_val),
      fn_prefix      = sprintf("yeojohnson[lambda=%.2f]", lam_val),
      folder_tag     = f_tag,
      is_transformed = TRUE,
      is_log         = FALSE,
      lambda         = lam_val,
      scale          = 1.0,
      valid_fn       = function(y) !is.na(y) & is.finite(y),
      fwd            = function(y) {
        y <- as.numeric(y)
        z <- rep(NA_real_, length(y))
        pos <- is.finite(y) & (y >= 0)
        neg <- is.finite(y) & (y < 0)
        if (any(pos)) {
          yp <- y[pos]
          z[pos] <- if (abs(lam_val) < 1e-4) log(yp + 1.0) else (((yp + 1.0)^lam_val) - 1.0) / lam_val
        }
        if (any(neg)) {
          yn <- y[neg]
          z[neg] <- if (abs(2.0 - lam_val) < 1e-4) -log(-yn + 1.0) else -(((-yn + 1.0)^(2.0 - lam_val)) - 1.0) / (2.0 - lam_val)
        }
        z
      },
      inv            = function(z) {
        z <- as.numeric(z)
        y <- rep(NA_real_, length(z))
        pos <- is.finite(z) & (z >= 0)
        neg <- is.finite(z) & (z < 0)
        if (any(pos)) {
          zp <- z[pos]
          yp <- if (abs(lam_val) < 1e-4) exp(zp) - 1.0 else (pmax(1e-6, zp * lam_val + 1.0))^(1.0 / lam_val) - 1.0
          y[pos] <- if (lam_val < 0) pmin(1e6, yp) else yp
        }
        if (any(neg)) {
          zn <- z[neg]
          yn <- if (abs(2.0 - lam_val) < 1e-4) 1.0 - exp(-zn) else 1.0 - (pmax(1e-6, -zn * (2.0 - lam_val) + 1.0))^(1.0 / (2.0 - lam_val))
          y[neg] <- if (lam_val > 2.0) pmax(-1e6, yn) else yn
        }
        y
      }
    ))
  }

  if (ttype == "asin") {
    y_fin <- if (!is.null(y_obs)) y_obs[is.finite(y_obs)] else numeric(0)
    sc_val <- if (!is.null(scale) && is.finite(as.numeric(scale))) as.numeric(scale) else if (length(y_fin) > 0 && max(y_fin) > 1.05) 100.0 else 1.0
    return(list(
      type           = "asin",
      label          = sprintf("Arcsine-Square-Root (0..%g)", sc_val),
      short_label    = "asin_sqrt",
      fn_prefix      = "asin_sqrt",
      folder_tag     = "asin",
      is_transformed = TRUE,
      is_log         = FALSE,
      lambda         = NA_real_,
      scale          = sc_val,
      valid_fn       = function(y) !is.na(y) & y >= 0,
      fwd            = function(y) asin(sqrt(pmin(1.0, pmax(0.0, y / sc_val)))),
      inv            = function(z) sc_val * (sin(pmin(pi / 2.0, pmax(0.0, z))))^2
    ))
  }

  if (ttype == "logit") {
    y_fin <- if (!is.null(y_obs)) y_obs[is.finite(y_obs)] else numeric(0)
    sc_val <- if (!is.null(scale) && is.finite(as.numeric(scale))) as.numeric(scale) else if (length(y_fin) > 0 && max(y_fin) > 1.05) 100.0 else 1.0
    return(list(
      type           = "logit",
      label          = sprintf("Logit (0..%g)", sc_val),
      short_label    = "logit",
      fn_prefix      = "logit",
      folder_tag     = "logit",
      is_transformed = TRUE,
      is_log         = FALSE,
      lambda         = NA_real_,
      scale          = sc_val,
      valid_fn       = function(y) !is.na(y) & y >= 0,
      fwd            = function(y) {
        p <- pmin(1.0 - 1e-4, pmax(1e-4, y / sc_val))
        qlogis(p)
      },
      inv            = function(z) sc_val * plogis(z)
    ))
  }

  list(
    type           = "none",
    label          = "Untransformed (Raw / Identity)",
    short_label    = "raw",
    fn_prefix      = "",
    folder_tag     = "raw",
    is_transformed = FALSE,
    is_log         = FALSE,
    lambda         = NA_real_,
    scale          = 1.0,
    valid_fn       = function(y) !is.na(y),
    fwd            = function(y) y,
    inv            = function(z) z
  )
}

as_transformer <- function(obj, y_obs = NULL, lambda = NULL, scale = NULL) {
  if (is.list(obj) && !is.null(obj$fwd) && !is.null(obj$inv)) {
    return(obj)
  }
  build_transformer(as.character(obj), y_obs = y_obs, lambda = lambda, scale = scale)
}

calc_response_meas_variance <- function(raw_val, raw_unc, trans_obj) {
  n <- length(raw_val)
  v_resp <- numeric(n)
  if (is.null(raw_unc) || all(is.na(raw_unc))) return(v_resp)

  u_clean <- ifelse(is.na(raw_unc) | raw_unc < 0, 0, as.numeric(raw_unc))
  x_raw   <- ifelse(is.na(raw_val), 1.0, as.numeric(raw_val))
  x_pos   <- pmax(1e-6, x_raw)

  tt <- tolower(trans_obj$type)
  if (tt == "log") {
    v_resp <- (u_clean / x_pos)^2
  } else if (tt == "log1p") {
    v_resp <- (u_clean / pmax(1e-6, 1.0 + x_raw))^2
  } else if (tt == "log10") {
    v_resp <- (u_clean / (x_pos * log(10)))^2
  } else if (tt %in% c("raw", "identity", "none")) {
    v_resp <- u_clean^2
  } else if (tt %in% c("sqrt", "square_root")) {
    v_resp <- (u_clean / (2.0 * sqrt(x_pos)))^2
  } else if (tt %in% c("cbrt", "cube_root")) {
    v_resp <- (u_clean / (3.0 * (pmax(1e-6, abs(x_raw))^(2.0 / 3.0))))^2
  } else if (tt %in% c("fourthrt", "qtrt")) {
    v_resp <- (u_clean / (4.0 * (x_pos^0.75)))^2
  } else if (tt %in% c("reciprocal", "inv", "inverse")) {
    v_resp <- (u_clean / (x_pos^2))^2
  } else if (tt %in% c("asin", "arcsin", "arcsine")) {
    sc <- if (!is.null(trans_obj$scale) && is.finite(trans_obj$scale) && trans_obj$scale > 0) trans_obj$scale else 1.0
    p_cl <- pmin(1.0 - 1e-6, pmax(1e-6, x_raw / sc))
    v_resp <- (u_clean / (2.0 * sc * sqrt(p_cl * (1.0 - p_cl))))^2
  } else if (tt == "logit") {
    sc <- if (!is.null(trans_obj$scale) && is.finite(trans_obj$scale) && trans_obj$scale > 0) trans_obj$scale else 1.0
    p_cl <- pmin(1.0 - 1e-6, pmax(1e-6, x_raw / sc))
    v_resp <- (u_clean / (sc * p_cl * (1.0 - p_cl)))^2
  } else if (tt %in% c("boxcox", "bc")) {
    lam <- if (!is.na(trans_obj$lambda)) trans_obj$lambda else 0.0
    v_resp <- (x_pos^(lam - 1.0) * u_clean)^2
  } else if (tt %in% c("yeojohnson", "yj")) {
    lam <- if (!is.na(trans_obj$lambda)) trans_obj$lambda else 0.0
    deriv <- ifelse(x_raw >= 0, (x_raw + 1.0)^(lam - 1.0), (-x_raw + 1.0)^(1.0 - lam))
    v_resp <- (deriv * u_clean)^2
  } else {
    v_resp <- u_clean^2
  }
  v_resp
}

format_transform_folder_tag <- function(trans_obj, default_type = "log") {
  if (is.null(trans_obj)) return("")
  tr <- as_transformer(trans_obj)
  if (tr$type == default_type) return("")
  paste0("trans_", tr$folder_tag)
}

print_available_transforms <- function() {
  cat("Available Data Transformations (--transform=<type> [--lambda=<val>]):\n\n")
  cat(sprintf("%-12s | %-16s | %-18s | %-16s | %s\n",
              "Type", "Short Aliases", "Formula z = f(y)", "Valid Domain", "Best Applied When"))
  cat(strrep("-", 125), "\n")
  cat(sprintf("%-12s | %-16s | %-18s | %-16s | %s\n",
              "log", "ln", "ln(y)", "y > 0", "Default for water chemistry & LCMS. Multiplicative concentrations spanning orders of magnitude."))
  cat(sprintf("%-12s | %-16s | %-18s | %-16s | %s\n",
              "sqrt", "square_root", "y^(1/2)", "y >= 0", "Strongly recommended for moderately skewed concentrations (e.g. UK Nitrate, skewness drop +1.00 -> +0.09) and taxa counts."))
  cat(sprintf("%-12s | %-16s | %-18s | %-16s | %s\n",
              "cbrt", "cube_root", "sign(y)*|y|^(1/3)", "All real y", "Wilson-Hilferty transformation for Gamma-distributed hydrological variables. Compresses high values without crushing near-zero observations."))
  cat(sprintf("%-12s | %-16s | %-18s | %-16s | %s\n",
              "log1p", "log_1p", "ln(1 + y)", "y >= 0", "Accommodates true zeros (y=0) without arbitrary LOD substitution. Gentle log compression for small values."))
  cat(sprintf("%-12s | %-16s | %-18s | %-16s | %s\n",
              "log10", "log_10", "log10(y)", "y > 0", "Order-of-magnitude chemical concentrations (e.g. pathogen counts, electrical conductivity)."))
  cat(sprintf("%-12s | %-16s | %-18s | %-16s | %s\n",
              "fourthrt", "qtrt", "y^(1/4)", "y >= 0", "Intermediate between square root and log; classic ecological standard for heavy-tailed abundance and biomass."))
  cat(sprintf("%-12s | %-16s | %-18s | %-16s | %s\n",
              "boxcox", "bc, box_cox", "(y^lambda-1)/lam", "y > 0", "Analytically optimizes residual normality via Profile Maximum Likelihood lambda in [-2, 2] (or custom --lambda)."))
  cat(sprintf("%-12s | %-16s | %-18s | %-16s | %s\n",
              "yeojohnson", "yj", "Yeo-Johnson(lam)", "All real y", "Extension of Box-Cox supporting zero and negative values. Automatic Profile MLE estimation of lambda."))
  cat(sprintf("%-12s | %-16s | %-18s | %-16s | %s\n",
              "reciprocal", "inv, inverse", "-1 / y", "y > 0", "Monotonic rate-dilution processes. Preserves positive spatial gradient direction (higher y -> higher z)."))
  cat(sprintf("%-12s | %-16s | %-18s | %-16s | %s\n",
              "asin", "arcsin, arcsine", "asin(sqrt(y/sc))", "0 <= y <= scale", "Angular transformation for bounded proportions and percentages (e.g. silt cover, diatom valences)."))
  cat(sprintf("%-12s | %-16s | %-18s | %-16s | %s\n",
              "logit", "logit", "ln(p / (1 - p))", "0 < y < scale", "S-shaped logistic transformation for bounded proportions and sediment percentages."))
  cat(sprintf("%-12s | %-16s | %-18s | %-16s | %s\n",
              "none", "raw, identity", "y", "All real y", "Default for ecological indices (WHPT ASPT, LIFE index, RHS scores) that are already normally distributed."))
}

# ==============================================================================
# NORMALISATION (canonical schema: pred_model_sc, obs_model_sc, pred_tu_sc, pred_td_sc)
# ==============================================================================

compute_norm_ref_stats <- function(ref_df) {
  if (is.null(ref_df) || !("pred_model_sc" %in% names(ref_df))) return(NULL)
  ref_x <- ref_df$pred_model_sc
  ref_x <- ref_x[is.finite(ref_x)]
  if (length(ref_x) == 0L) return(NULL)

  mu_x <- mean(ref_x, na.rm = TRUE)
  sd_x <- sd(ref_x, na.rm = TRUE)
  if (!is.finite(sd_x) || sd_x < 1e-12) sd_x <- 1.0

  med_x <- median(ref_x, na.rm = TRUE)
  mad_x <- mad(ref_x, center = med_x, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(mad_x) || mad_x < 1e-12) {
    iqr_x <- IQR(ref_x, na.rm = TRUE) / 1.349
    mad_x <- if (is.finite(iqr_x) && iqr_x >= 1e-12) iqr_x else sd_x
  }

  min_x <- min(ref_x, na.rm = TRUE)
  rng_x <- max(ref_x, na.rm = TRUE) - min_x
  if (!is.finite(rng_x) || rng_x < 1e-12) rng_x <- 1.0

  list(
    mu_x = mu_x, sd_x = sd_x, med_x = med_x, mad_x = mad_x,
    min_x = min_x, rng_x = rng_x, pct_fn = ecdf(ref_x)
  )
}

add_normalised_columns <- function(df, ref_df = df, ref_stats = NULL) {
  if (!("pred_model_sc" %in% names(df))) return(df)
  if (is.null(ref_stats)) {
    ref_stats <- compute_norm_ref_stats(ref_df)
    if (is.null(ref_stats)) return(df)
  }

  mu_x   <- ref_stats$mu_x
  sd_x   <- ref_stats$sd_x
  med_x  <- ref_stats$med_x
  mad_x  <- ref_stats$mad_x
  min_x  <- ref_stats$min_x
  rng_x  <- ref_stats$rng_x
  pct_fn <- ref_stats$pct_fn

  df$pred_zscore        <- round((df$pred_model_sc - mu_x) / sd_x, 5)
  df$pred_robust_zscore <- round((df$pred_model_sc - med_x) / mad_x, 5)
  df$pred_minmax        <- round((df$pred_model_sc - min_x) / rng_x, 5)
  df$pred_percentile    <- round(100.0 * pct_fn(df$pred_model_sc), 2)

  if ("pred_tu_sc" %in% names(df)) {
    df$pred_tu_zscore        <- round((df$pred_tu_sc - mu_x) / sd_x, 5)
    df$pred_tu_robust_zscore <- round((df$pred_tu_sc - med_x) / mad_x, 5)
    df$pred_tu_minmax        <- round((df$pred_tu_sc - min_x) / rng_x, 5)
    df$pred_tu_percentile    <- round(100.0 * pct_fn(df$pred_tu_sc), 2)
  }
  if ("pred_td_sc" %in% names(df)) {
    df$pred_td_zscore        <- round((df$pred_td_sc - mu_x) / sd_x, 5)
    df$pred_td_robust_zscore <- round((df$pred_td_sc - med_x) / mad_x, 5)
    df$pred_td_minmax        <- round((df$pred_td_sc - min_x) / rng_x, 5)
    df$pred_td_percentile    <- round(100.0 * pct_fn(df$pred_td_sc), 2)
  }
  if ("obs_model_sc" %in% names(df)) {
    has_o <- !is.na(df$obs_model_sc) & is.finite(df$obs_model_sc)
    df$obs_zscore        <- ifelse(has_o, round((df$obs_model_sc - mu_x) / sd_x, 5), NA_real_)
    df$obs_robust_zscore <- ifelse(has_o, round((df$obs_model_sc - med_x) / mad_x, 5), NA_real_)
    df$obs_minmax        <- ifelse(has_o, round((df$obs_model_sc - min_x) / rng_x, 5), NA_real_)
    df$obs_percentile    <- ifelse(has_o, round(100.0 * pct_fn(df$obs_model_sc), 2), NA_real_)
  }
  df
}

# Robust Lognormal Regression on Order Statistics (ROS, Helsel & Cohn 1988) — Vectorized O(n log n)
impute_censored_ros <- function(vals, is_cens) {
  is_cens <- !is.na(is_cens) & as.logical(is_cens) & is.finite(vals)
  if (!any(is_cens)) return(vals)
  n <- length(vals)
  c_limits <- sort(unique(vals[is_cens]))
  J <- length(c_limits)
  min_thr <- min(0.0, min(vals, na.rm = TRUE) - 1.0, na.rm = TRUE)
  thresholds <- c(min_thr, c_limits, Inf)

  cens_mask <- is_cens
  unc_mask  <- !is_cens & is.finite(vals)
  unc_bin   <- integer(n)
  if (any(unc_mask)) {
    unc_bin[unc_mask] <- findInterval(vals[unc_mask], c_limits, left.open = FALSE) + 1L
  }
  cens_bin <- integer(n)
  cens_bin[cens_mask] <- match(vals[cens_mask], c_limits) + 1L

  A_counts <- tabulate(unc_bin[unc_mask], nbins = J + 1L)
  p <- numeric(J + 2L)
  p[J + 2L] <- 0.0
  for (j in seq(J + 1L, 2L, by = -1L)) {
    d_j <- thresholds[j]
    A_j <- A_counts[j]
    B_j <- sum((cens_mask & vals <= d_j) | (unc_mask & vals < d_j))
    p[j] <- if (A_j + B_j > 0L) p[j + 1L] + (1.0 - p[j + 1L]) * (A_j / (A_j + B_j)) else p[j + 1L]
  }
  p[1L] <- if (A_counts[1L] > 0L) 1.0 else p[2L]

  pp <- numeric(n)
  for (j in seq_len(J + 1L)) {
    idx_a <- which(unc_mask & unc_bin == j)
    A_j   <- length(idx_a)
    if (A_j > 0L) {
      ranks_a <- rank(vals[idx_a], ties.method = "average")
      pp[idx_a] <- (1.0 - p[j]) + (p[j] - p[j + 1L]) * (ranks_a / (A_j + 1.0))
    }
  }
  for (j in 2:(J + 1L)) {
    idx_c <- which(cens_mask & cens_bin == j)
    C_j   <- length(idx_c)
    if (C_j > 0L) {
      pp[idx_c] <- (1.0 - p[j]) * (seq_len(C_j) / (C_j + 1.0))
    }
  }
  pp <- pmin(1.0 - 1e-6, pmax(1e-6, pp))

  unc_pos <- which(unc_mask & vals > 0)
  if (length(unc_pos) < 3L) {
    warning("Fewer than 3 positive uncensored samples for ROS; falling back to 0.5 * LOD.")
    return(ifelse(cens_mask, 0.5 * vals, vals))
  }
  z_unc <- stats::qnorm(pp[unc_pos])
  l_unc <- log(vals[unc_pos])
  slope <- max(0.05, stats::cov(z_unc, l_unc) / max(1e-9, stats::var(z_unc)))
  inter <- mean(l_unc) - slope * mean(z_unc)

  out <- vals
  c_idx <- which(cens_mask)
  out[c_idx] <- pmin(vals[c_idx], pmax(1e-6, exp(inter + slope * stats::qnorm(pp[c_idx]))))
  out
}
