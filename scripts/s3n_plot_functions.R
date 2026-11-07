# ==============================================================================
# s3n_plot_functions.R
#
# Modular base-R plotting library for Scalable Spatial Stream Network (S3N)
# models. Sourced by both:
#   - scripts/run_s3n_determinand.R
#   - scripts/03_generate_plots.R
#
# Provides:
#   - Color palettes & scale helpers: val_cols, uncert_cols, map_to_colors(),
#     draw_colorbar(), compute_color_scale(), compute_se_scale()
#   - Stream network renderer: plot_stream_network()
#   - Schema normalizers: normalize_pred_df(), normalize_cv_pts(), normalize_cv_tab()
#   - Publication figure generators:
#       * plot_fig1_england_map()
#       * plot_basin_map()
#       * plot_fig2_humber_map()
#       * plot_fig3_thames_severn_map()
#       * plot_fig4_diagnostics()
# ==============================================================================

# ==============================================================================
# RELATIVE PATH & CUSTOM DIRECTORY RESOLUTION HELPERS
# ==============================================================================

resolve_base_dir <- function() {
  if (file.exists("scripts/s3n_plot_functions.R") || dir.exists("configs/networks")) {
    return(".")
  }
  if (file.exists("s3n_plot_functions.R") && dir.exists("../configs/networks")) {
    return("..")
  }
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg  <- grep("^--file=", full_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_p  <- sub("^--file=", "", file_arg[1])
    cand_root <- dirname(dirname(script_p))
    if (nzchar(cand_root) && dir.exists(cand_root)) {
      return(cand_root)
    }
  }
  "."
}

rel_path <- function(base_dir, ...) {
  parts <- c(list(...))
  parts <- as.character(unlist(parts))
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) return(base_dir)
  p <- if (is.null(base_dir) || !nzchar(base_dir) || base_dir == ".") {
    do.call(file.path, as.list(parts))
  } else {
    do.call(file.path, as.list(c(base_dir, parts)))
  }
  sub("^\\./", "", p)
}

format_rel_display_path <- function(p, base_dir = ".") {
  p_clean <- sub("^\\./", "", p)
  if (!is.null(base_dir) && nzchar(base_dir) && base_dir != ".") {
    prefix <- paste0(gsub("/+$", "", base_dir), "/")
    if (startsWith(p_clean, prefix)) {
      p_clean <- substr(p_clean, nchar(prefix) + 1L, nchar(p_clean))
    }
  }
  p_clean
}

parse_custom_dir_args <- function(args = character(0), base_dir = ".") {
  gis_dir  <- Sys.getenv("S3N_GIS_DIR", unset = "")
  data_dir <- Sys.getenv("S3N_DATA_DIR", unset = "")
  out_dir  <- Sys.getenv("S3N_OUTPUT_DIR", unset = "")
  has_custom_out <- nzchar(out_dir)

  for (a in args) {
    if (grepl("^--(gis-dir|gis)=", a)) {
      gis_dir <- sub("^--(gis-dir|gis)=", "", a)
    } else if (grepl("^--(data-dir|data)=", a)) {
      data_dir <- sub("^--(data-dir|data)=", "", a)
    } else if (grepl("^--(out-dir|output-dir|out|output)=", a)) {
      out_dir <- sub("^--(out-dir|output-dir|out|output)=", "", a)
      has_custom_out <- TRUE
    }
  }

  if (!nzchar(gis_dir))  gis_dir  <- rel_path(base_dir, "gis")
  if (!nzchar(data_dir)) data_dir <- rel_path(base_dir, "data")
  if (!nzchar(out_dir))  out_dir  <- rel_path(base_dir, "output")

  cand_wq <- rel_path(data_dir, "NCEA_River_Surveillance_Network_ARD_Open")
  wq_dir <- if (dir.exists(cand_wq)) {
    cand_wq
  } else if (file.exists(rel_path(data_dir, "RSN Water Quality Samples.csv"))) {
    data_dir
  } else {
    cand_wq
  }

  cand_proc <- rel_path(data_dir, "processed")
  proc_root <- if (dir.exists(cand_proc)) {
    cand_proc
  } else if (file.exists(rel_path(data_dir, "england_rsn_reaches.csv"))) {
    data_dir
  } else {
    cand_proc
  }

  list(
    gis_dir        = sub("^\\./", "", gis_dir),
    data_dir       = sub("^\\./", "", data_dir),
    wq_dir         = sub("^\\./", "", wq_dir),
    proc_root      = sub("^\\./", "", proc_root),
    out_dir        = sub("^\\./", "", out_dir),
    has_custom_out = has_custom_out
  )
}

# ==============================================================================
# DATA SOURCE MANAGEMENT (RSN, SSN, RSS)
# ==============================================================================

resolve_source_tag <- function(source = "rsn") {
  if (is.null(source) || !nzchar(trimws(as.character(source)))) return("rsn")
  s <- tolower(trimws(as.character(source)))
  if (s %in% c("ssn", "small_stream", "small_streams", "small-streams", "small-stream")) {
    "ssn"
  } else if (s %in% c("rss", "combined", "both", "all", "rsn_ssn", "ssn_rsn", "rsn+ssn")) {
    "rss"
  } else {
    "rsn"
  }
}

format_source_abbr <- function(source = "rsn") {
  toupper(resolve_source_tag(source))
}

format_source_label <- function(source = "rsn") {
  s <- resolve_source_tag(source)
  switch(s,
    ssn = "Small Streams Network (SSN)",
    rss = "Combined River Surveillance & Small Streams (RSS)",
    "River Surveillance Network (RSN)"
  )
}

format_source_folder_tag <- function(source = "rsn") {
  s <- resolve_source_tag(source)
  if (s == "rsn") "" else paste0("src_", s)
}

parse_source_arg <- function(args) {
  src_arg <- grep("^--(source|src|data-source)=", args, value = TRUE)
  if (length(src_arg) > 0) {
    return(resolve_source_tag(sub("^--(source|src|data-source)=", "", src_arg[1])))
  }
  src_idx <- which(args %in% c("--source", "--src", "--data-source"))
  if (length(src_idx) > 0 && src_idx[1] < length(args)) {
    return(resolve_source_tag(args[src_idx[1] + 1]))
  }
  "rsn"
}

resolve_data_source_dirs <- function(data_dir, source = "rsn") {
  src_tag  <- resolve_source_tag(source)
  rsn_cand <- rel_path(data_dir, "NCEA_River_Surveillance_Network_ARD_Open")
  rsn_dir  <- if (dir.exists(rsn_cand)) rsn_cand else if (file.exists(rel_path(data_dir, "RSN Water Quality Samples.csv"))) data_dir else rsn_cand
  ssn_cand <- rel_path(data_dir, "NCEA_Small_Streams_Network_ARD_Open")
  ssn_dir  <- if (dir.exists(ssn_cand)) ssn_cand else if (file.exists(rel_path(data_dir, "SSN Water Quality Samples.csv"))) data_dir else ssn_cand

  list(
    source       = src_tag,
    source_abbr  = format_source_abbr(src_tag),
    source_label = format_source_label(src_tag),
    folder_tag   = format_source_folder_tag(src_tag),
    rsn_dir      = rsn_dir,
    ssn_dir      = ssn_dir,
    active_dirs  = switch(src_tag,
      rsn = list(rsn = rsn_dir),
      ssn = list(ssn = ssn_dir),
      rss = list(rsn = rsn_dir, ssn = ssn_dir)
    )
  )
}

print_available_sources <- function() {
  cat("Available Monitoring Data Sources:\n\n")
  cat(sprintf("  %-6s | %-48s | %-12s\n", "Source", "Description", "Default"))
  cat(strrep("-", 73), "\n")
  cat(sprintf("  %-6s | %-48s | %-12s\n", "rsn", "River Surveillance Network (1,382 WQ stations)", "YES (Default)"))
  cat(sprintf("  %-6s | %-48s | %-12s\n", "ssn", "Small Streams Network (520 headwater stations)", "No"))
  cat(sprintf("  %-6s | %-48s | %-12s\n", "rss", "Combined River Surveillance & Small Streams (RSS)", "No"))
  cat("\nUsage:\n")
  cat("  Rscript scripts/run_s3n_determinand.R 0117                  # Default RSN\n")
  cat("  Rscript scripts/run_s3n_determinand.R 0117 --source=ssn      # Small Streams Network only\n")
  cat("  Rscript scripts/run_s3n_determinand.R 0117 --source=rss      # Combined RSN + SSN\n")
  cat("  python3 scripts/01_preprocess_orn.py --source=ssn           # Preprocess SSN reaches\n")
  cat("  python3 scripts/01_preprocess_orn.py --source=rss           # Preprocess combined RSS reaches\n\n")
}

# ==============================================================================
# DATA TRANSFORMATION ENGINE (--transform=<type>)
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
    jac <- sum(sign(y_val) * log(abs(y_val) + 1.0) * ifelse(pos, lam - 1.0, 1.0 - lam))
    -0.5 * n * log(var_z) + jac
  }

  opt <- tryCatch(
    optimize(ll_fn, interval = c(lower, upper), maximum = TRUE),
    error = function(e) list(maximum = 0.5)
  )
  round(opt$maximum, 4)
}

build_transformer <- function(transform_type = "log", y_obs = NULL, lambda = NULL, scale = NULL) {
  t_raw <- tolower(trimws(as.character(transform_type)))
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
  ttype <- if (t_raw %in% names(alias_map)) unname(alias_map[t_raw]) else t_raw

  # 1. Natural Logarithm (Default for wq / lcms)
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

  # 2. Base-10 Logarithm
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

  # 3. Shifted Logarithm (ln(1 + y), handles exact zero)
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

  # 4. Square Root
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

  # 5. Cube Root (Wilson-Hilferty Gamma normalizer)
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

  # 6. Fourth Root (Biomass / Abundance normalizer)
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

  # 7. Monotonic Negative Reciprocal (-1/y)
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

  # 8. Box-Cox Power Transformation
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
          inner <- pmax(1e-9, lam_val * z + 1.0)
          pmax(0.0, inner^(1.0 / lam_val))
        }
      }
    ))
  }

  # 9. Yeo-Johnson Transformation (Handles zeros and negative values)
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
      valid_fn       = function(y) !is.na(y),
      fwd            = function(y) {
        pos <- (y >= 0)
        z <- numeric(length(y))
        if (any(pos)) {
          yp <- y[pos]
          z[pos] <- if (abs(lam_val) < 1e-4) log(yp + 1.0) else (((yp + 1.0)^lam_val) - 1.0) / lam_val
        }
        if (any(!pos)) {
          yn <- y[!pos]
          z[!pos] <- if (abs(2.0 - lam_val) < 1e-4) -log(-yn + 1.0) else -(((-yn + 1.0)^(2.0 - lam_val)) - 1.0) / (2.0 - lam_val)
        }
        z
      },
      inv            = function(z) {
        pos <- (z >= 0)
        y <- numeric(length(z))
        if (any(pos)) {
          zp <- z[pos]
          y[pos] <- if (abs(lam_val) < 1e-4) exp(zp) - 1.0 else (pmax(1e-9, zp * lam_val + 1.0))^(1.0 / lam_val) - 1.0
        }
        if (any(!pos)) {
          zn <- z[!pos]
          y[!pos] <- if (abs(2.0 - lam_val) < 1e-4) 1.0 - exp(-zn) else 1.0 - (pmax(1e-9, -zn * (2.0 - lam_val) + 1.0))^(1.0 / (2.0 - lam_val))
        }
        y
      }
    ))
  }

  # 10. Arcsine-Square-Root (Angular transform for bounded percentages / fractions)
  if (ttype == "asin") {
    sc_val <- if (!is.null(scale)) as.numeric(scale) else if (!is.null(y_obs) && max(y_obs, na.rm = TRUE) > 1.05) 100.0 else 1.0
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

  # 11. Logit (Percentages / Proportions)
  if (ttype == "logit") {
    sc_val <- if (!is.null(scale)) as.numeric(scale) else if (!is.null(y_obs) && max(y_obs, na.rm = TRUE) > 1.05) 100.0 else 1.0
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

  # 12. None / Raw (Identity transform - Default for inverts / diatoms / macrophytes / rhs)
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

as_transformer <- function(obj, y_obs = NULL, lambda = NULL) {
  if (is.list(obj) && !is.null(obj$fwd) && !is.null(obj$inv)) {
    return(obj)
  }
  if (is.logical(obj)) {
    return(if (isTRUE(obj)) build_transformer("log") else build_transformer("none"))
  }
  build_transformer(as.character(obj), y_obs = y_obs, lambda = lambda)
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

val_cols     <- colorRampPalette(c("#1b4f72", "#2874a6", "#1abc9c", "#f4d03f", "#e67e22", "#c0392b", "#641e16"))(256)
val_cols_rev <- rev(val_cols)
uncert_cols  <- colorRampPalette(c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"))(256)

get_val_pal <- function(reverse_color = FALSE) {
  if (isTRUE(reverse_color)) val_cols_rev else val_cols
}

map_to_colors <- function(vals, vmin, vmax, pal) {
  u <- pmin(1, pmax(0, (vals - vmin) / max(1e-9, vmax - vmin)))
  idx <- pmax(1, pmin(length(pal), round(u * (length(pal) - 1)) + 1))
  pal[idx]
}

draw_colorbar <- function(x0, y0, x1, y1, vmin, vmax, pal, ticks, tick_labels, title_txt) {
  n <- length(pal)
  ys <- seq(y0, y1, length.out = n + 1)
  for (i in seq_len(n)) {
    rect(x0, ys[i], x1, ys[i + 1], col = pal[i], border = NA)
  }
  rect(x0, y0, x1, y1, border = "#333333", lwd = 0.8)
  for (j in seq_along(ticks)) {
    ty <- y0 + (ticks[j] - vmin) / max(1e-9, vmax - vmin) * (y1 - y0)
    if (ty >= y0 - 0.02 * (y1 - y0) && ty <= y1 + 0.02 * (y1 - y0)) {
      segments(x1, ty, x1 + (x1 - x0) * 0.35, ty, col = "#333333", lwd = 0.8)
      text(x1 + (x1 - x0) * 0.45, ty, labels = tick_labels[j], adj = c(0, 0.5), cex = 0.78, col = "#222222")
    }
  }
  text(x0, y1 + (y1 - y0) * 0.08, labels = title_txt, adj = c(0, 0), cex = 0.82, font = 2, col = "#111111")
}

compute_color_scale <- function(vals_model_sc, use_log, vmin_raw = NULL, vmax_raw = NULL) {
  tr <- as_transformer(use_log)
  qs <- quantile(vals_model_sc, probs = c(0.01, 0.99), na.rm = TRUE)
  vmin <- as.numeric(qs[1])
  vmax <- as.numeric(qs[2])
  if (!is.null(vmin_raw) && !is.na(vmin_raw)) {
    vmin <- tr$fwd(as.numeric(vmin_raw))
  }
  if (!is.null(vmax_raw) && !is.na(vmax_raw)) {
    vmax <- tr$fwd(as.numeric(vmax_raw))
  }
  if (!is.finite(vmin) || !is.finite(vmax) || vmax - vmin < 1e-4) {
    vmin <- ifelse(is.finite(vmin), vmin, 0.0) - 0.5
    vmax <- ifelse(is.finite(vmax), vmax, 0.0) + 0.5
  }
  if (tr$is_transformed) {
    raw_min <- tr$inv(vmin)
    raw_max <- tr$inv(vmax)
    if (raw_max < raw_min) {
      tmp <- raw_min; raw_min <- raw_max; raw_max <- tmp
    }
    if (tr$type %in% c("log", "log10", "reciprocal") ||
        (tr$type %in% c("boxcox", "yeojohnson") && !is.na(tr$lambda) && tr$lambda <= 0.15)) {
      cands <- c(0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5,
                 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000)
    } else {
      p_ticks <- pretty(c(raw_min, raw_max), n = 6)
      log_cands <- c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 15, 20, 30, 50, 100, 200, 500, 1000)
      cands <- sort(unique(c(p_ticks, log_cands)))
      cands <- cands[tr$valid_fn(cands)]
    }
    lo_b <- if (raw_min > 0) raw_min * 0.95 else raw_min - 0.02 * max(1e-6, raw_max - raw_min)
    hi_b <- if (raw_max > 0) raw_max * 1.05 else raw_max + 0.02 * max(1e-6, raw_max - raw_min)
    raw_in <- cands[cands >= lo_b & cands <= hi_b]
    if (length(raw_in) > 7) {
      raw_in <- raw_in[unique(round(seq(1, length(raw_in), length.out = 6)))]
    }
    if (length(raw_in) >= 3) {
      ticks <- tr$fwd(raw_in)
      labels <- formatC(raw_in, format = "fg", digits = 2)
      if (is.null(vmin_raw) || is.na(vmin_raw)) vmin <- min(vmin, min(ticks))
      if (is.null(vmax_raw) || is.na(vmax_raw)) vmax <- max(vmax, max(ticks))
    } else {
      ticks <- seq(vmin, vmax, length.out = 6)
      labels <- formatC(tr$inv(ticks), format = "fg", digits = 2)
    }
  } else {
    ticks <- pretty(c(vmin, vmax), n = 6)
    ticks <- ticks[ticks >= vmin & ticks <= vmax]
    if (length(ticks) < 3) ticks <- seq(vmin, vmax, length.out = 5)
    labels <- formatC(ticks, format = "fg", digits = 3)
  }
  list(vmin = vmin, vmax = vmax, ticks = ticks, labels = labels)
}

compute_se_scale <- function(se_vals, vmin_se = NULL, vmax_se = NULL) {
  qs <- quantile(se_vals, probs = c(0.02, 0.98), na.rm = TRUE)
  vmin <- max(0.01, floor(qs[1] * 20) / 20)
  vmax <- max(vmin + 0.15, ceiling(qs[2] * 20) / 20)
  if (!is.null(vmin_se) && !is.na(vmin_se)) vmin <- max(0.001, as.numeric(vmin_se))
  if (!is.null(vmax_se) && !is.na(vmax_se)) vmax <- max(vmin + 0.05, as.numeric(vmax_se))
  ticks <- pretty(c(vmin, vmax), n = 5)
  ticks <- ticks[ticks >= vmin & ticks <= vmax]
  if (length(ticks) < 3) ticks <- seq(vmin, vmax, length.out = 5)
  list(vmin = vmin, vmax = vmax, ticks = ticks, labels = sprintf("%.2f", ticks))
}

# Schema normalizers so both 02_fit_s3n_models.R outputs and run_s3n_determinand.R outputs work seamlessly
normalize_pred_df <- function(df) {
  if (!("pred_model_sc" %in% names(df)) && ("pred_log_combined" %in% names(df))) {
    df$pred_model_sc <- df$pred_log_combined
  }
  if (!("se_model_sc" %in% names(df)) && ("se_log_combined" %in% names(df))) {
    df$se_model_sc <- df$se_log_combined
  }
  if (!("obs_model_sc" %in% names(df)) && ("obs_log_nitrate" %in% names(df))) {
    df$obs_model_sc <- df$obs_log_nitrate
  }
  if (!("pred_tu_sc" %in% names(df)) && ("pred_log_tu" %in% names(df))) {
    df$pred_tu_sc <- df$pred_log_tu
  }
  if (!("pred_td_sc" %in% names(df)) && ("pred_log_td" %in% names(df))) {
    df$pred_td_sc <- df$pred_log_td
  }
  if (!("pred_combined" %in% names(df)) && ("pred_mgL_combined" %in% names(df))) {
    df$pred_combined <- df$pred_mgL_combined
  }
  df
}

normalize_cv_pts <- function(cv_df) {
  if (is.null(cv_df)) return(NULL)
  if (!("obs_y" %in% names(cv_df)) && ("obs_log" %in% names(cv_df))) {
    cv_df$obs_y <- cv_df$obs_log
  }
  cv_df
}

normalize_cv_tab <- function(cv_tab) {
  if (is.null(cv_tab)) return(NULL)
  if (!("R2" %in% names(cv_tab)) && ("R2_log" %in% names(cv_tab))) {
    cv_tab$R2 <- cv_tab$R2_log
  }
  if (!("RMSE" %in% names(cv_tab)) && ("RMSE_log" %in% names(cv_tab))) {
    cv_tab$RMSE <- cv_tab$RMSE_log
  }
  if (!("MAE" %in% names(cv_tab)) && ("MAE_log" %in% names(cv_tab))) {
    cv_tab$MAE <- cv_tab$MAE_log
  }
  cv_tab
}

plot_stream_network <- function(df, val_col, vmin, vmax, pal, main_title, sub_title,
                                cbar_ticks, cbar_labels, cbar_title,
                                show_obs = FALSE, obs_col = NULL, base_lwd = 0.35, max_lwd = 2.6) {
  ord <- order(df$Shreve)
  d   <- df[ord, ]

  xr <- range(c(d$start_X, d$end_X), na.rm = TRUE)
  yr <- range(c(d$start_Y, d$end_Y), na.rm = TRUE)
  dx <- diff(xr)
  dy <- diff(yr)

  plot(NA, xlim = c(xr[1] - 0.03 * dx, xr[2] + 0.16 * dx),
       ylim = c(yr[1] - 0.03 * dy, yr[2] + 0.06 * dy),
       asp = 1, axes = FALSE, xlab = "", ylab = "")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#f8f9fa", border = "#d0d4d9", lwd = 1)

  lwds <- base_lwd + (max_lwd - base_lwd) * pmin(1, log10(pmax(1, d$Shreve)) / log10(max(10, max(d$Shreve))))
  cols <- map_to_colors(d[[val_col]], vmin, vmax, pal)
  segments(d$start_X, d$start_Y, d$end_X, d$end_Y, col = cols, lwd = lwds)

  if (show_obs && !is.null(obs_col)) {
    obs_d <- d[d$has_obs == 1, ]
    obs_c <- map_to_colors(obs_d[[obs_col]], vmin, vmax, pal)
    points(obs_d$mid_X, obs_d$mid_Y, pch = 21, bg = obs_c, col = "#111111", cex = 0.85, lwd = 0.7)
  }

  # Scale bar
  sb_len <- if (dx > 250000) 100000 else if (dx > 80000) 25000 else 15000
  sb_lab <- paste0(sb_len / 1000, " km")
  sb_x0  <- if (dx > 250000) xr[2] - sb_len - 0.02 * dx else xr[1] + 0.02 * dx
  sb_y0  <- yr[1] + 0.015 * dy
  rect(sb_x0 - 0.01 * dx, sb_y0 - 0.015 * dy, sb_x0 + sb_len + 0.01 * dx, sb_y0 + 0.038 * dy,
       col = rgb(1, 1, 1, 0.85), border = NA)
  segments(sb_x0, sb_y0, sb_x0 + sb_len, sb_y0, lwd = 2.5, col = "#222222")
  segments(sb_x0, sb_y0 - 0.008 * dy, sb_x0, sb_y0 + 0.008 * dy, lwd = 1.5, col = "#222222")
  segments(sb_x0 + sb_len, sb_y0 - 0.008 * dy, sb_x0 + sb_len, sb_y0 + 0.008 * dy, lwd = 1.5, col = "#222222")
  text(sb_x0 + 0.5 * sb_len, sb_y0 + 0.020 * dy, labels = sb_lab, cex = 0.78, font = 2, col = "#222222")

  # Colorbar
  cb_x0 <- xr[2] + 0.02 * dx
  cb_x1 <- xr[2] + 0.055 * dx
  cb_y0 <- yr[1] + 0.12 * dy
  cb_y1 <- yr[1] + 0.72 * dy
  draw_colorbar(cb_x0, cb_y0, cb_x1, cb_y1, vmin, vmax, pal, cbar_ticks, cbar_labels, cbar_title)

  mtext(main_title, side = 3, line = 1.2, adj = 0, cex = 0.96, font = 2, col = "#111111")
  mtext(sub_title,  side = 3, line = 0.15, adj = 0, cex = 0.80, col = "#444444")
}

shorten_det_name <- function(det_name) {
  tag <- if (grepl("\\[.*\\]$", det_name)) sub("^.*(\\[.*\\])$", " \\1", det_name) else ""
  base <- sub(" \\[.*\\]$", "", det_name)
  if (grepl("\\{DOC\\}", base)) return(paste0("DOC", tag))
  if (grepl("^Oxygen, Dissolved", base)) return(paste0("Dissolved Oxygen", tag))
  if (grepl("^Nitrogen, Total Oxidised", base)) return(paste0("Total Oxidised N", tag))
  if (grepl("^Nitrogen, Total", base)) return(paste0("Total Nitrogen", tag))
  if (grepl("^Phosphorus, Total", base)) return(paste0("Total Phosphorus", tag))
  s <- sub(",.*| :-.*| to pH.*| at 25.*", "", base)
  paste0(substr(s, 1, 22), tag)
}

format_top_banner <- function(network_label = NULL, vars_label = NULL, source_label = NULL) {
  parts <- character(0)
  if (!is.null(network_label) && nzchar(network_label)) {
    parts <- c(parts, sprintf("Water Network: %s", network_label))
  }
  if (!is.null(source_label) && nzchar(source_label)) {
    parts <- c(parts, sprintf("Data Source: %s", source_label))
  }
  if (!is.null(vars_label) && nzchar(vars_label)) {
    parts <- c(parts, sprintf("Fixed-Effect Regression Covariates: %s", vars_label))
  }
  paste(parts, collapse = "   |   ")
}

# Figure 1: Full England Stream Network Map
plot_fig1_england_map <- function(eng_pred, det_code, det_name, det_unit, use_log, out_file,
                                  reverse_color = FALSE, sc_val = NULL, sc_se = NULL,
                                  vars_label = NULL, network_label = NULL, source_label = NULL) {
  tr       <- as_transformer(use_log)
  eng_pred <- normalize_pred_df(eng_pred)
  pal_use  <- get_val_pal(reverse_color)
  has_vars <- !is.null(vars_label) && nzchar(vars_label)
  has_net  <- !is.null(network_label) && nzchar(network_label)
  banner   <- format_top_banner(network_label, vars_label, source_label)
  png(out_file, width = 2600, height = 1450, res = 200, type = "cairo")
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, 2), mar = c(1.2, 1.2, 3.0, 1.2),
      oma = if (nzchar(banner)) c(0, 0, 1.8, 0) else c(0, 0, 0, 0), bg = "white")

  if (is.null(sc_val)) sc_val <- compute_color_scale(eng_pred$pred_model_sc, tr)
  if (is.null(sc_se))  sc_se  <- compute_se_scale(eng_pred$se_model_sc)
  n_obs_eng <- sum(eng_pred$has_obs == 1)
  n_tot_eng <- nrow(eng_pred)
  det_short <- shorten_det_name(det_name)
  det_cbar  <- sub(" \\[.*\\]$", "", det_short)
  cbar_lbl  <- sprintf("%s\n(%s)", det_cbar, det_unit)
  net_tag   <- if (has_net) sprintf(" [%s]", network_label) else ""
  src_tag   <- if (!is.null(source_label) && nzchar(source_label)) format_source_abbr(source_label) else "RSN"
  sub_a <- if (has_vars) {
    sprintf("S3N across %d catchments%s | Covariates: %s",
            length(unique(eng_pred$CatchID)), net_tag, vars_label)
  } else {
    sprintf("Scalable Spatial Stream Network (S3N) across %d catchments%s",
            length(unique(eng_pred$CatchID)), net_tag)
  }

  plot_stream_network(
    eng_pred, "pred_model_sc", vmin = sc_val$vmin, vmax = sc_val$vmax, pal = pal_use,
    main_title = sprintf("A. S3N Predicted %s Across England (%s Reaches)", det_short, format(n_tot_eng, big.mark = ",")),
    sub_title  = sub_a,
    cbar_ticks = sc_val$ticks, cbar_labels = sc_val$labels, cbar_title = cbar_lbl,
    show_obs = FALSE, base_lwd = 0.18, max_lwd = 1.6
  )

  se_scale_txt <- if (tr$type == "none") det_unit else sprintf("%s scale", tr$short_label)
  se_cbar_txt  <- if (tr$type == "none") sprintf("Kriging SE\n(%s)", det_unit) else sprintf("Kriging SE\n(%s scale)", tr$short_label)

  plot_stream_network(
    eng_pred, "se_model_sc", vmin = sc_se$vmin, vmax = sc_se$vmax, pal = uncert_cols,
    main_title = sprintf("B. S3N Prediction Uncertainty & %s Stations (n = %s)", src_tag, format(n_obs_eng, big.mark = ",")),
    sub_title  = sprintf("Kriging standard error (%s) with observed %s stations overlaid%s",
                         se_scale_txt, src_tag, net_tag),
    cbar_ticks = sc_se$ticks, cbar_labels = sc_se$labels,
    cbar_title = se_cbar_txt,
    show_obs = FALSE, base_lwd = 0.18, max_lwd = 1.6
  )
  obs_eng <- eng_pred[eng_pred$has_obs == 1, ]
  points(obs_eng$mid_X, obs_eng$mid_Y, pch = 21,
         bg = map_to_colors(obs_eng$obs_model_sc, sc_val$vmin, sc_val$vmax, pal_use),
         col = "#111111", cex = 0.65, lwd = 0.5)
  legend("topleft", legend = sprintf("Observed %s Station (colored by %s)", src_tag, det_short),
         pch = 21, pt.bg = "#e67e22", col = "#111111", pt.cex = 1.0, bty = "n", cex = 0.80)
  if (nzchar(banner)) {
    mtext(banner, side = 3, line = 0.3, outer = TRUE, cex = 0.88, font = 3, col = "#2c3e50")
  }
}

# Generic 2-Panel River Basin High-Resolution S3N Map (used for Figure 2 Humber & Figure 5+ Custom Basins)
plot_basin_map <- function(basin_pred, basin_title, domain_label, cv_tab, sc_val, det_code, det_name, det_unit, use_log, out_file,
                           reverse_color = FALSE, sc_se = NULL, vars_label = NULL, network_label = NULL, source_label = NULL) {
  tr         <- as_transformer(use_log)
  basin_pred <- normalize_pred_df(basin_pred)
  cv_tab     <- normalize_cv_tab(cv_tab)
  pal_use    <- get_val_pal(reverse_color)
  has_vars   <- !is.null(vars_label) && nzchar(vars_label)
  has_net    <- !is.null(network_label) && nzchar(network_label)
  banner     <- format_top_banner(network_label, vars_label, source_label)
  png(out_file, width = 2600, height = 1350, res = 200, type = "cairo")
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, 2), mar = c(1.2, 1.2, 3.0, 1.2),
      oma = if (nzchar(banner)) c(0, 0, 1.8, 0) else c(0, 0, 0, 0), bg = "white")

  if (is.null(sc_se)) sc_se <- compute_se_scale(basin_pred$se_model_sc)
  n_obs_basin <- sum(basin_pred$has_obs == 1)
  det_short   <- shorten_det_name(det_name)
  net_tag     <- if (has_net) sprintf(" [%s]", network_label) else ""
  src_tag     <- if (!is.null(source_label) && nzchar(source_label)) format_source_abbr(source_label) else "RSN"
  cv_row      <- subset(cv_tab, (Domain == domain_label | grepl(paste0("^", basin_title), Domain)) &
                                Model == "S3N Combined (Tail-Up + Tail-Down)")
  sub_a <- if (nrow(cv_row) > 0) {
    if (has_vars) {
      sprintf("5-Fold CV R2 = %.3f (RMSE = %.3f) | n = %d sites%s | Covariates: %s",
              cv_row$R2[1], cv_row$RMSE[1], n_obs_basin, net_tag, vars_label)
    } else {
      sprintf("5-Fold CV R2 = %.3f (RMSE = %.3f) | %d observed %s stations%s",
              cv_row$R2[1], cv_row$RMSE[1], n_obs_basin, src_tag, net_tag)
    }
  } else {
    if (has_vars) {
      sprintf("England S3N subset (n = %d sites)%s | Covariates: %s", n_obs_basin, net_tag, vars_label)
    } else {
      sprintf("England S3N model subset | %d observed %s stations in %s%s", n_obs_basin, src_tag, basin_title, net_tag)
    }
  }

  plot_stream_network(
    basin_pred, "pred_model_sc", vmin = sc_val$vmin, vmax = sc_val$vmax, pal = pal_use,
    main_title = sprintf("A. %s: S3N Predicted %s (%s Reaches)",
                         basin_title, det_short, format(nrow(basin_pred), big.mark = ",")),
    sub_title  = sub_a,
    cbar_ticks = sc_val$ticks, cbar_labels = sc_val$labels,
    cbar_title = sprintf("%s\n(%s)", sub(" \\[.*\\]$", "", det_short), det_unit),
    show_obs = TRUE, obs_col = "obs_model_sc", base_lwd = 0.35, max_lwd = 2.8
  )
  legend("topleft", legend = sprintf("Observed %s Station (same color scale)", src_tag),
         pch = 21, pt.bg = "#f4d03f", col = "#111111", pt.cex = 1.1, bty = "n", cex = 0.82)

  se_cbar_txt <- if (tr$type == "none") sprintf("Kriging SE\n(%s)", det_unit) else sprintf("Kriging SE\n(%s scale)", tr$short_label)

  plot_stream_network(
    basin_pred, "se_model_sc", vmin = sc_se$vmin, vmax = sc_se$vmax, pal = uncert_cols,
    main_title = sprintf("B. %s: S3N Kriging Standard Error (%s)", basin_title, det_short),
    sub_title  = sprintf("Lower along monitored flow-connected paths, higher in headwaters%s", net_tag),
    cbar_ticks = sc_se$ticks, cbar_labels = sc_se$labels,
    cbar_title = se_cbar_txt,
    show_obs = FALSE, base_lwd = 0.35, max_lwd = 2.8
  )
  obs_b <- basin_pred[basin_pred$has_obs == 1, ]
  if (nrow(obs_b) > 0) {
    points(obs_b$mid_X, obs_b$mid_Y, pch = 3, col = "#c0392b", cex = 0.55, lwd = 0.8)
  }
  legend("topleft", legend = sprintf("%s Monitoring Site Location (n = %d)", src_tag, n_obs_basin),
         pch = 3, col = "#c0392b", pt.cex = 0.9, bty = "n", cex = 0.82)
  if (nzchar(banner)) {
    mtext(banner, side = 3, line = 0.3, outer = TRUE, cex = 0.88, font = 3, col = "#2c3e50")
  }
}

# Figure 2: Humber River Basin High-Resolution S3N Map
plot_fig2_humber_map <- function(hum_pred, cv_tab, sc_val, det_code, det_name, det_unit, use_log, out_file,
                                 reverse_color = FALSE, sc_se = NULL, vars_label = NULL, network_label = NULL, source_label = NULL) {
  plot_basin_map(hum_pred, "Humber Basin", "Humber Basin (CatchID=80)",
                 cv_tab, sc_val, det_code, det_name, det_unit, use_log, out_file,
                 reverse_color = reverse_color, sc_se = sc_se,
                 vars_label = vars_label, network_label = network_label, source_label = source_label)
}

# Figure 3: Thames & Severn River Basins S3N Maps
plot_fig3_thames_severn_map <- function(tha_pred, sev_pred, sc_val, det_code, det_name, det_unit, use_log, out_file,
                                        reverse_color = FALSE, vars_label = NULL, network_label = NULL, source_label = NULL) {
  tha_pred <- normalize_pred_df(tha_pred)
  sev_pred <- normalize_pred_df(sev_pred)
  pal_use  <- get_val_pal(reverse_color)
  has_vars <- !is.null(vars_label) && nzchar(vars_label)
  has_net  <- !is.null(network_label) && nzchar(network_label)
  banner   <- format_top_banner(network_label, vars_label, source_label)
  png(out_file, width = 2600, height = 1250, res = 200, type = "cairo")
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, 2), mar = c(1.2, 1.2, 3.0, 1.2),
      oma = if (nzchar(banner)) c(0, 0, 1.8, 0) else c(0, 0, 0, 0), bg = "white")

  det_short <- shorten_det_name(det_name)
  det_cbar  <- sub(" \\[.*\\]$", "", det_short)
  cbar_lbl  <- sprintf("%s\n(%s)", det_cbar, det_unit)
  net_tag   <- if (has_net) sprintf(" [%s]", network_label) else ""
  src_tag   <- if (!is.null(source_label) && nzchar(source_label)) format_source_abbr(source_label) else "RSN"
  sub_tha <- if (has_vars) {
    sprintf("n = %d %s sites%s | S3N Tail-Up | Covariates: %s", sum(tha_pred$has_obs == 1), src_tag, net_tag, vars_label)
  } else {
    sprintf("n = %d observed %s stations%s | Flow-connected Tail-Up predictions",
            sum(tha_pred$has_obs == 1), src_tag, net_tag)
  }
  sub_sev <- if (has_vars) {
    sprintf("n = %d %s sites%s | S3N Combined | Covariates: %s", sum(sev_pred$has_obs == 1), src_tag, net_tag, vars_label)
  } else {
    sprintf("n = %d observed %s stations%s | Combined Tail-Up + Tail-Down predictions",
            sum(sev_pred$has_obs == 1), src_tag, net_tag)
  }

  plot_stream_network(
    tha_pred, "pred_tu_sc", vmin = sc_val$vmin, vmax = sc_val$vmax, pal = pal_use,
    main_title = sprintf("A. Thames Basin: S3N Tail-Up Predicted %s (%s Reaches)",
                         det_short, format(nrow(tha_pred), big.mark = ",")),
    sub_title  = sub_tha,
    cbar_ticks = sc_val$ticks, cbar_labels = sc_val$labels, cbar_title = cbar_lbl,
    show_obs = TRUE, obs_col = "obs_model_sc", base_lwd = 0.40, max_lwd = 2.8
  )

  plot_stream_network(
    sev_pred, "pred_model_sc", vmin = sc_val$vmin, vmax = sc_val$vmax, pal = pal_use,
    main_title = sprintf("B. Severn Basin: S3N Predicted %s (%s Reaches)",
                         det_short, format(nrow(sev_pred), big.mark = ",")),
    sub_title  = sub_sev,
    cbar_ticks = sc_val$ticks, cbar_labels = sc_val$labels, cbar_title = cbar_lbl,
    show_obs = TRUE, obs_col = "obs_model_sc", base_lwd = 0.40, max_lwd = 2.8
  )
  if (nzchar(banner)) {
    mtext(banner, side = 3, line = 0.3, outer = TRUE, cex = 0.88, font = 3, col = "#2c3e50")
  }
}

# Figure 4: 4-Panel Model Diagnostics, Cross-Validation & Parameter CIs
plot_fig4_diagnostics <- function(eng_pred, eng_cv, hum_cv, cv_tab, param_tab, det_code, det_name, det_unit, use_log, out_file,
                                  extra_cv_list = list(), lims_a = NULL, ylim_d = NULL,
                                  vars_label = NULL, network_label = NULL, source_label = NULL) {
  tr       <- as_transformer(use_log)
  eng_pred <- normalize_pred_df(eng_pred)
  eng_cv   <- normalize_cv_pts(eng_cv)
  hum_cv   <- normalize_cv_pts(hum_cv)
  cv_tab   <- normalize_cv_tab(cv_tab)
  banner   <- format_top_banner(network_label, vars_label, source_label)

  png(out_file, width = 2400, height = 1800, res = 200, type = "cairo")
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(2, 2), mar = c(4.5, 5.2, 3.2, 1.8),
      oma = if (nzchar(banner)) c(0, 0, 1.8, 0) else c(0, 0, 0, 0), bg = "white")

  det_short <- shorten_det_name(det_name)
  sc_lbl <- if (tr$type == "none") {
    sprintf("%s (%s)", det_short, det_unit)
  } else {
    sprintf("%s(%s, %s)", tr$fn_prefix, det_short, det_unit)
  }

  # Panel A: Observed vs Predicted (5-Fold CV) on Humber, Custom Basin(s) & All England
  if (is.null(lims_a)) {
    lims <- range(c(eng_cv$obs_y, eng_cv$cv_comb), na.rm = TRUE)
    pad  <- 0.08 * diff(lims)
    lims <- c(lims[1] - pad, lims[2] + pad)
  } else {
    lims <- lims_a
  }

  plot(eng_cv$obs_y, eng_cv$cv_comb,
       pch = 16, col = rgb(0.16, 0.45, 0.65, 0.38), cex = 0.85,
       xlim = lims, ylim = lims,
       xlab = paste("Observed", sc_lbl),
       ylab = paste("5-Fold CV Predicted", sc_lbl),
       main = "A. Out-of-Sample 5-Fold Cross-Validation Accuracy")
  grid(col = "#e0e0e0", lty = "dotted")
  if (!is.null(hum_cv) && nrow(hum_cv) > 0) {
    points(hum_cv$obs_y, hum_cv$cv_comb, pch = 21, bg = rgb(0.90, 0.49, 0.13, 0.75), col = "#222222", cex = 0.95, lwd = 0.6)
  }

  eng_dom <- if (any(grepl("All England", cv_tab$Domain))) {
    unique(cv_tab$Domain[grep("All England", cv_tab$Domain)])[1]
  } else {
    "All England RSN Network"
  }
  eng_row <- subset(cv_tab, Domain == eng_dom & Model == "S3N Combined (Tail-Up + Tail-Down)")
  if (nrow(eng_row) == 0) {
    eng_row <- subset(cv_tab, Model == "S3N Combined (Tail-Up + Tail-Down)")
  }
  hum_row <- subset(cv_tab, grepl("^Humber Basin", Domain) & Model == "S3N Combined (Tail-Up + Tail-Down)")

  leg_txt  <- sprintf("All England (n = %d): R2 = %.3f, RMSE = %.3f", nrow(eng_cv), eng_row$R2[1], eng_row$RMSE[1])
  leg_pch  <- 16
  leg_lty  <- NA
  leg_lwd  <- NA
  leg_ptbg <- NA
  leg_col  <- rgb(0.16, 0.45, 0.65, 0.8)

  if (nrow(hum_row) > 0) {
    leg_txt  <- c(leg_txt, sprintf("Humber Basin (n = %d): R2 = %.3f, RMSE = %.3f", hum_row$n_obs[1], hum_row$R2[1], hum_row$RMSE[1]))
    leg_pch  <- c(leg_pch, 21)
    leg_lty  <- c(leg_lty, NA)
    leg_lwd  <- c(leg_lwd, NA)
    leg_ptbg <- c(leg_ptbg, rgb(0.90, 0.49, 0.13, 0.85))
    leg_col  <- c(leg_col, "#222222")
  }

  extra_pchs <- c(24, 22, 23, 25)
  extra_bgs  <- c(rgb(0.15, 0.68, 0.38, 0.85), rgb(0.56, 0.27, 0.68, 0.85),
                  rgb(0.95, 0.61, 0.07, 0.85), rgb(0.09, 0.63, 0.52, 0.85))
  if (length(extra_cv_list) > 0) {
    for (k in seq_along(extra_cv_list)) {
      ecv <- extra_cv_list[[k]]
      pts_k <- normalize_cv_pts(ecv$cv_pts)
      if (!is.null(pts_k) && nrow(pts_k) > 0) {
        pk <- extra_pchs[(k - 1) %% length(extra_pchs) + 1]
        bk <- extra_bgs[(k - 1) %% length(extra_bgs) + 1]
        points(pts_k$obs_y, pts_k$cv_comb, pch = pk, bg = bk, col = "#111111", cex = 1.0, lwd = 0.7)
        erow <- subset(cv_tab, Domain == ecv$domain_label & Model == "S3N Combined (Tail-Up + Tail-Down)")
        if (nrow(erow) > 0) {
          leg_txt  <- c(leg_txt, sprintf("%s (n = %d): R2 = %.3f, RMSE = %.3f", ecv$short_title, erow$n_obs[1], erow$R2[1], erow$RMSE[1]))
          leg_pch  <- c(leg_pch, pk)
          leg_lty  <- c(leg_lty, NA)
          leg_lwd  <- c(leg_lwd, NA)
          leg_ptbg <- c(leg_ptbg, bk)
          leg_col  <- c(leg_col, "#111111")
        }
      }
    }
  }

  abline(0, 1, col = "#c0392b", lwd = 2, lty = 2)
  leg_txt  <- c(leg_txt, "1:1 Perfect Prediction Line")
  leg_pch  <- c(leg_pch, NA)
  leg_lty  <- c(leg_lty, 2)
  leg_lwd  <- c(leg_lwd, 2)
  leg_ptbg <- c(leg_ptbg, NA)
  leg_col  <- c(leg_col, "#c0392b")

  legend("topleft", legend = leg_txt, pch = leg_pch, lty = leg_lty, lwd = leg_lwd,
         pt.bg = leg_ptbg, col = leg_col, bty = "n", cex = 0.80)

  # Panel B: 5-Fold CV R2 Comparison Across Domains
  avail_doms <- unique(cv_tab$Domain)
  r2_mat <- sapply(avail_doms, function(d) {
    sub <- cv_tab[cv_tab$Domain == d, ]
    pmax(0, c(
      sub$R2[sub$Model == "OLS (Non-Spatial)"],
      sub$R2[sub$Model == "S3N Tail-Down (Stream Distance)"],
      sub$R2[sub$Model == "S3N Tail-Up (Flow-Connected + AFV)"],
      sub$R2[sub$Model == "S3N Combined (Tail-Up + Tail-Down)"]
    ))
  })
  short_map <- c(
    "Humber Basin (CatchID=80)"   = "Humber Basin",
    "Severn Basin (CatchID=2786)" = "Severn Basin",
    "Thames Basin (CatchID=3055)" = "Thames Basin"
  )
  short_map[eng_dom] <- sprintf("All England\n(n=%d)", eng_row$n_obs[1])
  colnames(r2_mat) <- ifelse(avail_doms %in% names(short_map), short_map[avail_doms],
                             sub(" \\((CatchID|RBD)=.*\\)$", "", avail_doms))
  ymax_b <- max(0.88, max(r2_mat, na.rm = TRUE) * 1.18)
  cex_names_b <- if (ncol(r2_mat) > 4) 0.78 else 0.90
  cex_txt_b   <- if (ncol(r2_mat) > 4) 0.64 else 0.72
  bp <- barplot(r2_mat, beside = TRUE, ylim = c(0, ymax_b),
                col = c("#95a5a6", "#3498db", "#2ecc71", "#e67e22"),
                border = "#222222", cex.names = cex_names_b, ylab = "5-Fold Cross-Validation R2",
                main = "B. Out-of-Sample Variance Explained (R2) by Model")
  grid(nx = NA, ny = NULL, col = "#e0e0e0", lty = "dotted")
  text(bp, r2_mat + 0.03 * ymax_b, labels = sprintf("%.2f", r2_mat), cex = cex_txt_b, font = 2)
  legend("topright",
         legend = c("OLS (Non-Spatial)", "S3N Tail-Down (Stream)", "S3N Tail-Up (FC + AFV)", "S3N Combined (TU + TD)"),
         fill = c("#95a5a6", "#3498db", "#2ecc71", "#e67e22"), border = "#222222", bty = "n", cex = 0.80)

  # Panel C: Forest Plot of Fixed-Effect Estimates & 95% Bootstrap CIs (All England S3N)
  eng_p_dom <- if (!is.null(param_tab) && any(grepl("All England", param_tab$Domain))) {
    unique(param_tab$Domain[grep("All England", param_tab$Domain)])[1]
  } else {
    "All England RSN Network"
  }
  eng_td_p <- subset(param_tab, Domain == eng_p_dom & Model == "S3N Tail-Down (Stream)" &
                     !(Parameter %in% c("Intercept", "sigma.sq", "tau.sq", "phi")))
  eng_tu_p <- subset(param_tab, Domain == eng_p_dom & Model == "S3N Tail-Up (FC + AFV)" &
                     !(Parameter %in% c("Intercept", "sigma.sq", "tau.sq", "phi")))
  if (nrow(eng_td_p) == 0) {
    eng_td_p <- subset(param_tab, Domain == eng_p_dom & Model == "S3N Tail-Down (Stream)" &
                       !(Parameter %in% c("sigma.sq", "tau.sq", "phi")))
    eng_tu_p <- subset(param_tab, Domain == eng_p_dom & Model == "S3N Tail-Up (FC + AFV)" &
                       !(Parameter %in% c("sigma.sq", "tau.sq", "phi")))
  }

  nice_labels <- c(
    Intercept          = "Intercept (Spatial Mean)",
    Elevation_z        = "Elevation (std)",
    Slope_z            = "Channel Slope (std)",
    log10_US_Accum_z   = "Upstream Length (std)",
    Sinuosity_z        = "Channel Sinuosity (std)",
    SrcElevation_z     = "Source Elevation (std)",
    Dist2Mth_z         = "Distance to Mouth (std)",
    Length_z           = "Reach Length (std)",
    log10_Shreve_z     = "Shreve Order (log10 std)",
    Strahler_z         = "Strahler Order (std)",
    Is_Chalk_Limestone = "Chalk / Limestone Aquifer",
    Is_Sandstone       = "Sandstone Aquifer",
    Is_Aquifer         = "Aquifer Bedrock",
    Is_Peat            = "Peat Superficial Geology",
    Easting_z          = "Easting (std)",
    Northing_z         = "Northing (std)"
  )
  np <- nrow(eng_td_p)
  y_pos <- rev(seq_len(np))
  par_MAR_old <- par("mar")
  par(mar = c(4.5, 10.2, 3.2, 1.5))
  xr_ci <- range(c(eng_tu_p$CI_Lower, eng_tu_p$CI_Upper, eng_td_p$CI_Lower, eng_td_p$CI_Upper, -0.2, 0.2), na.rm = TRUE)
  xr_pad <- 0.12 * diff(xr_ci)
  plot(NA, xlim = c(xr_ci[1] - xr_pad, xr_ci[2] + xr_pad), ylim = c(0.4, np + 0.6),
       yaxt = "n", xlab = sprintf("Effect Size on %s [95%% Bootstrap CI]", sc_lbl), ylab = "",
       main = sprintf("C. Environmental Drivers of %s (England S3N)", det_short))
  grid(nx = NULL, ny = NA, col = "#e0e0e0", lty = "dotted")
  abline(v = 0, col = "#c0392b", lty = 2, lwd = 1.5)
  lbls_c <- ifelse(eng_td_p$Parameter %in% names(nice_labels), nice_labels[eng_td_p$Parameter], eng_td_p$Parameter)
  axis(2, at = y_pos, labels = lbls_c, las = 1, cex.axis = 0.85)

  segments(eng_tu_p$CI_Lower, y_pos + 0.14, eng_tu_p$CI_Upper, y_pos + 0.14, col = "#27ae60", lwd = 2.2)
  points(eng_tu_p$Estimate, y_pos + 0.14, pch = 16, col = "#27ae60", cex = 1.1)

  segments(eng_td_p$CI_Lower, y_pos - 0.14, eng_td_p$CI_Upper, y_pos - 0.14, col = "#2980b9", lwd = 2.2)
  points(eng_td_p$Estimate, y_pos - 0.14, pch = 17, col = "#2980b9", cex = 1.1)

  legend("bottomright", legend = c("S3N Tail-Up (95% Bootstrap CI)", "S3N Tail-Down (95% Bootstrap CI)"),
         col = c("#27ae60", "#2980b9"), pch = c(16, 17), lwd = 2.2, bty = "n", cex = 0.80)
  par(mar = par_MAR_old)

  # Panel D: Distribution across All Reaches (Overall) and by Hydrogeology & Elevation Regime
  reg_str <- ifelse(eng_pred$Elevation > 150, "Upland\n(>150m)",
             ifelse(eng_pred$Is_Aquifer == 1, "Lowland Aquifer\n(Chalk/Limestone/Sand)", "Lowland\nNon-Aquifer"))
  box_list <- list(
    "All Reaches\n(Overall)"                  = eng_pred$pred_combined,
    "Upland\n(>150m)"                         = eng_pred$pred_combined[reg_str == "Upland\n(>150m)"],
    "Lowland\nNon-Aquifer"                    = eng_pred$pred_combined[reg_str == "Lowland\nNon-Aquifer"],
    "Lowland Aquifer\n(Chalk/Limestone/Sand)" = eng_pred$pred_combined[reg_str == "Lowland Aquifer\n(Chalk/Limestone/Sand)"]
  )

  bp_d <- boxplot(box_list, outline = FALSE, ylim = ylim_d,
                  col = c("#95a5a6", "#2874a6", "#1abc9c", "#f4d03f"), border = "#222222",
                  cex.axis = 0.84,
                  ylab = sprintf("S3N Predicted %s (%s)", det_short, det_unit),
                  xlab = "Overall Network vs. Hydrogeological & Topographic Regime",
                  main = sprintf("D. Predicted %s Across %s Reaches by Regime",
                                 det_short, format(nrow(eng_pred), big.mark = ",")))
  grid(nx = NA, ny = NULL, col = "#e0e0e0", lty = "dotted")
  abline(v = 1.5, col = "#7f8c8d", lty = 2, lwd = 1.3)
  medians <- sapply(box_list, median, na.rm = TRUE)
  yr_d <- par("usr")[4] - par("usr")[3]
  y_txt <- pmin(par("usr")[4] - 0.07 * yr_d, bp_d$stats[4, ] + 0.07 * yr_d)
  text(1:4, y_txt, labels = sprintf("Median:\n%s %s", formatC(medians, format = "fg", digits = 3), det_unit),
       font = 2, cex = 0.78, col = "#111111")
  if (nzchar(banner)) {
    mtext(banner, side = 3, line = 0.3, outer = TRUE, cex = 0.88, font = 3, col = "#2c3e50")
  }
}

# Figure 5 (Temporal Multi-Panel Grid): All England S3N Comparison Across Seasons, Years, or Months
plot_fig5_temporal_england_grid <- function(period_items, by_mode, sc_val, det_code, det_name, det_unit,
                                            use_log, out_file, reverse_color = FALSE,
                                            vars_label = NULL, network_label = NULL, source_label = NULL) {
  period_items <- Filter(Negate(is.null), period_items)
  n_pan <- length(period_items)
  if (n_pan == 0) return(invisible(NULL))

  # Order seasonal panels in canonical Spring -> Summer -> Autumn -> Winter sequence
  if (tolower(by_mode) == "season") {
    s_order <- c("spring", "summer", "autumn", "winter")
    keys <- sapply(period_items, function(x) {
      lbl <- tolower(x$label)
      hit <- which(sapply(s_order, function(s) grepl(s, lbl)))
      if (length(hit) > 0) hit[1] else 99L
    })
    period_items <- period_items[order(keys)]
  }

  if (n_pan <= 2) {
    nrow_g <- 1L; ncol_g <- 2L; w_px <- 2600L; h_px <- 1450L
  } else if (n_pan <= 4) {
    nrow_g <- 2L; ncol_g <- 2L; w_px <- 2600L; h_px <- 2750L
  } else if (n_pan <= 6) {
    nrow_g <- 2L; ncol_g <- 3L; w_px <- 3700L; h_px <- 2700L
  } else if (n_pan <= 9) {
    nrow_g <- 3L; ncol_g <- 3L; w_px <- 3700L; h_px <- 3800L
  } else {
    nrow_g <- 3L; ncol_g <- 4L; w_px <- 4600L; h_px <- 3700L
  }

  pal_use   <- get_val_pal(reverse_color)
  det_short <- sub(" \\[.*\\]$", "", shorten_det_name(det_name))
  cbar_lbl  <- sprintf("%s\n(%s)", det_short, det_unit)
  mode_cap  <- paste0(toupper(substr(by_mode, 1, 1)), tolower(substr(by_mode, 2, nchar(by_mode))))

  png(out_file, width = w_px, height = h_px, res = 200, type = "cairo")
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(nrow_g, ncol_g), mar = c(1.0, 1.0, 2.8, 1.0), oma = c(0.5, 0.5, 3.4, 0.5), bg = "white")

  month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

  for (i in seq_along(period_items)) {
    item   <- period_items[[i]]
    eng_df <- normalize_pred_df(item$eng_pred)
    n_obs_i <- sum(eng_df$has_obs == 1, na.rm = TRUE)
    med_raw <- median(eng_df$pred_combined, na.rm = TRUE)

    # Format clean panel title
    raw_lbl <- item$label
    if (grepl("^Season:?\\s*", raw_lbl, ignore.case = TRUE)) {
      s_nm <- tolower(trimws(sub("^Season:?\\s*", "", raw_lbl, ignore.case = TRUE)))
      s_pretty <- switch(s_nm,
        spring = "Spring (Mar-May)",
        summer = "Summer (Jun-Aug)",
        autumn = "Autumn (Sep-Nov)",
        winter = "Winter (Dec-Feb)",
        paste0(toupper(substr(s_nm, 1, 1)), substr(s_nm, 2, nchar(s_nm)))
      )
      pan_lbl <- s_pretty
    } else if (grepl("^Month\\s+([0-9]+)$", raw_lbl, ignore.case = TRUE)) {
      m_idx <- as.integer(sub("^Month\\s+([0-9]+)$", "\\1", raw_lbl, ignore.case = TRUE))
      pan_lbl <- if (!is.na(m_idx) && m_idx >= 1 && m_idx <= 12) {
        sprintf("%s (Month %02d)", month_names[m_idx], m_idx)
      } else {
        raw_lbl
      }
    } else {
      pan_lbl <- raw_lbl
    }

    r2_str <- ""
    if (!is.null(item$cv_tab)) {
      cvt <- normalize_cv_tab(item$cv_tab)
      erow <- subset(cvt, grepl("All England", Domain) & Model == "S3N Combined (Tail-Up + Tail-Down)")
      if (nrow(erow) > 0 && !is.na(erow$R2[1])) {
        r2_str <- sprintf(" | CV R2 = %.3f", erow$R2[1])
      }
    }

    main_t <- sprintf("%s. %s — Predicted %s", LETTERS[i], pan_lbl, det_short)
    sub_t  <- sprintf("n = %s obs reaches%s | Network Median = %s %s",
                      format(n_obs_i, big.mark = ","), r2_str,
                      formatC(med_raw, format = "fg", digits = 3), det_unit)

    plot_stream_network(
      eng_df, "pred_model_sc", vmin = sc_val$vmin, vmax = sc_val$vmax, pal = pal_use,
      main_title = main_t,
      sub_title  = sub_t,
      cbar_ticks = sc_val$ticks, cbar_labels = sc_val$labels, cbar_title = cbar_lbl,
      show_obs = FALSE, base_lwd = 0.18, max_lwd = 1.6
    )
  }

  n_reaches_tot <- nrow(period_items[[1]]$eng_pred)
  n_catch_tot   <- length(unique(period_items[[1]]$eng_pred$CatchID))
  has_vars      <- !is.null(vars_label) && nzchar(vars_label)
  has_net       <- !is.null(network_label) && nzchar(network_label)
  has_src       <- !is.null(source_label) && nzchar(source_label)
  src_tag       <- if (has_src) format_source_abbr(source_label) else ""
  sub_outer <- sprintf("Unified color scale across all %d %s panels (%s reaches across %d catchments)%s%s%s",
                       n_pan, tolower(by_mode), format(n_reaches_tot, big.mark = ","), n_catch_tot,
                       if (has_src) paste0(" | Source: ", src_tag) else "",
                       if (has_net) paste0(" | Network: ", network_label) else "",
                       if (has_vars) paste0(" | Covariates: ", vars_label) else "")
  mtext(sprintf("All England S3N Stream Network Comparison by %s: %s (%s)",
                mode_cap, det_short, det_unit),
        side = 3, line = 1.5, outer = TRUE, cex = 1.25, font = 2, col = "#111111")
  mtext(sub_outer, side = 3, line = 0.2, outer = TRUE, cex = 0.92, col = "#444444")
}

# ==============================================================================
# SHARED FIXED-EFFECT DESIGN MATRIX BUILDER & COVARIATE SELECTION (--vars=...)
# ==============================================================================

COVARIATE_CATALOG <- list(
  Elevation_z = list(
    aliases = c("elevation", "elev", "elevation_z", "os_t50_ele"),
    col     = "Elevation",
    desc    = "Mean channel elevation (m AOD, z-score standardized)",
    default = "All domains",
    eval_fn = function(df) as.numeric(scale(df$Elevation))
  ),
  Slope_z = list(
    aliases = c("slope", "slope_z"),
    col     = "Slope",
    desc    = "Mean channel gradient (%, z-score standardized)",
    default = "All domains",
    eval_fn = function(df) as.numeric(scale(df$Slope))
  ),
  log10_US_Accum_z = list(
    aliases = c("us_accum", "accum", "us_accum_km", "log10_us_accum", "log10_us_accum_z"),
    col     = "US_Accum_km",
    desc    = "log10 upstream accumulated channel length (km, z-score standardized)",
    default = "All domains",
    eval_fn = function(df) as.numeric(scale(log10(pmax(1e-4, df$US_Accum_km))))
  ),
  Is_Chalk_Limestone = list(
    aliases = c("chalk_limestone", "chalk", "limestone", "is_chalk_limestone"),
    col     = "Is_Chalk_Limestone",
    desc    = "Binary indicator (0/1) for BGS Chalk / Limestone / Dolostone bedrock",
    default = "All England (default)",
    eval_fn = function(df) as.numeric(df$Is_Chalk_Limestone)
  ),
  Is_Sandstone = list(
    aliases = c("sandstone", "sand", "is_sandstone"),
    col     = "Is_Sandstone",
    desc    = "Binary indicator (0/1) for BGS Sandstone aquifer bedrock",
    default = "All England (default)",
    eval_fn = function(df) as.numeric(df$Is_Sandstone)
  ),
  Is_Peat = list(
    aliases = c("peat", "is_peat"),
    col     = "Is_Peat",
    desc    = "Binary indicator (0/1) for BGS Peat superficial deposits",
    default = "All England (default)",
    eval_fn = function(df) as.numeric(df$Is_Peat)
  ),
  Is_Aquifer = list(
    aliases = c("aquifer", "is_aquifer"),
    col     = "Is_Aquifer",
    desc    = "Composite binary indicator (0/1) for Chalk/Limestone or Sandstone aquifer",
    default = "Single-basin fits (default)",
    eval_fn = function(df) as.numeric(df$Is_Aquifer)
  ),
  Easting_z = list(
    aliases = c("easting", "mid_x", "easting_z"),
    col     = "mid_X",
    desc    = "Reach midpoint BNG Easting coordinate (m, z-score standardized)",
    default = "All domains",
    eval_fn = function(df) as.numeric(scale(df$mid_X))
  ),
  Northing_z = list(
    aliases = c("northing", "mid_y", "northing_z"),
    col     = "mid_Y",
    desc    = "Reach midpoint BNG Northing coordinate (m, z-score standardized)",
    default = "All domains",
    eval_fn = function(df) as.numeric(scale(df$mid_Y))
  ),
  Sinuosity_z = list(
    aliases = c("sinuosity", "sinuosity_z"),
    col     = "Sinuosity",
    desc    = "Mean channel sinuosity ratio (z-score standardized)",
    default = "Optional (--vars=+sinuosity)",
    eval_fn = function(df) as.numeric(scale(df$Sinuosity))
  ),
  SrcElevation_z = list(
    aliases = c("src_elevation", "srcelevation", "src_elev", "src_elevation_z"),
    col     = "SrcElevation",
    desc    = "Mean headwater source elevation (m AOD, z-score standardized)",
    default = "Optional (--vars=+src_elevation)",
    eval_fn = function(df) as.numeric(scale(df$SrcElevation))
  ),
  Dist2Mth_z = list(
    aliases = c("dist2mth", "dist2mth_km", "dist2mth_z"),
    col     = "Dist2Mth_km",
    desc    = "Hydrologic distance from upstream node to river mouth (km, z-score)",
    default = "Optional (--vars=+dist2mth)",
    eval_fn = function(df) as.numeric(scale(df$Dist2Mth_km))
  ),
  Length_z = list(
    aliases = c("length", "length_km", "length_z"),
    col     = "length_km",
    desc    = "Reach polyline length (km, z-score standardized)",
    default = "Optional (--vars=+length)",
    eval_fn = function(df) as.numeric(scale(df$length_km))
  ),
  log10_Shreve_z = list(
    aliases = c("shreve", "log10_shreve", "log10_shreve_z"),
    col     = "Shreve",
    desc    = "log10 Shreve stream magnitude (z-score standardized)",
    default = "Optional (--vars=+shreve)",
    eval_fn = function(df) as.numeric(scale(log10(pmax(1, df$Shreve))))
  ),
  Strahler_z = list(
    aliases = c("strahler", "strahler_z"),
    col     = "Strahler",
    desc    = "Strahler stream order (1..7, z-score standardized)",
    default = "Optional (--vars=+strahler)",
    eval_fn = function(df) as.numeric(scale(df$Strahler))
  )
)

print_available_covariates <- function() {
  cat("Available Fixed-Effect Regression Covariates (--vars=<list>):\n\n")
  cat(sprintf("%-20s | %-28s | %-18s | %s\n",
              "Model Term", "Short CLI Aliases", "Default Status", "Description"))
  cat(strrep("-", 118), "\n")
  for (nm in names(COVARIATE_CATALOG)) {
    info <- COVARIATE_CATALOG[[nm]]
    al_str <- paste(head(info$aliases, 2), collapse = ", ")
    cat(sprintf("%-20s | %-28s | %-18s | %s\n",
                nm, al_str, substr(info$default, 1, 18), info$desc))
  }
  cat("\nSpecial --vars Modes:\n")
  cat("  --vars=default                     # Default covariates (Elevation, Slope, US_Accum, Geology, Easting, Northing)\n")
  cat("  --vars=elevation,slope,aquifer     # Explicit subset of covariates (Intercept is always included)\n")
  cat("  --vars=-easting,-northing          # Start from default covariates and remove Easting & Northing\n")
  cat("  --vars=+sinuosity,+strahler        # Start from default covariates and add Sinuosity & Strahler order\n")
  cat("  --vars=none                        # Intercept-only model (pure spatial stream network kriging)\n")
  cat("  --vars=all                         # All available ORN topographic, channel, geological, and coordinate covariates\n")
}

match_covariate_token <- function(tok) {
  tok_clean <- tolower(trimws(tok))
  for (nm in names(COVARIATE_CATALOG)) {
    if (tok_clean == tolower(nm) || tok_clean %in% COVARIATE_CATALOG[[nm]]$aliases) {
      return(nm)
    }
  }
  valid_aliases <- sapply(names(COVARIATE_CATALOG), function(nm) COVARIATE_CATALOG[[nm]]$aliases[1])
  stop(sprintf("Unknown fixed-effect covariate '%s'. Run with --list-vars to see valid options (%s).",
               tok, paste(valid_aliases, collapse = ", ")))
}

resolve_covariate_terms <- function(df, include_peat = TRUE, var_specs = NULL) {
  has_split_geol <- include_peat &&
    ("Is_Peat" %in% colnames(df)) && ("Is_Sandstone" %in% colnames(df)) &&
    sd(df$Is_Peat[df$has_obs == 1], na.rm = TRUE) > 0 &&
    sd(df$Is_Sandstone[df$has_obs == 1], na.rm = TRUE) > 0

  default_terms <- if (has_split_geol) {
    c("Elevation_z", "Slope_z", "log10_US_Accum_z",
      "Is_Chalk_Limestone", "Is_Sandstone", "Is_Peat",
      "Easting_z", "Northing_z")
  } else {
    c("Elevation_z", "Slope_z", "log10_US_Accum_z",
      "Is_Aquifer", "Easting_z", "Northing_z")
  }

  if (is.null(var_specs) || length(var_specs) == 0) {
    return(default_terms)
  }
  specs <- trimws(var_specs[nchar(trimws(var_specs)) > 0])
  if (length(specs) == 0 || (length(specs) == 1 && tolower(specs[1]) == "default")) {
    return(default_terms)
  }
  if (length(specs) == 1 && tolower(specs[1]) %in% c("none", "intercept", "null", "0")) {
    return(character(0))
  }
  if (length(specs) == 1 && tolower(specs[1]) == "all") {
    all_t <- c("Elevation_z", "Slope_z", "log10_US_Accum_z",
               "Sinuosity_z", "SrcElevation_z", "Dist2Mth_z", "Length_z",
               "log10_Shreve_z", "Strahler_z",
               if (has_split_geol) c("Is_Chalk_Limestone", "Is_Sandstone", "Is_Peat") else "Is_Aquifer",
               "Easting_z", "Northing_z")
    return(all_t[sapply(all_t, function(nm) COVARIATE_CATALOG[[nm]]$col %in% colnames(df))])
  }

  # Check if relative modifier mode (+var / -var) is used
  is_mod <- grepl("^[+-]", specs)
  if (all(is_mod)) {
    terms <- default_terms
    for (sp in specs) {
      op  <- substr(sp, 1, 1)
      nm  <- match_covariate_token(substr(sp, 2, nchar(sp)))
      if (op == "+") {
        if (!(nm %in% terms)) terms <- c(terms, nm)
      } else {
        # If removing 'aquifer', also remove split geology terms (and vice versa) if user wrote -aquifer or -geology
        if (nm == "Is_Aquifer") {
          terms <- setdiff(terms, c("Is_Aquifer", "Is_Chalk_Limestone", "Is_Sandstone"))
        } else {
          terms <- setdiff(terms, nm)
        }
      }
    }
    return(terms)
  }

  # Explicit list of covariates (also allowing + prefix)
  terms <- character(0)
  for (sp in specs) {
    sp_clean <- sub("^\\+", "", sp)
    if (tolower(sp_clean) %in% c("intercept", "none")) next
    nm <- match_covariate_token(sp_clean)
    if (!(nm %in% terms)) terms <- c(terms, nm)
  }
  terms
}

build_design_matrix <- function(df, include_peat = TRUE, var_specs = NULL) {
  terms <- resolve_covariate_terms(df, include_peat = include_peat, var_specs = var_specs)
  X <- matrix(1.0, nrow = nrow(df), ncol = 1L, dimnames = list(NULL, "Intercept"))
  if (length(terms) == 0L) {
    return(X)
  }
  col_list <- list(Intercept = rep(1.0, nrow(df)))
  for (nm in terms) {
    info <- COVARIATE_CATALOG[[nm]]
    if (!(info$col %in% colnames(df))) {
      warning(sprintf("Column '%s' for covariate '%s' not found in dataframe; skipping.", info$col, nm))
      next
    }
    vec <- info$eval_fn(df)
    if (all(!is.na(vec)) && sd(vec[df$has_obs == 1], na.rm = TRUE) > 0) {
      col_list[[nm]] <- vec
    }
  }
  do.call(cbind, col_list)
}

COVARIATE_SHORT_SLUG <- c(
  Elevation_z        = "elev",
  Slope_z            = "slope",
  log10_US_Accum_z   = "accum",
  Is_Chalk_Limestone = "chalk",
  Is_Sandstone       = "sand",
  Is_Peat            = "peat",
  Is_Aquifer         = "aquifer",
  Easting_z          = "east",
  Northing_z         = "north",
  Sinuosity_z        = "sinuosity",
  SrcElevation_z     = "srcelev",
  Dist2Mth_z         = "dist2mth",
  Length_z           = "length",
  log10_Shreve_z     = "shreve",
  Strahler_z         = "strahler"
)

COVARIATE_SHORT_LABEL <- c(
  Elevation_z        = "Elevation",
  Slope_z            = "Slope",
  log10_US_Accum_z   = "US_Accum",
  Is_Chalk_Limestone = "Chalk/Limestone",
  Is_Sandstone       = "Sandstone",
  Is_Peat            = "Peat",
  Is_Aquifer         = "Aquifer",
  Easting_z          = "Easting",
  Northing_z         = "Northing",
  Sinuosity_z        = "Sinuosity",
  SrcElevation_z     = "SrcElevation",
  Dist2Mth_z         = "Dist2Mth",
  Length_z           = "Length",
  log10_Shreve_z     = "Shreve",
  Strahler_z         = "Strahler"
)

format_covariates_folder_tag <- function(var_specs = NULL) {
  if (is.null(var_specs) || length(var_specs) == 0) return("")
  specs <- trimws(var_specs[nchar(trimws(var_specs)) > 0])
  if (length(specs) == 0 || (length(specs) == 1 && tolower(specs[1]) == "default")) return("")
  if (length(specs) == 1 && tolower(specs[1]) %in% c("none", "intercept", "null", "0")) return("vars_none")
  if (length(specs) == 1 && tolower(specs[1]) == "all") return("vars_all")

  is_mod <- grepl("^[+-]", specs)
  if (all(is_mod)) {
    parts <- sapply(specs, function(sp) {
      op <- substr(sp, 1, 1)
      nm <- match_covariate_token(substr(sp, 2, nchar(sp)))
      sl <- unname(COVARIATE_SHORT_SLUG[nm])
      paste0(if (op == "+") "add_" else "no_", sl)
    })
    return(paste0("vars_", paste(parts, collapse = "_")))
  }

  terms <- character(0)
  for (sp in specs) {
    sp_clean <- sub("^\\+", "", sp)
    if (tolower(sp_clean) %in% c("intercept", "none")) next
    nm <- match_covariate_token(sp_clean)
    if (!(nm %in% terms)) terms <- c(terms, nm)
  }
  if (length(terms) == 0) return("vars_none")
  slugs <- unname(COVARIATE_SHORT_SLUG[terms])
  if (length(slugs) > 5) {
    slugs <- c(slugs[1:4], paste0("plus", length(slugs) - 4))
  }
  paste0("vars_", paste(slugs, collapse = "_"))
}

format_transform_folder_tag <- function(trans_obj, default_type = "log") {
  if (is.null(trans_obj)) return("")
  tr <- as_transformer(trans_obj)
  if (tr$type == default_type) return("")
  paste0("trans_", tr$folder_tag)
}

format_covariates_subtitle <- function(var_specs = NULL, param_tab = NULL) {
  default_eng <- c("Elevation_z", "Slope_z", "log10_US_Accum_z",
                   "Is_Chalk_Limestone", "Is_Sandstone", "Is_Peat",
                   "Easting_z", "Northing_z")
  default_bas <- c("Elevation_z", "Slope_z", "log10_US_Accum_z",
                   "Is_Aquifer", "Easting_z", "Northing_z")

  has_cli_vars <- !is.null(var_specs) && length(var_specs) > 0 &&
                  !(length(var_specs) == 1 && tolower(trimws(var_specs[1])) == "default")

  if (!is.null(param_tab) && nrow(param_tab) > 0) {
    eng_p <- unique(as.character(param_tab$Parameter[
      grepl("All England", param_tab$Domain) &
      !(param_tab$Parameter %in% c("Intercept", "sigma.sq", "tau.sq", "phi"))
    ]))
    if (!has_cli_vars && (identical(eng_p, default_eng) || identical(eng_p, default_bas))) {
      return(NULL)
    }
    if (length(eng_p) == 0) {
      return("None (Intercept-only spatial model)")
    }
    lbls <- ifelse(eng_p %in% names(COVARIATE_SHORT_LABEL), COVARIATE_SHORT_LABEL[eng_p], eng_p)
    return(paste(lbls, collapse = ", "))
  }

  if (!has_cli_vars) return(NULL)
  specs <- trimws(var_specs[nchar(trimws(var_specs)) > 0])
  if (length(specs) == 1 && tolower(specs[1]) %in% c("none", "intercept", "null", "0")) {
    return("None (Intercept-only spatial model)")
  }
  if (length(specs) == 1 && tolower(specs[1]) == "all") {
    return(paste(COVARIATE_SHORT_LABEL[names(COVARIATE_CATALOG) != "Is_Aquifer"], collapse = ", "))
  }
  is_mod <- grepl("^[+-]", specs)
  terms <- default_eng
  if (all(is_mod)) {
    for (sp in specs) {
      op <- substr(sp, 1, 1)
      nm <- match_covariate_token(substr(sp, 2, nchar(sp)))
      if (op == "+") {
        if (!(nm %in% terms)) terms <- c(terms, nm)
      } else {
        if (nm == "Is_Aquifer") {
          terms <- setdiff(terms, c("Is_Aquifer", "Is_Chalk_Limestone", "Is_Sandstone"))
        } else {
          terms <- setdiff(terms, nm)
        }
      }
    }
  } else {
    terms <- character(0)
    for (sp in specs) {
      sp_clean <- sub("^\\+", "", sp)
      if (tolower(sp_clean) %in% c("intercept", "none")) next
      nm <- match_covariate_token(sp_clean)
      if (!(nm %in% terms)) terms <- c(terms, nm)
    }
  }
  if (length(terms) == 0) return("None (Intercept-only spatial model)")
  paste(unname(COVARIATE_SHORT_LABEL[terms]), collapse = ", ")
}

# ==============================================================================
# WATER NETWORK METADATA & CORE BASIN RESOLUTION HELPERS
# ==============================================================================

NETWORK_ALIAS_MAP <- c(
  orn            = "orn",
  orn_v2         = "orn",
  os_orn         = "orn",
  euhydro_thames = "euhydro_thames",
  thames         = "euhydro_thames",
  eu_hydro_thames = "euhydro_thames",
  euhydro_tweed  = "euhydro_tweed",
  tweed          = "euhydro_tweed",
  eu_hydro_tweed = "euhydro_tweed",
  euhydro        = "euhydro",
  eu_hydro       = "euhydro"
)

parse_json_string_field <- function(lines, key, default = "") {
  pat <- paste0('"', key, '"\\s*:\\s*"([^"]*)"')
  hit <- grep(pat, lines, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0('^.*', pat, '.*$'), "\\1", hit[1])
}

parse_json_int_field <- function(lines, key, default = NA_integer_) {
  pat <- paste0('"', key, '"\\s*:\\s*([0-9]+)')
  hit <- grep(pat, lines, value = TRUE)
  if (length(hit) == 0) return(default)
  as.integer(sub(paste0('^.*', pat, '.*$'), "\\1", hit[1]))
}

resolve_network_processed_dir <- function(base_dir, network_arg = NULL, data_dir = NULL, gis_dir = NULL, source = "rsn") {
  src_tag   <- resolve_source_tag(source)
  src_suff  <- if (src_tag == "rsn") "" else paste0("_", src_tag)
  proc_root <- if (!is.null(data_dir) && nzchar(data_dir)) {
    cand_p <- rel_path(data_dir, "processed")
    if (dir.exists(cand_p)) cand_p else if (file.exists(rel_path(data_dir, "england_rsn_reaches.csv"))) data_dir else cand_p
  } else {
    rel_path(base_dir, "data", "processed")
  }

  if (is.null(network_arg) || !nzchar(trimws(network_arg))) {
    # Check active network_metadata.json in proc_root, defaulting to 'orn'
    meta_p <- rel_path(proc_root, "network_metadata.json")
    default_abbr <- "orn"
    if (file.exists(meta_p)) {
      lines <- readLines(meta_p, warn = FALSE)
      default_abbr <- parse_json_string_field(lines, "abbr", "orn")
    }
    abbr <- default_abbr
    net_cli_arg <- abbr
  } else {
    raw_net_orig <- trimws(network_arg)
    raw_net <- tolower(raw_net_orig)
    if (raw_net %in% names(NETWORK_ALIAS_MAP)) {
      abbr <- unname(NETWORK_ALIAS_MAP[raw_net])
      net_cli_arg <- abbr
    } else if (grepl("\\.(gpkg|parquet|geoparquet|pq|json)$", raw_net)) {
      stem <- sub("\\.(gpkg|parquet|geoparquet|pq|json)$", "", basename(raw_net))
      stem <- sub("\\.geo$", "", stem)
      abbr <- tolower(stem)
      net_cli_arg <- raw_net_orig
    } else {
      abbr <- raw_net
      net_cli_arg <- raw_net_orig
    }
  }

  sub_d <- rel_path(proc_root, paste0(abbr, src_suff))
  # If default RSN and sub_d doesn't exist, check root proc_root fallback
  if (src_tag == "rsn" && !dir.exists(sub_d) && file.exists(rel_path(proc_root, "england_rsn_reaches.csv"))) {
    return(proc_root)
  }

  if (!file.exists(rel_path(sub_d, "england_rsn_reaches.csv"))) {
    cat(sprintf("Preprocessing water network '%s' (source='%s') via scripts/01_preprocess_orn.py...\n", net_cli_arg, src_tag))
    extra_py_args <- ""
    if (!is.null(gis_dir) && nzchar(gis_dir)) {
      extra_py_args <- paste0(extra_py_args, " --gis-dir=", shQuote(gis_dir))
    }
    if (!is.null(data_dir) && nzchar(data_dir)) {
      extra_py_args <- paste0(extra_py_args, " --data-dir=", shQuote(data_dir))
    }
    cmd <- sprintf("python3 %s --network=%s --source=%s%s",
                   shQuote(rel_path(base_dir, "scripts", "01_preprocess_orn.py")),
                   shQuote(net_cli_arg),
                   shQuote(src_tag),
                   extra_py_args)
    status <- system(cmd)
    if (status != 0L) {
      stop(sprintf("Failed to preprocess water network '%s' with source '%s'. Run `python3 scripts/01_preprocess_orn.py --list-networks`.", net_cli_arg, src_tag))
    }
  }
  if (file.exists(rel_path(sub_d, "england_rsn_reaches.csv"))) {
    return(sub_d)
  }
  if (src_tag == "rsn" && file.exists(rel_path(proc_root, "england_rsn_reaches.csv"))) {
    return(proc_root)
  }
  sub_d
}

load_network_metadata <- function(proc_dir) {
  meta_p <- file.path(proc_dir, "network_metadata.json")
  if (!file.exists(meta_p)) {
    return(list(
      abbr         = "orn",
      label        = "OS Open Rivers Network (ORN v2)",
      short_label  = "ORN v2",
      source       = "rsn",
      source_abbr  = "RSN",
      source_label = "River Surveillance Network (RSN)",
      humber_cid   = 80L,
      thames_cid   = 3055L,
      severn_cid   = 2786L
    ))
  }
  lines <- readLines(meta_p, warn = FALSE)
  src_key <- parse_json_string_field(lines, "source", "rsn")
  list(
    abbr         = parse_json_string_field(lines, "abbr", "orn"),
    label        = parse_json_string_field(lines, "label", "OS Open Rivers Network (ORN v2)"),
    short_label  = parse_json_string_field(lines, "short_label", "ORN v2"),
    source       = src_key,
    source_abbr  = parse_json_string_field(lines, "source_abbr", format_source_abbr(src_key)),
    source_label = parse_json_string_field(lines, "source_label", format_source_label(src_key)),
    humber_cid   = parse_json_int_field(lines, "humber", NA_integer_),
    thames_cid   = parse_json_int_field(lines, "thames", NA_integer_),
    severn_cid   = parse_json_int_field(lines, "severn", NA_integer_)
  )
}

resolve_core_basin <- function(df_reaches, catch_lookup = NULL, net_meta = NULL, target = c("humber", "thames", "severn")) {
  target <- match.arg(target)
  avail_cids <- unique(df_reaches$CatchID)

  # 1. Check explicit core_basins in network_metadata.json
  meta_cid <- if (!is.null(net_meta)) switch(target,
    humber = net_meta$humber_cid,
    thames = net_meta$thames_cid,
    severn = net_meta$severn_cid
  ) else NA_integer_
  if (!is.na(meta_cid) && (meta_cid %in% avail_cids)) {
    lbl <- switch(target,
      humber = sprintf("Humber Basin (CatchID=%d)", meta_cid),
      thames = sprintf("Thames Basin (CatchID=%d)", meta_cid),
      severn = sprintf("Severn Basin (CatchID=%d)", meta_cid)
    )
    return(list(catch_id = meta_cid, domain_label = lbl))
  }

  # 1b. Known CatchIDs across ORN and EU-Hydro
  known_cids <- switch(target,
    humber = c(80L, 75L),
    thames = c(3055L, 306L),
    severn = c(2786L, 243L)
  )
  hit_known <- known_cids[known_cids %in% avail_cids]
  if (length(hit_known) > 0) {
    cid <- hit_known[1]
    lbl <- switch(target,
      humber = sprintf("Humber Basin (CatchID=%d)", cid),
      thames = sprintf("Thames Basin (CatchID=%d)", cid),
      severn = sprintf("Severn Basin (CatchID=%d)", cid)
    )
    return(list(catch_id = cid, domain_label = lbl))
  }

  # 2. Search catch_lookup by Main_River or RBD_Name
  if (!is.null(catch_lookup) && nrow(catch_lookup) > 0) {
    pat <- switch(target,
      humber = "humber|trent|ouse",
      thames = "thames",
      severn = "severn"
    )
    hits <- which(grepl(pat, catch_lookup$Main_River, ignore.case = TRUE) |
                  grepl(pat, catch_lookup$RBD_Name, ignore.case = TRUE))
    if (length(hits) > 0) {
      sub_lu <- catch_lookup[hits, ]
      sub_lu <- sub_lu[order(-sub_lu$n_obs_reaches, -sub_lu$n_reaches), ]
      cid <- sub_lu$CatchID[1]
      lbl <- switch(target,
        humber = sprintf("Humber Basin (CatchID=%d)", cid),
        thames = sprintf("Thames Basin (CatchID=%d)", cid),
        severn = sprintf("Severn Basin (CatchID=%d)", cid)
      )
      return(list(catch_id = cid, domain_label = lbl))
    }
  }

  # 3. Fallback to largest monitored catchment
  tab_c <- sort(table(df_reaches$CatchID), decreasing = TRUE)
  cid <- as.integer(names(tab_c)[1])
  lbl <- sprintf("%s Basin (CatchID=%d)", paste0(toupper(substr(target, 1, 1)), substr(target, 2, nchar(target))), cid)
  list(catch_id = cid, domain_label = lbl)
}

extract_network_label_from_tables <- function(cv_tab = NULL, param_tab = NULL, dir_path = NULL) {
  if (!is.null(cv_tab) && "Network_Label" %in% colnames(cv_tab) && nrow(cv_tab) > 0) {
    val <- as.character(cv_tab$Network_Label[1])
    if (!is.na(val) && nzchar(val)) return(val)
  }
  if (!is.null(param_tab) && "Network_Label" %in% colnames(param_tab) && nrow(param_tab) > 0) {
    val <- as.character(param_tab$Network_Label[1])
    if (!is.na(val) && nzchar(val)) return(val)
  }
  if (!is.null(dir_path)) {
    bname <- basename(normalizePath(dir_path, mustWork = FALSE))
    if (grepl("euhydro_thames", bname)) return("EU-Hydro (Thames)")
    if (grepl("euhydro_tweed", bname))  return("EU-Hydro (Tweed)")
    if (grepl("euhydro", bname))        return("EU-Hydro v1.3")
    if (grepl("_orn", bname))           return("ORN v2")
  }
  "ORN v2"
}


