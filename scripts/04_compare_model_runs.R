#!/usr/bin/env Rscript
# ==============================================================================
# 04_compare_model_runs.R
#
# Standalone script and reusable R function (`compare_s3n_runs`) for comparing
# two or more S3N model output directories:
#   - Cross-network comparisons (e.g., ORN v2 vs. EU-Hydro Thames)
#   - Fixed-effect covariate comparisons (e.g., default vs. custom --vars=...)
#   - Temporal comparisons (e.g., different years, seasons, or months)
#   - Cross-determinand comparisons & reach-level relationship analyses
#     (e.g., Nitrate [0117] vs. Orthophosphate [0180] or RHS HQA [hqa_score])
#
# Performs standard statistical tests (paired Student's t-test / Diebold-Mariano
# loss differential test, paired Wilcoxon signed-rank test, Welch's two-sample
# t-test, Mann-Whitney U test, two-sample Kolmogorov-Smirnov test, Fisher's z
# test on CV R^2, Fligner-Killeen variance test, Two-Sample Wald Z-test on
# bootstrap parameter CIs, Chi-Squared 4-quadrant co-occurrence test, OLS
# reach-scaling regression, and One-Way ANOVA / Kruskal-Wallis tests),
# computes comprehensive descriptive and stratified output statistics,
# generates a Markdown report (`comparison_report.md`), exports summary CSV
# tables, and renders 5 publication-quality base-R figures:
#   - fig_comp1_cv_performance.png
#   - fig_comp2_predictions_and_uncertainty.png
#   - fig_comp3_parameter_forest.png
#   - fig_comp4_side_by_side_maps.png
#   - fig_comp5_reach_relationships.png
#
# Usage:
#   Rscript scripts/04_compare_model_runs.R --list
#   Rscript scripts/04_compare_model_runs.R output/det_0117 output/det_0180_orn
#   Rscript scripts/04_compare_model_runs.R output/det_0117 output/det_0117_euhydro_thames
#   Rscript scripts/04_compare_model_runs.R output/det_0117 output/det_0117_vars_elev_slope_aquifer
#   Rscript scripts/04_compare_model_runs.R output/det_0117 --by=season
# ==============================================================================

BASE_DIR <- "/home/ben/code/r/spatial_stream_model"
source(file.path(BASE_DIR, "scripts/s3n_plot_functions.R"))

# ------------------------------------------------------------------------------
# Helper: Format p-values and significance stars
# ------------------------------------------------------------------------------
fmt_pval <- function(p) {
  if (is.null(p) || length(p) == 0 || is.na(p)) return("NA")
  if (p < 1e-16) return("< 1e-16")
  if (p < 1e-4)  return(sprintf("%.2e", p))
  sprintf("%.4f", p)
}

sig_stars <- function(p) {
  if (is.null(p) || length(p) == 0 || is.na(p)) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  if (p < 0.10)  return(".")
  "ns"
}

safe_paired_ttest <- function(x, y) {
  d <- x - y
  d <- d[!is.na(d)]
  if (length(d) < 3 || sd(d) < 1e-12) {
    return(list(statistic = 0.0, p.value = 1.0))
  }
  suppressWarnings(t.test(x, y, paired = TRUE))
}

safe_paired_wilcox <- function(x, y) {
  d <- x - y
  d <- d[!is.na(d) & abs(d) > 1e-12]
  if (length(d) < 3) {
    return(list(statistic = 0.0, p.value = 1.0))
  }
  suppressWarnings(wilcox.test(x, y, paired = TRUE))
}

canonical_domain <- function(dom_vec) {
  d <- as.character(dom_vec)
  d <- sub(" \\(CatchID=[0-9]+\\)$", "", d)
  d <- sub("^River ", "", d)
  trimws(d)
}

canonical_model <- function(mod_vec) {
  m <- as.character(mod_vec)
  m[grepl("OLS", m, ignore.case = TRUE)]       <- "OLS (Non-Spatial)"
  m[grepl("Tail-Down", m, ignore.case = TRUE) & !grepl("Combined", m, ignore.case = TRUE)] <- "S3N Tail-Down"
  m[grepl("Tail-Up", m, ignore.case = TRUE)   & !grepl("Combined", m, ignore.case = TRUE)] <- "S3N Tail-Up"
  m[grepl("Combined", m, ignore.case = TRUE)]  <- "S3N Combined"
  m
}

find_first_file <- function(dir_path, pattern, exclude_pattern = NULL) {
  if (!dir.exists(dir_path)) return(NULL)
  f <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  if (!is.null(exclude_pattern) && length(f) > 0) {
    f <- f[!grepl(exclude_pattern, basename(f))]
  }
  if (length(f) == 0) return(NULL)
  f[1]
}

calc_sample_skewness <- function(x) {
  x <- x[!is.na(x) & is.finite(x)]
  n <- length(x)
  if (n < 3) return(NA_real_)
  s <- sd(x)
  if (s < 1e-12) return(0.0)
  mean(((x - mean(x)) / s)^3)
}

calc_sample_kurtosis <- function(x) {
  x <- x[!is.na(x) & is.finite(x)]
  n <- length(x)
  if (n < 4) return(NA_real_)
  s <- sd(x)
  if (s < 1e-12) return(0.0)
  mean(((x - mean(x)) / s)^4) - 3.0
}

calc_geom_mean <- function(x) {
  x <- x[!is.na(x) & is.finite(x)]
  if (length(x) == 0 || any(x <= 0)) return(NA_real_)
  exp(mean(log(x)))
}

summarize_numeric_vector <- function(x) {
  x <- x[!is.na(x) & is.finite(x)]
  n <- length(x)
  if (n == 0) {
    return(data.frame(
      N = 0L, Mean = NA_real_, SD = NA_real_, CV_Pct = NA_real_,
      Geom_Mean = NA_real_, Skewness = NA_real_, Excess_Kurtosis = NA_real_,
      Min = NA_real_, P05 = NA_real_, P10 = NA_real_, Q25 = NA_real_,
      Median = NA_real_, Q75 = NA_real_, P90 = NA_real_, P95 = NA_real_,
      Max = NA_real_, IQR = NA_real_, MAD = NA_real_, stringsAsFactors = FALSE
    ))
  }
  qs <- as.numeric(quantile(x, probs = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95), na.rm = TRUE))
  m_val <- mean(x)
  s_val <- if (n >= 2) sd(x) else 0.0
  data.frame(
    N               = n,
    Mean            = round(m_val, 5),
    SD              = round(s_val, 5),
    CV_Pct          = if (abs(m_val) > 1e-9) round(100.0 * s_val / abs(m_val), 2) else NA_real_,
    Geom_Mean       = round(calc_geom_mean(x), 5),
    Skewness        = round(calc_sample_skewness(x), 4),
    Excess_Kurtosis = round(calc_sample_kurtosis(x), 4),
    Min             = round(min(x), 5),
    P05             = round(qs[1], 5),
    P10             = round(qs[2], 5),
    Q25             = round(qs[3], 5),
    Median          = round(qs[4], 5),
    Q75             = round(qs[5], 5),
    P90             = round(qs[6], 5),
    P95             = round(qs[7], 5),
    Max             = round(max(x), 5),
    IQR             = round(qs[5] - qs[3], 5),
    MAD             = round(mad(x, na.rm = TRUE), 5),
    stringsAsFactors = FALSE
  )
}

read_basin_catch_id <- function(file_path, fallback_id = NA_integer_) {
  if (!is.null(file_path) && file.exists(file_path)) {
    r1 <- tryCatch(read.csv(file_path, nrows = 1, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(r1) && ("CatchID" %in% names(r1)) && nrow(r1) > 0) {
      return(as.integer(r1$CatchID[1]))
    }
  }
  fallback_id
}

# ------------------------------------------------------------------------------
# Helper: Load a single model run bundle from an output directory
# ------------------------------------------------------------------------------
load_model_run_bundle <- function(dir_path, custom_label = NULL) {
  dir_p <- if (grepl("^/", dir_path)) dir_path else file.path(BASE_DIR, dir_path)
  if (!dir.exists(dir_p)) {
    stop(sprintf("Model output directory does not exist: %s", dir_p))
  }

  # Support both standard output/ layout (predictions/, tables/) and det_* layout
  if (dir.exists(file.path(dir_p, "predictions")) && dir.exists(file.path(dir_p, "tables"))) {
    pred_dir  <- file.path(dir_p, "predictions")
    tab_dir   <- file.path(dir_p, "tables")
    eng_file  <- find_first_file(pred_dir, "^england_.*predictions\\.csv$", "cv_predictions")
    hum_file  <- find_first_file(pred_dir, "^humber_.*predictions\\.csv$", "cv_predictions")
    tha_file  <- find_first_file(pred_dir, "^thames_.*predictions\\.csv$", "cv_predictions")
    sev_file  <- find_first_file(pred_dir, "^severn_.*predictions\\.csv$", "cv_predictions")
    cv_t_file <- file.path(tab_dir, "cv_model_comparison.csv")
    par_file  <- file.path(tab_dir, "s3n_parameter_estimates.csv")
    cv_e_file <- find_first_file(pred_dir, "^england_.*cv_predictions\\.csv$")
    cv_h_file <- find_first_file(pred_dir, "^humber_.*cv_predictions\\.csv$")
  } else {
    eng_file  <- find_first_file(dir_p, "^predictions_.*\\.csv$",
                                 "^predictions_(humber|thames|severn|river_|catchid_).*\\.csv$")
    hum_file  <- find_first_file(dir_p, "^predictions_humber_.*\\.csv$")
    tha_file  <- find_first_file(dir_p, "^predictions_thames_.*\\.csv$")
    sev_file  <- find_first_file(dir_p, "^predictions_severn_.*\\.csv$")
    cv_t_file <- find_first_file(dir_p, "^cv_comparison_.*\\.csv$")
    par_file  <- find_first_file(dir_p, "^parameters_.*\\.csv$")
    cv_e_file <- find_first_file(dir_p, "^cv_predictions_england_.*\\.csv$")
    cv_h_file <- find_first_file(dir_p, "^cv_predictions_humber_.*\\.csv$")
  }

  if (is.null(eng_file) || is.null(cv_t_file) || !file.exists(cv_t_file)) {
    stop(sprintf("Directory '%s' does not contain required predictions_*.csv and cv_comparison_*.csv files.", dir_p))
  }

  eng_pred  <- normalize_pred_df(read.csv(eng_file, stringsAsFactors = FALSE))
  cv_tab    <- normalize_cv_tab(read.csv(cv_t_file, stringsAsFactors = FALSE))
  param_tab <- if (!is.null(par_file) && file.exists(par_file)) read.csv(par_file, stringsAsFactors = FALSE) else NULL
  eng_cv    <- if (!is.null(cv_e_file) && file.exists(cv_e_file)) normalize_cv_pts(read.csv(cv_e_file, stringsAsFactors = FALSE)) else NULL
  hum_cv    <- if (!is.null(cv_h_file) && file.exists(cv_h_file)) normalize_cv_pts(read.csv(cv_h_file, stringsAsFactors = FALSE)) else NULL

  hum_cid <- read_basin_catch_id(hum_file, 80L)
  tha_cid <- read_basin_catch_id(tha_file, 3055L)
  sev_cid <- read_basin_catch_id(sev_file, 2786L)

  cv_tab$Canon_Domain <- canonical_domain(cv_tab$Domain)
  cv_tab$Canon_Model  <- canonical_model(cv_tab$Model)

  if (!is.null(param_tab)) {
    param_tab$Canon_Domain <- canonical_domain(param_tab$Domain)
    param_tab$Canon_Model  <- canonical_model(param_tab$Model)
  }

  # Extract metadata
  net_lbl  <- extract_network_label_from_tables(cv_tab, param_tab, dir_p)
  net_abbr <- if ("Network_Abbr" %in% names(cv_tab)) as.character(cv_tab$Network_Abbr[1]) else {
    if (grepl("euhydro_thames", dir_p)) "euhydro_thames"
    else if (grepl("euhydro_tweed", dir_p)) "euhydro_tweed"
    else if (grepl("euhydro", dir_p)) "euhydro"
    else "orn"
  }
  period_lbl <- if ("Period" %in% names(cv_tab) && !is.na(cv_tab$Period[1])) as.character(cv_tab$Period[1]) else {
    bname <- basename(dir_p)
    if (grepl("^(year|season|month)_", bname)) bname else "all"
  }
  det_code <- if ("Determinand_Code" %in% names(cv_tab)) {
    c0 <- as.character(cv_tab$Determinand_Code[1])
    if (nchar(c0) < 4 && grepl("^[0-9]+$", c0)) sprintf("%04d", as.integer(c0)) else c0
  } else {
    sub("^predictions_(.*)\\.csv$", "\\1", basename(eng_file))
  }
  known_names <- c(
    "0117" = "Nitrate as N",
    "0180" = "Orthophosphate as P",
    "0111" = "Ammoniacal Nitrogen as N",
    "0116" = "Total Oxidised Nitrogen",
    "0118" = "Nitrite as N",
    "0301" = "Dissolved Organic Carbon",
    "0077" = "Conductivity",
    "0162" = "Alkalinity",
    "0061" = "pH",
    "0076" = "Water Temperature",
    "0085" = "Dissolved Oxygen (% Sat)",
    "9924" = "Dissolved Oxygen",
    "whpt_aspt" = "WHPT ASPT Score",
    "whpt_ntaxa" = "WHPT N-Taxa Richness",
    "life_family_index" = "LIFE Flow Index",
    "psi_family_score" = "PSI Sediment Score",
    "rmni" = "River Macrophyte Nutrient Index (RMNI)",
    "hqa_score" = "RHS Habitat Quality Score (HQA)"
  )
  det_name <- if ("Determinand_Name" %in% names(cv_tab)) {
    as.character(cv_tab$Determinand_Name[1])
  } else if (det_code %in% names(known_names)) {
    unname(known_names[det_code])
  } else {
    det_code
  }
  det_unit <- if ("Unit" %in% names(cv_tab)) as.character(cv_tab$Unit[1]) else "mg/L"
  use_log  <- if ("Scale" %in% names(cv_tab)) (tolower(as.character(cv_tab$Scale[1])) == "log") else {
    obs_pos <- eng_pred$obs_value[eng_pred$has_obs == 1 & !is.na(eng_pred$obs_value) & eng_pred$obs_value > 0]
    sc_pos  <- eng_pred$obs_model_sc[eng_pred$has_obs == 1 & !is.na(eng_pred$obs_value) & eng_pred$obs_value > 0]
    (length(obs_pos) > 0 && max(abs(log(obs_pos) - sc_pos), na.rm = TRUE) < 1e-2)
  }

  vars_sub <- format_covariates_subtitle(param_tab = param_tab)
  user_supplied_label <- (!is.null(custom_label) && nzchar(custom_label))

  # Construct concise human-readable run label if not provided
  if (!user_supplied_label) {
    bname <- basename(normalizePath(dir_p, mustWork = FALSE))
    if (grepl("^(season|year|month)_", bname)) {
      p_clean <- sub("^season_", "Season: ", sub("^year_", "Year ", sub("^month_", "Month ", bname)))
      custom_label <- sprintf("%s (%s)", p_clean, net_lbl)
    } else if (!is.null(vars_sub) && nzchar(vars_sub)) {
      custom_label <- sprintf("%s [Vars: %s]", net_lbl, vars_sub)
    } else {
      custom_label <- net_lbl
    }
  }

  # Attach coordinates to eng_cv via RivID lookup from eng_pred
  if (!is.null(eng_cv) && "RivID" %in% names(eng_cv)) {
    idx <- match(eng_cv$RivID, eng_pred$RivID)
    eng_cv$mid_X     <- eng_pred$mid_X[idx]
    eng_cv$mid_Y     <- eng_pred$mid_Y[idx]
    eng_cv$obs_value <- eng_pred$obs_value[idx]
  }

  list(
    dir_path            = dir_p,
    dir_name            = basename(normalizePath(dir_p, mustWork = FALSE)),
    label               = custom_label,
    user_supplied_label = user_supplied_label,
    net_abbr            = net_abbr,
    net_label           = net_lbl,
    period              = period_lbl,
    det_code            = det_code,
    det_name            = det_name,
    det_unit            = det_unit,
    use_log             = use_log,
    vars_label          = if (is.null(vars_sub)) "Default (All 8 standard)" else vars_sub,
    hum_cid             = hum_cid,
    tha_cid             = tha_cid,
    sev_cid             = sev_cid,
    eng_pred            = eng_pred,
    cv_tab              = cv_tab,
    param_tab           = param_tab,
    eng_cv              = eng_cv,
    hum_cv              = hum_cv
  )
}

# ------------------------------------------------------------------------------
# Helper: Match reaches or observed stations between two runs
# Supports both same-network (exact RivID match) and cross-network (2D BNG
# nearest-neighbor spatial hash matching within max_dist_m)
# ------------------------------------------------------------------------------
match_reaches_between_runs <- function(df1, df2, same_network = FALSE, max_dist_m = 250.0,
                                       require_obs_match = FALSE) {
  if (same_network && ("RivID" %in% names(df1)) && ("RivID" %in% names(df2))) {
    common_ids <- intersect(df1$RivID, df2$RivID)
    if (length(common_ids) > 0) {
      i1 <- match(common_ids, df1$RivID)
      i2 <- match(common_ids, df2$RivID)
      return(data.frame(idx1 = i1, idx2 = i2, dist_m = rep(0.0, length(common_ids))))
    }
  }

  # Spatial grid hash matching in British National Grid (EPSG:27700) meters
  x1 <- df1$mid_X; y1 <- df1$mid_Y
  x2 <- df2$mid_X; y2 <- df2$mid_Y
  valid1 <- which(!is.na(x1) & !is.na(y1))
  valid2 <- which(!is.na(x2) & !is.na(y2))
  if (length(valid1) == 0 || length(valid2) == 0) {
    return(data.frame(idx1 = integer(0), idx2 = integer(0), dist_m = numeric(0)))
  }

  cell_size <- max(250.0, max_dist_m)
  gx2 <- as.integer(floor(x2[valid2] / cell_size))
  gy2 <- as.integer(floor(y2[valid2] / cell_size))
  # Expand target points into 3x3 neighboring grid cells so lookup is a single split()
  offsets_x <- rep(-1L:1L, times = 3L)
  offsets_y <- rep(-1L:1L, each  = 3L)
  rep_v2 <- rep(valid2, times = 9L)
  rep_gx <- rep(gx2, times = 9L) + rep(offsets_x, each = length(valid2))
  rep_gy <- rep(gy2, times = 9L) + rep(offsets_y, each = length(valid2))
  keys2  <- paste(rep_gx, rep_gy, sep = "_")
  bucket2 <- split(rep_v2, keys2)

  gx1 <- as.integer(floor(x1[valid1] / cell_size))
  gy1 <- as.integer(floor(y1[valid1] / cell_size))
  keys1 <- paste(gx1, gy1, sep = "_")
  in_b  <- which(keys1 %in% names(bucket2))
  if (length(in_b) == 0L) {
    return(data.frame(idx1 = integer(0), idx2 = integer(0), dist_m = numeric(0)))
  }
  cands_list <- bucket2[keys1[in_b]]

  out_i1   <- integer(length(in_b))
  out_i2   <- integer(length(in_b))
  out_dist <- rep(NA_real_, length(in_b))
  n_keep   <- 0L
  max_d2   <- max_dist_m * max_dist_m

  has_obs_y <- require_obs_match && ("obs_y" %in% names(df1)) && ("obs_y" %in% names(df2))

  for (m in seq_along(in_b)) {
    k  <- in_b[m]
    i1 <- valid1[k]
    cands <- cands_list[[m]]
    dx <- x2[cands] - x1[i1]
    dy <- y2[cands] - y1[i1]
    d2 <- dx * dx + dy * dy

    if (has_obs_y) {
      dy_obs <- abs(df2$obs_y[cands] - df1$obs_y[i1])
      ok <- which(d2 <= max_d2 & dy_obs < 0.25)
      if (length(ok) == 0L) ok <- which(d2 <= max_d2)
      if (length(ok) == 0L) next
      best_sub <- ok[which.min(d2[ok] + (dy_obs[ok] * 200.0)^2)]
      best_i2  <- cands[best_sub]
      best_d2  <- d2[best_sub]
    } else {
      min_idx <- which.min(d2)
      if (d2[min_idx] > max_d2) next
      best_i2 <- cands[min_idx]
      best_d2 <- d2[min_idx]
    }

    n_keep <- n_keep + 1L
    out_i1[n_keep]   <- i1
    out_i2[n_keep]   <- best_i2
    out_dist[n_keep] <- sqrt(best_d2)
  }

  if (n_keep == 0L) {
    return(data.frame(idx1 = integer(0), idx2 = integer(0), dist_m = numeric(0)))
  }
  res <- data.frame(
    idx1   = out_i1[1:n_keep],
    idx2   = out_i2[1:n_keep],
    dist_m = out_dist[1:n_keep]
  )
  # Deduplicate idx2 by keeping closest match
  res <- res[order(res$dist_m), ]
  res <- res[!duplicated(res$idx2), ]
  res[order(res$idx1), ]
}

# ------------------------------------------------------------------------------
# Helper: List all available model output directories under output/
# ------------------------------------------------------------------------------
list_available_model_runs <- function(out_root = file.path(BASE_DIR, "output")) {
  cat("==============================================================================\n")
  cat("Available Model Run Output Directories for Comparison:\n")
  cat("==============================================================================\n")
  all_dirs <- list.dirs(out_root, recursive = TRUE, full.names = TRUE)
  found <- 0L
  for (d in all_dirs) {
    if (grepl("/comparisons(/|$)", d)) next
    has_det_files <- length(list.files(d, pattern = "^cv_comparison_.*\\.csv$")) > 0 &&
                     length(list.files(d, pattern = "^predictions_.*\\.csv$")) > 0
    has_std_files <- file.exists(file.path(d, "tables/cv_model_comparison.csv"))
    if (has_det_files || has_std_files) {
      found <- found + 1L
      rel_d <- sub(paste0("^", BASE_DIR, "/?"), "", d)
      b <- tryCatch(load_model_run_bundle(d), error = function(e) NULL)
      if (!is.null(b)) {
        r_all <- subset(b$cv_tab, Canon_Domain == "All England RSN Network" & Canon_Model == "S3N Combined")
        r2_str <- if (nrow(r_all) > 0) sprintf("CV R2 = %.3f, RMSE = %.3f", r_all$R2[1], r_all$RMSE[1]) else "N/A"
        cat(sprintf("  [%2d] %-45s | %-20s | %-22s | %6d reaches (%4d obs) | %s\n",
                    found, rel_d, b$net_label, substr(b$det_name, 1, 22),
                    nrow(b$eng_pred), sum(b$eng_pred$has_obs == 1), r2_str))
      } else {
        cat(sprintf("  [%2d] %s\n", found, rel_d))
      }
    }
  }
  cat("==============================================================================\n")
}

# ------------------------------------------------------------------------------
# Helper: Comprehensive Descriptive Statistics of Model Output Results
# ------------------------------------------------------------------------------
compute_detailed_descriptives <- function(bundles) {
  desc_rows <- list()
  rep_rows  <- list()

  for (i in seq_along(bundles)) {
    b  <- bundles[[i]]
    ep <- b$eng_pred
    obs_sub   <- ep[ep$has_obs == 1 & !is.na(ep$obs_value), ]
    unmon_sub <- ep[ep$has_obs == 0, ]

    vars_to_sum <- list(
      list(grp = "Observed Stations (Raw)",            vec = obs_sub$obs_value,       scale = "raw"),
      list(grp = "S3N Combined Predictions (Raw)",     vec = ep$pred_combined,        scale = "raw"),
      list(grp = "Unmonitored Reach Predictions (Raw)",vec = unmon_sub$pred_combined, scale = "raw"),
      list(grp = "S3N Tail-Up Predictions (Raw)",      vec = ep$pred_tu,              scale = "raw"),
      list(grp = "S3N Tail-Down Predictions (Raw)",    vec = ep$pred_td,              scale = "raw"),
      list(grp = "S3N Combined (Model Scale)",         vec = ep$pred_model_sc,        scale = if (b$use_log) "log" else "raw"),
      list(grp = "Kriging Standard Error (SE)",        vec = ep$se_model_sc,          scale = if (b$use_log) "log_se" else "raw_se")
    )

    for (v_info in vars_to_sum) {
      s_df <- summarize_numeric_vector(v_info$vec)
      meta_df <- data.frame(
        Run_Index        = i,
        Run_Label        = b$label,
        Determinand_Code = b$det_code,
        Determinand_Name = b$det_name,
        Unit             = b$det_unit,
        Scale            = v_info$scale,
        Variable_Group   = v_info$grp,
        stringsAsFactors = FALSE
      )
      desc_rows[[length(desc_rows) + 1L]] <- cbind(meta_df, s_df)
    }

    # Representativeness test: Observed stations vs. Unmonitored network predictions
    if (nrow(obs_sub) >= 5 && nrow(unmon_sub) >= 5) {
      mw_rep <- suppressWarnings(wilcox.test(ep$pred_combined, obs_sub$obs_value, paired = FALSE))
      ks_rep <- suppressWarnings(ks.test(ep$pred_combined, obs_sub$obs_value))
      obs_med  <- median(obs_sub$obs_value, na.rm = TRUE)
      pred_med <- median(ep$pred_combined, na.rm = TRUE)
      obs_mn   <- mean(obs_sub$obs_value, na.rm = TRUE)
      pred_mn  <- mean(ep$pred_combined, na.rm = TRUE)
      ratio_med <- if (abs(obs_med) > 1e-9) pred_med / obs_med else NA_real_
      ratio_mn  <- if (abs(obs_mn) > 1e-9)  pred_mn / obs_mn   else NA_real_

      dir_txt <- if (!is.na(ratio_med) && ratio_med < 0.95) {
        "Full network median is lower than monitored station median (unmonitored upland headwaters dilute network average)"
      } else if (!is.na(ratio_med) && ratio_med > 1.05) {
        "Full network median is higher than monitored station median"
      } else {
        "Full network median closely matches monitored station median"
      }

      rep_rows[[length(rep_rows) + 1L]] <- data.frame(
        Run_Index              = i,
        Run_Label              = b$label,
        Determinand_Code       = b$det_code,
        Unit                   = b$det_unit,
        n_Obs_Stations         = nrow(obs_sub),
        n_Total_Reaches        = nrow(ep),
        Obs_Mean_Raw           = round(obs_mn, 5),
        Pred_Mean_Raw          = round(pred_mn, 5),
        Pred_to_Obs_Mean_Ratio = round(ratio_mn, 4),
        Obs_Median_Raw         = round(obs_med, 5),
        Pred_Median_Raw        = round(pred_med, 5),
        Pred_to_Obs_Med_Ratio  = round(ratio_med, 4),
        MW_W_Stat              = as.numeric(mw_rep$statistic),
        MW_P_Value             = mw_rep$p.value,
        KS_D_Stat              = round(as.numeric(ks_rep$statistic), 4),
        KS_P_Value             = ks_rep$p.value,
        Significance           = sig_stars(mw_rep$p.value),
        Interpretation         = dir_txt,
        stringsAsFactors       = FALSE
      )
    }
  }

  list(
    desc_df = do.call(rbind, desc_rows),
    rep_df  = if (length(rep_rows) > 0) do.call(rbind, rep_rows) else NULL
  )
}

# ------------------------------------------------------------------------------
# Helper: Stratified Descriptive Statistics of Model Output Results Across
#         Topography, Hydrogeology, Landscape Regime, Stream Order & Basins
# ------------------------------------------------------------------------------
compute_stratified_descriptives <- function(bundles) {
  strat_rows <- list()

  for (i in seq_along(bundles)) {
    b  <- bundles[[i]]
    ep <- b$eng_pred
    n_tot <- nrow(ep)
    z_all <- (ep$pred_model_sc - mean(ep$pred_model_sc, na.rm = TRUE)) /
             max(1e-9, sd(ep$pred_model_sc, na.rm = TRUE))

    add_stratum <- function(dim_name, strat_label, mask, test_name = "", p_val = NA_real_, ratio_vs_ref = NA_real_) {
      idx <- which(mask & !is.na(ep$pred_combined))
      if (length(idx) == 0) return(NULL)
      sub_ep <- ep[idx, ]
      obs_v  <- sub_ep$obs_value[sub_ep$has_obs == 1 & !is.na(sub_ep$obs_value)]
      qs     <- as.numeric(quantile(sub_ep$pred_combined, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE))
      strat_rows[[length(strat_rows) + 1L]] <<- data.frame(
        Run_Index          = i,
        Run_Label          = b$label,
        Determinand_Code   = b$det_code,
        Unit               = b$det_unit,
        Dimension          = dim_name,
        Stratum            = strat_label,
        n_reaches          = length(idx),
        Pct_Network        = round(100.0 * length(idx) / n_tot, 2),
        n_obs              = length(obs_v),
        Obs_Median_Raw     = if (length(obs_v) > 0) round(median(obs_v, na.rm = TRUE), 5) else NA_real_,
        Pred_Mean_Raw      = round(mean(sub_ep$pred_combined, na.rm = TRUE), 5),
        Pred_SD_Raw        = round(sd(sub_ep$pred_combined, na.rm = TRUE), 5),
        Pred_GeomMean_Raw  = round(calc_geom_mean(sub_ep$pred_combined), 5),
        Pred_P05_Raw       = round(qs[1], 5),
        Pred_Q25_Raw       = round(qs[2], 5),
        Pred_Median_Raw    = round(qs[3], 5),
        Pred_Q75_Raw       = round(qs[4], 5),
        Pred_P95_Raw       = round(qs[5], 5),
        Pred_IQR_Raw       = round(qs[4] - qs[2], 5),
        Pred_Mean_Z        = round(mean(z_all[idx], na.rm = TRUE), 4),
        SE_Mean            = round(mean(sub_ep$se_model_sc, na.rm = TRUE), 4),
        SE_Median          = round(median(sub_ep$se_model_sc, na.rm = TRUE), 4),
        Ratio_vs_Ref       = round(ratio_vs_ref, 4),
        Within_Run_Test    = test_name,
        Within_Run_P       = p_val,
        Within_Run_Sig     = sig_stars(p_val),
        stringsAsFactors   = FALSE
      )
    }

    # 1. Topography: Lowland (<150m) vs Upland (>=150m)
    if ("Elevation" %in% names(ep)) {
      is_up  <- !is.na(ep$Elevation) & ep$Elevation >= 150.0
      is_low <- !is.na(ep$Elevation) & ep$Elevation <  150.0
      p_topo <- if (sum(is_up) >= 5 && sum(is_low) >= 5) {
        suppressWarnings(wilcox.test(ep$pred_model_sc[is_low], ep$pred_model_sc[is_up]))$p.value
      } else NA_real_
      med_up <- median(ep$pred_combined[is_up], na.rm = TRUE)
      med_lo <- median(ep$pred_combined[is_low], na.rm = TRUE)
      add_stratum("1. Topography", "Upland (>= 150 m)", is_up, "Mann-Whitney U (Lowland vs Upland)", p_topo, 1.0)
      add_stratum("1. Topography", "Lowland (< 150 m)", is_low, "Mann-Whitney U (Lowland vs Upland)", p_topo,
                  if (!is.na(med_up) && abs(med_up) > 1e-9) med_lo / med_up else NA_real_)
    }

    # 2. Hydrogeology: Non-Aquifer (0) vs Principal Aquifer (1)
    if ("Is_Aquifer" %in% names(ep)) {
      is_aq  <- !is.na(ep$Is_Aquifer) & ep$Is_Aquifer == 1
      is_naq <- !is.na(ep$Is_Aquifer) & ep$Is_Aquifer == 0
      p_aq <- if (sum(is_aq) >= 5 && sum(is_naq) >= 5) {
        suppressWarnings(wilcox.test(ep$pred_model_sc[is_aq], ep$pred_model_sc[is_naq]))$p.value
      } else NA_real_
      med_naq <- median(ep$pred_combined[is_naq], na.rm = TRUE)
      med_aq  <- median(ep$pred_combined[is_aq], na.rm = TRUE)
      add_stratum("2. Hydrogeology", "Non-Aquifer (Bedrock = 0)", is_naq, "Mann-Whitney U (Aquifer vs Non-Aq)", p_aq, 1.0)
      add_stratum("2. Hydrogeology", "Principal Aquifer (Bedrock = 1)", is_aq, "Mann-Whitney U (Aquifer vs Non-Aq)", p_aq,
                  if (!is.na(med_naq) && abs(med_naq) > 1e-9) med_aq / med_naq else NA_real_)
    }

    # 3. Combined Landscape Typology
    if (("Elevation" %in% names(ep)) && ("Is_Aquifer" %in% names(ep))) {
      m_up   <- !is.na(ep$Elevation) & ep$Elevation >= 150.0
      m_lnaq <- !is.na(ep$Elevation) & ep$Elevation < 150.0 & ep$Is_Aquifer == 0
      m_laq  <- !is.na(ep$Elevation) & ep$Elevation < 150.0 & ep$Is_Aquifer == 1
      grp_l  <- ifelse(m_up, "Upland", ifelse(m_laq, "Low_Aq", "Low_NonAq"))
      p_kw_l <- if (length(unique(grp_l)) >= 2) {
        suppressWarnings(kruskal.test(ep$pred_model_sc ~ as.factor(grp_l)))$p.value
      } else NA_real_
      med_ref <- median(ep$pred_combined[m_up], na.rm = TRUE)
      add_stratum("3. Landscape Regime", "Upland (>= 150 m)", m_up, "Kruskal-Wallis (3 Regimes)", p_kw_l, 1.0)
      add_stratum("3. Landscape Regime", "Lowland Non-Aquifer (< 150 m)", m_lnaq, "Kruskal-Wallis (3 Regimes)", p_kw_l,
                  if (!is.na(med_ref) && abs(med_ref) > 1e-9) median(ep$pred_combined[m_lnaq], na.rm = TRUE) / med_ref else NA_real_)
      add_stratum("3. Landscape Regime", "Lowland Aquifer (< 150 m, Aq = 1)", m_laq, "Kruskal-Wallis (3 Regimes)", p_kw_l,
                  if (!is.na(med_ref) && abs(med_ref) > 1e-9) median(ep$pred_combined[m_laq], na.rm = TRUE) / med_ref else NA_real_)
    }

    # 4. Strahler Stream Order Gradient
    if ("Strahler" %in% names(ep)) {
      m_hw  <- !is.na(ep$Strahler) & ep$Strahler <= 2
      m_mid <- !is.na(ep$Strahler) & ep$Strahler %in% c(3, 4)
      m_ms  <- !is.na(ep$Strahler) & ep$Strahler >= 5
      ord_g <- ifelse(m_hw, "1_HW", ifelse(m_mid, "2_Mid", "3_Main"))
      p_ord <- if (length(unique(ord_g)) >= 2) {
        suppressWarnings(kruskal.test(ep$pred_model_sc ~ as.factor(ord_g)))$p.value
      } else NA_real_
      rho_ord <- suppressWarnings(cor(ep$Strahler, ep$pred_model_sc, use = "complete.obs", method = "spearman"))
      med_hw  <- median(ep$pred_combined[m_hw], na.rm = TRUE)
      t_lbl   <- sprintf("Kruskal-Wallis (Spearman rho = %+.3f)", rho_ord)
      add_stratum("4. Stream Order", "Headwaters (Order 1-2)", m_hw, t_lbl, p_ord, 1.0)
      add_stratum("4. Stream Order", "Medium Tributaries (Order 3-4)", m_mid, t_lbl, p_ord,
                  if (!is.na(med_hw) && abs(med_hw) > 1e-9) median(ep$pred_combined[m_mid], na.rm = TRUE) / med_hw else NA_real_)
      add_stratum("4. Stream Order", "Mainstem Rivers (Order >= 5)", m_ms, t_lbl, p_ord,
                  if (!is.na(med_hw) && abs(med_hw) > 1e-9) median(ep$pred_combined[m_ms], na.rm = TRUE) / med_hw else NA_real_)
    }

    # 5. Major River Basins
    med_eng <- median(ep$pred_combined, na.rm = TRUE)
    add_stratum("5. River Basin", "All England Network", rep(TRUE, nrow(ep)), "Reference Network", NA_real_, 1.0)
    if (!is.na(b$hum_cid) && any(ep$CatchID == b$hum_cid)) {
      m_b <- ep$CatchID == b$hum_cid
      p_b <- suppressWarnings(wilcox.test(ep$pred_model_sc[m_b], ep$pred_model_sc[!m_b]))$p.value
      add_stratum("5. River Basin", "Humber Basin", m_b, "Mann-Whitney U vs Rest of Network", p_b,
                  if (abs(med_eng) > 1e-9) median(ep$pred_combined[m_b], na.rm = TRUE) / med_eng else NA_real_)
    }
    if (!is.na(b$tha_cid) && any(ep$CatchID == b$tha_cid)) {
      m_b <- ep$CatchID == b$tha_cid
      p_b <- suppressWarnings(wilcox.test(ep$pred_model_sc[m_b], ep$pred_model_sc[!m_b]))$p.value
      add_stratum("5. River Basin", "Thames Basin", m_b, "Mann-Whitney U vs Rest of Network", p_b,
                  if (abs(med_eng) > 1e-9) median(ep$pred_combined[m_b], na.rm = TRUE) / med_eng else NA_real_)
    }
    if (!is.na(b$sev_cid) && any(ep$CatchID == b$sev_cid)) {
      m_b <- ep$CatchID == b$sev_cid
      p_b <- suppressWarnings(wilcox.test(ep$pred_model_sc[m_b], ep$pred_model_sc[!m_b]))$p.value
      add_stratum("5. River Basin", "Severn Basin", m_b, "Mann-Whitney U vs Rest of Network", p_b,
                  if (abs(med_eng) > 1e-9) median(ep$pred_combined[m_b], na.rm = TRUE) / med_eng else NA_real_)
    }
  }

  do.call(rbind, strat_rows)
}

# ------------------------------------------------------------------------------
# Helper: Reach-by-Reach Relationship, Co-Occurrence & Stoichiometric Analysis
#         (Supports both Cross-Determinand and Same-Determinand Model Pairs)
# ------------------------------------------------------------------------------
analyze_reach_relationships <- function(bundles, pair12_pred_match) {
  if (is.null(pair12_pred_match) || nrow(pair12_pred_match) < 20) return(NULL)

  b1 <- bundles[[1]]
  b2 <- bundles[[2]]
  pm <- pair12_pred_match
  same_det <- (b1$det_code == b2$det_code)

  # Standardized Z-scores across matched reaches
  s1 <- max(1e-9, sd(pm$pred1_sc, na.rm = TRUE))
  s2 <- max(1e-9, sd(pm$pred2_sc, na.rm = TRUE))
  pm$z1      <- (pm$pred1_sc - mean(pm$pred1_sc, na.rm = TRUE)) / s1
  pm$z2      <- (pm$pred2_sc - mean(pm$pred2_sc, na.rm = TRUE)) / s2
  pm$z_diff  <- pm$z2 - pm$z1
  pm$z_joint <- 0.5 * (pm$z1 + pm$z2)

  # 4-Quadrant Co-occurrence Classification relative to network medians
  med1_sc <- median(pm$pred1_sc, na.rm = TRUE)
  med2_sc <- median(pm$pred2_sc, na.rm = TRUE)
  med1_raw <- median(pm$pred1_raw, na.rm = TRUE)
  med2_raw <- median(pm$pred2_raw, na.rm = TRUE)

  hi1 <- pm$pred1_sc > med1_sc
  hi2 <- pm$pred2_sc > med2_sc
  pm$quadrant <- ifelse(hi1 & hi2, "Q1_High_High",
                 ifelse(!hi1 & !hi2, "Q2_Low_Low",
                 ifelse(hi1 & !hi2, "Q3_High1_Low2", "Q4_Low1_High2")))

  n_q1 <- sum(pm$quadrant == "Q1_High_High",  na.rm = TRUE)
  n_q2 <- sum(pm$quadrant == "Q2_Low_Low",    na.rm = TRUE)
  n_q3 <- sum(pm$quadrant == "Q3_High1_Low2", na.rm = TRUE)
  n_q4 <- sum(pm$quadrant == "Q4_Low1_High2", na.rm = TRUE)
  n_tot <- n_q1 + n_q2 + n_q3 + n_q4

  tab_2x2 <- matrix(c(n_q1, n_q3, n_q4, n_q2), nrow = 2, byrow = TRUE)
  chisq_res  <- suppressWarnings(chisq.test(tab_2x2))
  odds_ratio <- (as.numeric(n_q1) * as.numeric(n_q2)) / max(1.0, as.numeric(n_q3) * as.numeric(n_q4))
  pct_conc   <- 100.0 * (n_q1 + n_q2) / max(1, n_tot)

  # Top-Quartile (>P75) & Top-Decile (>P90) Hotspot Overlap
  q75_1 <- quantile(pm$pred1_sc, 0.75, na.rm = TRUE)
  q75_2 <- quantile(pm$pred2_sc, 0.75, na.rm = TRUE)
  h75_1 <- pm$pred1_sc > q75_1
  h75_2 <- pm$pred2_sc > q75_2
  n_h75_both <- sum(h75_1 & h75_2, na.rm = TRUE)
  jaccard_75 <- n_h75_both / max(1, sum(h75_1 | h75_2, na.rm = TRUE))
  cond_75    <- 100.0 * n_h75_both / max(1, sum(h75_1, na.rm = TRUE))

  q90_1 <- quantile(pm$pred1_sc, 0.90, na.rm = TRUE)
  q90_2 <- quantile(pm$pred2_sc, 0.90, na.rm = TRUE)
  h90_1 <- pm$pred1_sc > q90_1
  h90_2 <- pm$pred2_sc > q90_2
  n_h90_both <- sum(h90_1 & h90_2, na.rm = TRUE)
  jaccard_90 <- n_h90_both / max(1, sum(h90_1 | h90_2, na.rm = TRUE))
  cond_90    <- 100.0 * n_h90_both / max(1, sum(h90_1, na.rm = TRUE))

  # Bivariate correlations & OLS regression (All Matched Reaches)
  ct_p_sc  <- suppressWarnings(cor.test(pm$pred1_sc,  pm$pred2_sc,  method = "pearson"))
  ct_s_sc  <- suppressWarnings(cor.test(pm$pred1_sc,  pm$pred2_sc,  method = "spearman"))
  ct_p_raw <- suppressWarnings(cor.test(pm$pred1_raw, pm$pred2_raw, method = "pearson"))
  ct_s_raw <- suppressWarnings(cor.test(pm$pred1_raw, pm$pred2_raw, method = "spearman"))
  r_tu_sc  <- suppressWarnings(cor(pm$pred1_tu_sc, pm$pred2_tu_sc, use = "complete.obs", method = "pearson"))
  r_td_sc  <- suppressWarnings(cor(pm$pred1_td_sc, pm$pred2_td_sc, use = "complete.obs", method = "pearson"))

  ols_fit <- lm(pred2_sc ~ pred1_sc, data = pm)
  ols_sum <- summary(ols_fit)
  ols_coef <- ols_sum$coefficients
  ols_alpha <- ols_coef[1, 1]
  ols_beta  <- if (nrow(ols_coef) >= 2) ols_coef[2, 1] else NA_real_
  ols_se    <- if (nrow(ols_coef) >= 2) ols_coef[2, 2] else NA_real_
  ols_p     <- if (nrow(ols_coef) >= 2) ols_coef[2, 4] else NA_real_
  ols_r2    <- ols_sum$r.squared

  # Co-located Observed Station Relationship (where both runs have direct observations)
  both_obs_idx <- which(pm$has_obs_1 == 1 & pm$has_obs_2 == 1 & !is.na(pm$obs1_sc) & !is.na(pm$obs2_sc))
  n_both_obs   <- length(both_obs_idx)
  r_obs_sc     <- NA_real_; p_obs_sc <- NA_real_; rho_obs_raw <- NA_real_; ols_obs_beta <- NA_real_
  if (n_both_obs >= 5) {
    ct_obs_p   <- suppressWarnings(cor.test(pm$obs1_sc[both_obs_idx], pm$obs2_sc[both_obs_idx], method = "pearson"))
    ct_obs_s   <- suppressWarnings(cor.test(pm$obs1_raw[both_obs_idx], pm$obs2_raw[both_obs_idx], method = "spearman"))
    r_obs_sc   <- as.numeric(ct_obs_p$estimate)
    p_obs_sc   <- ct_obs_p$p.value
    rho_obs_raw <- as.numeric(ct_obs_s$estimate)
    fit_obs    <- lm(obs2_sc ~ obs1_sc, data = pm[both_obs_idx, ])
    if (length(coef(fit_obs)) >= 2) ols_obs_beta <- as.numeric(coef(fit_obs)[2])
  }

  # Nutrient Stoichiometry Check (Nitrogen vs Phosphorus pair)
  n_codes <- c("0117", "0116", "0111", "0118")
  p_codes <- c("0180")
  is_np_pair <- (b1$det_code %in% n_codes && b2$det_code %in% p_codes) ||
                (b2$det_code %in% n_codes && b1$det_code %in% p_codes)

  pm$ratio_1_to_2 <- ifelse(!is.na(pm$pred2_raw) & abs(pm$pred2_raw) > 1e-6,
                            pm$pred1_raw / pm$pred2_raw, NA_real_)

  np_info <- NULL
  if (is_np_pair) {
    n_raw <- if (b1$det_code %in% n_codes) pm$pred1_raw else pm$pred2_raw
    p_raw <- if (b1$det_code %in% n_codes) pm$pred2_raw else pm$pred1_raw
    n_lbl <- if (b1$det_code %in% n_codes) b1$label else b2$label
    p_lbl <- if (b1$det_code %in% n_codes) b2$label else b1$label

    pm$np_mass_ratio  <- n_raw / pmax(1e-6, p_raw)
    # Atomic conversion: (N mg/L / 14.0067 g/mol) / (P mg/L / 30.97376 g/mol) = 2.21135 * (N / P)
    pm$np_molar_ratio <- 2.21135 * pm$np_mass_ratio
    pm$np_regime      <- ifelse(pm$np_molar_ratio < 10.0, "N-Limited (Molar N:P < 10)",
                         ifelse(pm$np_molar_ratio <= 20.0, "Co-Limited (Molar N:P 10-20)",
                                "P-Limited (Molar N:P > 20)"))

    obs_np_mass  <- NA_real_; obs_np_molar <- NA_real_
    if (n_both_obs >= 5) {
      n_obs_v <- if (b1$det_code %in% n_codes) pm$obs1_raw[both_obs_idx] else pm$obs2_raw[both_obs_idx]
      p_obs_v <- if (b1$det_code %in% n_codes) pm$obs2_raw[both_obs_idx] else pm$obs1_raw[both_obs_idx]
      obs_np_mass  <- median(n_obs_v / pmax(1e-6, p_obs_v), na.rm = TRUE)
      obs_np_molar <- 2.21135 * obs_np_mass
    }

    np_info <- list(
      n_label        = n_lbl,
      p_label        = p_lbl,
      mass_summary   = summarize_numeric_vector(pm$np_mass_ratio),
      molar_summary  = summarize_numeric_vector(pm$np_molar_ratio),
      obs_med_mass   = obs_np_mass,
      obs_med_molar  = obs_np_molar,
      n_n_limited    = sum(pm$np_molar_ratio < 10.0, na.rm = TRUE),
      pct_n_limited  = 100.0 * mean(pm$np_molar_ratio < 10.0, na.rm = TRUE),
      n_co_limited   = sum(pm$np_molar_ratio >= 10.0 & pm$np_molar_ratio <= 20.0, na.rm = TRUE),
      pct_co_limited = 100.0 * mean(pm$np_molar_ratio >= 10.0 & pm$np_molar_ratio <= 20.0, na.rm = TRUE),
      n_p_limited    = sum(pm$np_molar_ratio > 20.0, na.rm = TRUE),
      pct_p_limited  = 100.0 * mean(pm$np_molar_ratio > 20.0, na.rm = TRUE)
    )
  }

  # Build Stratified Reach-Relationship Summary Table across Basins, Regimes & Stream Orders
  strat_rel_rows <- list()
  add_rel_stratum <- function(s_group, s_name, mask) {
    idx <- which(mask & !is.na(pm$pred1_sc) & !is.na(pm$pred2_sc))
    if (length(idx) < 5) return(NULL)
    sub_pm <- pm[idx, ]
    ct_p <- suppressWarnings(cor.test(sub_pm$pred1_sc, sub_pm$pred2_sc, method = "pearson"))
    ct_s <- suppressWarnings(cor.test(sub_pm$pred1_sc, sub_pm$pred2_sc, method = "spearman"))
    fit_s <- lm(pred2_sc ~ pred1_sc, data = sub_pm)
    beta_s <- if (length(coef(fit_s)) >= 2) as.numeric(coef(fit_s)[2]) else NA_real_
    r2_s   <- summary(fit_s)$r.squared

    h1_s <- sub_pm$pred1_sc > q75_1
    h2_s <- sub_pm$pred2_sc > q75_2
    jac_s <- sum(h1_s & h2_s, na.rm = TRUE) / max(1, sum(h1_s | h2_s, na.rm = TRUE))

    row_df <- data.frame(
      Stratum_Group       = s_group,
      Stratum             = s_name,
      n_reaches           = length(idx),
      Run1_Median_Raw     = round(median(sub_pm$pred1_raw, na.rm = TRUE), 4),
      Run1_IQR_Raw        = round(IQR(sub_pm$pred1_raw, na.rm = TRUE), 4),
      Run2_Median_Raw     = round(median(sub_pm$pred2_raw, na.rm = TRUE), 4),
      Run2_IQR_Raw        = round(IQR(sub_pm$pred2_raw, na.rm = TRUE), 4),
      Pearson_r           = round(as.numeric(ct_p$estimate), 4),
      Pearson_P           = ct_p$p.value,
      Spearman_rho        = round(as.numeric(ct_s$estimate), 4),
      Spearman_P          = ct_s$p.value,
      OLS_Slope_Beta      = round(beta_s, 4),
      OLS_R2              = round(r2_s, 4),
      Pct_Q1_HighHigh     = round(100.0 * mean(sub_pm$quadrant == "Q1_High_High",  na.rm = TRUE), 2),
      Pct_Q2_LowLow       = round(100.0 * mean(sub_pm$quadrant == "Q2_Low_Low",    na.rm = TRUE), 2),
      Pct_Q3_High1_Low2   = round(100.0 * mean(sub_pm$quadrant == "Q3_High1_Low2", na.rm = TRUE), 2),
      Pct_Q4_Low1_High2   = round(100.0 * mean(sub_pm$quadrant == "Q4_Low1_High2", na.rm = TRUE), 2),
      Pct_Concordant      = round(100.0 * mean(sub_pm$quadrant %in% c("Q1_High_High", "Q2_Low_Low"), na.rm = TRUE), 2),
      Hotspot75_Jaccard   = round(jac_s, 4),
      Mean_Z_Diff_2_vs_1  = round(mean(sub_pm$z_diff, na.rm = TRUE), 4),
      Median_Ratio_1_to_2 = round(median(sub_pm$ratio_1_to_2, na.rm = TRUE), 4),
      Median_NP_Mass      = if (is_np_pair) round(median(sub_pm$np_mass_ratio, na.rm = TRUE), 3) else NA_real_,
      Median_NP_Molar     = if (is_np_pair) round(median(sub_pm$np_molar_ratio, na.rm = TRUE), 3) else NA_real_,
      Pct_P_Limited_gt20  = if (is_np_pair) round(100.0 * mean(sub_pm$np_molar_ratio > 20.0, na.rm = TRUE), 2) else NA_real_,
      Pct_Co_Limited_10_20= if (is_np_pair) round(100.0 * mean(sub_pm$np_molar_ratio >= 10.0 & sub_pm$np_molar_ratio <= 20.0, na.rm = TRUE), 2) else NA_real_,
      Pct_N_Limited_lt10  = if (is_np_pair) round(100.0 * mean(sub_pm$np_molar_ratio < 10.0, na.rm = TRUE), 2) else NA_real_,
      stringsAsFactors    = FALSE
    )
    strat_rel_rows[[length(strat_rel_rows) + 1L]] <<- row_df
  }

  add_rel_stratum("Network", "All Matched Reaches", rep(TRUE, nrow(pm)))
  if (n_both_obs >= 5) {
    m_obs <- rep(FALSE, nrow(pm)); m_obs[both_obs_idx] <- TRUE
    add_rel_stratum("Network", "Co-located Monitored Reaches", m_obs)
  }
  if (!is.na(b2$hum_cid) && any(pm$CatchID_2 == b2$hum_cid)) {
    add_rel_stratum("Basin", "Humber Basin", pm$CatchID_2 == b2$hum_cid)
  }
  if (!is.na(b2$tha_cid) && any(pm$CatchID_2 == b2$tha_cid)) {
    add_rel_stratum("Basin", "Thames Basin", pm$CatchID_2 == b2$tha_cid)
  }
  if (!is.na(b2$sev_cid) && any(pm$CatchID_2 == b2$sev_cid)) {
    add_rel_stratum("Basin", "Severn Basin", pm$CatchID_2 == b2$sev_cid)
  }
  if (("Elevation" %in% names(pm)) && ("Is_Aquifer" %in% names(pm))) {
    add_rel_stratum("Landscape", "Upland (>= 150 m)", !is.na(pm$Elevation) & pm$Elevation >= 150.0)
    add_rel_stratum("Landscape", "Lowland Non-Aquifer (< 150 m)", !is.na(pm$Elevation) & pm$Elevation < 150.0 & pm$Is_Aquifer == 0)
    add_rel_stratum("Landscape", "Lowland Aquifer (< 150 m, Aq = 1)", !is.na(pm$Elevation) & pm$Elevation < 150.0 & pm$Is_Aquifer == 1)
  }
  if ("Strahler" %in% names(pm)) {
    add_rel_stratum("Stream Order", "Headwaters (Order 1-2)", !is.na(pm$Strahler) & pm$Strahler <= 2)
    add_rel_stratum("Stream Order", "Medium Tributaries (Order 3-4)", !is.na(pm$Strahler) & pm$Strahler %in% c(3, 4))
    add_rel_stratum("Stream Order", "Mainstem Rivers (Order >= 5)", !is.na(pm$Strahler) & pm$Strahler >= 5)
  }

  list(
    pm             = pm,
    same_det       = same_det,
    is_np_pair     = is_np_pair,
    np_info        = np_info,
    n_matched      = nrow(pm),
    med1_sc        = med1_sc,
    med2_sc        = med2_sc,
    med1_raw       = med1_raw,
    med2_raw       = med2_raw,
    n_q1           = n_q1,
    n_q2           = n_q2,
    n_q3           = n_q3,
    n_q4           = n_q4,
    pct_q1         = 100.0 * n_q1 / max(1, n_tot),
    pct_q2         = 100.0 * n_q2 / max(1, n_tot),
    pct_q3         = 100.0 * n_q3 / max(1, n_tot),
    pct_q4         = 100.0 * n_q4 / max(1, n_tot),
    pct_concordant = pct_conc,
    chisq_stat     = as.numeric(chisq_res$statistic),
    chisq_p        = chisq_res$p.value,
    odds_ratio     = odds_ratio,
    n_h75_both     = n_h75_both,
    jaccard_75     = jaccard_75,
    cond_75        = cond_75,
    n_h90_both     = n_h90_both,
    jaccard_90     = jaccard_90,
    cond_90        = cond_90,
    r_sc           = as.numeric(ct_p_sc$estimate),
    p_sc           = ct_p_sc$p.value,
    rho_sc         = as.numeric(ct_s_sc$estimate),
    p_rho_sc       = ct_s_sc$p.value,
    r_raw          = as.numeric(ct_p_raw$estimate),
    p_raw          = ct_p_raw$p.value,
    rho_raw        = as.numeric(ct_s_raw$estimate),
    p_rho_raw      = ct_s_raw$p.value,
    r_tu_sc        = r_tu_sc,
    r_td_sc        = r_td_sc,
    ols_alpha      = ols_alpha,
    ols_beta       = ols_beta,
    ols_se         = ols_se,
    ols_p          = ols_p,
    ols_r2         = ols_r2,
    n_both_obs     = n_both_obs,
    r_obs_sc       = r_obs_sc,
    p_obs_sc       = p_obs_sc,
    rho_obs_raw    = rho_obs_raw,
    ols_obs_beta   = ols_obs_beta,
    stratum_df     = do.call(rbind, strat_rel_rows)
  )
}

# ------------------------------------------------------------------------------
# Core Comparison & Single-Run Descriptive Engine:
#   compare_s3n_runs() / describe_s3n_run()
# ------------------------------------------------------------------------------
describe_s3n_run <- function(run_dir, label = NULL, out_dir = NULL, quiet = FALSE) {
  compare_s3n_runs(run_dirs = run_dir, labels = label, out_dir = out_dir, quiet = quiet)
}

compare_s3n_runs <- function(run_dirs, labels = NULL, out_dir = NULL, quiet = FALSE) {
  if (length(run_dirs) < 1) {
    stop("Please provide at least 1 model run directory.")
  }

  bundles <- lapply(seq_along(run_dirs), function(i) {
    lbl_i <- if (!is.null(labels) && length(labels) >= i) labels[i] else NULL
    load_model_run_bundle(run_dirs[i], custom_label = lbl_i)
  })
  K <- length(bundles)
  same_det_all <- (length(unique(sapply(bundles, function(b) b$det_code))) == 1L)

  # When comparing different determinands, include determinand name & code in default labels
  if (!same_det_all) {
    for (i in seq_len(K)) {
      if (!isTRUE(bundles[[i]]$user_supplied_label)) {
        bundles[[i]]$label <- sprintf("%s [%s] (%s)",
                                      shorten_det_name(bundles[[i]]$det_name),
                                      bundles[[i]]$det_code,
                                      bundles[[i]]$net_label)
      }
    }
  }

  run_labels <- sapply(bundles, function(b) b$label)
  # Ensure unique labels
  if (any(duplicated(run_labels))) {
    run_labels <- sprintf("%s [%d]", run_labels, seq_len(K))
    for (i in seq_len(K)) bundles[[i]]$label <- run_labels[i]
  }

  if (is.null(out_dir) || !nzchar(out_dir)) {
    if (K == 1L) {
      slug <- sprintf("single_%s", bundles[[1]]$dir_name)
    } else {
      slug <- paste(sapply(bundles, function(b) b$dir_name), collapse = "_vs_")
      if (nchar(slug) > 85) {
        slug <- sprintf("compare_%d_runs_%s_to_%s", K, bundles[[1]]$dir_name, bundles[[K]]$dir_name)
      }
    }
    out_dir <- file.path(BASE_DIR, "output", "comparisons", slug)
  } else if (!grepl("^/", out_dir)) {
    out_dir <- file.path(BASE_DIR, out_dir)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (!quiet) {
    cat("==============================================================================\n")
    if (K == 1L) {
      cat(sprintf("Single-Run Descriptive & Diagnostic Analysis (%s [%s]):\n",
                  bundles[[1]]$det_name, bundles[[1]]$det_code))
    } else {
      cat(sprintf("Comparing %d S3N Model Runs (%s):\n",
                  K, if (same_det_all) sprintf("Same Determinand: %s", bundles[[1]]$det_name) else "Cross-Determinand Comparison"))
    }
    for (i in seq_len(K)) {
      cat(sprintf("  Run %d: %-42s (%s | %s reaches, %d obs)\n",
                  i, bundles[[i]]$label, bundles[[i]]$dir_name,
                  format(nrow(bundles[[i]]$eng_pred), big.mark = ","),
                  sum(bundles[[i]]$eng_pred$has_obs == 1)))
    }
    cat(sprintf("Saving report, descriptive tables & figures to:\n  %s\n", out_dir))
    cat("==============================================================================\n")
  }

  # ----------------------------------------------------------------------------
  # 1. Assemble Combined Cross-Validation Metrics Table
  # ----------------------------------------------------------------------------
  cv_list <- lapply(seq_len(K), function(i) {
    b <- bundles[[i]]
    df <- b$cv_tab
    # Compute Spatial Gain over OLS within each domain
    ols_rows <- subset(df, Canon_Model == "OLS (Non-Spatial)")
    ols_r2   <- setNames(ols_rows$R2, ols_rows$Canon_Domain)
    ols_rmse <- setNames(ols_rows$RMSE, ols_rows$Canon_Domain)

    df$Run_Index        <- i
    df$Run_Label        <- b$label
    df$Determinand_Code <- b$det_code
    df$Determinand_Name <- b$det_name
    df$Unit             <- b$det_unit
    df$Dir_Name         <- b$dir_name
    df$Network_Label    <- b$net_label
    df$Covariates       <- b$vars_label
    df$Delta_R2_vs_OLS  <- df$R2 - ols_r2[df$Canon_Domain]
    df$Pct_RMSE_Red_OLS <- 100.0 * (ols_rmse[df$Canon_Domain] - df$RMSE) / ols_rmse[df$Canon_Domain]
    df[, c("Run_Index", "Run_Label", "Determinand_Code", "Determinand_Name", "Unit",
           "Dir_Name", "Network_Label", "Covariates",
           "Canon_Domain", "Domain", "n_obs", "n_reaches", "Canon_Model", "Model",
           "R2", "RMSE", "MAE", "Delta_R2_vs_OLS", "Pct_RMSE_Red_OLS")]
  })
  cv_comp_df <- do.call(rbind, cv_list)
  write.csv(cv_comp_df, file.path(out_dir, "cv_metrics_comparison.csv"), row.names = FALSE)

  # ----------------------------------------------------------------------------
  # 2. Assemble Reach Prediction & Uncertainty Distribution Summary Table
  #    + Detailed Descriptive Statistics & Stratified Results Tables
  # ----------------------------------------------------------------------------
  pred_sum_list <- lapply(seq_len(K), function(i) {
    b  <- bundles[[i]]
    ep <- b$eng_pred
    obs_sub <- ep[ep$has_obs == 1, ]
    q_raw <- quantile(ep$pred_combined, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE)
    q_sc  <- quantile(ep$pred_model_sc, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE)
    q_se  <- quantile(ep$se_model_sc,   probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE)
    data.frame(
      Run_Index        = i,
      Run_Label        = b$label,
      Determinand_Code = b$det_code,
      Determinand_Name = b$det_name,
      Unit             = b$det_unit,
      Scale            = if (b$use_log) "log" else "raw",
      Network_Label    = b$net_label,
      Period           = b$period,
      Covariates       = b$vars_label,
      n_reaches        = nrow(ep),
      n_catchments     = length(unique(ep$CatchID)),
      n_obs_reaches    = nrow(obs_sub),
      Obs_Mean_Raw     = mean(obs_sub$obs_value, na.rm = TRUE),
      Obs_SD_Raw       = sd(obs_sub$obs_value, na.rm = TRUE),
      Obs_Median_Raw   = median(obs_sub$obs_value, na.rm = TRUE),
      Obs_IQR_Raw      = IQR(obs_sub$obs_value, na.rm = TRUE),
      Pred_Mean_Raw    = mean(ep$pred_combined, na.rm = TRUE),
      Pred_SD_Raw      = sd(ep$pred_combined, na.rm = TRUE),
      Pred_CV_Pct_Raw  = 100.0 * sd(ep$pred_combined, na.rm = TRUE) / max(1e-9, abs(mean(ep$pred_combined, na.rm = TRUE))),
      Pred_GeomMean_Raw= calc_geom_mean(ep$pred_combined),
      Pred_Skew_Raw    = calc_sample_skewness(ep$pred_combined),
      Pred_P05_Raw     = as.numeric(q_raw[1]),
      Pred_Q25_Raw     = as.numeric(q_raw[2]),
      Pred_Median_Raw  = as.numeric(q_raw[3]),
      Pred_Q75_Raw     = as.numeric(q_raw[4]),
      Pred_P95_Raw     = as.numeric(q_raw[5]),
      Pred_IQR_Raw     = as.numeric(q_raw[4] - q_raw[2]),
      Pred_Mean_ModSc  = mean(ep$pred_model_sc, na.rm = TRUE),
      Pred_SD_ModSc    = sd(ep$pred_model_sc, na.rm = TRUE),
      Pred_Med_ModSc   = as.numeric(q_sc[3]),
      SE_Mean          = mean(ep$se_model_sc, na.rm = TRUE),
      SE_SD            = sd(ep$se_model_sc, na.rm = TRUE),
      SE_Median        = as.numeric(q_se[3]),
      SE_IQR           = as.numeric(q_se[4] - q_se[2]),
      SE_P95           = as.numeric(q_se[5]),
      stringsAsFactors = FALSE
    )
  })
  pred_sum_df <- do.call(rbind, pred_sum_list)
  write.csv(pred_sum_df, file.path(out_dir, "prediction_distribution_summary.csv"), row.names = FALSE)

  desc_bundle <- compute_detailed_descriptives(bundles)
  desc_df     <- desc_bundle$desc_df
  rep_df      <- desc_bundle$rep_df
  write.csv(desc_df, file.path(out_dir, "descriptive_stats_summary.csv"), row.names = FALSE)

  strat_df <- compute_stratified_descriptives(bundles)
  write.csv(strat_df, file.path(out_dir, "stratified_results_summary.csv"), row.names = FALSE)

  # ----------------------------------------------------------------------------
  # 3. Pairwise, Omnibus & Within-Run Statistical Tests
  # ----------------------------------------------------------------------------
  stat_rows <- list()
  add_stat_row <- function(comp_lbl, category, test_name, metric, n1, n2, mean1, mean2, diff_2_minus_1,
                           stat_name, stat_val, p_val, interpretation) {
    stat_rows[[length(stat_rows) + 1L]] <<- data.frame(
      Comparison     = comp_lbl,
      Category       = category,
      Test_Name      = test_name,
      Metric         = metric,
      n_1            = n1,
      n_2            = n2,
      Value_Run1     = round(mean1, 5),
      Value_Run2     = round(mean2, 5),
      Diff_Run2_Run1 = round(diff_2_minus_1, 5),
      Statistic_Name = stat_name,
      Statistic_Val  = round(stat_val, 5),
      P_Value        = p_val,
      Significance   = sig_stars(p_val),
      Interpretation = interpretation,
      stringsAsFactors = FALSE
    )
  }

  # Store matched data for pair (1, 2) for plotting in Figure 1, Figure 2 & Figure 5
  pair12_cv_match   <- NULL
  pair12_pred_match <- NULL

  if (K == 1L) {
    b1 <- bundles[[1]]
    ep <- b1$eng_pred
    lbl1 <- b1$label

    # 1. Cross-Validation Model Comparisons within the Single Run (S3N Combined vs OLS, TU, TD)
    if (!is.null(b1$eng_cv) && nrow(b1$eng_cv) >= 5) {
      e_comb <- b1$eng_cv$obs_y - b1$eng_cv$cv_comb
      if ("cv_ols" %in% names(b1$eng_cv)) {
        e_ols <- b1$eng_cv$obs_y - b1$eng_cv$cv_ols
        tt_ols <- safe_paired_ttest(e_comb^2, e_ols^2)
        wt_ols <- safe_paired_wilcox(abs(e_comb), abs(e_ols))
        add_stat_row("S3N Combined vs. OLS", "Cross-Validation",
                     "Paired t-Test / Diebold-Mariano (Station Squared CV Error)",
                     "Station MSE", length(e_ols), length(e_comb),
                     mean(e_ols^2, na.rm = TRUE), mean(e_comb^2, na.rm = TRUE),
                     mean(e_comb^2, na.rm = TRUE) - mean(e_ols^2, na.rm = TRUE),
                     "t", as.numeric(tt_ols$statistic), tt_ols$p.value,
                     sprintf("S3N Combined reduces station MSE vs non-spatial OLS (%s)", sig_stars(tt_ols$p.value)))
        add_stat_row("S3N Combined vs. OLS", "Cross-Validation",
                     "Paired Wilcoxon Signed-Rank Test (Station |CV Error|)",
                     "Mean |CV Error|", length(e_ols), length(e_comb),
                     mean(abs(e_ols), na.rm = TRUE), mean(abs(e_comb), na.rm = TRUE),
                     mean(abs(e_comb), na.rm = TRUE) - mean(abs(e_ols), na.rm = TRUE),
                     "V", as.numeric(wt_ols$statistic), wt_ols$p.value,
                     sprintf("S3N Combined reduces absolute station CV error vs OLS (%s)", sig_stars(wt_ols$p.value)))
      }
      if ("cv_tu" %in% names(b1$eng_cv)) {
        e_tu <- b1$eng_cv$obs_y - b1$eng_cv$cv_tu
        tt_tu <- safe_paired_ttest(e_comb^2, e_tu^2)
        add_stat_row("S3N Combined vs. Tail-Up", "Cross-Validation",
                     "Paired t-Test (Station Squared CV Error: Combined vs Tail-Up)",
                     "Station MSE", length(e_tu), length(e_comb),
                     mean(e_tu^2, na.rm = TRUE), mean(e_comb^2, na.rm = TRUE),
                     mean(e_comb^2, na.rm = TRUE) - mean(e_tu^2, na.rm = TRUE),
                     "t", as.numeric(tt_tu$statistic), tt_tu$p.value,
                     sprintf("Variance-weighted combination vs pure Tail-Up MSE (%s)", sig_stars(tt_tu$p.value)))
      }
      if ("cv_td" %in% names(b1$eng_cv)) {
        e_td <- b1$eng_cv$obs_y - b1$eng_cv$cv_td
        tt_td <- safe_paired_ttest(e_comb^2, e_td^2)
        add_stat_row("S3N Combined vs. Tail-Down", "Cross-Validation",
                     "Paired t-Test (Station Squared CV Error: Combined vs Tail-Down)",
                     "Station MSE", length(e_td), length(e_comb),
                     mean(e_td^2, na.rm = TRUE), mean(e_comb^2, na.rm = TRUE),
                     mean(e_comb^2, na.rm = TRUE) - mean(e_td^2, na.rm = TRUE),
                     "t", as.numeric(tt_td$statistic), tt_td$p.value,
                     sprintf("Variance-weighted combination vs pure Tail-Down MSE (%s)", sig_stars(tt_td$p.value)))
      }
    }

    # 2. Monitoring Representativeness: Observed Stations vs. Full Network Predictions
    obs_v <- ep$obs_value[ep$has_obs == 1 & !is.na(ep$obs_value)]
    if (length(obs_v) >= 5) {
      mw_rep <- suppressWarnings(wilcox.test(ep$pred_combined, obs_v, paired = FALSE))
      ks_rep <- suppressWarnings(ks.test(ep$pred_combined, obs_v))
      add_stat_row("Network vs. Monitored Sites", "Monitoring Representativeness",
                   "Mann-Whitney U Test (Full Network Predictions vs Observed Stations)",
                   sprintf("Median (%s)", b1$det_unit), length(obs_v), nrow(ep),
                   median(obs_v, na.rm = TRUE), median(ep$pred_combined, na.rm = TRUE),
                   median(ep$pred_combined, na.rm = TRUE) - median(obs_v, na.rm = TRUE),
                   "W", as.numeric(mw_rep$statistic), mw_rep$p.value,
                   sprintf("Pred/Obs median ratio = %.2fx (%s)", median(ep$pred_combined, na.rm = TRUE) / max(1e-9, median(obs_v, na.rm = TRUE)), sig_stars(mw_rep$p.value)))
      add_stat_row("Network vs. Monitored Sites", "Monitoring Representativeness",
                   "Two-Sample Kolmogorov-Smirnov Test (Predictions vs Observed Stations)",
                   "Empirical CDF D-stat", length(obs_v), nrow(ep),
                   mean(obs_v, na.rm = TRUE), mean(ep$pred_combined, na.rm = TRUE),
                   mean(ep$pred_combined, na.rm = TRUE) - mean(obs_v, na.rm = TRUE),
                   "D", as.numeric(ks_rep$statistic), ks_rep$p.value,
                   sprintf("Distributional shift between gauged sites and full network (D = %.3f, %s)", ks_rep$statistic, sig_stars(ks_rep$p.value)))
    }

    # 3. S3N Tail-Up vs. Tail-Down Component Concordance across All Reaches
    if ("pred_tu" %in% names(ep) && "pred_td" %in% names(ep)) {
      tu_v <- if ("pred_tu_sc" %in% names(ep)) ep$pred_tu_sc else ep$pred_tu
      td_v <- if ("pred_td_sc" %in% names(ep)) ep$pred_td_sc else ep$pred_td
      ct_tud <- suppressWarnings(cor.test(tu_v, td_v, method = "pearson"))
      cs_tud <- suppressWarnings(cor.test(tu_v, td_v, method = "spearman"))
      add_stat_row("Tail-Up vs. Tail-Down", "S3N Components",
                   "Pearson Correlation Test (Tail-Up vs Tail-Down Reach Predictions)",
                   "Pearson r (Spearman rho)", nrow(ep), nrow(ep),
                   as.numeric(ct_tud$estimate), as.numeric(cs_tud$estimate),
                   as.numeric(ct_tud$estimate)^2,
                   "t", as.numeric(ct_tud$statistic), ct_tud$p.value,
                   sprintf("r = %.4f (Spearman rho = %.4f, R2 = %.4f) across %s reaches (%s)",
                           ct_tud$estimate, cs_tud$estimate, ct_tud$estimate^2,
                           format(nrow(ep), big.mark = ","), sig_stars(ct_tud$p.value)))
    }

    # 4. Kriging Uncertainty at Monitored vs. Unmonitored Reaches
    se_mon   <- ep$se_model_sc[ep$has_obs == 1]
    se_unmon <- ep$se_model_sc[ep$has_obs == 0]
    if (length(se_mon) >= 5 && length(se_unmon) >= 5) {
      mw_se <- suppressWarnings(wilcox.test(se_unmon, se_mon, paired = FALSE))
      add_stat_row("Unmonitored vs. Monitored SE", "Kriging Uncertainty",
                   "Mann-Whitney U Test (Kriging SE: Unmonitored vs Monitored Reaches)",
                   "Mean Kriging SE", length(se_mon), length(se_unmon),
                   mean(se_mon, na.rm = TRUE), mean(se_unmon, na.rm = TRUE),
                   mean(se_unmon, na.rm = TRUE) - mean(se_mon, na.rm = TRUE),
                   "W", as.numeric(mw_se$statistic), mw_se$p.value,
                   sprintf("Unmonitored reaches have %+.4f higher mean Kriging SE than monitored reaches (%s)",
                           mean(se_unmon, na.rm = TRUE) - mean(se_mon, na.rm = TRUE), sig_stars(mw_se$p.value)))
    }

    # 5. Environmental Gradient Hypothesis Tests (Topography, Hydrogeology, Stream Order, Basins)
    if ("Elevation" %in% names(ep)) {
      v_up  <- ep$pred_combined[!is.na(ep$Elevation) & ep$Elevation >= 150.0]
      v_low <- ep$pred_combined[!is.na(ep$Elevation) & ep$Elevation <  150.0]
      if (length(v_up) >= 5 && length(v_low) >= 5) {
        mw_topo <- suppressWarnings(wilcox.test(v_low, v_up, paired = FALSE))
        add_stat_row("Lowland vs. Upland", "Environmental Gradient",
                     "Mann-Whitney U Test (Lowland < 150 m vs Upland >= 150 m)",
                     sprintf("Median (%s)", b1$det_unit), length(v_up), length(v_low),
                     median(v_up, na.rm = TRUE), median(v_low, na.rm = TRUE),
                     median(v_low, na.rm = TRUE) - median(v_up, na.rm = TRUE),
                     "W", as.numeric(mw_topo$statistic), mw_topo$p.value,
                     sprintf("Lowland/Upland median ratio = %.2fx (%s)",
                             median(v_low, na.rm = TRUE) / max(1e-9, median(v_up, na.rm = TRUE)), sig_stars(mw_topo$p.value)))
      }
    }
    if ("Is_Aquifer" %in% names(ep)) {
      v_noaq <- ep$pred_combined[!is.na(ep$Is_Aquifer) & ep$Is_Aquifer == 0]
      v_aq   <- ep$pred_combined[!is.na(ep$Is_Aquifer) & ep$Is_Aquifer == 1]
      if (length(v_noaq) >= 5 && length(v_aq) >= 5) {
        mw_aq <- suppressWarnings(wilcox.test(v_aq, v_noaq, paired = FALSE))
        add_stat_row("Aquifer vs. Non-Aquifer", "Environmental Gradient",
                     "Mann-Whitney U Test (Chalk/Limestone Aquifer vs Non-Aquifer Bedrock)",
                     sprintf("Median (%s)", b1$det_unit), length(v_noaq), length(v_aq),
                     median(v_noaq, na.rm = TRUE), median(v_aq, na.rm = TRUE),
                     median(v_aq, na.rm = TRUE) - median(v_noaq, na.rm = TRUE),
                     "W", as.numeric(mw_aq$statistic), mw_aq$p.value,
                     sprintf("Aquifer/Non-Aquifer median ratio = %.2fx (%s)",
                             median(v_aq, na.rm = TRUE) / max(1e-9, median(v_noaq, na.rm = TRUE)), sig_stars(mw_aq$p.value)))
      }
    }
    if ("Strahler" %in% names(ep)) {
      ord_grp <- ifelse(ep$Strahler <= 2, "Headwaters", ifelse(ep$Strahler <= 4, "Mid-Order", "Mainstem"))
      v_hw <- ep$pred_combined[ord_grp == "Headwaters"]
      v_ms <- ep$pred_combined[ord_grp == "Mainstem"]
      if (length(unique(ord_grp[!is.na(ord_grp)])) >= 2) {
        kw_ord <- suppressWarnings(kruskal.test(ep$pred_combined ~ as.factor(ord_grp)))
        add_stat_row("Strahler Stream Orders", "Environmental Gradient",
                     "Kruskal-Wallis Test Across Stream Order Classes (1-2, 3-4, >=5)",
                     sprintf("Mainstem vs Headw. Med (%s)", b1$det_unit), length(v_hw), length(v_ms),
                     median(v_hw, na.rm = TRUE), median(v_ms, na.rm = TRUE),
                     median(v_ms, na.rm = TRUE) - median(v_hw, na.rm = TRUE),
                     "ChiSq", as.numeric(kw_ord$statistic), kw_ord$p.value,
                     sprintf("Mainstem/Headwater median ratio = %.2fx (%s)",
                             median(v_ms, na.rm = TRUE) / max(1e-9, median(v_hw, na.rm = TRUE)), sig_stars(kw_ord$p.value)))
      }
    }
  }

  # Pairwise comparisons against Reference Run 1 (and all pairs if K <= 4) when K >= 2
  pairs_to_test <- if (K < 2) list() else if (K <= 4) combn(seq_len(K), 2, simplify = FALSE) else lapply(2:K, function(j) c(1L, j))

  for (pr in pairs_to_test) {
    i <- pr[1]; j <- pr[2]
    b1 <- bundles[[i]]; b2 <- bundles[[j]]
    comp_lbl <- sprintf("%s vs. %s", b1$label, b2$label)
    same_net <- (b1$net_abbr == b2$net_abbr) && (nrow(b1$eng_pred) == nrow(b2$eng_pred))
    same_det <- (b1$det_code == b2$det_code)

    # A. Fisher's z-test on All England S3N Combined CV R2
    r1_row <- subset(b1$cv_tab, Canon_Domain == "All England RSN Network" & Canon_Model == "S3N Combined")
    r2_row <- subset(b2$cv_tab, Canon_Domain == "All England RSN Network" & Canon_Model == "S3N Combined")
    if (nrow(r1_row) > 0 && nrow(r2_row) > 0) {
      r2_1 <- r1_row$R2[1]; n1 <- r1_row$n_obs[1]
      r2_2 <- r2_row$R2[1]; n2 <- r2_row$n_obs[1]
      cor1 <- sqrt(max(0, min(0.9999, r2_1)))
      cor2 <- sqrt(max(0, min(0.9999, r2_2)))
      z_f  <- (atanh(cor2) - atanh(cor1)) / sqrt(1 / max(4, n1 - 3) + 1 / max(4, n2 - 3))
      p_f  <- 2 * pnorm(-abs(z_f))
      better_cv <- if (r2_2 > r2_1) b2$label else b1$label
      add_stat_row(comp_lbl, "Cross-Validation", "Fisher's r-to-z Test (All England S3N Combined R2)",
                   "5-Fold CV R2", n1, n2, r2_1, r2_2, r2_2 - r2_1, "z", z_f, p_f,
                   sprintf("%s has higher CV R2 (%.4f vs %.4f, %s)", better_cv, max(r2_1, r2_2), min(r2_1, r2_2), sig_stars(p_f)))
    }

    # B. Out-of-Sample Station CV Error Tests (Paired co-located & Unpaired)
    if (!is.null(b1$eng_cv) && !is.null(b2$eng_cv)) {
      e1_all <- b1$eng_cv$obs_y - b1$eng_cv$cv_comb
      e2_all <- b2$eng_cv$obs_y - b2$eng_cv$cv_comb

      # Unpaired Mann-Whitney U test on |e|
      wt_unp <- suppressWarnings(wilcox.test(abs(e1_all), abs(e2_all), paired = FALSE))
      add_stat_row(comp_lbl, "Cross-Validation", "Mann-Whitney U Test (Unpaired Absolute CV Error)",
                   "Mean |CV Error|", length(e1_all), length(e2_all),
                   mean(abs(e1_all), na.rm = TRUE), mean(abs(e2_all), na.rm = TRUE),
                   mean(abs(e2_all), na.rm = TRUE) - mean(abs(e1_all), na.rm = TRUE),
                   "W", as.numeric(wt_unp$statistic), wt_unp$p.value,
                   sprintf("Unpaired station |error| comparison (%s)", sig_stars(wt_unp$p.value)))

      # Unpaired Welch's t-test on squared CV error e^2
      tt_unp <- suppressWarnings(t.test(e2_all^2, e1_all^2, paired = FALSE))
      add_stat_row(comp_lbl, "Cross-Validation", "Welch's Two-Sample t-Test (Unpaired Squared CV Error)",
                   "Mean Squared CV Error (MSE)", length(e1_all), length(e2_all),
                   mean(e1_all^2, na.rm = TRUE), mean(e2_all^2, na.rm = TRUE),
                   mean(e2_all^2, na.rm = TRUE) - mean(e1_all^2, na.rm = TRUE),
                   "t", as.numeric(tt_unp$statistic), tt_unp$p.value,
                   sprintf("Unpaired station MSE comparison (%s)", sig_stars(tt_unp$p.value)))

      # Match co-located observed stations across Run i and Run j
      m_cv <- match_reaches_between_runs(b1$eng_cv, b2$eng_cv, same_network = same_net,
                                         max_dist_m = 500.0, require_obs_match = same_det)
      if (nrow(m_cv) >= 10) {
        e1_m <- e1_all[m_cv$idx1]
        e2_m <- e2_all[m_cv$idx2]
        if (i == 1L && j == 2L) {
          pair12_cv_match <- data.frame(
            RivID_1 = b1$eng_cv$RivID[m_cv$idx1],
            RivID_2 = b2$eng_cv$RivID[m_cv$idx2],
            dist_m  = m_cv$dist_m,
            obs_1   = b1$eng_cv$obs_y[m_cv$idx1],
            obs_2   = b2$eng_cv$obs_y[m_cv$idx2],
            cv_1    = b1$eng_cv$cv_comb[m_cv$idx1],
            cv_2    = b2$eng_cv$cv_comb[m_cv$idx2],
            abs_e1  = abs(e1_m),
            abs_e2  = abs(e2_m),
            sq_e1   = e1_m^2,
            sq_e2   = e2_m^2
          )
        }
        # Paired Diebold-Mariano / Paired t-test on squared errors
        tt_pair <- safe_paired_ttest(e2_m^2, e1_m^2)
        better_p <- if (mean(e2_m^2, na.rm = TRUE) < mean(e1_m^2, na.rm = TRUE)) b2$label else b1$label
        add_stat_row(comp_lbl, "Cross-Validation", "Paired t-Test / Diebold-Mariano (Co-located Squared CV Error)",
                     "Paired Station MSE", nrow(m_cv), nrow(m_cv),
                     mean(e1_m^2, na.rm = TRUE), mean(e2_m^2, na.rm = TRUE),
                     mean(e2_m^2, na.rm = TRUE) - mean(e1_m^2, na.rm = TRUE),
                     "t", as.numeric(tt_pair$statistic), tt_pair$p.value,
                     sprintf("Lower paired MSE at %d co-located stations in %s (%s)", nrow(m_cv), better_p, sig_stars(tt_pair$p.value)))

        # Paired Wilcoxon Signed-Rank Test on |e|
        wt_pair <- safe_paired_wilcox(abs(e2_m), abs(e1_m))
        better_w <- if (median(abs(e2_m), na.rm = TRUE) < median(abs(e1_m), na.rm = TRUE)) b2$label else b1$label
        add_stat_row(comp_lbl, "Cross-Validation", "Paired Wilcoxon Signed-Rank Test (Co-located |CV Error|)",
                     "Paired Station |Error|", nrow(m_cv), nrow(m_cv),
                     mean(abs(e1_m), na.rm = TRUE), mean(abs(e2_m), na.rm = TRUE),
                     mean(abs(e2_m), na.rm = TRUE) - mean(abs(e1_m), na.rm = TRUE),
                     "V", as.numeric(wt_pair$statistic), wt_pair$p.value,
                     sprintf("Lower paired |error| at %d co-located stations in %s (%s)", nrow(m_cv), better_w, sig_stars(wt_pair$p.value)))
      }
    }

    # C. Reach-Level Prediction & Uncertainty Distribution Tests (All Reaches)
    p1_sc <- b1$eng_pred$pred_model_sc; p2_sc <- b2$eng_pred$pred_model_sc
    se1   <- b1$eng_pred$se_model_sc;   se2   <- b2$eng_pred$se_model_sc

    tt_pred <- suppressWarnings(t.test(p2_sc, p1_sc, paired = FALSE))
    add_stat_row(comp_lbl, "Reach Predictions", "Welch's Two-Sample t-Test (Predicted Model-Scale Values)",
                 "Mean pred_model_sc", length(p1_sc), length(p2_sc),
                 mean(p1_sc, na.rm = TRUE), mean(p2_sc, na.rm = TRUE),
                 mean(p2_sc, na.rm = TRUE) - mean(p1_sc, na.rm = TRUE),
                 "t", as.numeric(tt_pred$statistic), tt_pred$p.value,
                 sprintf("Difference in network-wide mean prediction (%s)", sig_stars(tt_pred$p.value)))

    ks_pred <- suppressWarnings(ks.test(p1_sc, p2_sc))
    add_stat_row(comp_lbl, "Reach Predictions", "Two-Sample Kolmogorov-Smirnov Test (Prediction Distribution)",
                 "Empirical CDF D-stat", length(p1_sc), length(p2_sc),
                 median(p1_sc, na.rm = TRUE), median(p2_sc, na.rm = TRUE),
                 median(p2_sc, na.rm = TRUE) - median(p1_sc, na.rm = TRUE),
                 "D", as.numeric(ks_pred$statistic), ks_pred$p.value,
                 sprintf("Distributional shift in reach predictions (D = %.4f, %s)", ks_pred$statistic, sig_stars(ks_pred$p.value)))

    wt_se <- suppressWarnings(wilcox.test(se2, se1, paired = FALSE))
    lower_se_lbl <- if (mean(se2, na.rm = TRUE) < mean(se1, na.rm = TRUE)) b2$label else b1$label
    add_stat_row(comp_lbl, "Kriging Uncertainty", "Mann-Whitney U Test (Kriging Standard Error se_model_sc)",
                 "Mean Kriging SE", length(se1), length(se2),
                 mean(se1, na.rm = TRUE), mean(se2, na.rm = TRUE),
                 mean(se2, na.rm = TRUE) - mean(se1, na.rm = TRUE),
                 "W", as.numeric(wt_se$statistic), wt_se$p.value,
                 sprintf("Lower network-wide prediction uncertainty in %s (%s)", lower_se_lbl, sig_stars(wt_se$p.value)))

    fk_se <- suppressWarnings(fligner.test(list(se1, se2)))
    add_stat_row(comp_lbl, "Kriging Uncertainty", "Fligner-Killeen Variance Homogeneity Test (Kriging SE)",
                 "SD of Kriging SE", length(se1), length(se2),
                 sd(se1, na.rm = TRUE), sd(se2, na.rm = TRUE),
                 sd(se2, na.rm = TRUE) - sd(se1, na.rm = TRUE),
                 "ChiSq", as.numeric(fk_se$statistic), fk_se$p.value,
                 sprintf("Dispersion of kriging standard errors (%s)", sig_stars(fk_se$p.value)))

    # D. Co-located Reach Spatial Agreement & Relationship Tests (Matched by RivID or <=250m BNG proximity)
    m_pred <- match_reaches_between_runs(b1$eng_pred, b2$eng_pred, same_network = same_net, max_dist_m = 250.0)
    if (nrow(m_pred) >= 20) {
      mp1 <- p1_sc[m_pred$idx1]; mp2 <- p2_sc[m_pred$idx2]
      mse1 <- se1[m_pred$idx1];  mse2 <- se2[m_pred$idx2]
      if (i == 1L && j == 2L) {
        pair12_pred_match <- data.frame(
          RivID_1     = b1$eng_pred$RivID[m_pred$idx1],
          RivID_2     = b2$eng_pred$RivID[m_pred$idx2],
          CatchID_1   = b1$eng_pred$CatchID[m_pred$idx1],
          CatchID_2   = b2$eng_pred$CatchID[m_pred$idx2],
          mid_X       = b2$eng_pred$mid_X[m_pred$idx2],
          mid_Y       = b2$eng_pred$mid_Y[m_pred$idx2],
          start_X     = b2$eng_pred$start_X[m_pred$idx2],
          start_Y     = b2$eng_pred$start_Y[m_pred$idx2],
          end_X       = b2$eng_pred$end_X[m_pred$idx2],
          end_Y       = b2$eng_pred$end_Y[m_pred$idx2],
          Shreve      = b2$eng_pred$Shreve[m_pred$idx2],
          Strahler    = if ("Strahler" %in% names(b2$eng_pred)) b2$eng_pred$Strahler[m_pred$idx2] else NA_integer_,
          Elevation   = if ("Elevation" %in% names(b2$eng_pred)) b2$eng_pred$Elevation[m_pred$idx2] else NA_real_,
          Is_Aquifer  = if ("Is_Aquifer" %in% names(b2$eng_pred)) b2$eng_pred$Is_Aquifer[m_pred$idx2] else NA_integer_,
          has_obs_1   = b1$eng_pred$has_obs[m_pred$idx1],
          has_obs_2   = b2$eng_pred$has_obs[m_pred$idx2],
          obs1_raw    = b1$eng_pred$obs_value[m_pred$idx1],
          obs2_raw    = b2$eng_pred$obs_value[m_pred$idx2],
          obs1_sc     = b1$eng_pred$obs_model_sc[m_pred$idx1],
          obs2_sc     = b2$eng_pred$obs_model_sc[m_pred$idx2],
          dist_m      = m_pred$dist_m,
          pred1_tu_sc = if ("pred_tu_sc" %in% names(b1$eng_pred)) b1$eng_pred$pred_tu_sc[m_pred$idx1] else mp1,
          pred2_tu_sc = if ("pred_tu_sc" %in% names(b2$eng_pred)) b2$eng_pred$pred_tu_sc[m_pred$idx2] else mp2,
          pred1_td_sc = if ("pred_td_sc" %in% names(b1$eng_pred)) b1$eng_pred$pred_td_sc[m_pred$idx1] else mp1,
          pred2_td_sc = if ("pred_td_sc" %in% names(b2$eng_pred)) b2$eng_pred$pred_td_sc[m_pred$idx2] else mp2,
          pred1_sc    = mp1,
          pred2_sc    = mp2,
          diff_sc     = mp2 - mp1,
          pred1_raw   = b1$eng_pred$pred_combined[m_pred$idx1],
          pred2_raw   = b2$eng_pred$pred_combined[m_pred$idx2],
          diff_raw    = b2$eng_pred$pred_combined[m_pred$idx2] - b1$eng_pred$pred_combined[m_pred$idx1],
          se1         = mse1,
          se2         = mse2
        )
      }
      ct_pear  <- suppressWarnings(cor.test(mp1, mp2, method = "pearson"))
      ct_spear <- suppressWarnings(cor.test(mp1, mp2, method = "spearman"))
      tt_mp    <- suppressWarnings(t.test(mp2, mp1, paired = TRUE))
      rmsd_sc  <- sqrt(mean((mp2 - mp1)^2, na.rm = TRUE))

      add_stat_row(comp_lbl, "Spatial Co-location", "Pearson Correlation Test (Co-located Reach Predictions)",
                   "Pearson r (Spearman rho)", nrow(m_pred), nrow(m_pred),
                   as.numeric(ct_pear$estimate), as.numeric(ct_spear$estimate),
                   as.numeric(ct_pear$estimate)^2,
                   "t", as.numeric(ct_pear$statistic), ct_pear$p.value,
                   sprintf("r = %.4f (Spearman rho = %.4f, R2 = %.4f) across %d matched reaches (%s)",
                           ct_pear$estimate, ct_spear$estimate, ct_pear$estimate^2, nrow(m_pred), sig_stars(ct_pear$p.value)))

      # 4-Quadrant Chi-Squared Co-occurrence Test
      med1_m <- median(mp1, na.rm = TRUE); med2_m <- median(mp2, na.rm = TRUE)
      q_tab  <- matrix(c(sum(mp1 > med1_m & mp2 > med2_m, na.rm = TRUE),
                         sum(mp1 > med1_m & mp2 <= med2_m, na.rm = TRUE),
                         sum(mp1 <= med1_m & mp2 > med2_m, na.rm = TRUE),
                         sum(mp1 <= med1_m & mp2 <= med2_m, na.rm = TRUE)), nrow = 2, byrow = TRUE)
      csq_m  <- suppressWarnings(chisq.test(q_tab))
      conc_m <- 100.0 * (q_tab[1, 1] + q_tab[2, 2]) / max(1, sum(q_tab))
      or_m   <- (as.numeric(q_tab[1, 1]) * as.numeric(q_tab[2, 2])) / max(1.0, as.numeric(q_tab[1, 2]) * as.numeric(q_tab[2, 1]))
      add_stat_row(comp_lbl, "Reach Co-Occurrence", "Chi-Squared 4-Quadrant Concordance Test (Above/Below Median)",
                   "% Concordant (Odds Ratio)", nrow(m_pred), nrow(m_pred),
                   conc_m, or_m, conc_m - 50.0,
                   "ChiSq", as.numeric(csq_m$statistic), csq_m$p.value,
                   sprintf("%.1f%% concordant reaches (High/High + Low/Low), Odds Ratio = %.2f (%s)",
                           conc_m, or_m, sig_stars(csq_m$p.value)))

      if (same_det) {
        add_stat_row(comp_lbl, "Spatial Co-location", "Paired Student's t-Test (Co-located Reach Mean Bias)",
                     "Mean Prediction Bias", nrow(m_pred), nrow(m_pred),
                     mean(mp1, na.rm = TRUE), mean(mp2, na.rm = TRUE),
                     mean(mp2 - mp1, na.rm = TRUE),
                     "t", as.numeric(tt_mp$statistic), tt_mp$p.value,
                     sprintf("Mean bias (Run2 - Run1) = %+.4f across %d co-located reaches (%s)",
                             mean(mp2 - mp1, na.rm = TRUE), nrow(m_pred), sig_stars(tt_mp$p.value)))
      }
    }
  }

  # E. Omnibus Multi-Run Statistical Tests (when K >= 3, e.g. seasonal/yearly/monthly comparisons)
  if (K >= 3) {
    all_pred_vals <- unlist(lapply(bundles, function(b) b$eng_pred$pred_model_sc))
    all_pred_grp  <- factor(rep(run_labels, sapply(bundles, function(b) nrow(b$eng_pred))), levels = run_labels)
    aov_pred <- summary(aov(all_pred_vals ~ all_pred_grp))[[1]]
    kw_pred  <- suppressWarnings(kruskal.test(all_pred_vals ~ all_pred_grp))

    add_stat_row("All Runs (Omnibus)", "Multi-Run Omnibus", "One-Way ANOVA F-Test (Reach Predictions Across Runs)",
                 "Group Means", K, length(all_pred_vals),
                 min(tapply(all_pred_vals, all_pred_grp, mean, na.rm = TRUE)),
                 max(tapply(all_pred_vals, all_pred_grp, mean, na.rm = TRUE)),
                 diff(range(tapply(all_pred_vals, all_pred_grp, mean, na.rm = TRUE))),
                 "F", as.numeric(aov_pred$`F value`[1]), as.numeric(aov_pred$`Pr(>F)`[1]),
                 sprintf("Significant variation in mean predictions across %d runs (%s)", K, sig_stars(aov_pred$`Pr(>F)`[1])))

    add_stat_row("All Runs (Omnibus)", "Multi-Run Omnibus", "Kruskal-Wallis Rank-Sum Test (Reach Predictions Across Runs)",
                 "Group Medians", K, length(all_pred_vals),
                 min(tapply(all_pred_vals, all_pred_grp, median, na.rm = TRUE)),
                 max(tapply(all_pred_vals, all_pred_grp, median, na.rm = TRUE)),
                 diff(range(tapply(all_pred_vals, all_pred_grp, median, na.rm = TRUE))),
                 "ChiSq", as.numeric(kw_pred$statistic), kw_pred$p.value,
                 sprintf("Non-parametric median shift across %d runs (%s)", K, sig_stars(kw_pred$p.value)))

    all_obs_vals <- unlist(lapply(bundles, function(b) b$eng_pred$obs_model_sc[b$eng_pred$has_obs == 1]))
    all_obs_grp  <- factor(rep(run_labels, sapply(bundles, function(b) sum(b$eng_pred$has_obs == 1))), levels = run_labels)
    if (length(unique(all_obs_grp)) >= 2 && length(all_obs_vals) > K) {
      kw_obs <- suppressWarnings(kruskal.test(all_obs_vals ~ all_obs_grp))
      add_stat_row("All Runs (Omnibus)", "Multi-Run Omnibus", "Kruskal-Wallis Test (Observed Station Values Across Runs)",
                   "Observed Median", K, length(all_obs_vals),
                   min(tapply(all_obs_vals, all_obs_grp, median, na.rm = TRUE)),
                   max(tapply(all_obs_vals, all_obs_grp, median, na.rm = TRUE)),
                   diff(range(tapply(all_obs_vals, all_obs_grp, median, na.rm = TRUE))),
                   "ChiSq", as.numeric(kw_obs$statistic), kw_obs$p.value,
                   sprintf("Variation in observed station values across %d runs (%s)", K, sig_stars(kw_obs$p.value)))
    }
  }

  stat_tests_df <- if (length(stat_rows) > 0) do.call(rbind, stat_rows) else data.frame()
  write.csv(stat_tests_df, file.path(out_dir, "statistical_tests_summary.csv"), row.names = FALSE)

  # ----------------------------------------------------------------------------
  # 4. Reach-by-Reach Relationship & Cross-Determinand Co-Occurrence Analysis
  # ----------------------------------------------------------------------------
  rel_info <- analyze_reach_relationships(bundles, pair12_pred_match)
  if (!is.null(rel_info) && !is.null(rel_info$stratum_df)) {
    pair12_pred_match <- rel_info$pm
    write.csv(rel_info$stratum_df, file.path(out_dir, "reach_relationship_summary.csv"), row.names = FALSE)
  }

  # ----------------------------------------------------------------------------
  # 5. Parameter Estimates Comparison & Wald Z-Tests
  # ----------------------------------------------------------------------------
  param_comp_df <- NULL
  if (!any(sapply(bundles, function(b) is.null(b$param_tab)))) {
    p_list <- lapply(seq_len(K), function(i) {
      pt <- bundles[[i]]$param_tab
      pt$Run_Index <- i
      pt$Run_Label <- bundles[[i]]$label
      pt$Approx_SE <- pmax(1e-6, (pt$CI_Upper - pt$CI_Lower) / (2 * 1.96))
      pt[, c("Run_Index", "Run_Label", "Canon_Domain", "Canon_Model", "Parameter",
             "Estimate", "CI_Lower", "CI_Upper", "Approx_SE")]
    })
    param_long <- do.call(rbind, p_list)

    if (K == 1L) {
      # Single-run mode: One-sample Wald Z-test against H0: Parameter = 0
      param_long$Ref_Estimate <- 0.0
      param_long$Diff_vs_Run1 <- param_long$Estimate
      z_val <- param_long$Estimate / param_long$Approx_SE
      p_val <- 2 * pnorm(-abs(z_val))
    } else {
      # Multi-run mode: Two-sample Wald Z-test comparing against Run 1
      ref_p <- subset(param_long, Run_Index == 1L)
      ref_key <- paste(ref_p$Canon_Domain, ref_p$Canon_Model, ref_p$Parameter, sep = "||")
      all_key <- paste(param_long$Canon_Domain, param_long$Canon_Model, param_long$Parameter, sep = "||")
      m_idx   <- match(all_key, ref_key)

      param_long$Ref_Estimate <- ref_p$Estimate[m_idx]
      param_long$Diff_vs_Run1 <- ifelse(param_long$Run_Index == 1L, 0.0, param_long$Estimate - param_long$Ref_Estimate)
      z_val <- ifelse(param_long$Run_Index == 1L | is.na(m_idx), NA_real_,
                      (param_long$Estimate - ref_p$Estimate[m_idx]) /
                        sqrt(param_long$Approx_SE^2 + ref_p$Approx_SE[m_idx]^2))
      p_val <- ifelse(is.na(z_val), NA_real_, 2 * pnorm(-abs(z_val)))
    }
    param_long$Wald_Z       <- round(z_val, 4)
    param_long$Wald_P       <- p_val
    param_long$Significance <- sapply(p_val, sig_stars)
    param_comp_df <- param_long
    write.csv(param_comp_df, file.path(out_dir, "parameter_comparison.csv"), row.names = FALSE)
  }

  # ----------------------------------------------------------------------------
  # 6. Render Comparison / Single-Run Figures (Pure Base R Graphics)
  # ----------------------------------------------------------------------------
  run_cols <- c("#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e",
                "#17becf", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22",
                "#393b79", "#ad494a")[seq_len(max(4, K))]

  fig1_path <- file.path(out_dir, "fig_comp1_cv_performance.png")
  fig2_path <- file.path(out_dir, "fig_comp2_predictions_and_uncertainty.png")
  fig3_path <- file.path(out_dir, "fig_comp3_parameter_forest.png")
  fig4_path <- file.path(out_dir, "fig_comp4_side_by_side_maps.png")
  fig5_path <- file.path(out_dir, "fig_comp5_reach_relationships.png")

  plot_comp1_cv(bundles, cv_comp_df, pair12_cv_match, stat_tests_df, run_cols, fig1_path)
  plot_comp2_predictions(bundles, pred_sum_df, pair12_pred_match, stat_tests_df, run_cols, fig2_path)
  if (!is.null(param_comp_df)) {
    plot_comp3_parameters(bundles, param_comp_df, run_cols, fig3_path)
  }
  plot_comp4_side_by_side_maps(bundles, fig4_path)
  if (!is.null(rel_info)) {
    plot_comp5_reach_relationships(bundles, strat_df, rel_info, run_cols, fig5_path)
  } else if (K == 1L) {
    plot_comp5_single_run_descriptives(bundles[[1]], strat_df, fig5_path)
  }

  # ----------------------------------------------------------------------------
  # 7. Write Comprehensive Markdown Report
  # ----------------------------------------------------------------------------
  report_path <- file.path(out_dir, "comparison_report.md")
  write_comparison_markdown_report(bundles, cv_comp_df, pred_sum_df, desc_df, rep_df,
                                   strat_df, stat_tests_df, rel_info,
                                   param_comp_df, pair12_cv_match, pair12_pred_match,
                                   out_dir, report_path)

  if (!quiet) {
    cat(sprintf("\n%s complete! Generated files:\n", if (K == 1L) "Single-run descriptive analysis" else "Comparison"))
    cat(sprintf("  - Markdown Report:   %s\n", report_path))
    cat(sprintf("  - CV Table:          %s\n", file.path(out_dir, "cv_metrics_comparison.csv")))
    cat(sprintf("  - Descriptives CSV:  %s\n", file.path(out_dir, "descriptive_stats_summary.csv")))
    cat(sprintf("  - Stratified CSV:    %s\n", file.path(out_dir, "stratified_results_summary.csv")))
    if (!is.null(rel_info)) {
      cat(sprintf("  - Reach Rel CSV:     %s\n", file.path(out_dir, "reach_relationship_summary.csv")))
    }
    cat(sprintf("  - Stat Tests:        %s\n", file.path(out_dir, "statistical_tests_summary.csv")))
    cat(sprintf("  - Pred Summary:      %s\n", file.path(out_dir, "prediction_distribution_summary.csv")))
    if (!is.null(param_comp_df)) {
      cat(sprintf("  - Param Table:       %s\n", file.path(out_dir, "parameter_comparison.csv")))
    }
    cat(sprintf("  - Figure 1 (CV):     %s\n", fig1_path))
    cat(sprintf("  - Figure 2 (Pred):   %s\n", fig2_path))
    if (!is.null(param_comp_df)) {
      cat(sprintf("  - Figure 3 (Pars):   %s\n", fig3_path))
    }
    cat(sprintf("  - Figure 4 (Maps):   %s\n", fig4_path))
    if (!is.null(rel_info) || K == 1L) {
      cat(sprintf("  - Figure 5 (Strat):  %s\n", fig5_path))
    }
  }

  invisible(list(
    out_dir       = out_dir,
    report_path   = report_path,
    cv_comp_df    = cv_comp_df,
    pred_sum_df   = pred_sum_df,
    desc_df       = desc_df,
    strat_df      = strat_df,
    rel_info      = rel_info,
    stat_tests_df = stat_tests_df,
    param_comp_df = param_comp_df
  ))
}

# ==============================================================================
# Plotting Functions for Model Run Comparison (Pure Base R)
# ==============================================================================

# Figure 1: 4-Panel Cross-Validation Performance Comparison
plot_comp1_cv <- function(bundles, cv_comp_df, pair12_cv_match, stat_tests_df, run_cols, out_file) {
  K <- length(bundles)
  run_labels <- sapply(bundles, function(b) b$label)
  det_short  <- shorten_det_name(bundles[[1]]$det_name)

  png(out_file, width = 2600, height = 2000, res = 200, type = "cairo")
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(2, 2), mar = c(5.2, 4.8, 3.6, 1.8), oma = c(0.5, 0.5, 2.2, 0.5), bg = "white")

  # Core domains present across runs
  core_doms <- c("All England RSN Network", "Humber Basin", "Thames Basin", "Severn Basin")
  avail_doms <- unique(cv_comp_df$Canon_Domain)
  doms_use <- c(intersect(core_doms, avail_doms), setdiff(avail_doms, core_doms))
  dom_short <- sub(" RSN Network$", "", sub(" Basin$", "", doms_use))

  # Panel A: S3N Combined 5-Fold CV R2 by Domain and Run
  r2_mat <- matrix(NA_real_, nrow = K, ncol = length(doms_use),
                   dimnames = list(run_labels, dom_short))
  rmse_mat <- r2_mat
  gain_mat <- r2_mat
  for (i in seq_len(K)) {
    for (d_idx in seq_along(doms_use)) {
      sub_r <- subset(cv_comp_df, Run_Index == i & Canon_Domain == doms_use[d_idx] & Canon_Model == "S3N Combined")
      if (nrow(sub_r) > 0) {
        r2_mat[i, d_idx]   <- sub_r$R2[1]
        rmse_mat[i, d_idx] <- sub_r$RMSE[1]
        gain_mat[i, d_idx] <- sub_r$Delta_R2_vs_OLS[1]
      }
    }
  }

  r2_plot <- ifelse(is.na(r2_mat), 0, pmax(0, r2_mat))
  ymax_a  <- max(0.85, max(r2_plot, na.rm = TRUE) * 1.22)
  bp_a <- matrix(barplot(r2_plot, beside = TRUE, col = run_cols, border = "#222222",
                         ylim = c(0, ymax_a), ylab = "5-Fold Cross-Validation R2",
                         main = sprintf("A. Out-of-Sample S3N Combined CV R2 (%s)", det_short),
                         cex.main = 1.12, cex.lab = 1.0, cex.names = 0.95, las = 1),
                 nrow = K)
  grid(nx = NA, ny = NULL, col = "#e0e0e0", lty = "dotted")
  abline(h = 0, col = "#333333", lwd = 1)
  for (i in seq_len(K)) {
    for (d_idx in seq_along(doms_use)) {
      val <- r2_mat[i, d_idx]
      if (!is.na(val)) {
        y_pos <- max(0, val) + 0.025 * ymax_a
        text(bp_a[i, d_idx], y_pos, labels = sprintf("%.3f", val), cex = 0.72, font = 2, col = "#111111")
      }
    }
  }
  legend("topright", legend = run_labels, fill = run_cols, border = "#222222",
         bty = "n", cex = 0.82)

  # Panel B: S3N Combined 5-Fold CV RMSE by Domain and Run
  rmse_plot <- ifelse(is.na(rmse_mat), 0, rmse_mat)
  ymax_b <- max(rmse_plot, na.rm = TRUE) * 1.25
  bp_b <- matrix(barplot(rmse_plot, beside = TRUE, col = run_cols, border = "#222222",
                         ylim = c(0, ymax_b), ylab = "5-Fold CV RMSE (Model Scale)",
                         main = sprintf("B. Out-of-Sample S3N Combined CV RMSE (Lower is Better)"),
                         cex.main = 1.12, cex.lab = 1.0, cex.names = 0.95, las = 1),
                 nrow = K)
  grid(nx = NA, ny = NULL, col = "#e0e0e0", lty = "dotted")
  for (i in seq_len(K)) {
    for (d_idx in seq_along(doms_use)) {
      val <- rmse_mat[i, d_idx]
      if (!is.na(val)) {
        text(bp_b[i, d_idx], val + 0.025 * ymax_b, labels = sprintf("%.3f", val), cex = 0.72, font = 2, col = "#111111")
      }
    }
  }
  legend("topright", legend = run_labels, fill = run_cols, border = "#222222",
         bty = "n", cex = 0.82)

  # Panel C: Empirical CDF of Out-of-Sample Absolute CV Errors (All England)
  has_cv_pts <- sapply(bundles, function(b) !is.null(b$eng_cv) && nrow(b$eng_cv) > 0)
  if (K == 1L && has_cv_pts[1]) {
    cv1 <- bundles[[1]]$eng_cv
    mod_errs <- list(
      "S3N Combined"      = abs(cv1$obs_y - cv1$cv_comb),
      "S3N Tail-Up"       = if ("cv_tu"  %in% names(cv1)) abs(cv1$obs_y - cv1$cv_tu)  else NULL,
      "S3N Tail-Down"     = if ("cv_td"  %in% names(cv1)) abs(cv1$obs_y - cv1$cv_td)  else NULL,
      "OLS (Non-Spatial)" = if ("cv_ols" %in% names(cv1)) abs(cv1$obs_y - cv1$cv_ols) else NULL
    )
    mod_errs <- mod_errs[!sapply(mod_errs, is.null)]
    mod_cols <- c("#1f77b4", "#2ca02c", "#ff7f0e", "#d62728")[seq_along(mod_errs)]
    xmax_c <- quantile(unlist(mod_errs), probs = 0.98, na.rm = TRUE)
    plot(NULL, xlim = c(0, xmax_c), ylim = c(0, 1),
         xlab = "Absolute Out-of-Sample CV Error |y_obs - y_cv| (Model Scale)",
         ylab = "Cumulative Proportion of Stations",
         main = "C. Empirical CDF of Station CV Errors Across Models",
         cex.main = 1.12, cex.lab = 1.0, las = 1)
    grid(col = "#e0e0e0", lty = "dotted")
    xs <- seq(0, xmax_c, length.out = 250)
    for (m_i in seq_along(mod_errs)) {
      ec <- ecdf(mod_errs[[m_i]])
      lines(xs, ec(xs), col = mod_cols[m_i], lwd = 2.4)
    }
    leg_c <- sapply(seq_along(mod_errs), function(m_i) {
      sprintf("%s (Med |e| = %.3f)", names(mod_errs)[m_i], median(mod_errs[[m_i]], na.rm = TRUE))
    })
    legend("bottomright", legend = leg_c, col = mod_cols, lwd = 2.5, bty = "n", cex = 0.80)
    if ("OLS (Non-Spatial)" %in% names(mod_errs)) {
      wt_p <- safe_paired_wilcox(mod_errs[["S3N Combined"]], mod_errs[["OLS (Non-Spatial)"]])$p.value
      dm_p <- safe_paired_ttest(mod_errs[["S3N Combined"]]^2, mod_errs[["OLS (Non-Spatial)"]]^2)$p.value
      mtext(sprintf("S3N Combined vs OLS (n = %d): Paired Wilcoxon p = %s (%s) | Paired MSE t-test p = %s (%s)",
                    nrow(cv1), fmt_pval(wt_p), sig_stars(wt_p), fmt_pval(dm_p), sig_stars(dm_p)),
            side = 3, line = 0.25, cex = 0.76, col = "#333333")
    }
  } else if (any(has_cv_pts)) {
    err_list <- lapply(bundles, function(b) {
      if (!is.null(b$eng_cv)) abs(b$eng_cv$obs_y - b$eng_cv$cv_comb) else numeric(0)
    })
    xmax_c <- quantile(unlist(err_list), probs = 0.98, na.rm = TRUE)
    plot(NULL, xlim = c(0, xmax_c), ylim = c(0, 1),
         xlab = "Absolute Out-of-Sample CV Error |y_obs - y_cv| (Model Scale)",
         ylab = "Cumulative Proportion of Stations",
         main = "C. Empirical CDF of Station CV Errors (All England)",
         cex.main = 1.12, cex.lab = 1.0, las = 1)
    grid(col = "#e0e0e0", lty = "dotted")
    for (i in which(has_cv_pts)) {
      ec <- ecdf(err_list[[i]])
      xs <- seq(0, xmax_c, length.out = 250)
      lines(xs, ec(xs), col = run_cols[i], lwd = 2.4)
    }
    leg_c <- sapply(which(has_cv_pts), function(i) {
      sprintf("%s (Med |e| = %.3f, n = %d)", run_labels[i], median(err_list[[i]], na.rm = TRUE), length(err_list[[i]]))
    })
    legend("bottomright", legend = leg_c, col = run_cols[which(has_cv_pts)], lwd = 2.5, bty = "n", cex = 0.80)
    if (!is.null(pair12_cv_match) && nrow(pair12_cv_match) >= 10) {
      wt_p <- safe_paired_wilcox(pair12_cv_match$abs_e2, pair12_cv_match$abs_e1)$p.value
      dm_p <- safe_paired_ttest(pair12_cv_match$sq_e2, pair12_cv_match$sq_e1)$p.value
      mtext(sprintf("Co-located Stations (n = %d): Paired Wilcoxon p = %s (%s) | Paired MSE t-test p = %s (%s)",
                    nrow(pair12_cv_match), fmt_pval(wt_p), sig_stars(wt_p), fmt_pval(dm_p), sig_stars(dm_p)),
            side = 3, line = 0.25, cex = 0.78, col = "#333333")
    }
  } else {
    plot.new()
  }

  # Panel D: Spatial Gain over Non-Spatial OLS (Delta R2 = R2_S3N - R2_OLS)
  gain_plot <- ifelse(is.na(gain_mat), 0, gain_mat)
  yr_d <- range(c(0, gain_plot), na.rm = TRUE)
  pad_d <- 0.22 * max(0.05, diff(yr_d))
  bp_d <- matrix(barplot(gain_plot, beside = TRUE, col = run_cols, border = "#222222",
                         ylim = c(min(0, yr_d[1] - pad_d), yr_d[2] + pad_d),
                         ylab = "Delta CV R2 (S3N Combined minus OLS)",
                         main = "D. Spatial Network Gain Over Non-Spatial OLS (Delta R2)",
                         cex.main = 1.12, cex.lab = 1.0, cex.names = 0.95, las = 1),
                 nrow = K)
  grid(nx = NA, ny = NULL, col = "#e0e0e0", lty = "dotted")
  abline(h = 0, col = "#333333", lwd = 1)
  for (i in seq_len(K)) {
    for (d_idx in seq_along(doms_use)) {
      val <- gain_mat[i, d_idx]
      if (!is.na(val)) {
        y_pos <- if (val >= 0) val + 0.08 * pad_d else val - 0.12 * pad_d
        text(bp_d[i, d_idx], y_pos, labels = sprintf("%+.3f", val), cex = 0.72, font = 2, col = "#111111")
      }
    }
  }
  legend("topright", legend = run_labels, fill = run_cols, border = "#222222",
         bty = "n", cex = 0.82)

  mtext(sprintf("S3N Model Run Comparison: Cross-Validation Performance (%s)", det_short),
        side = 3, line = 0.4, outer = TRUE, cex = 1.15, font = 2, col = "#111111")
}

# Figure 2: 4-Panel Reach Prediction & Uncertainty Comparison
plot_comp2_predictions <- function(bundles, pred_sum_df, pair12_pred_match, stat_tests_df, run_cols, out_file) {
  K <- length(bundles)
  run_labels <- sapply(bundles, function(b) b$label)
  same_det   <- (length(unique(sapply(bundles, function(b) b$det_code))) == 1L)
  det_short  <- if (same_det) shorten_det_name(bundles[[1]]$det_name) else sprintf("%s vs. %s", shorten_det_name(bundles[[1]]$det_name), shorten_det_name(bundles[[2]]$det_name))
  det_unit   <- bundles[[1]]$det_unit
  use_log    <- bundles[[1]]$use_log

  png(out_file, width = 2600, height = 2000, res = 200, type = "cairo")
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(2, 2), mar = c(4.8, 4.8, 3.6, 1.8), oma = c(0.5, 0.5, 2.2, 0.5), bg = "white")

  # Panel A: Co-located Reach Scatter Plot (Run 1 vs Run 2)
  if (!is.null(pair12_pred_match) && nrow(pair12_pred_match) >= 10) {
    pm <- pair12_pred_match
    set.seed(42)
    idx_pt <- if (nrow(pm) > 8000) sample(nrow(pm), 8000) else seq_len(nrow(pm))
    if (same_det) {
      xlim_a <- quantile(c(pm$pred1_sc, pm$pred2_sc), probs = c(0.005, 0.995), na.rm = TRUE)
      ylim_a <- xlim_a
    } else {
      xlim_a <- quantile(pm$pred1_sc, probs = c(0.005, 0.995), na.rm = TRUE)
      ylim_a <- quantile(pm$pred2_sc, probs = c(0.005, 0.995), na.rm = TRUE)
    }
    u1_lbl <- if (bundles[[1]]$use_log) sprintf("log(%s)", bundles[[1]]$det_unit) else bundles[[1]]$det_unit
    u2_lbl <- if (bundles[[2]]$use_log) sprintf("log(%s)", bundles[[2]]$det_unit) else bundles[[2]]$det_unit
    plot(pm$pred1_sc[idx_pt], pm$pred2_sc[idx_pt],
         xlim = xlim_a, ylim = ylim_a, pch = 16, cex = 0.45, col = rgb(0.12, 0.47, 0.71, 0.25),
         xlab = sprintf("Run 1: %s [%s]", run_labels[1], u1_lbl),
         ylab = sprintf("Run 2: %s [%s]", run_labels[2], u2_lbl),
         main = sprintf("A. Co-located Reach Predictions (n = %s matched reaches)",
                        format(nrow(pm), big.mark = ",")),
         cex.main = 1.12, cex.lab = 0.95, las = 1)
    grid(col = "#e0e0e0", lty = "dotted")
    fit_a <- lm(pred2_sc ~ pred1_sc, data = pm)
    abline(fit_a, col = "#111111", lwd = 2.0, lty = 1)
    if (same_det) {
      abline(0, 1, col = "#c0392b", lwd = 2, lty = 2)
    }
    r_p   <- cor(pm$pred1_sc, pm$pred2_sc, use = "complete.obs", method = "pearson")
    rho_s <- cor(pm$pred1_sc, pm$pred2_sc, use = "complete.obs", method = "spearman")
    rmsd  <- sqrt(mean((pm$pred2_sc - pm$pred1_sc)^2, na.rm = TRUE))
    bias  <- mean(pm$pred2_sc - pm$pred1_sc, na.rm = TRUE)
    if (same_det) {
      legend("topleft",
             legend = c(
               sprintf("Pearson r = %.4f (Spearman rho = %.4f)", r_p, rho_s),
               sprintf("RMSD = %.4f | Mean Bias = %+.4f", rmsd, bias),
               "OLS Fit Line",
               "1:1 Agreement Line"
             ),
             col = c(NA, NA, "#111111", "#c0392b"), lty = c(NA, NA, 1, 2), lwd = c(NA, NA, 2, 2),
             bty = "n", cex = 0.82, text.col = "#111111")
    } else {
      legend("topleft",
             legend = c(
               sprintf("Pearson r = %.4f (Spearman rho = %.4f)", r_p, rho_s),
               sprintf("OLS Slope = %+.3f (R2 = %.3f)", coef(fit_a)[2], summary(fit_a)$r.squared),
               "OLS Regression Fit"
             ),
             col = c(NA, NA, "#111111"), lty = c(NA, NA, 1), lwd = c(NA, NA, 2),
             bty = "n", cex = 0.82, text.col = "#111111")
    }
  } else {
    plot.new()
    title("A. Co-located Reach Predictions (No overlapping reaches)")
  }

  # Panel B: Kernel Density Curves of Predicted Concentrations Across Runs
  if (same_det) {
    dens_list <- lapply(bundles, function(b) density(b$eng_pred$pred_model_sc, na.rm = TRUE, n = 512))
    xlab_b <- if (use_log) sprintf("Predicted log(%s) [%s]", det_short, det_unit) else sprintf("Predicted %s (%s)", det_short, det_unit)
    main_b <- "B. Distribution of Reach-Level Predictions Across Runs"
  } else {
    dens_list <- lapply(bundles, function(b) {
      z_v <- (b$eng_pred$pred_model_sc - mean(b$eng_pred$pred_model_sc, na.rm = TRUE)) /
             max(1e-9, sd(b$eng_pred$pred_model_sc, na.rm = TRUE))
      density(z_v, na.rm = TRUE, n = 512)
    })
    xlab_b <- "Standardized Reach Prediction Z-Score ((y - mean) / SD)"
    main_b <- "B. Standardized Distribution Shape Across Determinands"
  }
  xr_b <- range(unlist(lapply(dens_list, function(d) d$x)), na.rm = TRUE)
  yr_b <- c(0, max(unlist(lapply(dens_list, function(d) d$y)), na.rm = TRUE) * 1.18)
  plot(NULL, xlim = xr_b, ylim = yr_b,
       xlab = xlab_b,
       ylab = "Probability Density Across Reaches",
       main = main_b,
       cex.main = 1.12, cex.lab = 0.95, las = 1)
  grid(col = "#e0e0e0", lty = "dotted")
  for (i in seq_len(K)) {
    lines(dens_list[[i]]$x, dens_list[[i]]$y, col = run_cols[i], lwd = 2.5)
    med_v <- if (same_det) median(bundles[[i]]$eng_pred$pred_model_sc, na.rm = TRUE) else 0.0
    abline(v = med_v, col = run_cols[i], lty = 3, lwd = 1.2)
  }
  leg_b <- sapply(seq_len(K), function(i) {
    sprintf("%s (Med = %.3g %s, n = %s)",
            run_labels[i], pred_sum_df$Pred_Median_Raw[i], bundles[[i]]$det_unit,
            format(pred_sum_df$n_reaches[i], big.mark = ","))
  })
  legend("topright", legend = leg_b, col = run_cols, lwd = 2.5, bty = "n", cex = 0.78)
  ks_p <- suppressWarnings(ks.test(bundles[[1]]$eng_pred$pred_model_sc, bundles[[2]]$eng_pred$pred_model_sc))
  mtext(sprintf("Two-Sample KS Test (Run 1 vs Run 2): D = %.4f, p = %s (%s)",
                ks_p$statistic, fmt_pval(ks_p$p.value), sig_stars(ks_p$p.value)),
        side = 3, line = 0.25, cex = 0.78, col = "#333333")

  # Panel C: Boxplot Comparison of Kriging Standard Error (se_model_sc) Across Runs
  se_list <- lapply(bundles, function(b) b$eng_pred$se_model_sc)
  short_lbls <- substr(run_labels, 1, 34)
  names(se_list) <- short_lbls
  bp_tmp <- boxplot(se_list, plot = FALSE)
  ylim_c <- c(min(bp_tmp$stats[1, ], na.rm = TRUE) * 0.94,
              max(bp_tmp$stats[5, ], na.rm = TRUE) * 1.12)
  bp_c <- boxplot(se_list, col = adjustcolor(run_cols, alpha.f = 0.65), border = "#222222",
                  outline = FALSE, las = 1, ylim = ylim_c,
                  ylab = if (same_det && use_log) "Kriging Standard Error (log scale)" else "Kriging Standard Error (Model Scale)",
                  main = "C. Kriging Prediction Uncertainty (se_model_sc) Across Runs",
                  cex.main = 1.12, cex.lab = 0.95, cex.axis = 0.85)
  grid(nx = NA, ny = NULL, col = "#e0e0e0", lty = "dotted")
  se_means <- sapply(se_list, mean, na.rm = TRUE)
  points(seq_len(K), se_means, pch = 18, cex = 1.5, col = "#111111")
  text(seq_len(K), bp_c$stats[5, ], labels = sprintf("Mean = %.3f", se_means), pos = 3, cex = 0.80, font = 2)
  wt_se_p <- suppressWarnings(wilcox.test(se_list[[1]], se_list[[2]], paired = FALSE))$p.value
  mtext(sprintf("Mann-Whitney U Test on Kriging SE (Run 1 vs Run 2): p = %s (%s)",
                fmt_pval(wt_se_p), sig_stars(wt_se_p)),
        side = 3, line = 0.25, cex = 0.78, col = "#333333")

  # Panel D: Spatial Map of Co-located Prediction Differences (Run 2 minus Run 1 or Z2 - Z1)
  if (!is.null(pair12_pred_match) && nrow(pair12_pred_match) >= 20) {
    pm <- pair12_pred_match
    diff_vec <- if (same_det) pm$diff_sc else if ("z_diff" %in% names(pm)) pm$z_diff else {
      z1 <- (pm$pred1_sc - mean(pm$pred1_sc, na.rm = TRUE)) / max(1e-9, sd(pm$pred1_sc, na.rm = TRUE))
      z2 <- (pm$pred2_sc - mean(pm$pred2_sc, na.rm = TRUE)) / max(1e-9, sd(pm$pred2_sc, na.rm = TRUE))
      z2 - z1
    }
    q_d <- quantile(abs(diff_vec), probs = 0.98, na.rm = TRUE)
    max_abs <- max(0.1, as.numeric(q_d))
    div_pal <- colorRampPalette(c("#2166ac", "#67a9cf", "#d1e5f0", "#f7f7f7", "#fddbc7", "#ef8a62", "#b2182b"))(100)

    par(mar = c(1.5, 1.5, 3.6, 1.5))
    xr <- range(c(pm$start_X, pm$end_X), na.rm = TRUE)
    yr <- range(c(pm$start_Y, pm$end_Y), na.rm = TRUE)
    pad_x <- diff(xr) * 0.16
    main_d <- if (same_det) {
      sprintf("D. Spatial Difference Map: Run 2 minus Run 1 (%s)", if (use_log) "log scale" else det_unit)
    } else {
      "D. Standardized Spatial Divergence Map: Z(Run 2) minus Z(Run 1)"
    }
    plot(NULL, xlim = c(xr[1] - diff(xr) * 0.03, xr[2] + pad_x), ylim = yr,
         asp = 1, axes = FALSE, xlab = "", ylab = "",
         main = main_d, cex.main = 1.08)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#f8f9fa", border = "#dcdcdc")
    cols_d <- map_to_colors(diff_vec, -max_abs, max_abs, div_pal)
    lwds_d <- pmin(2.2, 0.35 + 0.35 * log10(pmax(1, pm$Shreve)))
    segments(pm$start_X, pm$start_Y, pm$mid_X, pm$mid_Y, col = cols_d, lwd = lwds_d)
    segments(pm$mid_X,   pm$mid_Y,   pm$end_X, pm$end_Y, col = cols_d, lwd = lwds_d)

    ticks_d <- seq(-max_abs, max_abs, length.out = 5)
    lbls_d  <- sprintf("%+.2f", ticks_d)
    draw_colorbar(xr[2] + diff(xr) * 0.02, yr[1] + diff(yr) * 0.18,
                  xr[2] + diff(xr) * 0.05, yr[1] + diff(yr) * 0.72,
                  vmin = -max_abs, vmax = max_abs, pal = div_pal,
                  ticks = ticks_d, tick_labels = lbls_d,
                  title_txt = if (same_det) "Run 2 - Run 1" else "Z(R2) - Z(R1)")
    mtext(sprintf("Red = Relatively Higher in %s | Blue = Relatively Higher in %s", run_labels[2], run_labels[1]),
          side = 3, line = 0.2, cex = 0.76, col = "#333333")
  } else {
    plot.new()
  }

  mtext(sprintf("S3N Model Run Comparison: Predictions & Uncertainty (%s)", det_short),
        side = 3, line = 0.4, outer = TRUE, cex = 1.15, font = 2, col = "#111111")
}

# Figure 3: 2-Panel Fixed-Effect & Spatial Covariance Parameter Comparison
plot_comp3_parameters <- function(bundles, param_comp_df, run_cols, out_file) {
  K <- length(bundles)
  run_labels <- sapply(bundles, function(b) b$label)
  same_det   <- (length(unique(sapply(bundles, function(b) b$det_code))) == 1L)
  det_short  <- if (same_det) shorten_det_name(bundles[[1]]$det_name) else sprintf("%s vs. %s", shorten_det_name(bundles[[1]]$det_name), shorten_det_name(bundles[[2]]$det_name))

  png(out_file, width = 2600, height = 1350, res = 200, type = "cairo")
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, 2), mar = c(5.0, 9.5, 3.8, 2.0), oma = c(0.5, 0.5, 2.2, 0.5), bg = "white")

  # Panel A: Forest Plot of All England Tail-Up Fixed-Effect Coefficients
  cov_names <- c("sigma.sq", "tau.sq", "phi")
  fe_sub <- subset(param_comp_df, Canon_Domain == "All England RSN Network" &
                                  Canon_Model == "S3N Tail-Up" &
                                  !(Parameter %in% c("Intercept", cov_names)))
  if (nrow(fe_sub) > 0) {
    params_u <- unique(fe_sub$Parameter)
    n_p <- length(params_u)
    xlim_a <- range(c(fe_sub$CI_Lower, fe_sub$CI_Upper, 0), na.rm = TRUE)
    pad_a  <- 0.15 * max(0.2, diff(xlim_a))
    plot(NULL, xlim = c(xlim_a[1] - pad_a, xlim_a[2] + pad_a), ylim = c(0.4, n_p + 1.45),
         yaxt = "n", xlab = "Standardized Regression Coefficient (Estimate & 95% Bootstrap CI)",
         ylab = "", main = "A. Fixed-Effect Drivers (All England S3N Tail-Up)",
         cex.main = 1.10, cex.lab = 0.95, las = 1)
    grid(nx = NULL, ny = NA, col = "#e0e0e0", lty = "dotted")
    abline(v = 0, col = "#c0392b", lwd = 1.5, lty = 2)
    axis(2, at = seq_len(n_p), labels = rev(params_u), las = 1, cex.axis = 0.88)

    offsets <- seq(-0.22, 0.22, length.out = K)
    if (K == 1) offsets <- 0
    for (p_idx in seq_len(n_p)) {
      p_name <- rev(params_u)[p_idx]
      for (i in seq_len(K)) {
        r_p <- subset(fe_sub, Parameter == p_name & Run_Index == i)
        if (nrow(r_p) > 0) {
          y_pos <- p_idx + offsets[i]
          segments(r_p$CI_Lower[1], y_pos, r_p$CI_Upper[1], y_pos, col = run_cols[i], lwd = 2.2)
          points(r_p$Estimate[1], y_pos, pch = 16, cex = 1.15, col = run_cols[i])
          if (i > 1 && !is.na(r_p$Wald_P[1]) && r_p$Wald_P[1] < 0.05) {
            text(r_p$CI_Upper[1], y_pos, labels = sprintf(" %s", r_p$Significance[1]),
                 adj = c(0, 0.5), cex = 0.85, font = 2, col = run_cols[i])
          }
        }
      }
    }
    legend("topleft", legend = c(run_labels, "* p < 0.05 vs Run 1 (Wald Z)"),
           col = c(run_cols, NA), pch = c(rep(16, K), NA), lwd = c(rep(2, K), NA),
           bty = "n", cex = 0.80)
  } else {
    plot.new()
  }

  # Panel B: Spatial Variance Partitioning (Partial Sill sigma.sq vs Nugget tau.sq)
  par(mar = c(5.0, 5.0, 3.8, 2.0))
  cov_sub <- subset(param_comp_df, Canon_Domain == "All England RSN Network" & Parameter %in% c("sigma.sq", "tau.sq"))
  if (nrow(cov_sub) > 0) {
    mods_u <- c("S3N Tail-Up", "S3N Tail-Down")
    bar_vals <- matrix(0, nrow = 2, ncol = K * length(mods_u))
    bar_names <- character(K * length(mods_u))
    pct_spatial <- numeric(K * length(mods_u))
    col_idx <- 0L
    for (m_name in mods_u) {
      m_short <- if (m_name == "S3N Tail-Up") "TU" else "TD"
      for (i in seq_len(K)) {
        col_idx <- col_idx + 1L
        s_row <- subset(cov_sub, Canon_Model == m_name & Run_Index == i & Parameter == "sigma.sq")
        t_row <- subset(cov_sub, Canon_Model == m_name & Run_Index == i & Parameter == "tau.sq")
        s_val <- if (nrow(s_row) > 0) s_row$Estimate[1] else 0
        t_val <- if (nrow(t_row) > 0) t_row$Estimate[1] else 0
        bar_vals[1, col_idx] <- s_val
        bar_vals[2, col_idx] <- t_val
        bar_names[col_idx]   <- sprintf("R%d (%s)", i, m_short)
        pct_spatial[col_idx] <- if (s_val + t_val > 0) 100 * s_val / (s_val + t_val) else 0
      }
    }
    ymax_b <- max(colSums(bar_vals), na.rm = TRUE) * 1.25
    bp_b <- barplot(bar_vals, names.arg = bar_names, col = c("#2980b9", "#f39c12"),
                    border = "#222222", ylim = c(0, ymax_b), las = 2,
                    ylab = "Variance Component Estimate (Model Scale)",
                    main = "B. Spatial Sill (sigma^2) vs. Nugget (tau^2) Partitioning",
                    cex.main = 1.10, cex.lab = 0.95, cex.names = 0.85)
    grid(nx = NA, ny = NULL, col = "#e0e0e0", lty = "dotted")
    tot_v <- colSums(bar_vals)
    text(bp_b, tot_v + 0.03 * ymax_b, labels = sprintf("%.0f%% Sp", pct_spatial),
         cex = 0.76, font = 2, col = "#111111")
    legend("topright",
           legend = c("Spatial Partial Sill (sigma^2)", "Non-Spatial Nugget (tau^2)"),
           fill = c("#2980b9", "#f39c12"), border = "#222222", bty = "n", cex = 0.85)
  } else {
    plot.new()
  }

  mtext(sprintf("S3N Model Run Comparison: Fixed-Effect & Covariance Parameters (%s)", det_short),
        side = 3, line = 0.4, outer = TRUE, cex = 1.15, font = 2, col = "#111111")
}

# Figure 4: Side-by-Side All England Stream Network Maps
plot_comp4_side_by_side_maps <- function(bundles, out_file) {
  K <- length(bundles)
  same_det  <- (length(unique(sapply(bundles, function(b) b$det_code))) == 1L)
  use_log   <- bundles[[1]]$use_log
  det_short <- shorten_det_name(bundles[[1]]$det_name)
  det_unit  <- bundles[[1]]$det_unit

  all_vals <- unlist(lapply(bundles, function(b) {
    c(b$eng_pred$pred_model_sc, b$eng_pred$obs_model_sc[b$eng_pred$has_obs == 1])
  }))
  shared_sc <- compute_color_scale(all_vals, use_log)
  pal_use   <- get_val_pal(FALSE)

  if (K <= 2) {
    nr <- 1; nc <- K; w_px <- 1350 * K; h_px <- 1400
  } else if (K <= 4) {
    nr <- 2; nc <- 2; w_px <- 2500; h_px <- 2500
  } else {
    nc <- 3; nr <- ceiling(K / 3); w_px <- 3200; h_px <- 1100 * nr
  }

  png(out_file, width = w_px, height = h_px, res = 200, type = "cairo")
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(nr, nc), mar = c(1.2, 1.2, 3.2, 1.2), oma = c(0, 0, 2.5, 0), bg = "white")

  for (i in seq_len(K)) {
    b <- bundles[[i]]
    sc_i <- if (same_det) {
      shared_sc
    } else {
      vals_i <- c(b$eng_pred$pred_model_sc, b$eng_pred$obs_model_sc[b$eng_pred$has_obs == 1])
      compute_color_scale(vals_i, b$use_log)
    }
    d_short_i <- shorten_det_name(b$det_name)
    r_all <- subset(b$cv_tab, Canon_Domain == "All England RSN Network" & Canon_Model == "S3N Combined")
    sub_txt <- if (nrow(r_all) > 0) {
      sprintf("5-Fold CV R2 = %.3f (RMSE = %.3f) | Med = %.3g %s | %s reaches",
              r_all$R2[1], r_all$RMSE[1], median(b$eng_pred$pred_combined, na.rm = TRUE),
              b$det_unit, format(nrow(b$eng_pred), big.mark = ","))
    } else {
      sprintf("%s reaches (%d observed stations)", format(nrow(b$eng_pred), big.mark = ","), sum(b$eng_pred$has_obs == 1))
    }

    plot_stream_network(
      b$eng_pred, "pred_model_sc", vmin = sc_i$vmin, vmax = sc_i$vmax, pal = pal_use,
      main_title = sprintf("Run %d: %s", i, b$label),
      sub_title  = sub_txt,
      cbar_ticks = sc_i$ticks, cbar_labels = sc_i$labels,
      cbar_title = sprintf("%s\n(%s)", sub(" \\[.*\\]$", "", d_short_i), b$det_unit),
      show_obs = FALSE, base_lwd = 0.25, max_lwd = 2.2
    )
  }

  outer_title <- if (same_det) {
    sprintf("Side-by-Side All England S3N Comparison: %s (%s)", bundles[[1]]$det_name, det_unit)
  } else {
    sprintf("Side-by-Side All England S3N Comparison: %s", paste(sapply(bundles, function(b) shorten_det_name(b$det_name)), collapse = " vs. "))
  }
  mtext(outer_title, side = 3, line = 0.5, outer = TRUE, cex = 1.20, font = 2, col = "#111111")
}

# Figure 5: 4-Panel Reach-by-Reach Relationship, Stratified Results & Co-Occurrence Map
plot_comp5_reach_relationships <- function(bundles, strat_df, rel_info, run_cols, out_file) {
  K  <- length(bundles)
  b1 <- bundles[[1]]
  b2 <- bundles[[2]]
  pm <- rel_info$pm
  run_labels <- sapply(bundles, function(b) b$label)

  png(out_file, width = 2600, height = 2000, res = 200, type = "cairo")
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(2, 2), mar = c(5.0, 4.8, 3.8, 1.8), oma = c(0.5, 0.5, 2.2, 0.5), bg = "white")

  q_cols <- c(
    "Q1_High_High"  = "#c0392b",
    "Q2_Low_Low"    = "#2980b9",
    "Q3_High1_Low2" = "#e67e22",
    "Q4_Low1_High2" = "#16a085"
  )

  # Panel A: Reach-by-Reach Bivariate Relationship Colored by 4 Co-Occurrence Quadrants
  set.seed(42)
  idx_pt <- if (nrow(pm) > 9000) sample(nrow(pm), 9000) else seq_len(nrow(pm))
  xlim_a <- quantile(pm$pred1_sc, probs = c(0.005, 0.995), na.rm = TRUE)
  ylim_a <- quantile(pm$pred2_sc, probs = c(0.005, 0.995), na.rm = TRUE)
  u1_lbl <- if (b1$use_log) sprintf("log(%s)", b1$det_unit) else b1$det_unit
  u2_lbl <- if (b2$use_log) sprintf("log(%s)", b2$det_unit) else b2$det_unit
  pt_cols <- adjustcolor(q_cols[pm$quadrant[idx_pt]], alpha.f = 0.35)

  plot(pm$pred1_sc[idx_pt], pm$pred2_sc[idx_pt],
       xlim = xlim_a, ylim = ylim_a, pch = 16, cex = 0.48, col = pt_cols,
       xlab = sprintf("Run 1: %s [%s]", b1$label, u1_lbl),
       ylab = sprintf("Run 2: %s [%s]", b2$label, u2_lbl),
       main = sprintf("A. Reach-by-Reach Relationship & Co-Occurrence Quadrants (n = %s)",
                      format(nrow(pm), big.mark = ",")),
       cex.main = 1.10, cex.lab = 0.95, las = 1)
  grid(col = "#e0e0e0", lty = "dotted")
  abline(v = rel_info$med1_sc, col = "#333333", lty = 2, lwd = 1.4)
  abline(h = rel_info$med2_sc, col = "#333333", lty = 2, lwd = 1.4)
  if (!is.na(rel_info$ols_beta)) {
    abline(rel_info$ols_alpha, rel_info$ols_beta, col = "#111111", lwd = 2.4)
  }
  legend("topleft",
         legend = c(
           sprintf("Q1 High/High: %.1f%% (n = %s)", rel_info$pct_q1, format(rel_info$n_q1, big.mark = ",")),
           sprintf("Q2 Low/Low:   %.1f%% (n = %s)", rel_info$pct_q2, format(rel_info$n_q2, big.mark = ",")),
           sprintf("Q3 High1/Low2: %.1f%% (n = %s)", rel_info$pct_q3, format(rel_info$n_q3, big.mark = ",")),
           sprintf("Q4 Low1/High2: %.1f%% (n = %s)", rel_info$pct_q4, format(rel_info$n_q4, big.mark = ",")),
           sprintf("OLS Slope = %+.3f (R2 = %.3f)", rel_info$ols_beta, rel_info$ols_r2)
         ),
         col = c(q_cols, "#111111"), pch = c(16, 16, 16, 16, NA), lwd = c(NA, NA, NA, NA, 2.4),
         bty = "n", cex = 0.78)
  mtext(sprintf("r = %.4f (Spearman rho = %.4f) | Concordance = %.1f%% (OR = %.2f, ChiSq p = %s) | P75 Jaccard = %.3f",
                rel_info$r_sc, rel_info$rho_sc, rel_info$pct_concordant, rel_info$odds_ratio,
                fmt_pval(rel_info$chisq_p), rel_info$jaccard_75),
        side = 3, line = 0.25, cex = 0.72, col = "#333333")

  # Panel B: Stratified Comparison Across Landscape Regimes & Strahler Stream Orders
  target_strata <- c(
    "Upland (>= 150 m)",
    "Lowland Non-Aquifer (< 150 m)",
    "Lowland Aquifer (< 150 m, Aq = 1)",
    "Headwaters (Order 1-2)",
    "Medium Tributaries (Order 3-4)",
    "Mainstem Rivers (Order >= 5)"
  )
  short_strata <- c("Upland\n(>=150m)", "Lowland\nNon-Aq", "Lowland\nAquifer", "Order 1-2\n(Headw.)", "Order 3-4\n(Medium)", "Order >=5\n(Mainstem)")
  z_mat   <- matrix(0, nrow = K, ncol = length(target_strata), dimnames = list(run_labels, short_strata))
  med_mat <- matrix(NA_real_, nrow = K, ncol = length(target_strata))

  for (i in seq_len(K)) {
    for (s_i in seq_along(target_strata)) {
      sub_s <- subset(strat_df, Run_Index == i & Stratum == target_strata[s_i] &
                                Dimension %in% c("3. Landscape Regime", "4. Stream Order"))
      if (nrow(sub_s) > 0) {
        z_mat[i, s_i]   <- sub_s$Pred_Mean_Z[1]
        med_mat[i, s_i] <- sub_s$Pred_Median_Raw[1]
      }
    }
  }
  yr_b  <- range(c(-0.85, 0.85, z_mat), na.rm = TRUE)
  pad_b <- 0.26 * diff(yr_b)
  bp_b  <- matrix(barplot(z_mat, beside = TRUE, col = run_cols, border = "#222222",
                          ylim = c(yr_b[1] - pad_b, yr_b[2] + pad_b),
                          ylab = "Mean Standardized Prediction Z-Score ((y - mean) / SD)",
                          main = "B. Stratified Output Profile by Landscape Regime & Stream Order",
                          cex.main = 1.10, cex.lab = 0.95, cex.names = 0.84, las = 1),
                  nrow = K)
  grid(nx = NA, ny = NULL, col = "#e0e0e0", lty = "dotted")
  abline(h = 0, col = "#333333", lwd = 1)
  abline(v = mean(c(bp_b[K, 3], bp_b[1, 4])), col = "#7f8c8d", lwd = 1.5, lty = 2)
  for (i in seq_len(K)) {
    for (s_i in seq_along(target_strata)) {
      zv <- z_mat[i, s_i]
      mv <- med_mat[i, s_i]
      if (!is.na(mv)) {
        y_pos <- if (zv >= 0) zv + 0.04 * diff(yr_b) else zv - 0.055 * diff(yr_b)
        text(bp_b[i, s_i], y_pos, labels = sprintf("%.2g", mv), cex = 0.68, font = 2, col = "#111111")
      }
    }
  }
  legend("topleft", legend = c(run_labels, "(Numbers on bars = Stratum Raw Median)"),
         fill = c(run_cols, NA), border = c(rep("#222222", K), NA), bty = "n", cex = 0.78)

  # Panel C: Spatial Map of 4-Quadrant Co-Occurrence Across England
  par(mar = c(1.5, 1.5, 3.8, 1.5))
  xr <- range(c(pm$start_X, pm$end_X), na.rm = TRUE)
  yr <- range(c(pm$start_Y, pm$end_Y), na.rm = TRUE)
  pad_x <- diff(xr) * 0.24
  plot(NULL, xlim = c(xr[1] - diff(xr) * 0.03, xr[2] + pad_x), ylim = yr,
       asp = 1, axes = FALSE, xlab = "", ylab = "",
       main = "C. Spatial Co-Occurrence Map Across Stream Reaches (4 Quadrants)",
       cex.main = 1.08)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#f8f9fa", border = "#dcdcdc")
  cols_c <- q_cols[pm$quadrant]
  lwds_c <- pmin(2.1, 0.32 + 0.34 * log10(pmax(1, pm$Shreve)))
  segments(pm$start_X, pm$start_Y, pm$mid_X, pm$mid_Y, col = cols_c, lwd = lwds_c)
  segments(pm$mid_X,   pm$mid_Y,   pm$end_X, pm$end_Y, col = cols_c, lwd = lwds_c)
  legend(xr[2] - diff(xr) * 0.02, yr[1] + diff(yr) * 0.68,
         title = "Reach Co-Occurrence",
         legend = c(
           sprintf("High R1 & High R2 (%.1f%%)", rel_info$pct_q1),
           sprintf("Low R1 & Low R2 (%.1f%%)",   rel_info$pct_q2),
           sprintf("High R1 & Low R2 (%.1f%%)",  rel_info$pct_q3),
           sprintf("Low R1 & High R2 (%.1f%%)",  rel_info$pct_q4)
         ),
         fill = q_cols, border = "#222222", bty = "n", cex = 0.78)
  mtext(sprintf("Concordant Reaches (High/High + Low/Low): %.1f%% | Top-Quartile Hotspot Overlap (Jaccard): %.3f",
                rel_info$pct_concordant, rel_info$jaccard_75),
        side = 3, line = 0.2, cex = 0.75, col = "#333333")

  # Panel D: Basin & Regime Distribution of Stoichiometric N:P Ratio or Z-Score Divergence
  par(mar = c(5.2, 4.8, 3.8, 1.8))
  grp_masks <- list(
    "All Eng"     = rep(TRUE, nrow(pm)),
    "Humber"      = if (!is.na(b2$hum_cid)) (pm$CatchID_2 == b2$hum_cid) else rep(FALSE, nrow(pm)),
    "Thames"      = if (!is.na(b2$tha_cid)) (pm$CatchID_2 == b2$tha_cid) else rep(FALSE, nrow(pm)),
    "Severn"      = if (!is.na(b2$sev_cid)) (pm$CatchID_2 == b2$sev_cid) else rep(FALSE, nrow(pm)),
    "Upland"      = !is.na(pm$Elevation) & pm$Elevation >= 150.0,
    "Low Non-Aq"  = !is.na(pm$Elevation) & pm$Elevation < 150.0 & pm$Is_Aquifer == 0,
    "Low Aquifer" = !is.na(pm$Elevation) & pm$Elevation < 150.0 & pm$Is_Aquifer == 1
  )
  grp_masks <- grp_masks[sapply(grp_masks, function(m) sum(m, na.rm = TRUE) >= 10)]
  pal_d <- c("#34495e", "#2980b9", "#d35400", "#27ae60", "#8e44ad", "#16a085", "#c0392b")[seq_along(grp_masks)]

  if (isTRUE(rel_info$is_np_pair)) {
    box_d <- lapply(grp_masks, function(m) log10(pmax(0.5, pm$np_molar_ratio[m])))
    bp_tmp <- boxplot(box_d, plot = FALSE)
    ylim_d <- c(min(bp_tmp$stats[1, ], log10(8), na.rm = TRUE) - 0.12,
                max(bp_tmp$stats[5, ], log10(50), na.rm = TRUE) + 0.22)
    bp_d <- boxplot(box_d, col = adjustcolor(pal_d, alpha.f = 0.65), border = "#222222",
                    outline = FALSE, yaxt = "n", las = 2, ylim = ylim_d,
                    ylab = "Atomic / Molar N:P Ratio (mol N : mol P, log scale)",
                    main = "D. Reach Stoichiometric N:P Molar Ratio Across Basins & Regimes",
                    cex.main = 1.08, cex.lab = 0.95, cex.axis = 0.85)
    grid(nx = NA, ny = NULL, col = "#e0e0e0", lty = "dotted")
    np_ticks <- c(2, 5, 10, 16, 20, 50, 100, 200, 500, 1000)
    axis(2, at = log10(np_ticks), labels = paste0(np_ticks, ":1"), las = 1, cex.axis = 0.82)
    abline(h = log10(16), col = "#c0392b", lwd = 1.8, lty = 2)
    abline(h = log10(20), col = "#2980b9", lwd = 1.4, lty = 3)
    abline(h = log10(10), col = "#27ae60", lwd = 1.4, lty = 3)
    meds_np <- sapply(grp_masks, function(m) median(pm$np_molar_ratio[m], na.rm = TRUE))
    text(seq_along(grp_masks), bp_d$stats[5, ], labels = sprintf("%.0f:1", meds_np), pos = 3, cex = 0.76, font = 2)
    mtext(sprintf("All England: %.1f%% P-Limited (>20:1), %.1f%% Co-Limited (10-20:1), %.1f%% N-Limited (<10:1) | Redfield = 16:1",
                  rel_info$np_info$pct_p_limited, rel_info$np_info$pct_co_limited, rel_info$np_info$pct_n_limited),
          side = 3, line = 0.25, cex = 0.75, col = "#333333")
  } else {
    box_d <- lapply(grp_masks, function(m) pm$z_diff[m])
    bp_tmp <- boxplot(box_d, plot = FALSE)
    ylim_d <- c(min(bp_tmp$stats[1, ], -0.5, na.rm = TRUE) - 0.2,
                max(bp_tmp$stats[5, ],  0.5, na.rm = TRUE) + 0.35)
    bp_d <- boxplot(box_d, col = adjustcolor(pal_d, alpha.f = 0.65), border = "#222222",
                    outline = FALSE, las = 2, ylim = ylim_d,
                    ylab = "Standardized Difference: Z(Run 2) minus Z(Run 1)",
                    main = "D. Relative Dominance (Z2 - Z1) Across Basins & Landscape Regimes",
                    cex.main = 1.08, cex.lab = 0.95, cex.axis = 0.85)
    grid(nx = NA, ny = NULL, col = "#e0e0e0", lty = "dotted")
    abline(h = 0, col = "#c0392b", lwd = 1.8, lty = 2)
    mns_z <- sapply(box_d, mean, na.rm = TRUE)
    points(seq_along(grp_masks), mns_z, pch = 18, cex = 1.3, col = "#111111")
    text(seq_along(grp_masks), bp_d$stats[5, ], labels = sprintf("%+.2f", mns_z), pos = 3, cex = 0.76, font = 2)
    mtext(sprintf("Positive (>0) = Relatively Higher in %s | Negative (<0) = Relatively Higher in %s",
                  b2$label, b1$label),
          side = 3, line = 0.25, cex = 0.75, col = "#333333")
  }

  mtext(sprintf("S3N Output Descriptive & Reach-Level Relationship Analysis: %s vs. %s",
                b1$label, b2$label),
        side = 3, line = 0.4, outer = TRUE, cex = 1.15, font = 2, col = "#111111")
}

# ==============================================================================
# Markdown Report Generator
# ==============================================================================
write_comparison_markdown_report <- function(bundles, cv_comp_df, pred_sum_df, desc_df, rep_df,
                                             strat_df, stat_tests_df, rel_info,
                                             param_comp_df, pair12_cv_match, pair12_pred_match,
                                             out_dir, report_path) {
  K <- length(bundles)
  same_det_all <- (length(unique(sapply(bundles, function(b) b$det_code))) == 1L)
  lines <- character(0)
  add_line <- function(fmt = "", ...) {
    if (length(list(...)) == 0L) {
      lines <<- c(lines, fmt)
    } else {
      lines <<- c(lines, sprintf(fmt, ...))
    }
  }

  add_line("# S3N Model Run & Output Comparison Report")
  add_line("")
  if (same_det_all) {
    add_line("- **Comparison Mode**: Same Determinand (`%s` — %s, `%s`, %s scale)",
             bundles[[1]]$det_code, bundles[[1]]$det_name, bundles[[1]]$det_unit,
             if (bundles[[1]]$use_log) "log" else "raw")
  } else {
    det_desc <- paste(sapply(seq_len(K), function(i) {
      sprintf("Run %d: `%s` (%s, `%s`, %s scale)", i, bundles[[i]]$det_code,
              bundles[[i]]$det_name, bundles[[i]]$det_unit,
              if (bundles[[i]]$use_log) "log" else "raw")
    }), collapse = " vs. ")
    add_line("- **Comparison Mode**: Cross-Determinand Comparison (%s)", det_desc)
  }
  add_line("- **Number of Model Runs Compared**: `%d`", K)
  add_line("- **Output Directory**: `%s`", sub(paste0("^", BASE_DIR, "/?"), "", out_dir))
  add_line("")
  add_line("## 1. Overview of Compared Model Runs")
  add_line("")
  add_line("| Run | Run Label | Determinand | Unit (Scale) | Water Network | Period | Covariates | Reaches | Catchments | Observed Sites |")
  add_line("| :---: | :--- | :--- | :---: | :--- | :---: | :--- | ---: | ---: | ---: |")
  for (i in seq_len(K)) {
    b <- bundles[[i]]
    add_line("| **Run %d** | **%s** | %s (`%s`) | `%s` (%s) | %s | `%s` | %s | %s | %d | %d |",
             i, b$label, b$det_name, b$det_code, b$det_unit,
             if (b$use_log) "log" else "raw",
             b$net_label, b$period, b$vars_label,
             format(nrow(b$eng_pred), big.mark = ","),
             length(unique(b$eng_pred$CatchID)),
             sum(b$eng_pred$has_obs == 1))
  }
  add_line("")

  # ----------------------------------------------------------------------------
  # Executive Summary & Winner Determination
  # ----------------------------------------------------------------------------
  add_line("## 2. Executive Summary: Model Performance & Output Key Findings")
  add_line("")

  eng_s3n <- subset(cv_comp_df, Canon_Domain == "All England RSN Network" & Canon_Model == "S3N Combined")
  if (nrow(eng_s3n) >= 2) {
    best_r2_idx   <- which.max(eng_s3n$R2)
    best_rmse_idx <- which.min(eng_s3n$RMSE)
    best_se_idx   <- which.min(pred_sum_df$SE_Mean)

    w_r2   <- eng_s3n[best_r2_idx, ]
    w_rmse <- eng_s3n[best_rmse_idx, ]

    add_line("### 2.1 Overall Predictive Performance (All England 5-Fold Cross-Validation)")
    if (same_det_all && w_r2$Run_Index == w_rmse$Run_Index) {
      add_line("- **Overall Winner**: **%s (Run %d)** achieved both the highest out-of-sample 5-fold cross-validation **R² = %.4f** and the lowest **RMSE = %.4f** (MAE = %.4f) across the All England stream network.",
               w_r2$Run_Label, w_r2$Run_Index, w_r2$R2, w_r2$RMSE, w_r2$MAE)
    } else {
      add_line("- **Highest Variance Explained (CV R²)**: **%s (Run %d)** achieved the highest out-of-sample 5-fold CV **R² = %.4f** (RMSE = %.4f, MAE = %.4f, ΔR² over non-spatial OLS = `%+.4f`).",
               w_r2$Run_Label, w_r2$Run_Index, w_r2$R2, w_r2$RMSE, w_r2$MAE, w_r2$Delta_R2_vs_OLS)
      if (same_det_all) {
        add_line("- **Lowest Absolute Prediction Error (CV RMSE)**: **%s (Run %d)** achieved the lowest out-of-sample 5-fold CV **RMSE = %.4f** (MAE = %.4f, R² = %.4f) — note that when comparing different water networks, `RMSE` and `R²` can rank differently if one network snaps a slightly different subset of headwater vs. mainstem stations (`n = %d` vs `n = %d` observed reaches).",
                 w_rmse$Run_Label, w_rmse$Run_Index, w_rmse$RMSE, w_rmse$MAE, w_rmse$R2,
                 w_rmse$n_obs, w_r2$n_obs)
      } else {
        for (r_idx in seq_len(nrow(eng_s3n))) {
          er <- eng_s3n[r_idx, ]
          add_line("- **Run %d (%s)**: 5-Fold CV **R² = %.4f**, **RMSE = %.4f**, **MAE = %.4f**, Spatial Gain over OLS **ΔR² = %+.4f** (`%+.1f%%` RMSE reduction).",
                   er$Run_Index, er$Run_Label, er$R2, er$RMSE, er$MAE, er$Delta_R2_vs_OLS, er$Pct_RMSE_Red_OLS)
        }
      }
    }

    # Co-located station head-to-head result (when comparing same determinand)
    if (same_det_all && !is.null(pair12_cv_match) && nrow(pair12_cv_match) >= 10) {
      mse1 <- mean(pair12_cv_match$sq_e1, na.rm = TRUE)
      mse2 <- mean(pair12_cv_match$sq_e2, na.rm = TRUE)
      mae1 <- mean(pair12_cv_match$abs_e1, na.rm = TRUE)
      mae2 <- mean(pair12_cv_match$abs_e2, na.rm = TRUE)
      r2_m1 <- 1 - sum(pair12_cv_match$sq_e1) / sum((pair12_cv_match$obs_1 - mean(pair12_cv_match$obs_1))^2)
      r2_m2 <- 1 - sum(pair12_cv_match$sq_e2) / sum((pair12_cv_match$obs_2 - mean(pair12_cv_match$obs_2))^2)
      dm_test <- safe_paired_ttest(pair12_cv_match$sq_e2, pair12_cv_match$sq_e1)
      wx_test <- safe_paired_wilcox(pair12_cv_match$abs_e2, pair12_cv_match$abs_e1)
      coloc_winner <- if (mse1 <= mse2) bundles[[1]]$label else bundles[[2]]$label
      add_line("- **Head-to-Head at Identical Co-located Monitoring Stations (`n = %d` matched sites)**:", nrow(pair12_cv_match))
      add_line("  - **%s**: Matched CV R² = `%.4f`, Matched RMSE = `%.4f`, Matched MAE = `%.4f`",
               bundles[[1]]$label, r2_m1, sqrt(mse1), mae1)
      add_line("  - **%s**: Matched CV R² = `%.4f`, Matched RMSE = `%.4f`, Matched MAE = `%.4f`",
               bundles[[2]]$label, r2_m2, sqrt(mse2), mae2)
      add_line("  - **Matched-Site Winner**: **%s** (`Paired t-test on squared error p = %s [%s]`, `Paired Wilcoxon test on |error| p = %s [%s]`).",
               coloc_winner, fmt_pval(dm_test$p.value), sig_stars(dm_test$p.value),
               fmt_pval(wx_test$p.value), sig_stars(wx_test$p.value))
    }

    add_line("- **Lowest Network-Wide Kriging Uncertainty**: **%s (Run %d)** achieved the lowest mean kriging standard error (`SE = %.4f`, median `%.4f`).",
             pred_sum_df$Run_Label[best_se_idx], best_se_idx,
             pred_sum_df$SE_Mean[best_se_idx], pred_sum_df$SE_Median[best_se_idx])

    if (!is.null(rel_info)) {
      add_line("- **Reach-by-Reach Relationship (`n = %s` matched reaches)**:",
               format(rel_info$n_matched, big.mark = ","))
      add_line("  - Model-scale Pearson correlation **r = %.4f** (`p = %s`), Spearman rank correlation **ρ = %.4f** (`p = %s`), OLS scaling slope **β = %+.4f ± %.4f** (`R² = %.4f`).",
               rel_info$r_sc, fmt_pval(rel_info$p_sc), rel_info$rho_sc, fmt_pval(rel_info$p_rho_sc),
               rel_info$ols_beta, rel_info$ols_se, rel_info$ols_r2)
      add_line("  - **4-Quadrant Co-Occurrence Concordance**: **%.1f%%** of reaches are concordant (`%.1f%%` Co-High `Q1` + `%.1f%%` Co-Low `Q2` vs. `%.1f%%` discordant; `Odds Ratio = %.2f`, `Chi-squared p = %s`), with a top-quartile (`> P75`) hotspot Jaccard overlap of **`%.3f`** (`%.1f%%` conditional hotspot probability).",
               rel_info$pct_concordant, rel_info$pct_q1, rel_info$pct_q2,
               rel_info$pct_q3 + rel_info$pct_q4, rel_info$odds_ratio, fmt_pval(rel_info$chisq_p),
               rel_info$jaccard_75, rel_info$cond_75)
      if (isTRUE(rel_info$is_np_pair)) {
        add_line("  - **Nutrient Stoichiometry (`N:P` Ratio Across Reaches)**: Network median mass ratio is **`%.1f` mg N / mg P** and median atomic/molar ratio is **`%.1f : 1` mol N / mol P** (IQR `%.1f`–`%.1f`). Across England, **`%.1f%%`** of reaches are **P-limited (`Molar N:P > 20:1`)**, **`%.1f%%`** are **Co-limited (`10:1–20:1`)**, and **`%.1f%%`** are **N-limited (`< 10:1`)**.",
                 rel_info$np_info$mass_summary$Median, rel_info$np_info$molar_summary$Median,
                 rel_info$np_info$molar_summary$Q25, rel_info$np_info$molar_summary$Q75,
                 rel_info$np_info$pct_p_limited, rel_info$np_info$pct_co_limited, rel_info$np_info$pct_n_limited)
      }
    }
    add_line("")

    # Basin-by-basin winner summary
    add_line("### 2.2 Winner by Modelling Domain (S3N Combined)")
    add_line("")
    add_line("| Domain | Best Model Run (Highest CV R²) | Best CV R² | Best Model Run (Lowest CV RMSE) | Best CV RMSE | S3N Gain Over OLS (ΔR²) |")
    add_line("| :--- | :--- | ---: | :--- | ---: | ---: |")
    for (dom in unique(cv_comp_df$Canon_Domain)) {
      sub_d <- subset(cv_comp_df, Canon_Domain == dom & Canon_Model == "S3N Combined")
      if (nrow(sub_d) > 0) {
        ir2   <- which.max(sub_d$R2)
        irmse <- which.min(sub_d$RMSE)
        add_line("| **%s** | %s | `%.4f` | %s | `%.4f` | `%+.4f` |",
                 dom, sub_d$Run_Label[ir2], sub_d$R2[ir2],
                 sub_d$Run_Label[irmse], sub_d$RMSE[irmse], sub_d$Delta_R2_vs_OLS[ir2])
      }
    }
    add_line("")
  }

  # ----------------------------------------------------------------------------
  # 3. Full Cross-Validation Comparison Table
  # ----------------------------------------------------------------------------
  add_line("## 3. Cross-Validation Metrics Across All Domains & Models")
  add_line("")
  add_line("| Domain | Model | Run Label | n_obs | n_reaches | CV R² | CV RMSE | CV MAE | ΔR² vs OLS | % RMSE Red. vs OLS |")
  add_line("| :--- | :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
  for (r_i in seq_len(nrow(cv_comp_df))) {
    r <- cv_comp_df[r_i, ]
    add_line("| %s | %s | **%s** | %d | %s | `%.4f` | `%.4f` | `%.4f` | `%+.4f` | `%+.2f%%` |",
             r$Canon_Domain, r$Canon_Model, r$Run_Label, r$n_obs,
             format(r$n_reaches, big.mark = ","), r$R2, r$RMSE, r$MAE,
             r$Delta_R2_vs_OLS, r$Pct_RMSE_Red_OLS)
  }
  add_line("")

  # ----------------------------------------------------------------------------
  # 4. Statistical Hypothesis Tests Table
  # ----------------------------------------------------------------------------
  add_line("## 4. Statistical Hypothesis Tests")
  add_line("")
  add_line("Significance codes: `***` p < 0.001, `**` p < 0.01, `*` p < 0.05, `.` p < 0.10, `ns` not significant (p ≥ 0.10).")
  add_line("")
  add_line("| Comparison | Category | Statistical Test | Metric | Run 1 Val | Run 2 Val | Diff / R² | Stat | p-value | Sig | Finding |")
  add_line("| :--- | :--- | :--- | :--- | ---: | ---: | ---: | :--- | ---: | :---: | :--- |")
  for (r_i in seq_len(nrow(stat_tests_df))) {
    s <- stat_tests_df[r_i, ]
    add_line("| %s | %s | %s | %s | `%.4f` | `%.4f` | `%+.4f` | `%s = %.3f` | `%s` | `%s` | %s |",
             s$Comparison, s$Category, s$Test_Name, s$Metric,
             s$Value_Run1, s$Value_Run2, s$Diff_Run2_Run1,
             s$Statistic_Name, s$Statistic_Val, fmt_pval(s$P_Value),
             s$Significance, s$Interpretation)
  }
  add_line("")

  # ----------------------------------------------------------------------------
  # 5. Comprehensive Descriptive Statistics & Stratified Comparison of Outputs
  # ----------------------------------------------------------------------------
  add_line("## 5. Descriptive Statistics & Stratified Comparison of Model Output Results")
  add_line("")
  add_line("### 5.1 Full Network & Monitoring Station Descriptive Statistics")
  add_line("")
  add_line("| Run Label | Variable Group | Unit | N | Mean ± SD | CV (%) | Geom Mean | Skew | Kurt | Min | P05 | Q25 | Median | Q75 | P95 | Max | IQR |")
  add_line("| :--- | :--- | :---: | ---: | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
  for (r_i in seq_len(nrow(desc_df))) {
    d <- desc_df[r_i, ]
    gm_str <- if (is.na(d$Geom_Mean)) "—" else sprintf("%.4g", d$Geom_Mean)
    cv_str <- if (is.na(d$CV_Pct))    "—" else sprintf("%.1f%%", d$CV_Pct)
    add_line("| **%s** | %s | `%s` | %s | `%.4g ± %.4g` | `%s` | `%s` | `%+.2f` | `%+.2f` | `%.4g` | `%.4g` | `%.4g` | **`%.4g`** | `%.4g` | `%.4g` | `%.4g` | `%.4g` |",
             d$Run_Label, d$Variable_Group, d$Unit, format(d$N, big.mark = ","),
             d$Mean, d$SD, cv_str, gm_str, d$Skewness, d$Excess_Kurtosis,
             d$Min, d$P05, d$Q25, d$Median, d$Q75, d$P95, d$Max, d$IQR)
  }
  add_line("")

  if (!is.null(rep_df) && nrow(rep_df) > 0) {
    add_line("### 5.2 Monitoring Network Representativeness (Observed Stations vs. Full River Network)")
    add_line("")
    add_line("| Run Label | Monitored Sites | Network Reaches | Obs Median (Mean) | Network Pred Median (Mean) | Pred / Obs Median Ratio | Mann–Whitney p | KS D-stat (p) | Interpretation |")
    add_line("| :--- | ---: | ---: | :--- | :--- | ---: | ---: | :--- | :--- |")
    for (r_i in seq_len(nrow(rep_df))) {
      rp <- rep_df[r_i, ]
      add_line("| **%s** | %d | %s | `%.4g` (`%.4g`) | `%.4g` (`%.4g`) | **`%.3fx`** | `%s` (`%s`) | `%.3f` (`%s`) | %s |",
               rp$Run_Label, rp$n_Obs_Stations, format(rp$n_Total_Reaches, big.mark = ","),
               rp$Obs_Median_Raw, rp$Obs_Mean_Raw, rp$Pred_Median_Raw, rp$Pred_Mean_Raw,
               rp$Pred_to_Obs_Med_Ratio, fmt_pval(rp$MW_P_Value), rp$Significance,
               rp$KS_D_Stat, fmt_pval(rp$KS_P_Value), rp$Interpretation)
    }
    add_line("")
  }

  if (!is.null(strat_df) && nrow(strat_df) > 0) {
    add_line("### 5.3 Stratified Model Output Comparison (Topography, Hydrogeology, Stream Order & Major Basins)")
    add_line("")
    add_line("| Dimension | Stratum | Run Label | Reaches (% Net) | Obs Sites (Med) | Pred Mean ± SD | Pred Median (IQR) | Pred 5th–95th Pct | Mean Z-Score | Fold Ratio vs Ref | Within-Run Test (p-val) |")
    add_line("| :--- | :--- | :--- | ---: | ---: | :--- | :--- | :--- | ---: | ---: | :--- |")
    for (r_i in seq_len(nrow(strat_df))) {
      st <- strat_df[r_i, ]
      obs_s <- if (is.na(st$Obs_Median_Raw)) sprintf("%d (—)", st$n_obs) else sprintf("%d (`%.3g`)", st$n_obs, st$Obs_Median_Raw)
      rat_s <- if (is.na(st$Ratio_vs_Ref)) "—" else sprintf("%.2fx", st$Ratio_vs_Ref)
      tst_s <- if (is.na(st$Within_Run_P)) st$Within_Run_Test else sprintf("%s: p = %s (%s)", st$Within_Run_Test, fmt_pval(st$Within_Run_P), st$Within_Run_Sig)
      add_line("| %s | **%s** | **%s** | %s (`%.1f%%`) | %s | `%.4g ± %.4g` | **`%.4g`** (`%.4g`) | `%.4g` – `%.4g` | `%+.2f` | `%s` | %s |",
               st$Dimension, st$Stratum, st$Run_Label, format(st$n_reaches, big.mark = ","), st$Pct_Network,
               obs_s, st$Pred_Mean_Raw, st$Pred_SD_Raw, st$Pred_Median_Raw, st$Pred_IQR_Raw,
               st$Pred_P05_Raw, st$Pred_P95_Raw, st$Pred_Mean_Z, rat_s, tst_s)
    }
    add_line("")
  }

  # ----------------------------------------------------------------------------
  # 6. Reach-by-Reach Relationship & Cross-Determinand Co-Occurrence Analysis
  # ----------------------------------------------------------------------------
  if (!is.null(rel_info)) {
    add_line("## 6. Reach-by-Reach Relationship, Co-Occurrence & Stoichiometric Analysis")
    add_line("")
    add_line("Analysis comparing **Run 1 (%s)** and **Run 2 (%s)** across **`%s` spatially matched stream reaches** (`fig_comp5_reach_relationships.png`).",
             bundles[[1]]$label, bundles[[2]]$label, format(rel_info$n_matched, big.mark = ","))
    add_line("")

    add_line("### 6.1 Bivariate Correlation & Scaling Summary (Monitored Stations vs. Full Predicted Network)")
    add_line("")
    add_line("| Metric / Relationship Test | Value | Statistical Test / Precision | Interpretation |")
    add_line("| :--- | :--- | :--- | :--- |")
    add_line("| **Model-Scale Pearson Correlation (`r`)** | **`%.4f`** (`R² = %.4f`) | `p = %s` (`%s`) | Linear association on model scale across `%s` reaches |",
             rel_info$r_sc, rel_info$r_sc^2, fmt_pval(rel_info$p_sc), sig_stars(rel_info$p_sc), format(rel_info$n_matched, big.mark = ","))
    add_line("| **Rank-Order Spearman Correlation (`ρ`)** | **`%.4f`** | `p = %s` (`%s`) | Monotonic rank concordance across `%s` reaches |",
             rel_info$rho_sc, fmt_pval(rel_info$p_rho_sc), sig_stars(rel_info$p_rho_sc), format(rel_info$n_matched, big.mark = ","))
    add_line("| **Raw-Scale Pearson / Spearman Correlation** | `r = %.4f` / `ρ = %.4f` | `p = %s` | Correlation in untransformed units (`%s` vs `%s`) |",
             rel_info$r_raw, rel_info$rho_raw, fmt_pval(rel_info$p_raw), bundles[[1]]$det_unit, bundles[[2]]$det_unit)
    add_line("| **S3N Component Correlations (Tail-Up / Tail-Down)** | `r_TU = %.4f` / `r_TD = %.4f` | — | Flow-connected (`TU`) vs. hydrological distance (`TD`) concordance |",
             rel_info$r_tu_sc, rel_info$r_td_sc)
    add_line("| **OLS Scaling Regression (`Pred_2 = α + β * Pred_1`)** | `β = %+.4f ± %.4f`, `α = %+.4f` | `R² = %.4f`, `p = %s` | %s |",
             rel_info$ols_beta, rel_info$ols_se, rel_info$ols_alpha, rel_info$ols_r2, fmt_pval(rel_info$ols_p),
             if (bundles[[1]]$use_log && bundles[[2]]$use_log) sprintf("Power-law scaling exponent: Det2 scales as Det1^(%.3f)", rel_info$ols_beta) else "Linear scaling slope of Run 2 per unit Run 1")
    if (rel_info$n_both_obs >= 5) {
      add_line("| **Co-located Monitored Stations (`n = %d` sites)** | `r_obs = %.4f` (`ρ_obs = %.4f`) | `OLS β_obs = %+.4f`, `p = %s` | Empirical station correlation vs. S3N network-wide `r = %.4f` |",
               rel_info$n_both_obs, rel_info$r_obs_sc, rel_info$rho_obs_raw, rel_info$ols_obs_beta, fmt_pval(rel_info$p_obs_sc), rel_info$r_sc)
    }
    add_line("")

    add_line("### 6.2 Four-Quadrant Reach Co-Occurrence & Hotspot Concordance")
    add_line("")
    add_line("Reaches are classified relative to each run's network median (`Med_1 = %.4g %s`, `Med_2 = %.4g %s`) and upper-quartile (`> P75`) / top-decile (`> P90`) hotspot thresholds:",
             rel_info$med1_raw, bundles[[1]]$det_unit, rel_info$med2_raw, bundles[[2]]$det_unit)
    add_line("")
    add_line("| Co-Occurrence Quadrant / Hotspot Metric | Reach Count | % of Network | Statistical Concordance & Interpretation |")
    add_line("| :--- | ---: | ---: | :--- |")
    add_line("| **Q1: High Run 1 & High Run 2 (`> Med_1` & `> Med_2`)** | `%s` | **`%.2f%%`** | Co-elevated / Dual-pressure reaches (`> %.3g %s` and `> %.3g %s`) |",
             format(rel_info$n_q1, big.mark = ","), rel_info$pct_q1, rel_info$med1_raw, bundles[[1]]$det_unit, rel_info$med2_raw, bundles[[2]]$det_unit)
    add_line("| **Q2: Low Run 1 & Low Run 2 (`<= Med_1` & `<= Med_2`)** | `%s` | **`%.2f%%`** | Co-low / Reference reaches (`<= %.3g %s` and `<= %.3g %s`) |",
             format(rel_info$n_q2, big.mark = ","), rel_info$pct_q2, rel_info$med1_raw, bundles[[1]]$det_unit, rel_info$med2_raw, bundles[[2]]$det_unit)
    add_line("| **Q3: High Run 1 & Low Run 2 (`> Med_1` & `<= Med_2`)** | `%s` | `%.2f%%` | Decoupled reaches: Elevated in %s but low in %s |",
             format(rel_info$n_q3, big.mark = ","), rel_info$pct_q3, bundles[[1]]$label, bundles[[2]]$label)
    add_line("| **Q4: Low Run 1 & High Run 2 (`<= Med_1` & `> Med_2`)** | `%s` | `%.2f%%` | Decoupled reaches: Low in %s but elevated in %s |",
             format(rel_info$n_q4, big.mark = ","), rel_info$pct_q4, bundles[[1]]$label, bundles[[2]]$label)
    add_line("| **Overall Concordance (`Q1 + Q2`)** | `%s` | **`%.2f%%`** | `ChiSq = %.1f` (`p = %s`), **Odds Ratio = `%.2f`** (vs. `50.0%%` under independence) |",
             format(rel_info$n_q1 + rel_info$n_q2, big.mark = ","), rel_info$pct_concordant, rel_info$chisq_stat, fmt_pval(rel_info$chisq_p), rel_info$odds_ratio)
    add_line("| **Upper-Quartile Dual Hotspots (`Both > P75`)** | `%s` | `%.2f%%` | **Jaccard Index = `%.3f`**, Conditional `P(R2 > P75 | R1 > P75) = %.1f%%` (vs. `25.0%%` random) |",
             format(rel_info$n_h75_both, big.mark = ","), 100.0 * rel_info$n_h75_both / rel_info$n_matched, rel_info$jaccard_75, rel_info$cond_75)
    add_line("| **Top-Decile Extreme Dual Hotspots (`Both > P90`)** | `%s` | `%.2f%%` | **Jaccard Index = `%.3f`**, Conditional `P(R2 > P90 | R1 > P90) = %.1f%%` (vs. `10.0%%` random) |",
             format(rel_info$n_h90_both, big.mark = ","), 100.0 * rel_info$n_h90_both / rel_info$n_matched, rel_info$jaccard_90, rel_info$cond_90)
    add_line("")

    if (isTRUE(rel_info$is_np_pair) && !is.null(rel_info$np_info)) {
      np <- rel_info$np_info
      ms <- np$mass_summary
      mo <- np$molar_summary
      add_line("### 6.3 Nutrient Stoichiometry (`N:P` Mass & Molar Ratios and Ecological Limitation Regimes)")
      add_line("")
      add_line("Comparing **%s** and **%s** across `%s` reaches (`Atomic/Molar N:P = 2.21135 * Mass N:P`):",
               np$n_label, np$p_label, format(rel_info$n_matched, big.mark = ","))
      add_line("")
      add_line("- **Mass `N:P` Ratio (`mg N / mg P`)**: Mean = `%.1f ± %.1f`, Geometric Mean = `%.1f`, Median = **`%.1f`** (IQR `%.1f`–`%.1f`, 5th–95th percentile `%.1f`–`%.1f`; Monitored station median = `%.1f`).",
               ms$Mean, ms$SD, ms$Geom_Mean, ms$Median, ms$Q25, ms$Q75, ms$P05, ms$P95, np$obs_med_mass)
      add_line("- **Atomic / Molar `N:P` Ratio (`mol N : mol P`)**: Mean = `%.1f ± %.1f`, Geometric Mean = `%.1f : 1`, Median = **`%.1f : 1`** (IQR `%.1f : 1` to `%.1f : 1`, 5th–95th percentile `%.1f : 1` to `%.1f : 1`; Monitored station median = `%.1f : 1`).",
               mo$Mean, mo$SD, mo$Geom_Mean, mo$Median, mo$Q25, mo$Q75, mo$P05, mo$P95, np$obs_med_molar)
      add_line("- **Ecological Nutrient Limitation Classification (relative to Redfield `16:1` / Guildford & Hecky thresholds)**:")
      add_line("  - **Phosphorus-Limited (`Molar N:P > 20:1`)**: **`%s` reaches (`%.2f%%`)** — excess nitrogen relative to bioavailable orthophosphate; primary production is P-controlled.",
               format(np$n_p_limited, big.mark = ","), np$pct_p_limited)
      add_line("  - **Co-Limited / Balanced (`Molar N:P 10:1 to 20:1`)**: **`%s` reaches (`%.2f%%`)** — near-Redfield stoichiometric balance.",
               format(np$n_co_limited, big.mark = ","), np$pct_co_limited)
      add_line("  - **Nitrogen-Limited (`Molar N:P < 10:1`)**: **`%s` reaches (`%.2f%%`)** — typically upland/peat headwaters with very low nitrate or reaches with elevated point-source phosphate relative to nitrate.",
               format(np$n_n_limited, big.mark = ","), np$pct_n_limited)
      add_line("")
    }

    if (!is.null(rel_info$stratum_df) && nrow(rel_info$stratum_df) > 0) {
      add_line("### 6.4 Reach-by-Reach Relationship Breakdown Across Basins, Landscape Regimes & Stream Orders")
      add_line("")
      if (isTRUE(rel_info$is_np_pair)) {
        add_line("| Group | Stratum | Reaches | Run 1 Med (IQR) | Run 2 Med (IQR) | Pearson r | Spearman ρ | OLS β | % High/High | % Low/Low | % Concordant | P75 Jaccard | Med Molar N:P | % P-Limited (>20:1) |")
        add_line("| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
        for (r_i in seq_len(nrow(rel_info$stratum_df))) {
          sr <- rel_info$stratum_df[r_i, ]
          add_line("| %s | **%s** | %s | `%.3g` (`%.3g`) | `%.3g` (`%.3g`) | `%.3f` | `%.3f` | `%+.3f` | `%.1f%%` | `%.1f%%` | **`%.1f%%`** | `%.3f` | **`%.1f:1`** | `%.1f%%` |",
                   sr$Stratum_Group, sr$Stratum, format(sr$n_reaches, big.mark = ","),
                   sr$Run1_Median_Raw, sr$Run1_IQR_Raw, sr$Run2_Median_Raw, sr$Run2_IQR_Raw,
                   sr$Pearson_r, sr$Spearman_rho, sr$OLS_Slope_Beta,
                   sr$Pct_Q1_HighHigh, sr$Pct_Q2_LowLow, sr$Pct_Concordant,
                   sr$Hotspot75_Jaccard, sr$Median_NP_Molar, sr$Pct_P_Limited_gt20)
        }
      } else {
        add_line("| Group | Stratum | Reaches | Run 1 Med (IQR) | Run 2 Med (IQR) | Pearson r | Spearman ρ | OLS β | % High/High | % Low/Low | % Concordant | P75 Jaccard | Mean Z2 - Z1 | Med Ratio (R1/R2) |")
        add_line("| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
        for (r_i in seq_len(nrow(rel_info$stratum_df))) {
          sr <- rel_info$stratum_df[r_i, ]
          add_line("| %s | **%s** | %s | `%.3g` (`%.3g`) | `%.3g` (`%.3g`) | `%.3f` | `%.3f` | `%+.3f` | `%.1f%%` | `%.1f%%` | **`%.1f%%`** | `%.3f` | `%+.2f` | `%.3g` |",
                   sr$Stratum_Group, sr$Stratum, format(sr$n_reaches, big.mark = ","),
                   sr$Run1_Median_Raw, sr$Run1_IQR_Raw, sr$Run2_Median_Raw, sr$Run2_IQR_Raw,
                   sr$Pearson_r, sr$Spearman_rho, sr$OLS_Slope_Beta,
                   sr$Pct_Q1_HighHigh, sr$Pct_Q2_LowLow, sr$Pct_Concordant,
                   sr$Hotspot75_Jaccard, sr$Mean_Z_Diff_2_vs_1, sr$Median_Ratio_1_to_2)
        }
      }
      add_line("")
    }
  }

  # ----------------------------------------------------------------------------
  # 7. Fixed-Effect & Covariance Parameter Comparison (All England)
  # ----------------------------------------------------------------------------
  if (!is.null(param_comp_df)) {
    add_line("## 7. Fixed-Effect & Spatial Covariance Parameter Comparison (All England)")
    add_line("")
    add_line("Two-sample Wald Z-tests compare each parameter against **Run 1** using bootstrap standard errors (`SE = (CI_Upper - CI_Lower) / 3.92`).")
    add_line("")
    add_line("| Model | Parameter | Run Label | Estimate | 95% Bootstrap CI | Diff vs Run 1 | Wald Z | p-value | Sig |")
    add_line("| :--- | :--- | :--- | ---: | :---: | ---: | ---: | ---: | :---: |")
    eng_p <- subset(param_comp_df, Canon_Domain == "All England RSN Network")
    for (r_i in seq_len(nrow(eng_p))) {
      p <- eng_p[r_i, ]
      z_str <- if (is.na(p$Wald_Z)) "—" else sprintf("%.2f", p$Wald_Z)
      p_str <- if (is.na(p$Wald_P)) "—" else fmt_pval(p$Wald_P)
      s_str <- if (is.na(p$Wald_P)) "—" else p$Significance
      add_line("| %s | `%s` | **%s** | `%.4f` | `[%.4f, %.4f]` | `%+.4f` | `%s` | `%s` | `%s` |",
               p$Canon_Model, p$Parameter, p$Run_Label, p$Estimate,
               p$CI_Lower, p$CI_Upper, p$Diff_vs_Run1, z_str, p_str, s_str)
    }
    add_line("")
  }

  add_line("## 8. Generated Comparison Figures & Tables")
  add_line("")
  add_line("### Figures (Pure Base R Graphics)")
  add_line("1. `fig_comp1_cv_performance.png` — 4-panel comparison of 5-fold CV R², CV RMSE, empirical CDF of station CV errors, and spatial gain over non-spatial OLS.")
  add_line("2. `fig_comp2_predictions_and_uncertainty.png` — 4-panel comparison of co-located reach predictions, kernel density distributions, kriging standard error boxplots, and a spatial difference / Z-divergence map.")
  add_line("3. `fig_comp3_parameter_forest.png` — Multi-run forest plot of fixed-effect regression coefficients (with Wald Z-test significance flags) and spatial variance partitioning (`sigma.sq` vs `tau.sq`).")
  add_line("4. `fig_comp4_side_by_side_maps.png` — Side-by-side All England stream network prediction maps.")
  add_line("5. `fig_comp5_reach_relationships.png` — 4-panel reach-by-reach relationship scatter (colored by 4 co-occurrence quadrants), stratified environmental/topological output profile, England 4-quadrant co-occurrence map, and stoichiometric `N:P` ratio / Z-divergence boxplots across basins and regimes.")
  add_line("")
  add_line("### Exported Summary CSV Tables")
  add_line("- `cv_metrics_comparison.csv` — Cross-validation metrics (`R2`, `RMSE`, `MAE`, `Delta_R2_vs_OLS`) across all runs, domains, and models.")
  add_line("- `descriptive_stats_summary.csv` — Full descriptive statistics (`Mean`, `SD`, `CV_Pct`, `Geom_Mean`, `Skewness`, `Excess_Kurtosis`, `Min`..`Max`, `IQR`, `MAD`) for observed stations, S3N predictions, and kriging uncertainty.")
  add_line("- `stratified_results_summary.csv` — Stratified descriptive statistics across topography (`Upland` vs `Lowland`), hydrogeology (`Aquifer` vs `Non-Aquifer`), Strahler stream orders (`1-2`, `3-4`, `>=5`), and major river basins.")
  add_line("- `reach_relationship_summary.csv` — Reach-by-reach bivariate correlations, 4-quadrant co-occurrence percentages, hotspot Jaccard overlap, and stoichiometric `N:P` / Z-difference metrics by basin and regime.")
  add_line("- `statistical_tests_summary.csv` — All formal hypothesis test statistics, p-values, and interpretations.")
  add_line("- `parameter_comparison.csv` — Fixed-effect and spatial covariance parameter estimates with two-sample Wald Z-tests.")
  add_line("")

  writeLines(lines, report_path)
}

# ==============================================================================
# Command-Line Interface (when invoked directly via Rscript)
# ==============================================================================
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
    cat("Usage:\n")
    cat("  Rscript scripts/04_compare_model_runs.R --list\n")
    cat("  Rscript scripts/04_compare_model_runs.R <dir1> <dir2> [dir3 ...] [--labels=\"L1,L2,...\"] [--out=<dir>]\n")
    cat("  Rscript scripts/04_compare_model_runs.R <parent_dir> --by=season|year|month [--out=<dir>]\n")
    cat("  Rscript scripts/04_compare_model_runs.R --dirs=output/det_0117,output/det_0117_euhydro_thames\n")
    quit(save = "no", status = 0)
  }

  if (any(args %in% c("--list", "--list-runs"))) {
    list_available_model_runs()
    quit(save = "no", status = 0)
  }

  pos_dirs <- character(0)
  labels_cli <- NULL
  out_cli    <- NULL
  by_mode    <- NULL

  for (a in args) {
    if (grepl("^--dirs=", a)) {
      d_parts <- trimws(unlist(strsplit(sub("^--dirs=", "", a), ",")))
      pos_dirs <- c(pos_dirs, d_parts[nzchar(d_parts)])
    } else if (grepl("^--labels=", a)) {
      labels_cli <- trimws(unlist(strsplit(sub("^--labels=", "", a), ",")))
    } else if (grepl("^--out=", a) || grepl("^--output=", a)) {
      out_cli <- sub("^--(out|output)=", "", a)
    } else if (grepl("^--by=", a)) {
      by_mode <- tolower(sub("^--by=", "", a))
    } else if (a %in% c("--by-season", "--seasons")) {
      by_mode <- "season"
    } else if (a %in% c("--by-year", "--years")) {
      by_mode <- "year"
    } else if (a %in% c("--by-month", "--months")) {
      by_mode <- "month"
    } else if (!grepl("^--", a)) {
      pos_dirs <- c(pos_dirs, a)
    }
  }

  # If a single directory is passed (or --by=season|year|month is requested),
  # automatically discover temporal subdirectories inside it
  if (length(pos_dirs) == 1 && !is.null(by_mode)) {
    p_dir <- if (grepl("^/", pos_dirs[1])) pos_dirs[1] else file.path(BASE_DIR, pos_dirs[1])
    subs <- list.dirs(p_dir, recursive = FALSE, full.names = TRUE)
    if (by_mode == "season") {
      order_s <- c("season_winter", "season_spring", "season_summer", "season_autumn")
      matched <- subs[basename(subs) %in% order_s]
      pos_dirs <- matched[order(match(basename(matched), order_s))]
    } else if (by_mode == "year") {
      matched <- subs[grepl("^year_[0-9]{4}$", basename(subs))]
      pos_dirs <- matched[order(basename(matched))]
    } else if (by_mode == "month") {
      matched <- subs[grepl("^month_[0-9]{2}$", basename(subs))]
      pos_dirs <- matched[order(basename(matched))]
    }
  } else if (length(pos_dirs) == 1 && is.null(by_mode)) {
    # Check if single directory contains temporal subdirectories automatically
    p_dir <- if (grepl("^/", pos_dirs[1])) pos_dirs[1] else file.path(BASE_DIR, pos_dirs[1])
    subs <- list.dirs(p_dir, recursive = FALSE, full.names = TRUE)
    for (m_type in c("season", "year", "month")) {
      m_subs <- subs[grepl(paste0("^", m_type, "_"), basename(subs))]
      if (length(m_subs) >= 2) {
        pos_dirs <- sort(m_subs)
        break
      }
    }
  }

  compare_s3n_runs(run_dirs = pos_dirs, labels = labels_cli, out_dir = out_cli)
}

