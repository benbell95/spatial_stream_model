#!/usr/bin/env Rscript
# ==============================================================================
# 05_descriptive_stats.R
#
# Standalone script and reusable R functions for computing comprehensive
# descriptive statistics, monitoring network representativeness tests,
# stratified environmental/topological breakdowns, and base-R diagnostic
# plots for a single S3N model output directory.
#
# Also sourced by `scripts/04_compare_model_runs.R` so multi-run comparisons
# call the exact same descriptive statistics functions (`compute_detailed_descriptives`
# and `compute_stratified_descriptives`).
#
# Usage (Single Model Output):
#   Rscript scripts/05_descriptive_stats.R --list
#   Rscript scripts/05_descriptive_stats.R output/det_0117
#   Rscript scripts/05_descriptive_stats.R output/det_0180_orn
#   Rscript scripts/05_descriptive_stats.R output/det_0117_euhydro_thames --out=output/descriptives/thames_euhydro
# ==============================================================================

BASE_DIR <- if (file.exists("scripts/s3n_plot_functions.R")) "." else if (file.exists("s3n_plot_functions.R")) ".." else "."
source(file.path(BASE_DIR, "scripts/s3n_plot_functions.R"))
BASE_DIR <- resolve_base_dir()

# ------------------------------------------------------------------------------
# Helper: Formatting and Statistical Moment Utilities
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
load_model_run_bundle <- function(dir_path, custom_label = NULL, out_root = rel_path(BASE_DIR, "output")) {
  dir_p <- if (grepl("^/", dir_path)) dir_path else rel_path(BASE_DIR, dir_path)
  if (!dir.exists(dir_p) && !grepl("^/", dir_path)) {
    cand_in_out <- rel_path(out_root, dir_path)
    if (dir.exists(cand_in_out)) dir_p <- cand_in_out
  }
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

  trans_type <- if ("Transform_Type" %in% names(cv_tab) && nzchar(as.character(cv_tab$Transform_Type[1]))) {
    as.character(cv_tab$Transform_Type[1])
  } else if (use_log) {
    "log"
  } else {
    "none"
  }
  trans_lambda <- if ("Transform_Lambda" %in% names(cv_tab) && nzchar(as.character(cv_tab$Transform_Lambda[1]))) {
    as.numeric(cv_tab$Transform_Lambda[1])
  } else {
    NA_real_
  }
  trans_label <- if ("Transform_Label" %in% names(cv_tab) && nzchar(as.character(cv_tab$Transform_Label[1]))) {
    as.character(cv_tab$Transform_Label[1])
  } else if (use_log) {
    "Natural Logarithm ln(y)"
  } else {
    "Untransformed (raw)"
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
    } else if (trans_type != "log" && trans_type != "none") {
      custom_label <- sprintf("%s [%s]", net_lbl, trans_label)
    } else {
      custom_label <- net_lbl
    }
    if (grepl("_trans_", bname) && !grepl(trans_type, custom_label, ignore.case = TRUE)) {
      custom_label <- sprintf("%s [%s]", custom_label, trans_label)
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
    trans_type          = trans_type,
    trans_lambda        = trans_lambda,
    trans_label         = trans_label,
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
# Helper: List all available model output directories under output/
# ------------------------------------------------------------------------------
list_available_model_runs <- function(out_root = rel_path(BASE_DIR, "output")) {
  cat("==============================================================================\n")
  cat(sprintf("Available Model Run Output Directories (under '%s'):\n", out_root))
  cat("==============================================================================\n")
  all_dirs <- list.dirs(out_root, recursive = TRUE, full.names = TRUE)
  found <- 0L
  for (d in all_dirs) {
    if (grepl("/(comparisons|descriptives)(/|$)", d)) next
    has_det_files <- length(list.files(d, pattern = "^cv_comparison_.*\\.csv$")) > 0 &&
                     length(list.files(d, pattern = "^predictions_.*\\.csv$")) > 0
    has_std_files <- file.exists(file.path(d, "tables/cv_model_comparison.csv"))
    if (has_det_files || has_std_files) {
      found <- found + 1L
      rel_d <- format_rel_display_path(d, BASE_DIR)
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
# Core Function 1: Comprehensive Descriptive Statistics of Model Output Results
# Accepts a single bundle or a list of bundles.
# ------------------------------------------------------------------------------
compute_detailed_descriptives <- function(bundles) {
  if (!is.null(bundles$eng_pred)) {
    bundles <- list(bundles)
  }
  desc_rows <- list()
  rep_rows  <- list()

  for (i in seq_along(bundles)) {
    b  <- bundles[[i]]
    ep <- b$eng_pred
    obs_sub   <- ep[ep$has_obs == 1 & !is.na(ep$obs_value), ]
    unmon_sub <- ep[ep$has_obs == 0, ]

    vars_to_sum <- list(
      list(grp = "Observed Stations (Raw)",             vec = obs_sub$obs_value,       scale = "raw"),
      list(grp = "S3N Combined Predictions (Raw)",      vec = ep$pred_combined,        scale = "raw"),
      list(grp = "Unmonitored Reach Predictions (Raw)", vec = unmon_sub$pred_combined, scale = "raw"),
      list(grp = "S3N Tail-Up Predictions (Raw)",       vec = ep$pred_tu,              scale = "raw"),
      list(grp = "S3N Tail-Down Predictions (Raw)",     vec = ep$pred_td,              scale = "raw"),
      list(grp = "S3N Combined (Model Scale)",          vec = ep$pred_model_sc,        scale = if (b$use_log) "log" else "raw"),
      list(grp = "Kriging Standard Error (SE)",         vec = ep$se_model_sc,          scale = if (b$use_log) "log_se" else "raw_se")
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

    # Representativeness test: Observed stations vs. Full network predictions
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
# Core Function 2: Stratified Descriptive Statistics of Model Output Results
# Across Topography, Hydrogeology, Landscape Regime, Stream Order & Basins.
# Accepts a single bundle or a list of bundles.
# ------------------------------------------------------------------------------
compute_stratified_descriptives <- function(bundles) {
  if (!is.null(bundles$eng_pred)) {
    bundles <- list(bundles)
  }
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
# Core Function 3: Within-Run Statistical Hypothesis Tests for a Single Output
# ------------------------------------------------------------------------------
compute_single_run_stat_tests <- function(b1) {
  ep <- b1$eng_pred
  stat_rows <- list()
  add_stat_row <- function(comp_lbl, category, test_name, metric, n1, n2, val1, val2, diff_val,
                           stat_name, stat_val, p_val, interpretation) {
    stat_rows[[length(stat_rows) + 1L]] <<- data.frame(
      Comparison     = comp_lbl,
      Category       = category,
      Test_Name      = test_name,
      Metric         = metric,
      n_1            = n1,
      n_2            = n2,
      Value_Group1   = round(val1, 5),
      Value_Group2   = round(val2, 5),
      Diff_G2_G1     = round(diff_val, 5),
      Statistic_Name = stat_name,
      Statistic_Val  = round(stat_val, 5),
      P_Value        = p_val,
      Significance   = sig_stars(p_val),
      Interpretation = interpretation,
      stringsAsFactors = FALSE
    )
  }

  # 1. Cross-Validation Model Comparisons (S3N Combined vs OLS, Tail-Up, Tail-Down)
  if (!is.null(b1$eng_cv) && nrow(b1$eng_cv) >= 5) {
    e_comb <- b1$eng_cv$obs_y - b1$eng_cv$cv_comb
    if ("cv_ols" %in% names(b1$eng_cv)) {
      e_ols  <- b1$eng_cv$obs_y - b1$eng_cv$cv_ols
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
      e_tu  <- b1$eng_cv$obs_y - b1$eng_cv$cv_tu
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
      e_td  <- b1$eng_cv$obs_y - b1$eng_cv$cv_td
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
                 sprintf("Pred/Obs median ratio = %.2fx (%s)",
                         median(ep$pred_combined, na.rm = TRUE) / max(1e-9, median(obs_v, na.rm = TRUE)),
                         sig_stars(mw_rep$p.value)))
    add_stat_row("Network vs. Monitored Sites", "Monitoring Representativeness",
                 "Two-Sample Kolmogorov-Smirnov Test (Predictions vs Observed Stations)",
                 "Empirical CDF D-stat", length(obs_v), nrow(ep),
                 mean(obs_v, na.rm = TRUE), mean(ep$pred_combined, na.rm = TRUE),
                 mean(ep$pred_combined, na.rm = TRUE) - mean(obs_v, na.rm = TRUE),
                 "D", as.numeric(ks_rep$statistic), ks_rep$p.value,
                 sprintf("Distributional shift between gauged sites and full network (D = %.3f, %s)",
                         ks_rep$statistic, sig_stars(ks_rep$p.value)))
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

  # 5. Environmental Gradient Hypothesis Tests (Topography, Hydrogeology, Stream Order)
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
                           median(v_low, na.rm = TRUE) / max(1e-9, median(v_up, na.rm = TRUE)),
                           sig_stars(mw_topo$p.value)))
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
                           median(v_aq, na.rm = TRUE) / max(1e-9, median(v_noaq, na.rm = TRUE)),
                           sig_stars(mw_aq$p.value)))
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
                           median(v_ms, na.rm = TRUE) / max(1e-9, median(v_hw, na.rm = TRUE)),
                           sig_stars(kw_ord$p.value)))
    }
  }

  do.call(rbind, stat_rows)
}

# ------------------------------------------------------------------------------
# Single-Run Base-R Plotting Functions
# ------------------------------------------------------------------------------
plot_single_run_descriptives <- function(b1, desc_df, rep_df, strat_df, out_dir) {
  ep        <- b1$eng_pred
  det_short <- shorten_det_name(b1$det_name)
  det_unit  <- b1$det_unit
  use_log   <- b1$use_log
  u_lbl     <- if (use_log) sprintf("log(%s)", det_unit) else det_unit

  # Figure 1: 4-Panel Distributions, Representativeness, Tail-Up vs Tail-Down & Uncertainty
  fig1_file <- file.path(out_dir, "fig_desc1_distributions_and_representativeness.png")
  png(fig1_file, width = 2600, height = 2000, res = 200, type = "cairo")
  par(mfrow = c(2, 2), mar = c(4.8, 4.8, 3.6, 1.8), oma = c(0.5, 0.5, 2.2, 0.5), bg = "white")

  # Panel A: Kernel Density of Full Network Predictions vs. Observed Monitoring Stations
  obs_sc   <- ep$obs_model_sc[ep$has_obs == 1 & !is.na(ep$obs_model_sc)]
  unmon_sc <- ep$pred_model_sc[ep$has_obs == 0]
  d_all    <- density(ep$pred_model_sc, na.rm = TRUE, n = 512)
  d_unmon  <- if (length(unmon_sc) >= 10) density(unmon_sc, na.rm = TRUE, n = 512) else d_all
  d_obs    <- if (length(obs_sc) >= 10)   density(obs_sc,   na.rm = TRUE, n = 512) else d_all

  xr_a <- range(c(d_all$x, d_obs$x), na.rm = TRUE)
  yr_a <- c(0, max(c(d_all$y, d_obs$y), na.rm = TRUE) * 1.22)
  plot(NULL, xlim = xr_a, ylim = yr_a,
       xlab = sprintf("%s [%s]", det_short, u_lbl),
       ylab = "Probability Density",
       main = "A. Network Predictions vs. Observed Monitoring Stations",
       cex.main = 1.10, cex.lab = 0.95, las = 1)
  grid(col = "#e0e0e0", lty = "dotted")
  polygon(c(d_all$x, rev(d_all$x)), c(d_all$y, rep(0, length(d_all$y))),
          col = adjustcolor("#1f77b4", alpha.f = 0.22), border = NA)
  lines(d_all$x,   d_all$y,   col = "#1f77b4", lwd = 2.6, lty = 1)
  lines(d_unmon$x, d_unmon$y, col = "#2ca02c", lwd = 2.0, lty = 2)
  lines(d_obs$x,   d_obs$y,   col = "#d62728", lwd = 2.6, lty = 1)
  abline(v = median(ep$pred_model_sc, na.rm = TRUE), col = "#1f77b4", lty = 3, lwd = 1.4)
  abline(v = median(obs_sc, na.rm = TRUE),           col = "#d62728", lty = 3, lwd = 1.4)
  leg_pos_a <- if (median(ep$pred_model_sc, na.rm = TRUE) > mean(xr_a)) "topleft" else "topright"
  legend(leg_pos_a,
         legend = c(
           sprintf("All Reaches (Med = %.3g %s, n = %s)", median(ep$pred_combined, na.rm = TRUE), det_unit, format(nrow(ep), big.mark = ",")),
           sprintf("Unmonitored Reaches (Med = %.3g %s)", median(ep$pred_combined[ep$has_obs == 0], na.rm = TRUE), det_unit),
           sprintf("Observed Stations (Med = %.3g %s, n = %d)", median(ep$obs_value[ep$has_obs == 1], na.rm = TRUE), det_unit, length(obs_sc))
         ),
         col = c("#1f77b4", "#2ca02c", "#d62728"), lty = c(1, 2, 1), lwd = c(2.6, 2.0, 2.6),
         bty = "n", cex = 0.78)
  if (!is.null(rep_df) && nrow(rep_df) > 0) {
    mtext(sprintf("Pred/Obs Median Ratio = %.2fx | Mann-Whitney p = %s (%s) | KS D = %.3f (p = %s)",
                  rep_df$Pred_to_Obs_Med_Ratio[1], fmt_pval(rep_df$MW_P_Value[1]),
                  rep_df$Significance[1], rep_df$KS_D_Stat[1], fmt_pval(rep_df$KS_P_Value[1])),
          side = 3, line = 0.25, cex = 0.76, col = "#333333")
  }

  # Panel B: Empirical CDF Comparing S3N Combined, Tail-Up, Tail-Down & Observed Stations
  tu_sc <- if ("pred_tu_sc" %in% names(ep)) ep$pred_tu_sc else if (use_log) log(pmax(1e-6, ep$pred_tu)) else ep$pred_tu
  td_sc <- if ("pred_td_sc" %in% names(ep)) ep$pred_td_sc else if (use_log) log(pmax(1e-6, ep$pred_td)) else ep$pred_td
  xr_b  <- quantile(c(ep$pred_model_sc, obs_sc), probs = c(0.005, 0.995), na.rm = TRUE)
  xs_b  <- seq(xr_b[1], xr_b[2], length.out = 300)
  plot(NULL, xlim = xr_b, ylim = c(0, 1),
       xlab = sprintf("%s [%s]", det_short, u_lbl),
       ylab = "Cumulative Proportion",
       main = "B. Empirical CDF: S3N Combined, Tail-Up, Tail-Down & Observed",
       cex.main = 1.10, cex.lab = 0.95, las = 1)
  grid(col = "#e0e0e0", lty = "dotted")
  lines(xs_b, ecdf(ep$pred_model_sc)(xs_b), col = "#1f77b4", lwd = 2.5, lty = 1)
  lines(xs_b, ecdf(tu_sc)(xs_b),            col = "#2ca02c", lwd = 2.0, lty = 2)
  lines(xs_b, ecdf(td_sc)(xs_b),            col = "#9467bd", lwd = 2.0, lty = 4)
  if (length(obs_sc) >= 5) {
    lines(xs_b, ecdf(obs_sc)(xs_b),         col = "#d62728", lwd = 2.4, lty = 1)
  }
  legend("bottomright",
         legend = c(
           sprintf("S3N Combined (Med = %.3g %s)", median(ep$pred_combined, na.rm = TRUE), det_unit),
           sprintf("S3N Tail-Up (Med = %.3g %s)",  median(ep$pred_tu, na.rm = TRUE), det_unit),
           sprintf("S3N Tail-Down (Med = %.3g %s)",median(ep$pred_td, na.rm = TRUE), det_unit),
           sprintf("Observed Stations (Med = %.3g %s)", median(ep$obs_value[ep$has_obs == 1], na.rm = TRUE), det_unit)
         ),
         col = c("#1f77b4", "#2ca02c", "#9467bd", "#d62728"),
         lty = c(1, 2, 4, 1), lwd = c(2.5, 2.0, 2.0, 2.4), bty = "n", cex = 0.78)

  # Panel C: Reach-by-Reach Scatter of S3N Tail-Up vs. S3N Tail-Down Predictions
  set.seed(42)
  idx_pt <- if (nrow(ep) > 8500) sample(nrow(ep), 8500) else seq_len(nrow(ep))
  lim_c  <- quantile(c(tu_sc, td_sc), probs = c(0.005, 0.995), na.rm = TRUE)
  plot(tu_sc[idx_pt], td_sc[idx_pt],
       xlim = lim_c, ylim = lim_c, pch = 16, cex = 0.45, col = rgb(0.16, 0.50, 0.73, 0.25),
       xlab = sprintf("S3N Tail-Up Prediction [%s]", u_lbl),
       ylab = sprintf("S3N Tail-Down Prediction [%s]", u_lbl),
       main = sprintf("C. S3N Tail-Up vs. Tail-Down Component Agreement (n = %s)",
                      format(nrow(ep), big.mark = ",")),
       cex.main = 1.10, cex.lab = 0.95, las = 1)
  grid(col = "#e0e0e0", lty = "dotted")
  abline(0, 1, col = "#c0392b", lwd = 2.0, lty = 2)
  fit_c <- lm(td_sc ~ tu_sc)
  abline(fit_c, col = "#111111", lwd = 2.0, lty = 1)
  r_tud   <- cor(tu_sc, td_sc, use = "complete.obs", method = "pearson")
  rho_tud <- cor(tu_sc, td_sc, use = "complete.obs", method = "spearman")
  rmsd_c  <- sqrt(mean((td_sc - tu_sc)^2, na.rm = TRUE))
  legend("topleft",
         legend = c(
           sprintf("Pearson r = %.4f (Spearman rho = %.4f)", r_tud, rho_tud),
           sprintf("RMSD = %.4f | OLS Slope = %.3f", rmsd_c, coef(fit_c)[2]),
           "OLS Regression Fit",
           "1:1 Agreement Line"
         ),
         col = c(NA, NA, "#111111", "#c0392b"), lty = c(NA, NA, 1, 2), lwd = c(NA, NA, 2, 2),
         bty = "n", cex = 0.80)

  # Panel D: Kriging Standard Error (se_model_sc) by Monitoring Status & Stream Order
  se_groups <- list(
    "Monitored\nReaches"   = ep$se_model_sc[ep$has_obs == 1],
    "Unmonitored\nReaches" = ep$se_model_sc[ep$has_obs == 0],
    "Order 1-2\n(Headw.)"  = if ("Strahler" %in% names(ep)) ep$se_model_sc[!is.na(ep$Strahler) & ep$Strahler <= 2] else numeric(0),
    "Order 3-4\n(Medium)"  = if ("Strahler" %in% names(ep)) ep$se_model_sc[!is.na(ep$Strahler) & ep$Strahler %in% c(3, 4)] else numeric(0),
    "Order >=5\n(Mainstem)"= if ("Strahler" %in% names(ep)) ep$se_model_sc[!is.na(ep$Strahler) & ep$Strahler >= 5] else numeric(0)
  )
  se_groups <- se_groups[sapply(se_groups, length) >= 5]
  cols_d <- c("#27ae60", "#e67e22", "#3498db", "#2980b9", "#1f618d")[seq_along(se_groups)]
  bp_tmp <- boxplot(se_groups, plot = FALSE)
  ylim_d <- c(min(bp_tmp$stats[1, ], na.rm = TRUE) * 0.92,
              max(bp_tmp$stats[5, ], na.rm = TRUE) * 1.14)
  bp_d <- boxplot(se_groups, col = adjustcolor(cols_d, alpha.f = 0.65), border = "#222222",
                  outline = FALSE, las = 1, ylim = ylim_d,
                  ylab = sprintf("Kriging Standard Error (%s)", if (use_log) "log scale" else det_unit),
                  main = "D. Kriging Prediction Uncertainty (se_model_sc) Breakdown",
                  cex.main = 1.10, cex.lab = 0.95, cex.axis = 0.84)
  grid(nx = NA, ny = NULL, col = "#e0e0e0", lty = "dotted")
  se_mns <- sapply(se_groups, mean, na.rm = TRUE)
  points(seq_along(se_groups), se_mns, pch = 18, cex = 1.4, col = "#111111")
  text(seq_along(se_groups), bp_d$stats[5, ], labels = sprintf("Mean=%.3f", se_mns), pos = 3, cex = 0.76, font = 2)

  mtext(sprintf("S3N Single-Run Descriptive Diagnostics: %s (%s — %s)",
                b1$det_name, b1$label, b1$dir_name),
        side = 3, line = 0.4, outer = TRUE, cex = 1.15, font = 2, col = "#111111")
  dev.off()

  # Figure 2: 4-Panel Stratified Environmental & Topological Gradients
  fig2_file <- file.path(out_dir, "fig_desc2_stratified_gradients.png")
  png(fig2_file, width = 2600, height = 2000, res = 200, type = "cairo")
  par(mfrow = c(2, 2), mar = c(5.0, 4.8, 3.6, 1.8), oma = c(0.5, 0.5, 2.2, 0.5), bg = "white")

  # Panel A: Boxplot by Combined Landscape Regime (Upland, Lowland Non-Aquifer, Lowland Aquifer)
  if (("Elevation" %in% names(ep)) && ("Is_Aquifer" %in% names(ep))) {
    reg_list <- list(
      "Upland\n(>= 150 m)"         = ep$pred_model_sc[!is.na(ep$Elevation) & ep$Elevation >= 150.0],
      "Lowland\nNon-Aquifer"       = ep$pred_model_sc[!is.na(ep$Elevation) & ep$Elevation < 150.0 & ep$Is_Aquifer == 0],
      "Lowland\nAquifer (Chalk/Lm)"= ep$pred_model_sc[!is.na(ep$Elevation) & ep$Elevation < 150.0 & ep$Is_Aquifer == 1]
    )
    raw_meds_a <- c(
      median(ep$pred_combined[!is.na(ep$Elevation) & ep$Elevation >= 150.0], na.rm = TRUE),
      median(ep$pred_combined[!is.na(ep$Elevation) & ep$Elevation < 150.0 & ep$Is_Aquifer == 0], na.rm = TRUE),
      median(ep$pred_combined[!is.na(ep$Elevation) & ep$Elevation < 150.0 & ep$Is_Aquifer == 1], na.rm = TRUE)
    )
    bp_t <- boxplot(reg_list, plot = FALSE)
    ylim_a <- c(min(bp_t$stats[1, ], na.rm = TRUE) - 0.1 * diff(range(bp_t$stats, na.rm = TRUE)),
                max(bp_t$stats[5, ], na.rm = TRUE) + 0.18 * diff(range(bp_t$stats, na.rm = TRUE)))
    bp_a <- boxplot(reg_list, col = adjustcolor(c("#2980b9", "#27ae60", "#e67e22"), alpha.f = 0.65),
                    border = "#222222", outline = FALSE, las = 1, ylim = ylim_a,
                    ylab = sprintf("Predicted %s [%s]", det_short, u_lbl),
                    main = "A. Predictions Across Landscape & Hydrogeological Regimes",
                    cex.main = 1.10, cex.lab = 0.95, cex.axis = 0.88)
    grid(nx = NA, ny = NULL, col = "#e0e0e0", lty = "dotted")
    text(1:3, bp_a$stats[5, ], labels = sprintf("Med = %.3g %s\n(n = %s)", raw_meds_a, det_unit,
                                                format(sapply(reg_list, length), big.mark = ",")),
         pos = 3, cex = 0.76, font = 2)
  } else {
    plot.new()
  }

  # Panel B: Boxplot Across Individual Strahler Stream Orders
  if ("Strahler" %in% names(ep)) {
    ord_capped <- pmin(6L, pmax(1L, as.integer(ep$Strahler)))
    ord_levels <- sort(unique(ord_capped[!is.na(ord_capped)]))
    ord_list   <- lapply(ord_levels, function(o) ep$pred_model_sc[ord_capped == o])
    names(ord_list) <- ifelse(ord_levels == 6L, "Order >=6", paste0("Order ", ord_levels))
    raw_meds_b <- sapply(ord_levels, function(o) median(ep$pred_combined[ord_capped == o], na.rm = TRUE))
    bp_t <- boxplot(ord_list, plot = FALSE)
    ylim_b <- c(min(bp_t$stats[1, ], na.rm = TRUE) - 0.08 * diff(range(bp_t$stats, na.rm = TRUE)),
                max(bp_t$stats[5, ], na.rm = TRUE) + 0.16 * diff(range(bp_t$stats, na.rm = TRUE)))
    pal_b <- colorRampPalette(c("#85c1e9", "#2e86c1", "#1b4f72"))(length(ord_list))
    bp_b <- boxplot(ord_list, col = adjustcolor(pal_b, alpha.f = 0.70), border = "#222222",
                    outline = FALSE, las = 1, ylim = ylim_b,
                    ylab = sprintf("Predicted %s [%s]", det_short, u_lbl),
                    main = "B. Longitudinal Gradient Across Strahler Stream Orders",
                    cex.main = 1.10, cex.lab = 0.95, cex.axis = 0.86)
    grid(nx = NA, ny = NULL, col = "#e0e0e0", lty = "dotted")
    text(seq_along(ord_list), bp_b$stats[5, ], labels = sprintf("%.3g", raw_meds_b),
         pos = 3, cex = 0.76, font = 2)
  } else {
    plot.new()
  }

  # Panel C: Boxplot Across Major River Basins
  bas_masks <- list(
    "All England" = rep(TRUE, nrow(ep)),
    "Humber"      = if (!is.na(b1$hum_cid)) (ep$CatchID == b1$hum_cid) else rep(FALSE, nrow(ep)),
    "Thames"      = if (!is.na(b1$tha_cid)) (ep$CatchID == b1$tha_cid) else rep(FALSE, nrow(ep)),
    "Severn"      = if (!is.na(b1$sev_cid)) (ep$CatchID == b1$sev_cid) else rep(FALSE, nrow(ep))
  )
  bas_masks <- bas_masks[sapply(bas_masks, sum) >= 5]
  bas_list  <- lapply(bas_masks, function(m) ep$pred_model_sc[m])
  raw_meds_c <- sapply(bas_masks, function(m) median(ep$pred_combined[m], na.rm = TRUE))
  bp_t <- boxplot(bas_list, plot = FALSE)
  ylim_c <- c(min(bp_t$stats[1, ], na.rm = TRUE) - 0.08 * diff(range(bp_t$stats, na.rm = TRUE)),
              max(bp_t$stats[5, ], na.rm = TRUE) + 0.16 * diff(range(bp_t$stats, na.rm = TRUE)))
  cols_c <- c("#34495e", "#2980b9", "#d35400", "#27ae60")[seq_along(bas_list)]
  bp_c <- boxplot(bas_list, col = adjustcolor(cols_c, alpha.f = 0.65), border = "#222222",
                  outline = FALSE, las = 1, ylim = ylim_c,
                  ylab = sprintf("Predicted %s [%s]", det_short, u_lbl),
                  main = "C. Regional Comparison Across Major River Basins",
                  cex.main = 1.10, cex.lab = 0.95, cex.axis = 0.88)
  grid(nx = NA, ny = NULL, col = "#e0e0e0", lty = "dotted")
  text(seq_along(bas_list), bp_c$stats[5, ],
       labels = sprintf("Med = %.3g %s\n(n = %s)", raw_meds_c, det_unit,
                        format(sapply(bas_list, length), big.mark = ",")),
       pos = 3, cex = 0.75, font = 2)

  # Panel D: Binned Decile Profile of Predicted Value vs. Elevation (m)
  if ("Elevation" %in% names(ep) && sum(!is.na(ep$Elevation)) >= 50) {
    q_el <- unique(quantile(ep$Elevation, probs = seq(0, 1, by = 0.1), na.rm = TRUE))
    if (length(q_el) >= 4) {
      el_bin <- cut(ep$Elevation, breaks = q_el, include.lowest = TRUE)
      el_mid <- tapply(ep$Elevation,     el_bin, median, na.rm = TRUE)
      y_med  <- tapply(ep$pred_combined, el_bin, median, na.rm = TRUE)
      y_q25  <- tapply(ep$pred_combined, el_bin, quantile, probs = 0.25, na.rm = TRUE)
      y_q75  <- tapply(ep$pred_combined, el_bin, quantile, probs = 0.75, na.rm = TRUE)
      ok <- !is.na(el_mid) & !is.na(y_med)
      el_mid <- as.numeric(el_mid[ok]); y_med <- as.numeric(y_med[ok])
      y_q25  <- as.numeric(y_q25[ok]);  y_q75 <- as.numeric(y_q75[ok])

      plot(el_mid, y_med, type = "n",
           xlim = range(el_mid, na.rm = TRUE),
           ylim = c(0, max(y_q75, na.rm = TRUE) * 1.15),
           xlab = "Reach Elevation Decile Midpoint (m AOD)",
           ylab = sprintf("Predicted %s (%s, Median & IQR)", det_short, det_unit),
           main = "D. Altitudinal Gradient (Decile Median & IQR Profile)",
           cex.main = 1.10, cex.lab = 0.95, las = 1)
      grid(col = "#e0e0e0", lty = "dotted")
      polygon(c(el_mid, rev(el_mid)), c(y_q25, rev(y_q75)),
              col = adjustcolor("#1f77b4", alpha.f = 0.22), border = NA)
      lines(el_mid, y_med, col = "#1f77b4", lwd = 2.6, type = "b", pch = 16, cex = 1.1)
      abline(v = 150.0, col = "#c0392b", lwd = 1.6, lty = 2)
      legend("topright",
             legend = c("Decile Median Prediction", "Interquartile Range (Q25 - Q75)", "150 m Upland Threshold"),
             col = c("#1f77b4", adjustcolor("#1f77b4", alpha.f = 0.35), "#c0392b"),
             lwd = c(2.6, 8, 1.6), lty = c(1, 1, 2), pch = c(16, NA, NA), bty = "n", cex = 0.80)
    } else {
      plot.new()
    }
  } else {
    plot.new()
  }

  mtext(sprintf("S3N Stratified Environmental & Topological Gradients: %s (%s)",
                b1$det_name, b1$label),
        side = 3, line = 0.4, outer = TRUE, cex = 1.15, font = 2, col = "#111111")
  dev.off()

  invisible(c(fig1_file, fig2_file))
}

# ------------------------------------------------------------------------------
# Single-Run Markdown Report Writer
# ------------------------------------------------------------------------------
write_descriptive_markdown_report <- function(b1, desc_df, rep_df, strat_df, stat_df, out_dir, report_path) {
  lines <- character(0)
  add_line <- function(fmt = "", ...) {
    if (length(list(...)) == 0L) {
      lines <<- c(lines, fmt)
    } else {
      lines <<- c(lines, sprintf(fmt, ...))
    }
  }

  add_line("# S3N Single-Run Descriptive Statistics & Diagnostic Report")
  add_line("")
  add_line("- **Determinand**: `%s` — **%s** (`%s`, %s scale)",
           b1$det_code, b1$det_name, b1$det_unit, if (b1$use_log) "log" else "raw")
  add_line("- **Data Transformation**: **%s** (`%s`%s)",
           b1$trans_label, b1$trans_type,
           if (!is.na(b1$trans_lambda)) sprintf(", λ = %.4f", b1$trans_lambda) else "")
  add_line("- **Water Network**: %s (`%s` reaches across `%d` catchments, `%d` observed monitoring sites)",
           b1$net_label, format(nrow(b1$eng_pred), big.mark = ","),
           length(unique(b1$eng_pred$CatchID)), sum(b1$eng_pred$has_obs == 1))
  add_line("- **Model Run Directory**: `%s`", b1$dir_name)
  add_line("- **Fixed-Effect Covariates**: %s", b1$vars_label)
  add_line("- **Output Directory**: `%s`", format_rel_display_path(out_dir, BASE_DIR))
  add_line("")

  add_line("## 1. Full Network & Monitoring Station Descriptive Statistics")
  add_line("")
  add_line("| Variable Group | Unit | N | Mean ± SD | CV (%) | Geom Mean | Skew | Kurt | Min | P05 | Q25 | Median | Q75 | P95 | Max | IQR | MAD |")
  add_line("| :--- | :---: | ---: | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
  for (r_i in seq_len(nrow(desc_df))) {
    d <- desc_df[r_i, ]
    gm_str <- if (is.na(d$Geom_Mean)) "—" else sprintf("%.4g", d$Geom_Mean)
    cv_str <- if (is.na(d$CV_Pct))    "—" else sprintf("%.1f%%", d$CV_Pct)
    add_line("| **%s** | `%s` | %s | `%.4g ± %.4g` | `%s` | `%s` | `%+.2f` | `%+.2f` | `%.4g` | `%.4g` | `%.4g` | **`%.4g`** | `%.4g` | `%.4g` | `%.4g` | `%.4g` | `%.4g` |",
             d$Variable_Group, d$Unit, format(d$N, big.mark = ","),
             d$Mean, d$SD, cv_str, gm_str, d$Skewness, d$Excess_Kurtosis,
             d$Min, d$P05, d$Q25, d$Median, d$Q75, d$P95, d$Max, d$IQR, d$MAD)
  }
  add_line("")

  if (!is.null(rep_df) && nrow(rep_df) > 0) {
    add_line("## 2. Monitoring Network Representativeness (Observed Stations vs. Full River Network)")
    add_line("")
    add_line("| Monitored Sites | Network Reaches | Obs Median (Mean) | Network Pred Median (Mean) | Pred / Obs Median Ratio | Mann–Whitney p | KS D-stat (p) | Interpretation |")
    add_line("| ---: | ---: | :--- | :--- | ---: | ---: | ---: | :--- |")
    for (r_i in seq_len(nrow(rep_df))) {
      rp <- rep_df[r_i, ]
      add_line("| %d | %s | `%.4g` (`%.4g`) | `%.4g` (`%.4g`) | **`%.3fx`** | `%s` (`%s`) | `%.3f` (`%s`) | %s |",
               rp$n_Obs_Stations, format(rp$n_Total_Reaches, big.mark = ","),
               rp$Obs_Median_Raw, rp$Obs_Mean_Raw, rp$Pred_Median_Raw, rp$Pred_Mean_Raw,
               rp$Pred_to_Obs_Med_Ratio, fmt_pval(rp$MW_P_Value), rp$Significance,
               rp$KS_D_Stat, fmt_pval(rp$KS_P_Value), rp$Interpretation)
    }
    add_line("")
  }

  if (!is.null(strat_df) && nrow(strat_df) > 0) {
    add_line("## 3. Stratified Model Output Breakdown (Topography, Hydrogeology, Stream Order & Major Basins)")
    add_line("")
    add_line("| Dimension | Stratum | Reaches (% Net) | Obs Sites (Med) | Pred Mean ± SD | Pred Median (IQR) | Pred 5th–95th Pct | Mean Z-Score | Kriging SE Med | Fold Ratio vs Ref | Within-Run Test (p-val) |")
    add_line("| :--- | :--- | ---: | ---: | :--- | :--- | :--- | ---: | ---: | ---: | :--- |")
    for (r_i in seq_len(nrow(strat_df))) {
      st <- strat_df[r_i, ]
      obs_s <- if (is.na(st$Obs_Median_Raw)) sprintf("%d (—)", st$n_obs) else sprintf("%d (`%.3g`)", st$n_obs, st$Obs_Median_Raw)
      rat_s <- if (is.na(st$Ratio_vs_Ref)) "—" else sprintf("%.2fx", st$Ratio_vs_Ref)
      tst_s <- if (is.na(st$Within_Run_P)) st$Within_Run_Test else sprintf("%s: p = %s (%s)", st$Within_Run_Test, fmt_pval(st$Within_Run_P), st$Within_Run_Sig)
      add_line("| %s | **%s** | %s (`%.1f%%`) | %s | `%.4g ± %.4g` | **`%.4g`** (`%.4g`) | `%.4g` – `%.4g` | `%+.2f` | `%.4f` | `%s` | %s |",
               st$Dimension, st$Stratum, format(st$n_reaches, big.mark = ","), st$Pct_Network,
               obs_s, st$Pred_Mean_Raw, st$Pred_SD_Raw, st$Pred_Median_Raw, st$Pred_IQR_Raw,
               st$Pred_P05_Raw, st$Pred_P95_Raw, st$Pred_Mean_Z, st$SE_Median, rat_s, tst_s)
    }
    add_line("")
  }

  if (!is.null(stat_df) && nrow(stat_df) > 0) {
    add_line("## 4. Within-Run Statistical Hypothesis Tests")
    add_line("")
    add_line("| Comparison | Category | Statistical Test | Metric | Group 1 Val | Group 2 Val | Diff / R² | Stat | p-value | Sig | Finding |")
    add_line("| :--- | :--- | :--- | :--- | ---: | ---: | ---: | :--- | ---: | :---: | :--- |")
    for (r_i in seq_len(nrow(stat_df))) {
      s <- stat_df[r_i, ]
      add_line("| %s | %s | %s | %s | `%.4f` | `%.4f` | `%+.4f` | `%s = %.3f` | `%s` | `%s` | %s |",
               s$Comparison, s$Category, s$Test_Name, s$Metric,
               s$Value_Group1, s$Value_Group2, s$Diff_G2_G1,
               s$Statistic_Name, s$Statistic_Val, fmt_pval(s$P_Value),
               s$Significance, s$Interpretation)
    }
    add_line("")
  }

  add_line("## 5. Cross-Validation Performance Summary")
  add_line("")
  add_line("| Domain | Model | n_obs | n_reaches | CV R² | CV RMSE | CV MAE |")
  add_line("| :--- | :--- | ---: | ---: | ---: | ---: | ---: |")
  for (r_i in seq_len(nrow(b1$cv_tab))) {
    r <- b1$cv_tab[r_i, ]
    add_line("| %s | %s | %d | %s | `%.4f` | `%.4f` | `%.4f` |",
             r$Canon_Domain, r$Canon_Model, r$n_obs,
             format(r$n_reaches, big.mark = ","), r$R2, r$RMSE, r$MAE)
  }
  add_line("")

  add_line("## 6. Generated Files")
  add_line("")
  add_line("- `descriptive_report.md` — Full single-run descriptive & diagnostic Markdown report.")
  add_line("- `descriptive_stats_summary.csv` — 18-column distributional summary across all model output and observation variables.")
  add_line("- `monitoring_representativeness.csv` — Mann–Whitney U and KS test comparing gauged stations against the full stream network.")
  add_line("- `stratified_results_summary.csv` — Stratified descriptive statistics across topography, hydrogeology, stream order, and basins.")
  add_line("- `within_run_statistical_tests.csv` — Formal statistical tests for CV model improvements, component agreement, and environmental gradients.")
  add_line("- `fig_desc1_distributions_and_representativeness.png` — 4-panel base-R plot of distributions, ECDF, Tail-Up vs Tail-Down agreement, and Kriging SE.")
  add_line("- `fig_desc2_stratified_gradients.png` — 4-panel base-R plot of landscape regimes, Strahler stream orders, major river basins, and altitudinal profile.")
  add_line("")

  writeLines(lines, report_path)
}

# ------------------------------------------------------------------------------
# Main Single-Run Entry Point: run_descriptive_stats()
# ------------------------------------------------------------------------------
run_descriptive_stats <- function(run_dir, label = NULL, out_dir = NULL, out_root = rel_path(BASE_DIR, "output"), quiet = FALSE) {
  b1 <- load_model_run_bundle(run_dir, custom_label = label, out_root = out_root)

  if (is.null(out_dir) || !nzchar(out_dir)) {
    out_dir <- rel_path(b1$dir_path, "descriptives")
  } else if (!grepl("^/", out_dir)) {
    out_dir <- rel_path(BASE_DIR, out_dir)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (!quiet) {
    cat("==============================================================================\n")
    cat(sprintf("Running Single-Output Descriptive Statistics & Diagnostics:\n"))
    cat(sprintf("  Run Label:    %s (%s)\n", b1$label, b1$dir_name))
    cat(sprintf("  Determinand:  %s [%s] (%s, %s scale)\n",
                b1$det_name, b1$det_code, b1$det_unit, if (b1$use_log) "log" else "raw"))
    cat(sprintf("  Network:      %s (%s reaches, %d observed stations)\n",
                b1$net_label, format(nrow(b1$eng_pred), big.mark = ","), sum(b1$eng_pred$has_obs == 1)))
    cat(sprintf("  Output Dir:   %s\n", out_dir))
    cat("==============================================================================\n")
  }

  desc_bundle <- compute_detailed_descriptives(list(b1))
  desc_df     <- desc_bundle$desc_df
  rep_df      <- desc_bundle$rep_df
  strat_df    <- compute_stratified_descriptives(list(b1))
  stat_df     <- compute_single_run_stat_tests(b1)

  desc_csv  <- file.path(out_dir, "descriptive_stats_summary.csv")
  rep_csv   <- file.path(out_dir, "monitoring_representativeness.csv")
  strat_csv <- file.path(out_dir, "stratified_results_summary.csv")
  stat_csv  <- file.path(out_dir, "within_run_statistical_tests.csv")
  rep_md    <- file.path(out_dir, "descriptive_report.md")

  write.csv(desc_df,  desc_csv,  row.names = FALSE)
  if (!is.null(rep_df)) write.csv(rep_df, rep_csv, row.names = FALSE)
  write.csv(strat_df, strat_csv, row.names = FALSE)
  write.csv(stat_df,  stat_csv,  row.names = FALSE)

  fig_paths <- plot_single_run_descriptives(b1, desc_df, rep_df, strat_df, out_dir)
  write_descriptive_markdown_report(b1, desc_df, rep_df, strat_df, stat_df, out_dir, rep_md)

  if (!quiet) {
    cat("\nDescriptive analysis complete! Generated files:\n")
    cat(sprintf("  - Markdown Report:      %s\n", rep_md))
    cat(sprintf("  - Descriptives CSV:     %s\n", desc_csv))
    if (!is.null(rep_df)) {
      cat(sprintf("  - Representativeness:   %s\n", rep_csv))
    }
    cat(sprintf("  - Stratified CSV:       %s\n", strat_csv))
    cat(sprintf("  - Statistical Tests:    %s\n", stat_csv))
    cat(sprintf("  - Figure 1 (Dist/SE):   %s\n", fig_paths[1]))
    cat(sprintf("  - Figure 2 (Gradients): %s\n", fig_paths[2]))
  }

  invisible(list(
    bundle      = b1,
    out_dir     = out_dir,
    report_path = rep_md,
    desc_df     = desc_df,
    rep_df      = rep_df,
    strat_df    = strat_df,
    stat_df     = stat_df,
    fig_paths   = fig_paths
  ))
}

# ==============================================================================
# Command-Line Interface (only runs when 05_descriptive_stats.R is invoked directly)
# ==============================================================================
if (sys.nframe() == 0L) {
  args    <- commandArgs(trailingOnly = TRUE)
  dir_cfg <- parse_custom_dir_args(args, base_dir = BASE_DIR)

  if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
    cat("Usage:\n")
    cat("  Rscript scripts/05_descriptive_stats.R --list [--output-dir=<dir>]\n")
    cat("  Rscript scripts/05_descriptive_stats.R <output_dir> [--label=\"Custom Label\"] [--out=<dir>] [--output-dir=<dir>]\n")
    cat("Examples:\n")
    cat("  Rscript scripts/05_descriptive_stats.R output/det_0117\n")
    cat("  Rscript scripts/05_descriptive_stats.R output/det_0180_orn\n")
    cat("  Rscript scripts/05_descriptive_stats.R output/det_0117_euhydro_thames\n")
    quit(save = "no", status = 0)
  }

  if (any(args %in% c("--list", "--list-runs"))) {
    list_available_model_runs(out_root = dir_cfg$out_dir)
    quit(save = "no", status = 0)
  }

  pos_dirs  <- character(0)
  label_cli <- NULL
  out_cli   <- NULL

  for (a in args) {
    if (grepl("^--dir=", a)) {
      pos_dirs <- c(pos_dirs, sub("^--dir=", "", a))
    } else if (grepl("^--label=", a)) {
      label_cli <- sub("^--label=", "", a)
    } else if (grepl("^--(out|output)=", a)) {
      out_cli <- sub("^--(out|output)=", "", a)
    } else if (grepl("^--(output-dir|out-dir|data-dir|data|gis-dir|gis)=", a)) {
      # Handled via parse_custom_dir_args()
    } else if (!grepl("^--", a)) {
      pos_dirs <- c(pos_dirs, a)
    }
  }

  if (length(pos_dirs) == 0) {
    stop("Please provide a model output directory (e.g., Rscript scripts/05_descriptive_stats.R output/det_0117).")
  }

  for (d_i in seq_along(pos_dirs)) {
    out_i <- if (!is.null(out_cli) && length(pos_dirs) > 1) {
      rel_path(out_cli, basename(pos_dirs[d_i]))
    } else {
      out_cli
    }
    run_descriptive_stats(run_dir = pos_dirs[d_i], label = label_cli, out_dir = out_i, out_root = dir_cfg$out_dir)
  }
}
