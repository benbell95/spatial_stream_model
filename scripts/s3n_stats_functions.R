# ==============================================================================
# s3n_stats_functions.R
#
# Shared model run bundle loader, path resolution, run discovery, statistical
# moment / hypothesis-test utilities, descriptive & stratified summary engines,
# and buffered Markdown writer shared by 04_compare_model_runs.R and
# 05_descriptive_stats.R.
# ==============================================================================

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
  d[grepl("^All England", d)] <- "All England Stream Network"
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

resolve_run_dir_path <- function(dir_path, out_root = NULL, base_dir = BASE_DIR) {
  default_out <- load_directory_config(base_dir)$output_dir
  if (is.null(out_root) || !nzchar(out_root)) {
    out_root <- default_out
  }
  if (grepl("^/", dir_path)) return(dir_path)

  out_roots <- unique(c(out_root, default_out))
  stripped <- sub("^(output|results)/+", "", dir_path)

  cands <- rel_path(base_dir, dir_path)
  for (root in out_roots) {
    cands <- c(cands,
               rel_path(root, "results", dir_path),
               rel_path(root, dir_path))
    if (!identical(stripped, dir_path)) {
      cands <- c(cands,
                 rel_path(root, "results", stripped),
                 rel_path(root, stripped))
    }
  }
  for (cp in unique(cands)) {
    if (dir.exists(cp)) return(cp)
  }
  cands[1]
}

# ------------------------------------------------------------------------------
# Helper: Load a single model run bundle from an output directory
# ------------------------------------------------------------------------------
load_model_run_bundle <- function(dir_path, custom_label = NULL, out_root = NULL) {
  if (is.null(out_root) || !nzchar(out_root)) {
    out_root <- load_directory_config(BASE_DIR)$output_dir
  }
  dir_p <- resolve_run_dir_path(dir_path, out_root = out_root, base_dir = BASE_DIR)
  if (!dir.exists(dir_p)) {
    stop(sprintf("Model output directory does not exist: %s", dir_p))
  }

  eng_file  <- find_first_file(dir_p, "^predictions_.*\\.csv$",
                               "^predictions_(humber|thames|severn|basin_|river_|catchid_).*\\.csv$")
  hum_file  <- find_first_file(dir_p, "^predictions_humber_.*\\.csv$")
  tha_file  <- find_first_file(dir_p, "^predictions_thames_.*\\.csv$")
  sev_file  <- find_first_file(dir_p, "^predictions_severn_.*\\.csv$")
  cv_t_file <- find_first_file(dir_p, "^cv_comparison_.*\\.csv$")
  par_file  <- find_first_file(dir_p, "^parameters_.*\\.csv$")
  cv_e_file <- find_first_file(dir_p, "^cv_predictions_england_.*\\.csv$")
  cv_h_file <- find_first_file(dir_p, "^cv_predictions_humber_.*\\.csv$")

  if (is.null(eng_file) || is.null(cv_t_file) || !file.exists(cv_t_file)) {
    stop(sprintf("Directory '%s' does not contain required predictions_*.csv and cv_comparison_*.csv files.", dir_p))
  }

  eng_pred  <- read.csv(eng_file, stringsAsFactors = FALSE)
  if (!all(c("pred_zscore", "pred_robust_zscore", "pred_minmax", "pred_percentile") %in% names(eng_pred))) {
    eng_pred <- add_normalised_columns(eng_pred, ref_df = eng_pred)
  }
  cv_tab    <- read.csv(cv_t_file, stringsAsFactors = FALSE)
  param_tab <- if (!is.null(par_file) && file.exists(par_file)) read.csv(par_file, stringsAsFactors = FALSE) else NULL
  eng_cv    <- if (!is.null(cv_e_file) && file.exists(cv_e_file)) read.csv(cv_e_file, stringsAsFactors = FALSE) else NULL
  hum_cv    <- if (!is.null(cv_h_file) && file.exists(cv_h_file)) read.csv(cv_h_file, stringsAsFactors = FALSE) else NULL

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
  det_meta_row <- get_determinand_meta(det_code, BASE_DIR, default_unit = "mg/L")
  det_name <- if ("Determinand_Name" %in% names(cv_tab)) {
    as.character(cv_tab$Determinand_Name[1])
  } else {
    unname(det_meta_row["name"])
  }
  det_unit <- if ("Unit" %in% names(cv_tab)) {
    as.character(cv_tab$Unit[1])
  } else if (nzchar(det_meta_row["unit"])) {
    unname(det_meta_row["unit"])
  } else {
    "mg/L"
  }

  trans_type <- if ("Transform_Type" %in% names(cv_tab) && nzchar(as.character(cv_tab$Transform_Type[1]))) {
    as.character(cv_tab$Transform_Type[1])
  } else {
    "log"
  }
  trans_lambda <- if ("Transform_Lambda" %in% names(cv_tab) && nzchar(as.character(cv_tab$Transform_Lambda[1]))) {
    as.numeric(cv_tab$Transform_Lambda[1])
  } else {
    NA_real_
  }
  trans_scale <- if ("Transform_Scale" %in% names(cv_tab) && nzchar(as.character(cv_tab$Transform_Scale[1]))) {
    as.numeric(cv_tab$Transform_Scale[1])
  } else {
    NA_real_
  }
  trans_label <- if ("Transform_Label" %in% names(cv_tab) && nzchar(as.character(cv_tab$Transform_Label[1]))) {
    as.character(cv_tab$Transform_Label[1])
  } else if (trans_type == "log") {
    "Natural Logarithm ln(y)"
  } else {
    build_transformer(trans_type, lambda = if (is.na(trans_lambda)) NULL else trans_lambda,
                      scale = if (is.na(trans_scale)) NULL else trans_scale)$label
  }

  source_abbr <- if ("Source" %in% names(cv_tab) && nzchar(as.character(cv_tab$Source[1]))) {
    as.character(cv_tab$Source[1])
  } else {
    bname <- basename(dir_p)
    b_no_t <- sub("_(trans|vars)_.*$", "", bname)
    reg_s <- load_source_registry(BASE_DIR)
    cand_s <- names(reg_s$sources)
    cand_s <- cand_s[order(-nchar(cand_s))]
    hit_s <- NA_character_
    for (sk in cand_s) {
      if (grepl(paste0("_", sk, "($|_)"), b_no_t)) {
        hit_s <- sk
        break
      }
    }
    if (!is.na(hit_s)) format_source_abbr(hit_s, base_dir = BASE_DIR) else "RSN"
  }
  source_label <- if ("Source_Desc" %in% names(cv_tab) && nzchar(as.character(cv_tab$Source_Desc[1]))) {
    as.character(cv_tab$Source_Desc[1])
  } else {
    format_source_label(source_abbr, base_dir = BASE_DIR)
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
    if (grepl("_trans_", bname, fixed = TRUE) && !grepl(tolower(trans_type), tolower(custom_label), fixed = TRUE)) {
      custom_label <- sprintf("%s [%s]", custom_label, trans_label)
    }
    if (source_abbr != "RSN" && !grepl(tolower(source_abbr), tolower(custom_label), fixed = TRUE)) {
      custom_label <- sprintf("%s [Source: %s]", custom_label, source_abbr)
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
    source              = source_abbr,
    source_label        = source_label,
    period              = period_lbl,
    det_code            = det_code,
    det_name            = det_name,
    det_unit            = det_unit,
    trans_type          = trans_type,
    trans_lambda        = trans_lambda,
    trans_scale         = trans_scale,
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
list_available_model_runs <- function(out_root = NULL) {
  if (is.null(out_root) || !nzchar(out_root)) {
    out_root <- load_directory_config(BASE_DIR)$output_dir
  }
  cat("==============================================================================\n")
  cat(sprintf("Available Model Run Output Directories (under '%s'):\n", out_root))
  cat("==============================================================================\n")
  all_dirs <- if (dir.exists(out_root)) list.dirs(out_root, recursive = TRUE, full.names = TRUE) else character(0)
  out_norm <- normalizePath(out_root, mustWork = FALSE)
  found <- 0L
  for (d in all_dirs) {
    if (normalizePath(d, mustWork = FALSE) == out_norm) next
    if (grepl("/(comparisons|descriptives|station_snapping)(/|$)", d)) next
    cv_det_files <- list.files(d, pattern = "^cv_comparison_.*\\.csv$", full.names = TRUE)
    has_det_files <- length(cv_det_files) > 0 && length(list.files(d, pattern = "^predictions_.*\\.csv$")) > 0
    if (has_det_files) {
      found <- found + 1L
      rel_d <- format_rel_display_path(d, BASE_DIR)
      cv_f  <- cv_det_files[1]
      info <- tryCatch({
        cv_t <- read.csv(cv_f, stringsAsFactors = FALSE)
        cv_t$Canon_Domain <- canonical_domain(cv_t$Domain)
        cv_t$Canon_Model  <- canonical_model(cv_t$Model)
        r_all <- subset(cv_t, grepl("All England", Canon_Domain) & Canon_Model == "S3N Combined")
        if (nrow(r_all) == 0L) r_all <- cv_t[1, , drop = FALSE]
        net_lbl <- if ("Network_Label" %in% names(cv_t) && nzchar(as.character(cv_t$Network_Label[1]))) {
          as.character(cv_t$Network_Label[1])
        } else {
          extract_network_label_from_tables(cv_t, NULL, d)
        }
        d_name <- if ("Determinand_Name" %in% names(cv_t) && nzchar(as.character(cv_t$Determinand_Name[1]))) {
          as.character(cv_t$Determinand_Name[1])
        } else {
          basename(d)
        }
        r2_str <- if (nrow(r_all) > 0) sprintf("CV R2 = %.3f, RMSE = %.3f", r_all$R2[1], r_all$RMSE[1]) else "N/A"
        n_rch  <- if (nrow(r_all) > 0 && "n_reaches" %in% names(r_all)) as.integer(r_all$n_reaches[1]) else 0L
        n_ob   <- if (nrow(r_all) > 0 && "n_obs" %in% names(r_all)) as.integer(r_all$n_obs[1]) else 0L
        list(net_label = net_lbl, det_name = d_name, n_reaches = n_rch, n_obs = n_ob, r2_str = r2_str)
      }, error = function(e) NULL)
      if (!is.null(info)) {
        cat(sprintf("  [%2d] %-45s | %-20s | %-22s | %6d reaches (%4d obs) | %s\n",
                    found, rel_d, substr(info$net_label, 1, 20), substr(info$det_name, 1, 22),
                    info$n_reaches, info$n_obs, info$r2_str))
      } else {
        cat(sprintf("  [%2d] %s\n", found, rel_d))
      }
    }
  }
  if (found == 0L) {
    cat(sprintf("  (No fitted model runs found under '%s'.)\n", format_rel_display_path(out_root, BASE_DIR)))
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
      list(grp = "S3N Combined (Model Scale)",          vec = ep$pred_model_sc,        scale = b$trans_type),
      list(grp = "Kriging Standard Error (SE)",         vec = ep$se_model_sc,          scale = paste0(b$trans_type, "_se"))
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
    sd_sc <- sd(ep$pred_model_sc, na.rm = TRUE)
    if (!is.finite(sd_sc) || sd_sc < 1e-9) sd_sc <- 1e-9
    z_all <- (ep$pred_model_sc - mean(ep$pred_model_sc, na.rm = TRUE)) / sd_sc

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
      p_b <- if (sum(m_b) >= 5L && sum(!m_b) >= 5L) suppressWarnings(wilcox.test(ep$pred_model_sc[m_b], ep$pred_model_sc[!m_b]))$p.value else NA_real_
      add_stratum("5. River Basin", "Humber Basin", m_b, "Mann-Whitney U vs Rest of Network", p_b,
                  if (!is.na(med_eng) && abs(med_eng) > 1e-9) median(ep$pred_combined[m_b], na.rm = TRUE) / med_eng else NA_real_)
    }
    if (!is.na(b$tha_cid) && any(ep$CatchID == b$tha_cid)) {
      m_b <- ep$CatchID == b$tha_cid
      p_b <- if (sum(m_b) >= 5L && sum(!m_b) >= 5L) suppressWarnings(wilcox.test(ep$pred_model_sc[m_b], ep$pred_model_sc[!m_b]))$p.value else NA_real_
      add_stratum("5. River Basin", "Thames Basin", m_b, "Mann-Whitney U vs Rest of Network", p_b,
                  if (!is.na(med_eng) && abs(med_eng) > 1e-9) median(ep$pred_combined[m_b], na.rm = TRUE) / med_eng else NA_real_)
    }
    if (!is.na(b$sev_cid) && any(ep$CatchID == b$sev_cid)) {
      m_b <- ep$CatchID == b$sev_cid
      p_b <- if (sum(m_b) >= 5L && sum(!m_b) >= 5L) suppressWarnings(wilcox.test(ep$pred_model_sc[m_b], ep$pred_model_sc[!m_b]))$p.value else NA_real_
      add_stratum("5. River Basin", "Severn Basin", m_b, "Mann-Whitney U vs Rest of Network", p_b,
                  if (!is.na(med_eng) && abs(med_eng) > 1e-9) median(ep$pred_combined[m_b], na.rm = TRUE) / med_eng else NA_real_)
    }
  }

  do.call(rbind, strat_rows)
}


# ------------------------------------------------------------------------------
# Helper: Buffered Markdown Report Writer & Prediction Distribution Summarizer
# ------------------------------------------------------------------------------
new_md_writer <- function(init_cap = 512L) {
  env <- new.env(parent = emptyenv())
  env$buf <- character(init_cap)
  env$n   <- 0L
  add_line <- function(fmt = "", ...) {
    line_str <- if (length(list(...)) == 0L) fmt else sprintf(fmt, ...)
    env$n <- env$n + 1L
    if (env$n > length(env$buf)) {
      length(env$buf) <- length(env$buf) * 2L
    }
    env$buf[env$n] <- line_str
  }
  get_lines <- function() {
    if (env$n == 0L) character(0) else env$buf[seq_len(env$n)]
  }
  list(add_line = add_line, get_lines = get_lines)
}

summarize_prediction_bundles <- function(bundles) {
  if (!is.null(bundles$eng_pred)) {
    bundles <- list(bundles)
  }
  pred_sum_list <- lapply(seq_along(bundles), function(i) {
    b  <- bundles[[i]]
    ep <- b$eng_pred
    obs_sub <- ep[ep$has_obs == 1, ]
    q_raw <- quantile(ep$pred_combined, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE)
    q_sc  <- quantile(ep$pred_model_sc, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE)
    q_se  <- quantile(ep$se_model_sc,   probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE)
    data.frame(
      Run_Index         = i,
      Run_Label         = b$label,
      Determinand_Code  = b$det_code,
      Determinand_Name  = b$det_name,
      Unit              = b$det_unit,
      Scale             = b$trans_type,
      Transform_Type    = b$trans_type,
      Transform_Lambda  = ifelse(is.na(b$trans_lambda), "", as.character(round(b$trans_lambda, 4))),
      Transform_Label   = b$trans_label,
      Network_Label     = b$net_label,
      Period            = b$period,
      Covariates        = b$vars_label,
      n_reaches         = nrow(ep),
      n_catchments      = length(unique(ep$CatchID)),
      n_obs_reaches     = nrow(obs_sub),
      Obs_Mean_Raw      = mean(obs_sub$obs_value, na.rm = TRUE),
      Obs_SD_Raw        = sd(obs_sub$obs_value, na.rm = TRUE),
      Obs_Median_Raw    = median(obs_sub$obs_value, na.rm = TRUE),
      Obs_IQR_Raw       = IQR(obs_sub$obs_value, na.rm = TRUE),
      Pred_Mean_Raw     = mean(ep$pred_combined, na.rm = TRUE),
      Pred_SD_Raw       = sd(ep$pred_combined, na.rm = TRUE),
      Pred_CV_Pct_Raw   = 100.0 * sd(ep$pred_combined, na.rm = TRUE) / max(1e-9, abs(mean(ep$pred_combined, na.rm = TRUE))),
      Pred_GeomMean_Raw = calc_geom_mean(ep$pred_combined),
      Pred_Skew_Raw     = calc_sample_skewness(ep$pred_combined),
      Pred_P05_Raw      = as.numeric(q_raw[1]),
      Pred_Q25_Raw      = as.numeric(q_raw[2]),
      Pred_Median_Raw   = as.numeric(q_raw[3]),
      Pred_Q75_Raw      = as.numeric(q_raw[4]),
      Pred_P95_Raw      = as.numeric(q_raw[5]),
      Pred_IQR_Raw      = as.numeric(q_raw[4] - q_raw[2]),
      Pred_Mean_ModSc   = mean(ep$pred_model_sc, na.rm = TRUE),
      Pred_SD_ModSc     = sd(ep$pred_model_sc, na.rm = TRUE),
      Pred_Med_ModSc    = as.numeric(q_sc[3]),
      SE_Mean           = mean(ep$se_model_sc, na.rm = TRUE),
      SE_SD             = sd(ep$se_model_sc, na.rm = TRUE),
      SE_Median         = as.numeric(q_se[3]),
      SE_IQR            = as.numeric(q_se[4] - q_se[2]),
      SE_P95            = as.numeric(q_se[5]),
      stringsAsFactors  = FALSE
    )
  })
  do.call(rbind, pred_sum_list)
}
