#!/usr/bin/env Rscript
# ==============================================================================
# 03_generate_plots.R
#
# Standalone plotting script that sources the shared plotting module
# (scripts/s3n_plot_functions.R) to generate publication-quality S3N maps and
# diagnostic figures from saved prediction and cross-validation tables:
#   - Figure 1: Full England Stream Network Map (Predictions & Kriging SE)
#   - Figure 2: Humber River Basin High-Resolution S3N Map
#   - Figure 3: Thames & Severn River Basins S3N Maps
#   - Figure 4: 4-Panel Model Diagnostics, Cross-Validation & Parameter CIs
#
# Supports both directory layouts and automatically synchronizes color scale
# bars across seasonal, yearly, or monthly subdirectories!
#
# Usage:
#   Rscript scripts/03_generate_plots.R
#   Rscript scripts/03_generate_plots.R --dir=output/det_0117
#   Rscript scripts/03_generate_plots.R --dir=output/rhs_hqa_score --reverse-color
#   Rscript scripts/03_generate_plots.R --dir=output/det_0117 --vmin=0.2 --vmax=20
# ==============================================================================

BASE_DIR <- if (file.exists("scripts/s3n_plot_functions.R")) "." else if (file.exists("s3n_plot_functions.R")) ".." else "."
source(file.path(BASE_DIR, "scripts/s3n_plot_functions.R"))
BASE_DIR <- resolve_base_dir()

args    <- commandArgs(trailingOnly = TRUE)
dir_cfg <- parse_custom_dir_args(args, base_dir = BASE_DIR)

target_dir    <- dir_cfg$out_dir
det_code_cli  <- NULL
det_name_cli  <- NULL
det_unit_cli  <- NULL
use_log_cli   <- NULL
reverse_color <- FALSE
vmin_raw      <- NULL
vmax_raw      <- NULL
custom_basins <- character(0)
norm_cli      <- parse_norm_mode(args, default = "none")
norm_meta     <- get_norm_meta(norm_cli)

for (a in args) {
  if (grepl("^--(dir|out-dir|output-dir|out|output)=", a)) {
    d_val <- sub("^--(dir|out-dir|output-dir|out|output)=", "", a)
    cand_direct <- if (grepl("^/", d_val)) d_val else rel_path(BASE_DIR, d_val)
    if (!dir.exists(cand_direct) && !grepl("^/", d_val)) {
      cand_res <- rel_path(dir_cfg$out_dir, "results", d_val)
      cand_out <- rel_path(dir_cfg$out_dir, d_val)
      if (dir.exists(cand_res)) {
        cand_direct <- cand_res
      } else if (dir.exists(cand_out)) {
        cand_direct <- cand_out
      }
    }
    target_dir <- cand_direct
  } else if (grepl("^--code=", a)) {
    det_code_cli <- sub("^--code=", "", a)
  } else if (grepl("^--name=", a)) {
    det_name_cli <- sub("^--name=", "", a)
  } else if (grepl("^--unit=", a)) {
    det_unit_cli <- sub("^--unit=", "", a)
  } else if (a == "--no-log") {
    use_log_cli <- FALSE
  } else if (a == "--log") {
    use_log_cli <- TRUE
  } else if (a %in% c("--reverse-color", "--reverse-colour", "--rev-color", "--rev-colour", "--reverse-scale")) {
    reverse_color <- TRUE
  } else if (a %in% c("--normal-color", "--no-reverse-color")) {
    reverse_color <- FALSE
  } else if (grepl("^--vmin=", a)) {
    vmin_raw <- as.numeric(sub("^--vmin=", "", a))
  } else if (grepl("^--vmax=", a)) {
    vmax_raw <- as.numeric(sub("^--vmax=", "", a))
  } else if (grepl("^--basin=", a) || grepl("^--catchid=", a)) {
    b_val <- sub("^--(basin|catchid)=", "", a)
    b_parts <- trimws(unlist(strsplit(b_val, ",")))
    custom_basins <- c(custom_basins, b_parts[nzchar(b_parts)])
  }
}

find_first_csv <- function(dir_path, pattern, exclude_pattern = NULL, required = TRUE) {
  files <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  if (!is.null(exclude_pattern)) {
    files <- files[!grepl(exclude_pattern, basename(files))]
  }
  if (length(files) == 0) {
    if (required) stop("Could not find file matching '", pattern, "' in ", dir_path)
    return(NULL)
  }
  files[1]
}

known_meta <- list(
  nitrate           = c(name = "Nitrate as N",         unit = "mg/L"),
  "0117"            = c(name = "Nitrate as N",         unit = "mg/L"),
  whpt_aspt         = c(name = "WHPT ASPT",            unit = "ASPT"),
  whpt_ntaxa        = c(name = "WHPT N-Taxa",          unit = "taxa"),
  life_family_index = c(name = "LIFE Flow Index",      unit = "LIFE"),
  psi_family_score  = c(name = "PSI Sediment Score",   unit = "%"),
  epsi_family_score = c(name = "EPSI Sediment Score",  unit = "%"),
  rmni              = c(name = "RMNI Nutrient Index",  unit = "RMNI"),
  hqa_score         = c(name = "RHS HQA Score",        unit = "HQA"),
  hms_score         = c(name = "RHS HMS Modification", unit = "HMS")
)

load_plot_bundle <- function(dir_p) {
  if (dir.exists(file.path(dir_p, "predictions"))) {
    PRED_DIR <- file.path(dir_p, "predictions")
    TAB_DIR  <- file.path(dir_p, "tables")
    FIG_DIR  <- file.path(dir_p, "figures")
    dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

    eng_pred_file <- find_first_csv(PRED_DIR, "^england_.*_predictions\\.csv$", "_cv_predictions\\.csv$")
    hum_pred_file <- find_first_csv(PRED_DIR, "^humber_.*_predictions\\.csv$",  "_cv_predictions\\.csv$")
    tha_pred_file <- find_first_csv(PRED_DIR, "^thames_.*_predictions\\.csv$",  "_cv_predictions\\.csv$")
    sev_pred_file <- find_first_csv(PRED_DIR, "^severn_.*_predictions\\.csv$",  "_cv_predictions\\.csv$")

    eng_cv_file   <- find_first_csv(PRED_DIR, "^england_.*cv_predictions\\.csv$")
    hum_cv_file   <- find_first_csv(PRED_DIR, "^humber_.*cv_predictions\\.csv$", required = FALSE)

    eng_pred  <- normalize_pred_df(read.csv(eng_pred_file))
    hum_pred  <- normalize_pred_df(read.csv(hum_pred_file))
    tha_pred  <- normalize_pred_df(read.csv(tha_pred_file))
    sev_pred  <- normalize_pred_df(read.csv(sev_pred_file))

    cv_tab    <- normalize_cv_tab(read.csv(file.path(TAB_DIR, "cv_model_comparison.csv")))
    param_tab <- read.csv(file.path(TAB_DIR, "s3n_parameter_estimates.csv"))
    eng_cv    <- normalize_cv_pts(read.csv(eng_cv_file))
    hum_cv    <- if (!is.null(hum_cv_file)) normalize_cv_pts(read.csv(hum_cv_file)) else NULL
    is_default_out <- TRUE
  } else {
    FIG_DIR <- dir_p
    eng_pred_file <- find_first_csv(dir_p, "^predictions_.*\\.csv$",
                                    "^predictions_(humber|thames|severn|river_|catchid_).*\\.csv$")
    hum_pred_file <- find_first_csv(dir_p, "^predictions_humber_.*\\.csv$", required = FALSE)
    tha_pred_file <- find_first_csv(dir_p, "^predictions_thames_.*\\.csv$", required = FALSE)
    sev_pred_file <- find_first_csv(dir_p, "^predictions_severn_.*\\.csv$", required = FALSE)

    cv_tab_file   <- find_first_csv(dir_p, "^cv_comparison_.*\\.csv$")
    param_file    <- find_first_csv(dir_p, "^parameters_.*\\.csv$")
    eng_cv_file   <- find_first_csv(dir_p, "^cv_predictions_england_.*\\.csv$", required = FALSE)
    hum_cv_file   <- find_first_csv(dir_p, "^cv_predictions_humber_.*\\.csv$",  required = FALSE)

    eng_pred  <- normalize_pred_df(read.csv(eng_pred_file))
    hum_cid   <- resolve_core_basin(eng_pred, target = "humber")$catch_id
    tha_cid   <- resolve_core_basin(eng_pred, target = "thames")$catch_id
    sev_cid   <- resolve_core_basin(eng_pred, target = "severn")$catch_id
    hum_pred  <- if (!is.null(hum_pred_file)) normalize_pred_df(read.csv(hum_pred_file)) else subset(eng_pred, CatchID == hum_cid)
    tha_pred  <- if (!is.null(tha_pred_file)) normalize_pred_df(read.csv(tha_pred_file)) else subset(eng_pred, CatchID == tha_cid)
    sev_pred  <- if (!is.null(sev_pred_file)) normalize_pred_df(read.csv(sev_pred_file)) else subset(eng_pred, CatchID == sev_cid)

    cv_tab    <- normalize_cv_tab(read.csv(cv_tab_file))
    param_tab <- read.csv(param_file)

    if (!is.null(eng_cv_file)) {
      eng_cv <- normalize_cv_pts(read.csv(eng_cv_file))
    } else {
      obs_e  <- subset(eng_pred, has_obs == 1)
      eng_cv <- data.frame(RivID = obs_e$RivID, obs_y = obs_e$obs_model_sc, cv_comb = obs_e$pred_model_sc)
    }
    if (!is.null(hum_cv_file)) {
      hum_cv <- normalize_cv_pts(read.csv(hum_cv_file))
    } else {
      obs_h  <- subset(hum_pred, has_obs == 1)
      hum_cv <- if (nrow(obs_h) > 0) data.frame(RivID = obs_h$RivID, obs_y = obs_h$obs_model_sc, cv_comb = obs_h$pred_model_sc) else NULL
    }
    is_default_out <- FALSE
  }

  det_code <- det_code_cli
  if (is.null(det_code)) {
    if ("Determinand_Code" %in% names(cv_tab)) {
      det_code <- as.character(cv_tab$Determinand_Code[1])
      if (nchar(det_code) < 4 && grepl("^[0-9]+$", det_code)) {
        det_code <- sprintf("%04d", as.integer(det_code))
      }
    } else if (is_default_out) {
      det_code <- sub("^england_(.*)_predictions\\.csv$", "\\1", basename(eng_pred_file))
    } else {
      det_code <- sub("^predictions_(.*)\\.csv$", "\\1", basename(eng_pred_file))
    }
  }

  det_name <- det_name_cli
  if (is.null(det_name)) {
    if ("Determinand_Name" %in% names(cv_tab)) {
      det_name <- as.character(cv_tab$Determinand_Name[1])
    } else if (det_code %in% names(known_meta)) {
      det_name <- unname(known_meta[[det_code]]["name"])
    } else {
      det_name <- det_code
    }
  }

  det_unit <- det_unit_cli
  if (is.null(det_unit)) {
    if ("Unit" %in% names(cv_tab)) {
      det_unit <- as.character(cv_tab$Unit[1])
    } else if (det_code %in% names(known_meta)) {
      det_unit <- unname(known_meta[[det_code]]["unit"])
    } else {
      det_unit <- "mg/L"
    }
  }

  use_log <- use_log_cli
  if (is.null(use_log)) {
    if ("Scale" %in% names(cv_tab)) {
      use_log <- (as.character(cv_tab$Scale[1]) == "log")
    } else {
      use_log <- ("pred_log_combined" %in% names(eng_pred)) ||
                 any(abs(exp(eng_pred$pred_model_sc[1:min(20, nrow(eng_pred))]) -
                         eng_pred$pred_combined[1:min(20, nrow(eng_pred))]) < 1e-4, na.rm = TRUE)
    }
  }

  if (norm_meta$mode != "none") {
    if (!(norm_meta$pred_col %in% names(eng_pred))) eng_pred <- add_normalised_columns(eng_pred, ref_df = eng_pred)
    if (!(norm_meta$pred_col %in% names(hum_pred))) hum_pred <- add_normalised_columns(hum_pred, ref_df = eng_pred)
    if (!(norm_meta$pred_col %in% names(tha_pred))) tha_pred <- add_normalised_columns(tha_pred, ref_df = eng_pred)
    if (!(norm_meta$pred_col %in% names(sev_pred))) sev_pred <- add_normalised_columns(sev_pred, ref_df = eng_pred)
  }

  list(
    FIG_DIR        = FIG_DIR,
    is_default_out = is_default_out,
    det_code       = det_code,
    det_name       = det_name,
    det_unit       = det_unit,
    use_log        = use_log,
    eng_pred       = eng_pred,
    hum_pred       = hum_pred,
    tha_pred       = tha_pred,
    sev_pred       = sev_pred,
    cv_tab         = cv_tab,
    param_tab      = param_tab,
    eng_cv         = eng_cv,
    hum_cv         = hum_cv
  )
}

# Discover whether target_dir has temporal subdirectories (year_*, season_*, month_*)
sub_dirs <- list.dirs(target_dir, recursive = FALSE, full.names = TRUE)
period_dirs <- sub_dirs[grepl("^(year_|season_|month_)", basename(sub_dirs))]
has_root_preds <- dir.exists(file.path(target_dir, "predictions")) ||
                  length(list.files(target_dir, pattern = "^predictions_.*\\.csv$")) > 0

dirs_to_plot <- character(0)
if (has_root_preds) dirs_to_plot <- c(dirs_to_plot, target_dir)
if (length(period_dirs) > 0) dirs_to_plot <- c(dirs_to_plot, period_dirs)

if (length(dirs_to_plot) == 0) {
  stop("No prediction files or period subdirectories found in ", target_dir)
}

cat(sprintf("Loading S3N predictions across %d directory(s) under %s...\n", length(dirs_to_plot), target_dir))
bundles <- lapply(dirs_to_plot, load_plot_bundle)

# Also check if target_dir is a single period subfolder (e.g. output/det_0117/season_summer)
# and pool scale ranges with any sibling period folders in dirname(target_dir)
scale_vals <- unlist(lapply(bundles, function(b) b$eng_pred[[norm_meta$pred_col]]))
scale_se_e <- unlist(lapply(bundles, function(b) b$eng_pred$se_model_sc))
scale_se_h <- unlist(lapply(bundles, function(b) b$hum_pred$se_model_sc))

if (length(bundles) == 1 && grepl("^(year_|season_|month_)", basename(target_dir))) {
  parent_d <- dirname(target_dir)
  sib_csvs <- list.files(parent_d, pattern = "^predictions_.*\\.csv$", recursive = TRUE, full.names = TRUE)
  sib_csvs <- sib_csvs[!grepl("^predictions_(humber|thames|severn|river_|catchid_)", basename(sib_csvs))]
  if (length(sib_csvs) > 1) {
    for (sc in sib_csvs) {
      sdf <- tryCatch(normalize_pred_df(read.csv(sc)), error = function(e) NULL)
      if (!is.null(sdf)) {
        if (norm_meta$mode != "none" && !(norm_meta$pred_col %in% names(sdf))) {
          sdf <- add_normalised_columns(sdf)
        }
        scale_vals <- c(scale_vals, sdf[[norm_meta$pred_col]])
        scale_se_e <- c(scale_se_e, sdf$se_model_sc)
      }
    }
  }
}

use_log_shared <- bundles[[1]]$use_log
shared_sc_val  <- if (norm_meta$mode == "none") {
  compute_color_scale(scale_vals, use_log_shared, vmin_raw = vmin_raw, vmax_raw = vmax_raw)
} else {
  compute_norm_color_scale(scale_vals, norm_meta$mode, vmin_raw = vmin_raw, vmax_raw = vmax_raw)
}
shared_sc_se_e <- compute_se_scale(scale_se_e)
shared_sc_se_h <- compute_se_scale(scale_se_h)

all_cv_vals   <- unlist(lapply(bundles, function(b) c(b$eng_cv$obs_y, b$eng_cv$cv_comb)))
lims_a_shared <- range(all_cv_vals, na.rm = TRUE)
pad_a         <- 0.08 * diff(lims_a_shared)
lims_a_shared <- c(lims_a_shared[1] - pad_a, lims_a_shared[2] + pad_a)

all_raw_preds <- unlist(lapply(bundles, function(b) quantile(b$eng_pred$pred_combined, probs = c(0.01, 0.96), na.rm = TRUE)))
yr_d          <- range(all_raw_preds, na.rm = TRUE)
ylim_d_shared <- if (length(bundles) > 1) c(max(0, yr_d[1] - 0.05 * diff(yr_d)), yr_d[2] + 0.15 * diff(yr_d)) else NULL

for (b in bundles) {
  ftag    <- norm_meta$file_tag
  f1_name <- if (b$is_default_out) sprintf("fig1_england_stream_nitrate_map%s.png", ftag) else sprintf("fig1_england_stream_%s_map%s.png", b$det_code, ftag)
  f2_name <- if (b$is_default_out) sprintf("fig2_humber_basin_s3n_map%s.png", ftag)       else sprintf("fig2_humber_basin_%s_map%s.png", b$det_code, ftag)
  f3_name <- if (b$is_default_out) sprintf("fig3_thames_severn_s3n_maps%s.png", ftag)     else sprintf("fig3_thames_severn_%s_maps%s.png", b$det_code, ftag)
  f4_name <- if (b$is_default_out) "fig4_cv_and_diagnostics.png"         else sprintf("fig4_cv_and_diagnostics_%s.png", b$det_code)

  sc_lo_disp <- if (norm_meta$mode == "none") (if (b$use_log) exp(shared_sc_val$vmin) else shared_sc_val$vmin) else shared_sc_val$vmin
  sc_hi_disp <- if (norm_meta$mode == "none") (if (b$use_log) exp(shared_sc_val$vmax) else shared_sc_val$vmax) else shared_sc_val$vmax
  cat(sprintf("Rendering 4 figures in %s (shared scale [%s, %s]%s)...\n",
              b$FIG_DIR,
              formatC(sc_lo_disp, format = "fg", digits = 3),
              formatC(sc_hi_disp, format = "fg", digits = 3),
              if (norm_meta$mode == "none") "" else paste0(", norm = ", norm_meta$mode)))

  vars_sub <- format_covariates_subtitle(param_tab = b$param_tab)
  net_lbl  <- extract_network_label_from_tables(b$cv_tab, b$param_tab, b$FIG_DIR)

  plot_fig1_england_map(b$eng_pred, b$det_code, b$det_name, b$det_unit, b$use_log,
                        file.path(b$FIG_DIR, f1_name), reverse_color = reverse_color,
                        sc_val = shared_sc_val, sc_se = shared_sc_se_e,
                        vars_label = vars_sub, network_label = net_lbl, norm_mode = norm_cli)

  plot_fig2_humber_map(b$hum_pred, b$cv_tab, shared_sc_val, b$det_code, b$det_name, b$det_unit, b$use_log,
                       file.path(b$FIG_DIR, f2_name), reverse_color = reverse_color,
                       sc_se = shared_sc_se_h, vars_label = vars_sub, network_label = net_lbl, norm_mode = norm_cli)

  plot_fig3_thames_severn_map(b$tha_pred, b$sev_pred, shared_sc_val, b$det_code, b$det_name, b$det_unit, b$use_log,
                              file.path(b$FIG_DIR, f3_name), reverse_color = reverse_color,
                              vars_label = vars_sub, network_label = net_lbl, norm_mode = norm_cli)

  plot_fig4_diagnostics(b$eng_pred, b$eng_cv, b$hum_cv, b$cv_tab, b$param_tab,
                        b$det_code, b$det_name, b$det_unit, b$use_log,
                        file.path(b$FIG_DIR, f4_name),
                        lims_a = if (length(bundles) > 1) lims_a_shared else NULL,
                        ylim_d = ylim_d_shared,
                        vars_label = vars_sub, network_label = net_lbl)

  if (length(custom_basins) > 0) {
    for (bi in seq_along(custom_basins)) {
      b_spec <- custom_basins[bi]
      if (tolower(b_spec) %in% c("thames", "humber", "severn")) {
        b_res <- resolve_core_basin(b$eng_pred, target = tolower(b_spec))
        b_cid <- b_res$catch_id
        b_title <- b_res$domain_label
      } else if (grepl("^[0-9]+$", b_spec)) {
        b_cid <- as.integer(b_spec)
        b_title <- sprintf("CatchID %d Basin", b_cid)
      } else {
        b_cid <- NA_integer_
        b_title <- b_spec
      }
      b_pred <- subset(b$eng_pred, CatchID == b_cid)
      if (nrow(b_pred) > 0) {
        b_slug <- gsub("[^a-z0-9]+", "_", tolower(b_title))
        b_slug <- gsub("^_|_$", "", b_slug)
        f5_b_name <- file.path(b$FIG_DIR, sprintf("fig%d_%s_%s_map%s.png", 4 + bi, b_slug, b$det_code, ftag))
        cat(sprintf("Rendering custom basin figure: %s (%d reaches)...\n", basename(f5_b_name), nrow(b_pred)))
        plot_basin_map(b_pred, b_title, b_title, b$cv_tab, shared_sc_val, b$det_code, b$det_name, b$det_unit, b$use_log,
                       f5_b_name, reverse_color = reverse_color, vars_label = vars_sub, network_label = net_lbl, norm_mode = norm_cli)
      }
    }
  }
}

# If temporal subdirectories (year_*, season_*, month_*) exist under target_dir,
# also generate Figure 5 (multi-panel All England comparison grid) for each temporal group
for (m_type in c("season", "year", "month")) {
  m_idx <- which(grepl(paste0("^", m_type, "_"), basename(dirs_to_plot)))
  if (length(m_idx) >= 2) {
    m_bundles <- bundles[m_idx]
    grid_items <- lapply(m_bundles, function(b) {
      d_base <- basename(b$FIG_DIR)
      lbl <- if (m_type == "season") {
        paste0("Season: ", sub("^season_", "", d_base))
      } else if (m_type == "year") {
        paste0("Year ", sub("^year_", "", d_base))
      } else {
        sprintf("Month %02d", as.integer(sub("^month_", "", d_base)))
      }
      list(label = lbl, eng_pred = b$eng_pred, cv_tab = b$cv_tab)
    })
    b0 <- m_bundles[[1]]
    vars_sub_g <- format_covariates_subtitle(param_tab = b0$param_tab)
    net_lbl_g  <- extract_network_label_from_tables(b0$cv_tab, b0$param_tab, target_dir)
    fig5_path <- file.path(target_dir, sprintf("fig5_england_by_%s_%s_map%s.png", m_type, b0$det_code, norm_meta$file_tag))
    cat(sprintf("Rendering Figure 5 (%d-panel All England %s comparison grid): %s...\n",
                length(grid_items), m_type, fig5_path))
    plot_fig5_temporal_england_grid(grid_items, m_type, shared_sc_val,
                                    b0$det_code, b0$det_name, b0$det_unit, b0$use_log,
                                    fig5_path, reverse_color = reverse_color,
                                    vars_label = vars_sub_g, network_label = net_lbl_g, norm_mode = norm_cli)
  }
}

cat("Done!\n")
