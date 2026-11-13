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
trans_cli     <- NULL
lambda_cli    <- NULL
scale_cli     <- NULL
reverse_color <- FALSE
vmin_raw      <- NULL
vmax_raw      <- NULL
custom_basins <- character(0)
norm_cli      <- parse_norm_mode(args, default = "none")
norm_meta     <- get_norm_meta(norm_cli)

for (a in args) {
  if (grepl("^--(dir|out-dir|output-dir|out|output)=", a)) {
    d_val <- sub("^--(dir|out-dir|output-dir|out|output)=", "", a)
    cand_direct <- if (is_abs_path(d_val)) d_val else rel_path(BASE_DIR, d_val)
    if (!dir.exists(cand_direct) && !is_abs_path(d_val)) {
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
  } else if (grepl("^--(transform|trans|transformation)=", a)) {
    trans_cli <- tolower(sub("^--(transform|trans|transformation)=", "", a))
  } else if (grepl("^--(lambda|transform-lambda)=", a)) {
    lambda_cli <- suppressWarnings(as.numeric(sub("^--(lambda|transform-lambda)=", "", a)))
  } else if (grepl("^--(scale|transform-scale)=", a)) {
    scale_cli <- suppressWarnings(as.numeric(sub("^--(scale|transform-scale)=", "", a)))
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

# Determinand metadata resolved dynamically via get_determinand_meta() in configs/determinands.json

load_plot_bundle <- function(dir_p) {
  FIG_DIR <- dir_p
  eng_pred_file <- find_first_csv(dir_p, "^predictions_.*\\.csv$",
                                  "^predictions_(humber|thames|severn|basin_|river_|catchid_).*\\.csv$")
  hum_pred_file <- find_first_csv(dir_p, "^predictions_humber_.*\\.csv$", required = FALSE)
  tha_pred_file <- find_first_csv(dir_p, "^predictions_thames_.*\\.csv$", required = FALSE)
  sev_pred_file <- find_first_csv(dir_p, "^predictions_severn_.*\\.csv$", required = FALSE)

  cv_tab_file   <- find_first_csv(dir_p, "^cv_comparison_.*\\.csv$")
  param_file    <- find_first_csv(dir_p, "^parameters_.*\\.csv$")
  eng_cv_file   <- find_first_csv(dir_p, "^cv_predictions_england_.*\\.csv$", required = FALSE)
  hum_cv_file   <- find_first_csv(dir_p, "^cv_predictions_humber_.*\\.csv$",  required = FALSE)

  eng_pred  <- read.csv(eng_pred_file, stringsAsFactors = FALSE)
  hum_cid   <- resolve_core_basin(eng_pred, target = "humber")$catch_id
  tha_cid   <- resolve_core_basin(eng_pred, target = "thames")$catch_id
  sev_cid   <- resolve_core_basin(eng_pred, target = "severn")$catch_id
  hum_pred  <- if (!is.null(hum_pred_file)) read.csv(hum_pred_file, stringsAsFactors = FALSE) else subset(eng_pred, CatchID == hum_cid)
  tha_pred  <- if (!is.null(tha_pred_file)) read.csv(tha_pred_file, stringsAsFactors = FALSE) else subset(eng_pred, CatchID == tha_cid)
  sev_pred  <- if (!is.null(sev_pred_file)) read.csv(sev_pred_file, stringsAsFactors = FALSE) else subset(eng_pred, CatchID == sev_cid)

  cv_tab    <- read.csv(cv_tab_file, stringsAsFactors = FALSE)
  param_tab <- read.csv(param_file, stringsAsFactors = FALSE)
  eng_cv    <- if (!is.null(eng_cv_file)) read.csv(eng_cv_file, stringsAsFactors = FALSE) else NULL
  hum_cv    <- if (!is.null(hum_cv_file)) read.csv(hum_cv_file, stringsAsFactors = FALSE) else NULL

  det_code <- det_code_cli
  if (is.null(det_code)) {
    if ("Determinand_Code" %in% names(cv_tab)) {
      det_code <- as.character(cv_tab$Determinand_Code[1])
      if (nchar(det_code) < 4 && grepl("^[0-9]+$", det_code)) {
        det_code <- sprintf("%04d", as.integer(det_code))
      }
    } else {
      det_code <- sub("^predictions_(.*)\\.csv$", "\\1", basename(eng_pred_file))
    }
  }

  det_meta_row <- get_determinand_meta(det_code, BASE_DIR, default_unit = "mg/L")

  det_name <- det_name_cli
  if (is.null(det_name)) {
    if ("Determinand_Name" %in% names(cv_tab)) {
      det_name <- as.character(cv_tab$Determinand_Name[1])
    } else {
      det_name <- unname(det_meta_row["name"])
    }
  }

  det_unit <- det_unit_cli
  if (is.null(det_unit)) {
    if ("Unit" %in% names(cv_tab)) {
      det_unit <- as.character(cv_tab$Unit[1])
    } else {
      det_unit <- if (nzchar(det_meta_row["unit"])) unname(det_meta_row["unit"]) else "mg/L"
    }
  }

  if (!is.null(trans_cli) && nzchar(trans_cli)) {
    trans_obj <- build_transformer(trans_cli, lambda = lambda_cli, scale = scale_cli)
  } else {
    t_type  <- if ("Transform_Type" %in% names(cv_tab) && nzchar(as.character(cv_tab$Transform_Type[1]))) {
      tolower(trimws(as.character(cv_tab$Transform_Type[1])))
    } else {
      "log"
    }
    t_lam   <- if ("Transform_Lambda" %in% names(cv_tab)) suppressWarnings(as.numeric(as.character(cv_tab$Transform_Lambda[1]))) else NA_real_
    t_scale <- if ("Transform_Scale"  %in% names(cv_tab)) suppressWarnings(as.numeric(as.character(cv_tab$Transform_Scale[1])))  else NA_real_
    eff_lam   <- if (!is.null(lambda_cli) && is.finite(lambda_cli)) lambda_cli else if (is.finite(t_lam)) t_lam else NULL
    eff_scale <- if (!is.null(scale_cli)  && is.finite(scale_cli))  scale_cli  else if (is.finite(t_scale)) t_scale else NULL
    trans_obj <- build_transformer(t_type, lambda = eff_lam, scale = eff_scale)
  }

  source_lbl <- if ("Source_Desc" %in% names(cv_tab) && nzchar(as.character(cv_tab$Source_Desc[1]))) {
    as.character(cv_tab$Source_Desc[1])
  } else if ("Source" %in% names(cv_tab) && nzchar(as.character(cv_tab$Source[1]))) {
    as.character(cv_tab$Source[1])
  } else {
    NULL
  }

  if (norm_meta$mode != "none") {
    if (!(norm_meta$pred_col %in% names(eng_pred))) eng_pred <- add_normalised_columns(eng_pred, ref_df = eng_pred)
    if (!(norm_meta$pred_col %in% names(hum_pred))) hum_pred <- add_normalised_columns(hum_pred, ref_df = eng_pred)
    if (!(norm_meta$pred_col %in% names(tha_pred))) tha_pred <- add_normalised_columns(tha_pred, ref_df = eng_pred)
    if (!(norm_meta$pred_col %in% names(sev_pred))) sev_pred <- add_normalised_columns(sev_pred, ref_df = eng_pred)
  }

  list(
    bundle_dir     = dir_p,
    eng_pred_file  = eng_pred_file,
    FIG_DIR        = FIG_DIR,
    det_code       = det_code,
    det_name       = det_name,
    det_unit       = det_unit,
    trans_obj      = trans_obj,
    source_label   = source_lbl,
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
sub_dirs <- list.dirs(target_dir, recursive = TRUE, full.names = TRUE)
period_dirs <- sub_dirs[grepl("^(year_|season_|month_)", basename(sub_dirs))]
has_root_preds <- length(list.files(target_dir, pattern = "^predictions_.*\\.csv$")) > 0

dirs_to_plot <- character(0)
if (has_root_preds) dirs_to_plot <- c(dirs_to_plot, target_dir)
if (length(period_dirs) > 0) dirs_to_plot <- c(dirs_to_plot, period_dirs)

if (length(dirs_to_plot) == 0) {
  stop("No prediction files or period subdirectories found in ", target_dir)
}

cat(sprintf("Loading S3N predictions across %d directory(s) under %s...\n", length(dirs_to_plot), target_dir))
bundles <- lapply(dirs_to_plot, load_plot_bundle)

# Also check if target_dir is a single period subfolder
# and pool scale ranges with any sibling period folders in dirname(target_dir)
scale_vals <- unlist(lapply(bundles, function(b) b$eng_pred[[norm_meta$pred_col]]))
scale_se_e <- unlist(lapply(bundles, function(b) b$eng_pred$se_model_sc))
scale_se_h <- unlist(lapply(bundles, function(b) b$hum_pred$se_model_sc))

if (length(bundles) == 1 && grepl("^(year_|season_|month_)", basename(target_dir))) {
  parent_d <- dirname(target_dir)
  self_eng_p <- normalizePath(bundles[[1]]$eng_pred_file, mustWork = FALSE)
  sib_csvs <- list.files(parent_d, pattern = "^predictions_.*\\.csv$", recursive = TRUE, full.names = TRUE)
  sib_csvs <- sib_csvs[!grepl("^predictions_(humber|thames|severn|basin_|river_|catchid_)", basename(sib_csvs))]
  sib_csvs <- sib_csvs[normalizePath(sib_csvs, mustWork = FALSE) != self_eng_p]
  if (length(sib_csvs) > 0) {
    for (sc in sib_csvs) {
      sdf <- tryCatch(read.csv(sc, stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(sdf)) {
        if (norm_meta$mode != "none" && !(norm_meta$pred_col %in% names(sdf))) {
          sdf <- add_normalised_columns(sdf)
        }
        scale_vals <- c(scale_vals, sdf[[norm_meta$pred_col]])
        scale_se_e <- c(scale_se_e, sdf$se_model_sc)
      }
      hf <- file.path(dirname(sc), paste0("predictions_humber_", bundles[[1]]$det_code, ".csv"))
      if (file.exists(hf)) {
        hdf <- tryCatch(read.csv(hf, stringsAsFactors = FALSE), error = function(e) NULL)
        if (!is.null(hdf) && "se_model_sc" %in% names(hdf)) {
          scale_se_h <- c(scale_se_h, hdf$se_model_sc)
        }
      }
    }
  }
}

trans_shared   <- bundles[[1]]$trans_obj
shared_sc_val  <- if (norm_meta$mode == "none") {
  compute_color_scale(scale_vals, trans_shared, vmin_raw = vmin_raw, vmax_raw = vmax_raw)
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

catch_lu_plot <- NULL
if (length(custom_basins) > 0) {
  b1_net <- if ("Network_Abbr" %in% names(bundles[[1]]$cv_tab)) tolower(as.character(bundles[[1]]$cv_tab$Network_Abbr[1])) else "orn"
  proc_roots <- unique(c(dir_cfg$processed_dir, rel_path(dir_cfg$data_dir, "processed")))
  cands_lu <- unlist(lapply(proc_roots[dir.exists(proc_roots)], function(pr) {
    list.files(pr, pattern = "_catchments_lookup\\.csv$", recursive = TRUE, full.names = TRUE)
  }))
  if (length(cands_lu) > 0) {
    net_hits <- cands_lu[grepl(b1_net, cands_lu, fixed = TRUE)]
    lu_path <- if (length(net_hits) > 0) net_hits[1] else cands_lu[1]
    catch_lu_plot <- tryCatch(read.csv(lu_path, stringsAsFactors = FALSE), error = function(e) NULL)
  }
}

for (b in bundles) {
  ftag    <- norm_meta$file_tag
  f1_name <- sprintf("fig1_england_stream_%s_map%s.png", b$det_code, ftag)
  f2_name <- sprintf("fig2_humber_basin_%s_map%s.png", b$det_code, ftag)
  f3_name <- sprintf("fig3_thames_severn_%s_maps%s.png", b$det_code, ftag)
  f4_name <- sprintf("fig4_cv_and_diagnostics_%s.png", b$det_code)

  sc_lo_disp <- if (norm_meta$mode == "none") b$trans_obj$inv(shared_sc_val$vmin) else shared_sc_val$vmin
  sc_hi_disp <- if (norm_meta$mode == "none") b$trans_obj$inv(shared_sc_val$vmax) else shared_sc_val$vmax
  cat(sprintf("Rendering 4 figures in %s (shared scale [%s, %s]%s)...\n",
              b$FIG_DIR,
              formatC(sc_lo_disp, format = "fg", digits = 3),
              formatC(sc_hi_disp, format = "fg", digits = 3),
              if (norm_meta$mode == "none") "" else paste0(", norm = ", norm_meta$mode)))

  vars_sub <- format_covariates_subtitle(param_tab = b$param_tab)
  net_lbl  <- extract_network_label_from_tables(b$cv_tab, b$param_tab, b$FIG_DIR)

  plot_fig1_england_map(b$eng_pred, b$det_code, b$det_name, b$det_unit, b$trans_obj,
                        file.path(b$FIG_DIR, f1_name), reverse_color = reverse_color,
                        sc_val = shared_sc_val, sc_se = shared_sc_se_e,
                        vars_label = vars_sub, network_label = net_lbl,
                        source_label = b$source_label, norm_mode = norm_cli)

  plot_fig2_humber_map(b$hum_pred, b$cv_tab, shared_sc_val, b$det_code, b$det_name, b$det_unit, b$trans_obj,
                       file.path(b$FIG_DIR, f2_name), reverse_color = reverse_color,
                       sc_se = shared_sc_se_h, vars_label = vars_sub, network_label = net_lbl,
                       source_label = b$source_label, norm_mode = norm_cli)

  plot_fig3_thames_severn_map(b$tha_pred, b$sev_pred, shared_sc_val, b$det_code, b$det_name, b$det_unit, b$trans_obj,
                              file.path(b$FIG_DIR, f3_name), reverse_color = reverse_color,
                              vars_label = vars_sub, network_label = net_lbl,
                              source_label = b$source_label, norm_mode = norm_cli)

  plot_fig4_diagnostics(b$eng_pred, b$eng_cv, b$hum_cv, b$cv_tab, b$param_tab,
                        b$det_code, b$det_name, b$det_unit, b$trans_obj,
                        file.path(b$FIG_DIR, f4_name),
                        lims_a = if (length(bundles) > 1) lims_a_shared else NULL,
                        ylim_d = ylim_d_shared,
                        vars_label = vars_sub, network_label = net_lbl,
                        source_label = b$source_label)

  if (length(custom_basins) > 0) {
    for (bi in seq_along(custom_basins)) {
      b_spec <- custom_basins[bi]
      b_cids <- integer(0)
      b_title <- b_spec
      if (tolower(b_spec) %in% c("thames", "humber", "severn")) {
        b_res <- resolve_core_basin(b$eng_pred, catch_lu_plot, target = tolower(b_spec))
        b_cids <- b_res$catch_id
        b_title <- b_res$domain_label
      } else if (grepl("^[0-9]+$", b_spec)) {
        cid_num <- as.integer(b_spec)
        b_cids <- cid_num
        if (!is.null(catch_lu_plot) && any(catch_lu_plot$CatchID == cid_num)) {
          b_title <- sprintf("%s (CatchID=%d)", catch_lu_plot$Main_River[catch_lu_plot$CatchID == cid_num][1], cid_num)
        } else {
          b_title <- sprintf("CatchID %d Basin", cid_num)
        }
      } else if (!is.null(catch_lu_plot)) {
        if (grepl("^rbd:", tolower(b_spec))) {
          rid <- suppressWarnings(as.integer(sub("^rbd:", "", tolower(b_spec))))
          sub_r <- catch_lu_plot[catch_lu_plot$RBD_ID == rid, ]
          if (nrow(sub_r) > 0) {
            b_cids <- unique(sub_r$CatchID)
            b_title <- sprintf("%s RBD (RBD=%d)", sub_r$RBD_Name[1], rid)
          }
        } else {
          m_idx_lu <- grep(b_spec, catch_lu_plot$Main_River, ignore.case = TRUE)
          if (length(m_idx_lu) > 0) {
            sub_m <- catch_lu_plot[m_idx_lu, ]
            sub_m <- sub_m[order(-sub_m$n_obs_reaches, -sub_m$n_reaches), ]
            b_cids <- sub_m$CatchID[1]
            b_title <- sprintf("%s (CatchID=%d)", sub_m$Main_River[1], b_cids)
          } else {
            r_idx_lu <- grep(b_spec, catch_lu_plot$RBD_Name, ignore.case = TRUE)
            if (length(r_idx_lu) > 0) {
              rid <- catch_lu_plot$RBD_ID[r_idx_lu[1]]
              b_cids <- unique(catch_lu_plot$CatchID[catch_lu_plot$RBD_ID == rid])
              b_title <- sprintf("%s RBD (RBD=%d)", catch_lu_plot$RBD_Name[r_idx_lu[1]], rid)
            }
          }
        }
      }
      b_pred <- subset(b$eng_pred, CatchID %in% b_cids)
      if (nrow(b_pred) > 0) {
        b_slug <- gsub("[^a-z0-9]+", "_", tolower(b_title))
        b_slug <- gsub("^_|_$", "", b_slug)
        f5_b_name <- file.path(b$FIG_DIR, sprintf("fig%d_%s_%s_map%s.png", 4 + bi, b_slug, b$det_code, ftag))
        cat(sprintf("Rendering custom basin figure: %s (%d reaches)...\n", basename(f5_b_name), nrow(b_pred)))
        plot_basin_map(b_pred, b_title, b_title, b$cv_tab, shared_sc_val, b$det_code, b$det_name, b$det_unit, b$trans_obj,
                       f5_b_name, reverse_color = reverse_color, vars_label = vars_sub, network_label = net_lbl,
                       source_label = b$source_label, norm_mode = norm_cli)
      }
    }
  }

  if ("compliance_tier" %in% names(b$eng_pred) && "exceed_prob" %in% names(b$eng_pred)) {
    ex_t_val <- NA_real_
    ex_t_dir <- "above"
    ex_sum_csv <- find_first_csv(b$bundle_dir, "^exceedance_summary_.*\\.csv$", required = FALSE)
    if (!is.null(ex_sum_csv) && file.exists(ex_sum_csv)) {
      ex_df <- read.csv(ex_sum_csv, stringsAsFactors = FALSE)
      if (nrow(ex_df) > 0L && "Threshold_Raw" %in% names(ex_df)) {
        ex_t_val <- as.numeric(ex_df$Threshold_Raw[1])
        if ("Direction" %in% names(ex_df)) ex_t_dir <- as.character(ex_df$Direction[1])
      }
    }
    ex_eng_path <- file.path(b$FIG_DIR, sprintf("fig_exceedance_england_%s.png", b$det_code))
    ex_bas_path <- file.path(b$FIG_DIR, sprintf("fig_exceedance_basins_%s.png", b$det_code))
    plot_fig_exceedance_england(
      b$eng_pred, b$det_code, b$det_name, b$det_unit, ex_eng_path,
      threshold_val = ex_t_val, threshold_dir = ex_t_dir,
      vars_label = vars_sub, network_label = net_lbl, source_label = b$source_label
    )
    plot_fig_exceedance_basins(
      b$hum_pred, b$tha_pred, b$sev_pred, b$det_code, b$det_name, b$det_unit, ex_bas_path,
      threshold_val = ex_t_val, threshold_dir = ex_t_dir,
      vars_label = vars_sub, network_label = net_lbl, source_label = b$source_label
    )
    cat(sprintf("  + Saved regulatory exceedance & management tier maps: %s, %s\n",
                basename(ex_eng_path), basename(ex_bas_path)))
  }

  des_csv <- find_first_csv(b$bundle_dir, "^optimal_monitoring_design_.*\\.csv$", required = FALSE)
  if (!is.null(des_csv) && file.exists(des_csv)) {
    des_df <- tryCatch(read.csv(des_csv, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(des_df) && nrow(des_df) > 0L && "Domain" %in% names(des_df)) {
      des_dom_specs <- list(
        list(pred = b$eng_pred, pattern = "^All England", slug = "england", sc = shared_sc_se_e),
        list(pred = b$hum_pred, pattern = "^Humber",      slug = "humber",  sc = shared_sc_se_h),
        list(pred = b$tha_pred, pattern = "^Thames",      slug = "thames",  sc = NULL),
        list(pred = b$sev_pred, pattern = "^Severn",      slug = "severn",  sc = NULL)
      )
      for (ds in des_dom_specs) {
        sub_d <- des_df[grepl(ds$pattern, des_df$Domain, ignore.case = TRUE), , drop = FALSE]
        if (nrow(sub_d) > 0L && !is.null(ds$pred) && nrow(ds$pred) > 0L) {
          des_out <- file.path(b$FIG_DIR, sprintf("fig_optimal_design_%s_%s.png", ds$slug, b$det_code))
          plot_fig_optimal_design(
            ds$pred, sub_d, domain_label = sub_d$Domain[1],
            b$det_code, b$det_name, b$det_unit, b$trans_obj, des_out,
            sc_se = ds$sc, vars_label = vars_sub, network_label = net_lbl, source_label = b$source_label
          )
          cat(sprintf("  + Saved optimal monitoring design map: %s\n", basename(des_out)))
        }
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
      d_base <- basename(b$bundle_dir)
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
                                    b0$det_code, b0$det_name, b0$det_unit, b0$trans_obj,
                                    fig5_path, reverse_color = reverse_color,
                                    vars_label = vars_sub_g, network_label = net_lbl_g,
                                    source_label = b0$source_label, norm_mode = norm_cli)
  }
}

cat("Done!\n")
