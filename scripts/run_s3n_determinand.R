#!/usr/bin/env Rscript
# ==============================================================================
# run_s3n_determinand.R
#
# Reusable runner to fit the Scalable Spatial Stream Network (S3N) model
# (Tail-Up, Tail-Down, and Combined), run 5-fold out-of-sample cross-validation,
# predict across all 90,876 stream reaches (England, Humber, Thames, Severn),
# and generate the 4 publication-quality base-R figures for ANY water chemistry
# variable in RSN Water Quality Samples.csv.
#
# Usage:
#   Rscript scripts/run_s3n_determinand.R --list
#   Rscript scripts/run_s3n_determinand.R 0180          # Orthophosphate as P
#   Rscript scripts/run_s3n_determinand.R 0301          # Dissolved Organic Carbon (DOC)
#   Rscript scripts/run_s3n_determinand.R 0077          # Conductivity at 25 C
#   Rscript scripts/run_s3n_determinand.R 0061 --raw    # pH (modeled on raw scale)
# ==============================================================================

args     <- commandArgs(trailingOnly = TRUE)
BASE_DIR <- if (file.exists("scripts/s3n_plot_functions.R")) "." else if (file.exists("s3n_plot_functions.R")) ".." else "."
source(file.path(BASE_DIR, "scripts/s3n_plot_functions.R"))
source(file.path(BASE_DIR, "scripts/s3n_model_functions.R"))
BASE_DIR <- resolve_base_dir()

dir_cfg        <- parse_custom_dir_args(args, base_dir = BASE_DIR)
GIS_DIR        <- dir_cfg$gis_dir
DATA_DIR       <- dir_cfg$data_dir
OUT_BASE       <- dir_cfg$out_dir

# Custom sample / site arguments
custom_samples_arg <- grep("^--custom-samples=", args, value = TRUE)
custom_samples_cli <- if (length(custom_samples_arg) > 0) sub("^--custom-samples=", "", custom_samples_arg[1]) else NULL
custom_sites_arg   <- grep("^--custom-sites=", args, value = TRUE)
custom_sites_cli   <- if (length(custom_sites_arg) > 0) sub("^--custom-sites=", "", custom_sites_arg[1]) else NULL
custom_dir_arg     <- grep("^--custom-dir=", args, value = TRUE)
custom_dir_cli     <- if (length(custom_dir_arg) > 0) sub("^--custom-dir=", "", custom_dir_arg[1]) else NULL

src_cli        <- parse_source_arg(args)
source_cli     <- resolve_source_tag(src_cli, base_dir = BASE_DIR)
source_abbr    <- format_source_abbr(src_cli, base_dir = BASE_DIR)
source_label   <- format_source_label(src_cli, base_dir = BASE_DIR)
source_tag     <- format_source_folder_tag(source_cli)

src_dirs       <- resolve_data_source_dirs(DATA_DIR, source = src_cli, base_dir = BASE_DIR,
                                          custom_dir = custom_dir_cli,
                                          custom_samples = custom_samples_cli,
                                          custom_sites = custom_sites_cli)
RSN_DIR        <- src_dirs$rsn_dir
SSN_DIR        <- src_dirs$ssn_dir
n_omp_cli      <- parse_n_omp_arg(args, default = 1L)

if ("--list-sources" %in% args) {
  print_available_sources(base_dir = BASE_DIR)
  quit(save = "no", status = 0)
}

DATASET_FILES <- list(
  wq          = list(file = "Water Quality Samples.csv",      dt_col = "sample_datetime", label = "Water Quality Samples"),
  inverts     = list(file = "Macroinvertebrate Samples.csv",  dt_col = "sample_datetime", label = "Macroinvertebrate Samples"),
  diatoms     = list(file = "Diatom Samples.csv",             dt_col = "sample_datetime", label = "Diatom Samples"),
  macrophytes = list(file = "Macrophyte Samples.csv",         dt_col = "sample_datetime", label = "Macrophyte Samples"),
  rhs         = list(file = "River Habitat Surveys.csv",      dt_col = "survey_date",     label = "River Habitat Surveys"),
  lcms        = list(file = "LCMS GCMS Samples.csv",          dt_col = "sample_datetime", label = "LCMS/GCMS Organic Pollutants")
)

WIDE_VAR_META <- load_determinand_catalog(BASE_DIR)

ds_arg      <- grep("^--dataset=", args, value = TRUE)
ds_mode     <- if (length(ds_arg) > 0) tolower(sub("^--dataset=", "", ds_arg[1])) else NA_character_
net_arg_vec <- grep("^--network=", args, value = TRUE)
network_cli <- if (length(net_arg_vec) > 0) sub("^--network=", "", net_arg_vec[1]) else "orn"

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

if ("--list-transforms" %in% args || "--list-transformations" %in% args) {
  print_available_transforms()
  quit(save = "no", status = 0)
}

PROC_DIR       <- resolve_network_processed_dir(BASE_DIR, network_cli, data_dir = DATA_DIR, gis_dir = GIS_DIR, source = source_cli, processed_dir = dir_cfg$processed_dir)
net_meta       <- load_network_metadata(PROC_DIR)
REACHES_CSV    <- rel_path(PROC_DIR, "england_rsn_reaches.csv")
SITES_CSV      <- rel_path(PROC_DIR, "rsn_sites.csv")
CATCHMENTS_CSV <- rel_path(PROC_DIR, "rsn_catchments_lookup.csv")

if ("--list-basins" %in% args) {
  cat(sprintf("Available River Basin Districts (RBD) and Top Monitored Catchments [Network: %s (%s), Source: %s]:\n\n",
              net_meta$label, net_meta$abbr, source_abbr))
  c_lu <- read.csv(CATCHMENTS_CSV, stringsAsFactors = FALSE)
  rbd_tab <- aggregate(cbind(n_reaches, n_obs_reaches) ~ RBD_ID + RBD_Name, data = c_lu, FUN = sum)
  rbd_tab <- rbd_tab[order(-rbd_tab$n_obs_reaches), ]
  cat("=== 1. River Basin Districts (specify via --basin=\"<RBD_Name>\" or --basin=rbd:<RBD_ID>) ===\n")
  cat(sprintf("%-8s | %-22s | %-12s | %-12s\n", "RBD_ID", "RBD_Name", "Reaches", "Obs Reaches"))
  cat(strrep("-", 62), "\n")
  for (i in seq_len(nrow(rbd_tab))) {
    cat(sprintf("%-8d | %-22s | %-12d | %-12d\n",
                rbd_tab$RBD_ID[i], rbd_tab$RBD_Name[i], rbd_tab$n_reaches[i], rbd_tab$n_obs_reaches[i]))
  }
  cat("\n=== 2. Top 30 Monitored Catchments (specify via --basin=<CatchID> or --basin=\"<Main_River>\") ===\n")
  cat(sprintf("%-8s | %-28s | %-18s | %-10s | %-11s\n",
              "CatchID", "Main_River", "RBD_Name", "Reaches", "Obs Reaches"))
  cat(strrep("-", 86), "\n")
  top_c <- head(c_lu[order(-c_lu$n_obs_reaches, -c_lu$n_reaches), ], 30)
  for (i in seq_len(nrow(top_c))) {
    cat(sprintf("%-8d | %-28s | %-18s | %-10d | %-11d\n",
                top_c$CatchID[i], substr(top_c$Main_River[i], 1, 28),
                substr(top_c$RBD_Name[i], 1, 18), top_c$n_reaches[i], top_c$n_obs_reaches[i]))
  }
  quit(save = "no", status = 0)
}

if (length(args) == 0 || "--list" %in% args || "-h" %in% args || "--help" %in% args) {
  list_ds <- if (!is.na(ds_mode)) ds_mode else "wq"
  if (!(list_ds %in% names(DATASET_FILES))) {
    stop(sprintf("Unknown --dataset='%s'. Choose from: %s", list_ds, paste(names(DATASET_FILES), collapse = ", ")))
  }
  ds_info <- DATASET_FILES[[list_ds]]
  targets <- src_dirs$active_dirs
  cat(sprintf("Scanning %s (--dataset=%s, --source=%s [%s]) for available variables...\n\n",
              ds_info$label, list_ds, source_cli, source_label))

  if (list_ds == "wq") {
    dfs <- list()
    for (s_name in names(targets)) {
      f_path <- file.path(targets[[s_name]], paste0(toupper(s_name), " Water Quality Samples.csv"))
      if (file.exists(f_path)) {
        raw <- read.csv(f_path, stringsAsFactors = FALSE, colClasses = c(det_code = "character"), fileEncoding = "UTF-8-BOM")
        dfs[[s_name]] <- aggregate(network_id ~ det_code + det_desc + unit_short_desc, data = raw, FUN = length)
        rm(raw)
      }
    }
    tab_all <- do.call(rbind, dfs)
    tab <- aggregate(network_id ~ det_code + det_desc + unit_short_desc, data = tab_all, FUN = sum)
    tab <- tab[order(-tab$network_id), ]
    tab <- head(tab, 25)
    cat(sprintf("%-8s | %-10s | %-44s | %-10s\n", "det_code", "Samples", "Determinand Description", "Units"))
    cat(strrep("-", 82), "\n")
    for (i in seq_len(nrow(tab))) {
      cat(sprintf("%-8s | %-10d | %-44s | %-10s\n",
                  tab$det_code[i], tab$network_id[i], substr(tab$det_desc[i], 1, 44), tab$unit_short_desc[i]))
    }
  } else if (list_ds == "lcms") {
    if (source_cli == "ssn") {
      cat("LCMS/GCMS Organic Pollutants are only monitored on RSN and are not available in Small Streams Network (SSN).\n")
      quit(save = "no", status = 0)
    }
    csv_p <- file.path(RSN_DIR, "RSN LCMS GCMS Samples.csv")
    raw <- read.csv(csv_p, stringsAsFactors = FALSE, colClasses = c(cas_number = "character"), fileEncoding = "UTF-8-BOM")
    tab <- aggregate(network_id ~ cas_number + compound_name + unit_short_desc, data = raw, FUN = length)
    tab <- tab[order(-tab$network_id), ]
    tab <- head(tab, 25)
    cat(sprintf("%-12s | %-8s | %-42s | %-8s\n", "cas_number", "Samples", "Compound Name", "Units"))
    cat(strrep("-", 80), "\n")
    for (i in seq_len(nrow(tab))) {
      cat(sprintf("%-12s | %-8d | %-42s | %-8s\n",
                  tab$cas_number[i], tab$network_id[i], substr(tab$compound_name[i], 1, 42), tab$unit_short_desc[i]))
    }
  } else {
    dfs <- list()
    for (s_name in names(targets)) {
      f_path <- file.path(targets[[s_name]], paste0(toupper(s_name), " ", ds_info$file))
      if (file.exists(f_path)) {
        raw <- read.csv(f_path, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
        skip_c <- c("network_id", "biosys_id", "ecosys_id", "sample_id", "survey_id",
                    "sample_version", "analysis_id", "survey_format")
        cand_c <- setdiff(colnames(raw), skip_c)
        counts <- sapply(cand_c, function(cn) sum(!is.na(suppressWarnings(as.numeric(raw[[cn]])))))
        dfs[[s_name]] <- data.frame(Variable = names(counts), Count = as.integer(counts), stringsAsFactors = FALSE)
        rm(raw)
      }
    }
    tab_all <- do.call(rbind, dfs)
    tab <- aggregate(Count ~ Variable, data = tab_all, FUN = sum)
    tab <- tab[tab$Count >= 100, ]
    tab <- tab[order(-tab$Count), ]
    cat(sprintf("%-26s | %-8s | %-40s | %-8s\n", "Variable Name", "Samples", "Description", "Units"))
    cat(strrep("-", 90), "\n")
    for (i in seq_len(nrow(tab))) {
      cn <- tab$Variable[i]
      m_info <- if (cn %in% names(WIDE_VAR_META)) WIDE_VAR_META[[cn]] else c(desc = cn, unit = "index")
      cat(sprintf("%-26s | %-8d | %-40s | %-8s\n",
                  cn, tab$Count[i], substr(m_info["desc"], 1, 40), m_info["unit"]))
    }
  }
cat("Usage:\n  Rscript scripts/run_s3n_determinand.R <code_or_var> [--dataset=wq|inverts|diatoms|macrophytes|rhs|lcms]\n")
  cat("        [--source=rsn|ssn|rss] [--list-sources] [--transform=<type>] [--lambda=<val>] [--scale=<val>]\n")
  cat("        [--reverse-color] [--censor=ros|half_lod] [--boot=50] [--basin=<CatchID_or_Name>] [--list-basins]\n")
  cat("        [--network=<name>] [--list-networks] [--gis-dir=<dir>] [--data-dir=<dir>] [--output-dir=<dir>]\n")
  cat("        [--year=YYYY] [--season=winter|spring|summer|autumn] [--month=1..12] [--by=year|season|month]\n")
  cat("Examples:\n")
  cat("  Rscript scripts/run_s3n_determinand.R 0117                   # Nitrate on RSN (default)\n")
  cat("  Rscript scripts/run_s3n_determinand.R 0117 --source=ssn      # Nitrate on Small Streams Network (SSN)\n")
  cat("  Rscript scripts/run_s3n_determinand.R 0117 --source=rss      # Nitrate on Combined RSN + SSN (RSS)\n")
  cat("  Rscript scripts/run_s3n_determinand.R whpt_aspt --reverse-color # Macroinvertebrates ASPT\n")
  cat("  Rscript scripts/run_s3n_determinand.R 0117 --basin=Eden        # Nitrate + River Eden (CatchID=3125)\n")
  quit(save = "no", status = 0)
}

pos_args <- args[!grepl("^--", args)]
if (length(pos_args) == 0) {
  stop("Please provide a determinand code or variable name (e.g. 0117 or whpt_aspt). Run with --list to see options.")
}
det_code <- pos_args[1]

# Auto-detect dataset if --dataset was not explicitly provided
if (is.na(ds_mode)) {
  if (grepl("^[0-9]{4}$", det_code)) {
    ds_mode <- "wq"
  } else {
    found_ds <- NA_character_
    for (cand_ds in c("inverts", "diatoms", "macrophytes", "rhs")) {
      for (s_name in names(src_dirs$active_dirs)) {
        s_dir <- src_dirs$active_dirs[[s_name]]
        cand_p <- file.path(s_dir, paste0(toupper(s_name), " ", DATASET_FILES[[cand_ds]]$file))
        if (!file.exists(cand_p)) cand_p <- file.path(s_dir, DATASET_FILES[[cand_ds]]$file)
        if (file.exists(cand_p)) {
          hdr <- colnames(read.csv(cand_p, nrows = 1, fileEncoding = "UTF-8-BOM"))
          if (det_code %in% hdr) {
            found_ds <- cand_ds
            break
          }
        }
      }
      if (!is.na(found_ds)) break
    }
    ds_mode <- if (!is.na(found_ds)) found_ds else "wq"
  }
}
if (!(ds_mode %in% names(DATASET_FILES))) {
  stop(sprintf("Invalid --dataset='%s'. Choose from: %s", ds_mode, paste(names(DATASET_FILES), collapse = ", ")))
}

# Default scale & transformation:
# --transform=<type> (log, sqrt, cbrt, fourthrt, log1p, log10, reciprocal, boxcox, yeojohnson, asin, logit, none/raw)
# --lambda=<val> (optional Box-Cox / Yeo-Johnson power parameter override; auto-estimated by default)
trans_arg  <- grep("^--(transform|trans|transformation)=", args, value = TRUE)
trans_cli  <- if (length(trans_arg) > 0) tolower(sub("^--(transform|trans|transformation)=", "", trans_arg[1])) else NA_character_
lambda_arg <- grep("^--(lambda|transform-lambda)=", args, value = TRUE)
lambda_cli <- if (length(lambda_arg) > 0) as.numeric(sub("^--(lambda|transform-lambda)=", "", lambda_arg[1])) else NULL
scale_arg  <- grep("^--(scale|transform-scale)=", args, value = TRUE)
scale_cli  <- if (length(scale_arg) > 0) as.numeric(sub("^--(scale|transform-scale)=", "", scale_arg[1])) else NULL

if (is.na(trans_cli)) {
  trans_cli <- if (ds_mode %in% c("wq", "lcms")) "log" else "none"
}

trans_obj <- build_transformer(trans_cli, lambda = lambda_cli, scale = scale_cli)

reverse_color <- any(args %in% c("--reverse-color", "--reverse-colour", "--rev-color", "--rev-colour", "--reverse-scale"))
vmin_arg   <- grep("^--vmin=", args, value = TRUE)
vmin_raw   <- if (length(vmin_arg) > 0) as.numeric(sub("^--vmin=", "", vmin_arg[1])) else NULL
vmax_arg   <- grep("^--vmax=", args, value = TRUE)
vmax_raw   <- if (length(vmax_arg) > 0) as.numeric(sub("^--vmax=", "", vmax_arg[1])) else NULL
boot_arg   <- grep("^--boot=", args, value = TRUE)
n_boot     <- if (length(boot_arg) > 0) as.integer(sub("^--boot=", "", boot_arg[1])) else 50L
censor_arg <- grep("^--censor=", args, value = TRUE)
censor_m   <- if (length(censor_arg) > 0) sub("^--censor=", "", censor_arg[1]) else "ros"
year_arg   <- grep("^--year=", args, value = TRUE)
year_f     <- if (length(year_arg) > 0) as.integer(sub("^--year=", "", year_arg[1])) else NA_integer_
season_arg <- grep("^--season=", args, value = TRUE)
season_f   <- if (length(season_arg) > 0) tolower(sub("^--season=", "", season_arg[1])) else NA_character_
if (!is.na(season_f) && season_f == "fall") season_f <- "autumn"
month_arg  <- grep("^--month=", args, value = TRUE)
month_f    <- if (length(month_arg) > 0) as.integer(sub("^--month=", "", month_arg[1])) else NA_integer_
by_arg     <- grep("^--by=", args, value = TRUE)
by_mode    <- if (length(by_arg) > 0) tolower(sub("^--by=", "", by_arg[1])) else NA_character_
basin_args <- grep("^--basin=", args, value = TRUE)
basin_specs <- if (length(basin_args) > 0) {
  unlist(strsplit(sub("^--basin=", "", basin_args), ",", fixed = TRUE))
} else {
  character(0)
}
basin_specs <- trimws(basin_specs[nchar(trimws(basin_specs)) > 0])

vars_args <- grep("^--(vars|covariates)=", args, value = TRUE)
var_specs <- if (length(vars_args) > 0) {
  unlist(strsplit(sub("^--(vars|covariates)=", "", vars_args), ",", fixed = TRUE))
} else {
  character(0)
}
var_specs <- trimws(var_specs[nchar(trimws(var_specs)) > 0])

load_s3n_backend(BASE_DIR)
meas_error_arg <- any(c("--meas-error", "--gam") %in% args)
norm_cli       <- parse_norm_mode(args, default = "none")
norm_meta      <- get_norm_meta(norm_cli)

# Range mode handling and folder tagging
range_arg <- grep("^--range-mode=", args, value = TRUE)
range_cli <- if (length(range_arg) > 0) tolower(sub("^--range-mode=", "", range_arg[1])) else NA_character_

src_tokens_init <- unlist(strsplit(source_cli, "[+,]+"))
src_tokens_init <- trimws(src_tokens_init[nzchar(trimws(src_tokens_init))])

detected_range_mode <- NA_character_
if (!is.na(range_cli) && !(range_cli %in% c("auto", "none"))) {
  detected_range_mode <- range_cli
} else if (!is.na(range_cli) && range_cli == "none") {
  detected_range_mode <- NA_character_
} else {
  # 1. Suffix in source_cli or any '+' constituent token
  for (rm in c("tobit", "geometric", "midpoint", "lower", "upper")) {
    if (any(grepl(paste0("_", rm, "$"), src_tokens_init))) {
      detected_range_mode <- rm
      break
    }
  }
  # 2. Network metadata in PROC_DIR
  if (is.na(detected_range_mode) && !is.null(net_meta$range_mode) && nzchar(net_meta$range_mode)) {
    detected_range_mode <- net_meta$range_mode
  }
  # 3. Source config registry (single or composite constituent)
  if (is.na(detected_range_mode)) {
    reg <- load_source_registry(BASE_DIR)
    for (tok in src_tokens_init) {
      if (tok %in% names(reg$sources)) {
        rm_cfg <- reg$sources[[tok]]$range_mode
        if (!is.null(rm_cfg) && nzchar(rm_cfg)) {
          detected_range_mode <- rm_cfg
          break
        }
      }
    }
  }
  # 4. Check meas_method in samples.csv
  if (is.na(detected_range_mode)) {
    cand_smp <- custom_samples_cli
    if (is.null(cand_smp) && !is.null(custom_dir_cli)) {
      cand_smp <- file.path(custom_dir_cli, "samples.csv")
    }
    if (is.null(cand_smp)) {
      for (tok in names(src_dirs$active_dirs)) {
        s_d <- src_dirs$active_dirs[[tok]]
        cand_try <- file.path(s_d, "samples.csv")
        if (file.exists(cand_try)) {
          cand_smp <- cand_try
          break
        }
      }
    }
    if (!is.null(cand_smp) && file.exists(cand_smp)) {
      tryCatch({
        samp_hdr <- colnames(read.csv(cand_smp, nrows = 1, fileEncoding = "UTF-8-BOM"))
        if ("meas_method" %in% samp_hdr) {
          samp_preview <- read.csv(cand_smp, nrows = 200, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
          mm_vals <- tolower(unique(as.character(samp_preview$meas_method)))
          mm_vals <- mm_vals[!is.na(mm_vals) & nzchar(mm_vals)]
          if (any(grepl("tobit", mm_vals))) detected_range_mode <- "tobit"
          else if (any(grepl("geometric", mm_vals))) detected_range_mode <- "geometric"
          else if (any(grepl("midpoint|arithmetic", mm_vals))) detected_range_mode <- "midpoint"
          else if (any(grepl("lower", mm_vals))) detected_range_mode <- "lower"
          else if (any(grepl("upper", mm_vals))) detected_range_mode <- "upper"
        }
      }, error = function(e) NULL)
    }
  }
}

if (!is.na(detected_range_mode)) {
  has_token_mode <- any(grepl(paste0("_", detected_range_mode, "$"), src_tokens_init))
  if (!has_token_mode) {
    if (nzchar(source_tag) && !grepl(paste0("_", detected_range_mode, "$"), source_tag)) {
      source_tag <- paste0(source_tag, "_", detected_range_mode)
    } else if (!nzchar(source_tag)) {
      source_tag <- detected_range_mode
    }
    if (!grepl("[+,]", source_cli) && !grepl(paste0("_", detected_range_mode, "$"), source_cli)) {
      source_cli <- paste0(source_cli, "_", detected_range_mode)
    }
  }
  source_abbr  <- format_source_abbr(source_cli, base_dir = BASE_DIR)
  source_label <- format_source_label(source_cli, base_dir = BASE_DIR)
  PROC_DIR_NEW <- resolve_network_processed_dir(BASE_DIR, network_cli, data_dir = DATA_DIR, gis_dir = GIS_DIR, source = source_cli, processed_dir = dir_cfg$processed_dir)
  if (!identical(PROC_DIR_NEW, PROC_DIR)) {
    PROC_DIR       <- PROC_DIR_NEW
    net_meta       <- load_network_metadata(PROC_DIR)
    REACHES_CSV    <- rel_path(PROC_DIR, "england_rsn_reaches.csv")
    SITES_CSV      <- rel_path(PROC_DIR, "rsn_sites.csv")
    CATCHMENTS_CSV <- rel_path(PROC_DIR, "rsn_catchments_lookup.csv")
  }
}

if (!is.na(detected_range_mode)) {
  cat(sprintf("   [Range Handling] Active range mode: %s (Output folder tagged with _%s)\n",
              toupper(detected_range_mode), detected_range_mode))
}

# ==============================================================================
# DATA LOADING & MULTI-DOMAIN S3N ESTIMATION + 5-FOLD CV
# ==============================================================================

# Robust Lognormal Regression on Order Statistics (ROS, Helsel & Cohn 1988) in Base R
impute_censored_ros <- function(vals, is_cens) {
  if (!any(is_cens)) return(vals)
  n <- length(vals)
  c_limits <- sort(unique(vals[is_cens]))
  J <- length(c_limits)
  min_thr <- min(0.0, min(vals, na.rm = TRUE) - 1.0, na.rm = TRUE)
  thresholds <- c(min_thr, c_limits, Inf)

  A <- vector("list", J + 1)
  C <- vector("list", J + 1)
  for (i in seq_len(n)) {
    if (is_cens[i]) {
      j <- match(vals[i], c_limits) + 1L
      C[[j]] <- c(C[[j]], i)
    } else {
      j_find <- max(which(thresholds[1:(J + 1)] <= vals[i]))
      A[[j_find]] <- c(A[[j_find]], i)
    }
  }
  for (j in seq_len(J + 1)) {
    if (length(A[[j]]) > 1) A[[j]] <- A[[j]][order(vals[A[[j]]])]
  }

  p <- numeric(J + 2)
  p[J + 2] <- 0.0
  for (j in seq(J + 1, 2, by = -1)) {
    d_j <- thresholds[j]
    A_j <- length(A[[j]])
    B_j <- sum((is_cens & vals <= d_j) | (!is_cens & vals < d_j))
    p[j] <- if (A_j + B_j > 0) p[j + 1] + (1.0 - p[j + 1]) * (A_j / (A_j + B_j)) else p[j + 1]
  }
  p[1] <- if (length(A[[1]]) > 0) 1.0 else p[2]

  pp <- numeric(n)
  for (j in seq_len(J + 1)) {
    A_j <- length(A[[j]])
    if (A_j > 0) {
      # Average plotting positions across exact ties so tie order never affects ROS fit
      ranks_a <- ave(seq_len(A_j), vals[A[[j]]], FUN = mean)
      pp[A[[j]]] <- (1.0 - p[j]) + (p[j] - p[j + 1]) * (ranks_a / (A_j + 1.0))
    }
  }
  for (j in 2:(J + 1)) {
    C_j <- length(C[[j]])
    if (C_j > 0) {
      pp[C[[j]]] <- (1.0 - p[j]) * ((seq_len(C_j)) / (C_j + 1.0))
    }
  }
  pp <- pmin(1.0 - 1e-6, pmax(1e-6, pp))

  unc_idx <- which(!is_cens & vals > 0)
  if (length(unc_idx) < 3) {
    warning("Fewer than 3 positive uncensored samples for ROS; falling back to 0.5 * LOD.")
    return(ifelse(is_cens, 0.5 * vals, vals))
  }
  z_unc <- qnorm(pp[unc_idx])
  l_unc <- log(vals[unc_idx])
  slope <- max(0.05, cov(z_unc, l_unc) / max(1e-9, var(z_unc)))
  inter <- mean(l_unc) - slope * mean(z_unc)

  out <- vals
  c_idx <- which(is_cens)
  out[c_idx] <- pmin(vals[c_idx], pmax(1e-6, exp(inter + slope * qnorm(pp[c_idx]))))
  out
}

cat(sprintf("======================================================================\n"))
cat(sprintf("  Running S3N Pipeline: dataset = '%s', variable = '%s' (transform = %s, censor = %s)\n",
            ds_mode, det_code, trans_obj$label, censor_m))
cat(sprintf("======================================================================\n"))

ds_info <- DATASET_FILES[[ds_mode]]

find_first_col <- function(cols, candidates) {
  cols_low <- tolower(trimws(cols))
  for (cand in candidates) {
    hit <- which(cols_low == tolower(trimws(cand)))
    if (length(hit) > 0) return(cols[hit[1]])
  }
  NULL
}

read_wq_samples_file <- function(f_path, source_label_tag, det_code, allow_substr_det = FALSE) {
  if (!file.exists(f_path)) return(NULL)
  raw <- read.csv(f_path, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  fnames <- colnames(raw)
  id_c   <- find_first_col(fnames, c("network_id", "site_id", "station_id", "id", "location_id", "site"))
  dt_c   <- find_first_col(fnames, c("sample_datetime", "datetime", "date_time", "date", "survey_date", "timestamp"))
  val_c  <- find_first_col(fnames, c("meas_result", "result", "value", "val", "reading", "score", "concentration"))
  det_c  <- find_first_col(fnames, c("det_code", "determinand_code", "param_code", "cas_number", "determinand", "code", "metric_name"))
  desc_c <- find_first_col(fnames, c("det_desc", "det_name", "compound_name", "description", "determinand"))
  sign_c <- find_first_col(fnames, c("meas_sign", "less_than_flag", "censored", "is_censored", "qualifier"))
  unit_c <- find_first_col(fnames, c("unit_short_desc", "unit", "units"))

  if (is.null(id_c)) id_c <- fnames[1]
  keep_rows <- rep(FALSE, nrow(raw))
  if (!is.null(det_c)) {
    det_vals  <- as.character(raw[[det_c]])
    det_clean <- sub("^0+", "", det_vals)
    tgt_clean <- sub("^0+", "", det_code)
    keep_rows <- (det_vals == det_code) | (det_clean == tgt_clean) |
      grepl(paste0("^", det_code, "$"), det_vals, ignore.case = TRUE)
    if (allow_substr_det) {
      keep_rows <- keep_rows | grepl(det_code, det_vals, fixed = TRUE)
    }
  } else if (det_code %in% fnames) {
    val_c <- det_code
    keep_rows <- !is.na(raw[[det_code]])
  } else if (!is.null(val_c)) {
    keep_rows <- !is.na(raw[[val_c]])
  }

  sub_raw <- raw[keep_rows, , drop = FALSE]
  rm(raw)
  if (nrow(sub_raw) == 0) return(NULL)
  unc_c  <- find_first_col(fnames, c("meas_uncertainty", "uncertainty", "se", "std_err", "standard_error", "sigma_error"))
  meth_c <- find_first_col(fnames, c("meas_method", "method", "model_type", "type"))
  data.frame(
    network_id       = as.character(sub_raw[[id_c]]),
    sample_datetime  = if (!is.null(dt_c)) as.character(sub_raw[[dt_c]]) else rep("2023-06-01", nrow(sub_raw)),
    det_desc         = if (!is.null(desc_c)) as.character(sub_raw[[desc_c]]) else rep(det_code, nrow(sub_raw)),
    meas_sign        = if (!is.null(sign_c)) as.character(sub_raw[[sign_c]]) else rep("=", nrow(sub_raw)),
    meas_result      = as.character(sub_raw[[val_c]]),
    meas_uncertainty = if (!is.null(unc_c)) suppressWarnings(as.numeric(sub_raw[[unc_c]])) else rep(NA_real_, nrow(sub_raw)),
    meas_method      = if (!is.null(meth_c)) tolower(trimws(as.character(sub_raw[[meth_c]]))) else rep("", nrow(sub_raw)),
    unit_short_desc  = if (!is.null(unit_c)) as.character(sub_raw[[unit_c]]) else rep("mg/L", nrow(sub_raw)),
    data_source      = source_label_tag,
    stringsAsFactors = FALSE
  )
}

if (ds_mode == "wq") {
  sample_dfs <- list()
  for (s_name in names(src_dirs$active_dirs)) {
    s_dir  <- src_dirs$active_dirs[[s_name]]
    f_path <- file.path(s_dir, paste0(toupper(s_name), " Water Quality Samples.csv"))
    if (!file.exists(f_path)) f_path <- file.path(s_dir, "samples.csv")
    if (!file.exists(f_path)) f_path <- file.path(s_dir, "Water Quality Samples.csv")
    if (!file.exists(f_path) && !is.null(custom_samples_cli) && file.exists(custom_samples_cli)) {
      f_path <- custom_samples_cli
    }
    if (!file.exists(f_path) && dir.exists(s_dir)) {
      cands <- list.files(s_dir, pattern = "(sample|wq).*\\.csv$", ignore.case = TRUE, full.names = TRUE)
      if (length(cands) > 0) f_path <- cands[1]
    }
    sub_df <- read_wq_samples_file(f_path, toupper(s_name), det_code, allow_substr_det = FALSE)
    if (!is.null(sub_df)) sample_dfs[[s_name]] <- sub_df
  }

  # Direct custom samples injection if specified via CLI
  if (!is.null(custom_samples_cli) && file.exists(custom_samples_cli) && !"custom" %in% names(sample_dfs)) {
    c_sub <- read_wq_samples_file(custom_samples_cli, "CUSTOM", det_code, allow_substr_det = TRUE)
    if (!is.null(c_sub)) sample_dfs[["custom"]] <- c_sub
  }

  if (length(sample_dfs) == 0) {
    stop(sprintf("No water quality samples found for det_code = '%s' in %s. Run with --list to see valid codes.",
                 det_code, source_label))
  }
  det_samples <- do.call(rbind, sample_dfs)
  det_name <- det_samples$det_desc[1]
  det_unit <- det_samples$unit_short_desc[1]
  raw_vals <- suppressWarnings(as.numeric(det_samples$meas_result))
  is_cens  <- (det_samples$meas_sign %in% c("<", "yes", "true", "1")) & !grepl("tobit", det_samples$meas_method)
} else if (ds_mode == "lcms") {
  if (source_cli == "ssn") {
    stop("LCMS/GCMS Organic Pollutants are only monitored on RSN and are not available in Small Streams Network (SSN). Use --source=rsn or --source=rss.")
  }
  csv_p <- file.path(RSN_DIR, "RSN LCMS GCMS Samples.csv")
  all_samples <- read.csv(csv_p, stringsAsFactors = FALSE, colClasses = c(cas_number = "character"), fileEncoding = "UTF-8-BOM")
  keep_m      <- (all_samples$cas_number == det_code) | grepl(det_code, all_samples$compound_name, ignore.case = TRUE)
  det_samples <- all_samples[keep_m, c("network_id", "sample_datetime", "compound_name", "less_than_flag", "meas_result", "unit_short_desc")]
  rm(all_samples)
  if (nrow(det_samples) == 0) {
    stop(sprintf("No LCMS/GCMS samples found for '%s'. Run with --list --dataset=lcms to see valid compounds.", det_code))
  }
  det_samples$data_source <- "RSN"
  det_samples$meas_uncertainty <- NA_real_
  det_samples$meas_method <- ""
  det_name <- det_samples$compound_name[1]
  det_unit <- det_samples$unit_short_desc[1]
  raw_vals <- suppressWarnings(as.numeric(det_samples$meas_result))
  is_cens  <- (tolower(trimws(det_samples$less_than_flag)) == "yes")
} else {
  if (ds_mode == "rhs" && source_cli == "ssn") {
    stop("River Habitat Surveys (RHS) are only surveyed on RSN and are not available in Small Streams Network (SSN). Use --source=rsn or --source=rss.")
  }
  sample_dfs <- list()
  for (s_name in names(src_dirs$active_dirs)) {
    s_dir  <- src_dirs$active_dirs[[s_name]]
    f_path <- file.path(s_dir, paste0(toupper(s_name), " ", ds_info$file))
    if (!file.exists(f_path)) f_path <- file.path(s_dir, ds_info$file)
    if (!file.exists(f_path)) f_path <- file.path(s_dir, "samples.csv")
    if (!file.exists(f_path) && !is.null(custom_samples_cli) && file.exists(custom_samples_cli)) {
      f_path <- custom_samples_cli
    }

    if (file.exists(f_path)) {
      raw <- read.csv(f_path, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
      fnames <- colnames(raw)
      id_c   <- find_first_col(fnames, c("network_id", "site_id", "station_id", "id"))
      dt_c   <- find_first_col(fnames, c(ds_info$dt_col, "sample_datetime", "datetime", "date_time", "date", "survey_date"))

      if (is.null(id_c)) id_c <- fnames[1]

      if (det_code %in% fnames) {
        sub_df <- data.frame(
          network_id      = as.character(raw[[id_c]]),
          sample_datetime = if (!is.null(dt_c)) as.character(raw[[dt_c]]) else rep("2023-06-01", nrow(raw)),
          meas_result     = as.character(raw[[det_code]]),
          data_source     = toupper(s_name),
          stringsAsFactors = FALSE
        )
        sample_dfs[[s_name]] <- sub_df
      }
      rm(raw)
    }
  }
  if (length(sample_dfs) == 0) {
    stop(sprintf("Column '%s' not found in %s across %s. Run with --list --dataset=%s to see available variables.",
                 det_code, ds_info$file, source_label, ds_mode))
  }
  det_samples <- do.call(rbind, sample_dfs)
  m_info   <- if (det_code %in% names(WIDE_VAR_META)) WIDE_VAR_META[[det_code]] else c(short = det_code, desc = det_code, unit = "index")
  det_name <- unname(m_info["short"])
  det_unit <- unname(m_info["unit"])
  raw_vals <- suppressWarnings(as.numeric(det_samples$meas_result))
  is_cens  <- rep(FALSE, length(raw_vals))
}
if (!"meas_uncertainty" %in% names(det_samples)) {
  det_samples$meas_uncertainty <- NA_real_
}
if (!"meas_method" %in% names(det_samples)) {
  det_samples$meas_method <- rep("", nrow(det_samples))
} else {
  det_samples$meas_method <- ifelse(is.na(det_samples$meas_method), "", as.character(det_samples$meas_method))
}

valid_k     <- !is.na(raw_vals) & !is.na(det_samples$sample_datetime) & nchar(trimws(det_samples$sample_datetime)) >= 7
det_samples <- det_samples[valid_k, ]
raw_vals    <- raw_vals[valid_k]
is_cens     <- is_cens[valid_k]

if (any(is_cens) && censor_m == "ros") {
  det_samples$val <- impute_censored_ros(raw_vals, is_cens)
} else if (any(is_cens)) {
  det_samples$val <- ifelse(is_cens, 0.5 * raw_vals, raw_vals)
} else {
  det_samples$val <- raw_vals
}

valid_t     <- trans_obj$valid_fn(det_samples$val)
det_samples <- det_samples[valid_t, ]
is_cens     <- is_cens[valid_t]

# Optimize profile Maximum Likelihood lambda for Box-Cox / Yeo-Johnson or scale for asin / logit AFTER censoring imputation
if ((trans_obj$type %in% c("boxcox", "yeojohnson") && is.na(trans_obj$lambda)) ||
    (is.null(scale_cli) && trans_obj$type %in% c("asin", "logit"))) {
  trans_obj <- build_transformer(trans_cli, y_obs = det_samples$val, lambda = lambda_cli, scale = scale_cli)
  if (trans_obj$type %in% c("boxcox", "yeojohnson")) {
    cat(sprintf("  Estimated Profile Maximum Likelihood lambda = %.4f for %s (%s)\n",
                trans_obj$lambda, det_name, trans_obj$label))
  }
}

default_trans   <- if (ds_mode %in% c("wq", "lcms")) "log" else "none"
trans_tag       <- format_transform_folder_tag(trans_obj, default_type = default_trans)
vars_tag        <- format_covariates_folder_tag(var_specs)
out_folder_name <- if (ds_mode == "wq") paste0("det_", det_code) else paste0(ds_mode, "_", gsub("[^a-zA-Z0-9_]+", "_", det_code))
out_folder_name <- paste0(out_folder_name, "_", net_meta$abbr)
if (nzchar(source_tag)) {
  out_folder_name <- paste0(out_folder_name, "_", source_tag)
}
if (nzchar(trans_tag)) {
  out_folder_name <- paste0(out_folder_name, "_", trans_tag)
}
if (nzchar(vars_tag)) {
  out_folder_name <- paste0(out_folder_name, "_", vars_tag)
}
OUT_DIR <- rel_path(OUT_BASE, "results", out_folder_name)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

det_samples$year   <- as.integer(substr(det_samples$sample_datetime, 1, 4))
det_samples$month  <- as.integer(substr(det_samples$sample_datetime, 6, 7))
season_lu          <- c("winter", "winter", "spring", "spring", "spring",
                        "summer", "summer", "summer", "autumn", "autumn", "autumn", "winter")
det_samples$season <- season_lu[det_samples$month]

cat(sprintf("Extracted %d total valid samples (%d censored '<' [%.2f%%] via '%s') for '%s' (%s) across %d sites.\n",
            nrow(det_samples), sum(is_cens), 100 * mean(is_cens), censor_m,
            det_name, det_unit, length(unique(det_samples$network_id))))

sites_lu   <- read.csv(SITES_CSV, stringsAsFactors = FALSE)
sites_lu   <- sites_lu[sites_lu$included_250m == 1, c("site_id", "RivID", "obs_up_dist_km")]
df_reaches <- read.csv(REACHES_CSV, stringsAsFactors = FALSE)

# Fitting/CV engine (load_s3n_backend, parse_n_omp_arg, fit_s3n_domain) lives in scripts/s3n_model_functions.R

catch_lookup <- if (file.exists(CATCHMENTS_CSV)) {
  read.csv(CATCHMENTS_CSV, stringsAsFactors = FALSE)
} else {
  NULL
}

resolve_basin_spec <- function(spec, c_lu) {
  spec_trim <- trimws(spec)
  # 1. Explicit rbd:<id>
  if (grepl("^rbd:", tolower(spec_trim))) {
    rid <- suppressWarnings(as.integer(sub("^rbd:", "", tolower(spec_trim))))
    sub_r <- c_lu[c_lu$RBD_ID == rid, ]
    if (nrow(sub_r) == 0) stop(sprintf("No River Basin District found for '%s'. Run --list-basins.", spec))
    rname <- sub_r$RBD_Name[1]
    slug  <- gsub("[^a-z0-9]+", "_", tolower(paste0("rbd_", rname)))
    return(list(
      slug         = gsub("^_|_$", "", slug),
      short_title  = paste0(rname, " RBD"),
      domain_label = sprintf("%s RBD (RBD=%d)", rname, rid),
      catch_ids    = unique(sub_r$CatchID)
    ))
  }
  # 2. Integer CatchID
  if (grepl("^[0-9]+$", spec_trim)) {
    cid <- as.integer(spec_trim)
    sub_c <- if (!is.null(c_lu)) c_lu[c_lu$CatchID == cid, ] else data.frame()
    rname <- if (nrow(sub_c) > 0) sub_c$Main_River[1] else paste0("Catchment ", cid)
    slug  <- gsub("[^a-z0-9]+", "_", tolower(paste0(rname, "_", cid)))
    return(list(
      slug         = gsub("^_|_$", "", slug),
      short_title  = rname,
      domain_label = sprintf("%s (CatchID=%d)", rname, cid),
      catch_ids    = cid
    ))
  }
  if (is.null(c_lu)) stop("Catchment lookup CSV not found; specify numeric CatchID.")
  # 3. Exact or substring match on RBD_Name (if not matching a specific river first, or if exact RBD match)
  rbd_exact <- which(tolower(c_lu$RBD_Name) == tolower(spec_trim))
  if (length(rbd_exact) > 0 && !any(tolower(c_lu$Main_River) == tolower(spec_trim))) {
    rid   <- c_lu$RBD_ID[rbd_exact[1]]
    rname <- c_lu$RBD_Name[rbd_exact[1]]
    slug  <- gsub("[^a-z0-9]+", "_", tolower(paste0("rbd_", rname)))
    return(list(
      slug         = gsub("^_|_$", "", slug),
      short_title  = paste0(rname, " RBD"),
      domain_label = sprintf("%s RBD (RBD=%d)", rname, rid),
      catch_ids    = unique(c_lu$CatchID[c_lu$RBD_ID == rid])
    ))
  }
  # 4. Substring match on Main_River (preferring highest n_obs_reaches)
  m_idx <- grep(spec_trim, c_lu$Main_River, ignore.case = TRUE)
  if (length(m_idx) > 0) {
    sub_m <- c_lu[m_idx, ]
    sub_m <- sub_m[order(-sub_m$n_obs_reaches, -sub_m$n_reaches), ]
    if (nrow(sub_m) > 1) {
      cat(sprintf("  Note: '--basin=%s' matched %d catchments; selecting largest monitored: %s (CatchID=%d, %d obs reaches).\n",
                  spec_trim, nrow(sub_m), sub_m$Main_River[1], sub_m$CatchID[1], sub_m$n_obs_reaches[1]))
    }
    cid   <- sub_m$CatchID[1]
    rname <- sub_m$Main_River[1]
    slug  <- gsub("[^a-z0-9]+", "_", tolower(paste0(rname, "_", cid)))
    return(list(
      slug         = gsub("^_|_$", "", slug),
      short_title  = rname,
      domain_label = sprintf("%s (CatchID=%d)", rname, cid),
      catch_ids    = cid
    ))
  }
  # 5. Substring match on RBD_Name
  r_idx <- grep(spec_trim, c_lu$RBD_Name, ignore.case = TRUE)
  if (length(r_idx) > 0) {
    rid   <- c_lu$RBD_ID[r_idx[1]]
    rname <- c_lu$RBD_Name[r_idx[1]]
    slug  <- gsub("[^a-z0-9]+", "_", tolower(paste0("rbd_", rname)))
    return(list(
      slug         = gsub("^_|_$", "", slug),
      short_title  = paste0(rname, " RBD"),
      domain_label = sprintf("%s RBD (RBD=%d)", rname, rid),
      catch_ids    = unique(c_lu$CatchID[c_lu$RBD_ID == rid])
    ))
  }
  stop(sprintf("Could not match '--basin=%s' to any CatchID, Main_River, or RBD_Name. Run with --list-basins to see options.", spec))
}

extra_basins <- lapply(basin_specs, resolve_basin_spec, c_lu = catch_lookup)

slice_ctx <- list(
  sites_lu       = sites_lu,
  df_reaches     = df_reaches,
  trans_obj      = trans_obj,
  var_specs      = var_specs,
  meas_error_arg = meas_error_arg,
  net_meta       = net_meta,
  catch_lookup   = catch_lookup,
  source_cli     = source_cli,
  source_abbr    = source_abbr,
  source_label   = source_label,
  n_boot         = n_boot,
  n_omp          = n_omp_cli,
  extra_basins   = extra_basins,
  det_code       = det_code,
  det_name       = det_name,
  det_unit       = det_unit
)

fit_temporal_slice <- function(sub_samples, out_subdir, period_tag, ctx) {
  dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("\n======================================================================\n"))
  cat(sprintf("  Temporal Period: %s (%d samples across %d stations)\n",
              period_tag, nrow(sub_samples), length(unique(sub_samples$network_id))))
  cat(sprintf("======================================================================\n"))

  det_WithRiv <- merge(sub_samples, ctx$sites_lu, by.x = "network_id", by.y = "site_id")
  if (nrow(det_WithRiv) == 0) {
    cat(sprintf("  No snapped samples found for period '%s'; skipping.\n", period_tag))
    return(NULL)
  }
  det_WithRiv$y_sample <- ctx$trans_obj$fwd(det_WithRiv$val)
  riv_groups <- split(seq_len(nrow(det_WithRiv)), det_WithRiv$RivID)
  riv_ids_u  <- names(riv_groups)
  reach_agg <- data.frame(
    RivID          = as.integer(riv_ids_u),
    val            = vapply(riv_groups, function(idx) mean(det_WithRiv$val[idx]), numeric(1)),
    y_resp         = vapply(riv_groups, function(idx) mean(det_WithRiv$y_sample[idx]), numeric(1)),
    obs_up_dist_km = vapply(riv_groups, function(idx) mean(det_WithRiv$obs_up_dist_km[idx]), numeric(1)),
    stringsAsFactors = FALSE
  )

  has_samp_unc <- "meas_uncertainty" %in% names(det_WithRiv) && any(!is.na(det_WithRiv$meas_uncertainty))
  if (has_samp_unc) {
    reach_agg$meas_uncertainty <- vapply(riv_groups, function(idx) {
      u <- det_WithRiv$meas_uncertainty[idx]
      u_clean <- u[!is.na(u) & u >= 0]
      if (length(u_clean) == 0L) return(NA_real_)
      sqrt(sum(u_clean^2)) / length(u_clean)
    }, numeric(1))
  }
  has_samp_meth <- "meas_method" %in% names(det_WithRiv) && isTRUE(any(nzchar(det_WithRiv$meas_method[!is.na(det_WithRiv$meas_method)])))
  if (has_samp_meth) {
    reach_agg$meas_method <- vapply(riv_groups, function(idx) {
      m <- det_WithRiv$meas_method[idx]
      m_clean <- m[!is.na(m) & nzchar(m)]
      if (length(m_clean) == 0L) return("")
      m_clean[1]
    }, character(1))
  }

  df_full <- ctx$df_reaches
  idx_match   <- match(df_full$RivID, reach_agg$RivID)
  has_obs_idx <- which(!is.na(idx_match))
  matched_row <- idx_match[has_obs_idx]

  df_full$has_obs <- 0L
  df_full$has_obs[has_obs_idx] <- 1L
  df_full$obs_val <- NA_real_
  df_full$obs_val[has_obs_idx] <- reach_agg$val[matched_row]
  df_full$obs_y <- NA_real_
  df_full$obs_y[has_obs_idx] <- reach_agg$y_resp[matched_row]
  df_full$obs_up_dist_km <- df_full$up_dist_mid_km
  df_full$obs_up_dist_km[has_obs_idx] <- reach_agg$obs_up_dist_km[matched_row]
  if ("meas_uncertainty" %in% names(reach_agg)) {
    df_full$meas_uncertainty <- NA_real_
    df_full$meas_uncertainty[has_obs_idx] <- reach_agg$meas_uncertainty[matched_row]
  } else if (!"meas_uncertainty" %in% names(df_full)) {
    df_full$meas_uncertainty <- NA_real_
  }
  df_full$meas_method <- rep("", nrow(df_full))
  if ("meas_method" %in% names(reach_agg)) {
    mm <- as.character(reach_agg$meas_method[matched_row])
    mm[is.na(mm)] <- ""
    df_full$meas_method[has_obs_idx] <- mm
  }

  active_meas_error <- isTRUE(ctx$meas_error_arg) ||
    (has_samp_meth && isTRUE(any(tolower(det_WithRiv$meas_method) == "gam", na.rm = TRUE))) ||
    ("meas_method" %in% names(df_full) && isTRUE(any(tolower(df_full$meas_method[df_full$has_obs == 1L]) == "gam", na.rm = TRUE))) ||
    (!is.null(ctx$net_meta$meas_method) && identical(tolower(ctx$net_meta$meas_method), "gam"))
  active_meas_error <- isTRUE(active_meas_error)

  fit_domain <- function(df_dom, dom_id, dom_label) {
    fit_s3n_domain(
      df = df_dom,
      domain_id = dom_id,
      domain_label = dom_label,
      trans_obj = ctx$trans_obj,
      var_specs = ctx$var_specs,
      m_neighbors = 15L,
      n_boot = ctx$n_boot,
      n_omp = ctx$n_omp,
      use_meas_error = active_meas_error,
      net_meta = ctx$net_meta,
      obs_col = "obs_val"
    )
  }

  humber_info <- resolve_core_basin(df_full, ctx$catch_lookup, ctx$net_meta, "humber")
  thames_info <- resolve_core_basin(df_full, ctx$catch_lookup, ctx$net_meta, "thames")
  severn_info <- resolve_core_basin(df_full, ctx$catch_lookup, ctx$net_meta, "severn")

  eng_domain_label <- switch(ctx$source_cli,
    rsn = "All England RSN Network",
    ssn = "All England SSN Network",
    rss = "All England RSS Combined Network",
    sprintf("All England %s Network", ctx$source_abbr)
  )
  res_eng <- fit_domain(df_full, "england", eng_domain_label)
  if (is.null(res_eng)) return(NULL)
  res_hum <- fit_domain(subset(df_full, CatchID == humber_info$catch_id), "humber", humber_info$domain_label)
  res_tha <- fit_domain(subset(df_full, CatchID == thames_info$catch_id), "thames", thames_info$domain_label)
  res_sev <- fit_domain(subset(df_full, CatchID == severn_info$catch_id), "severn", severn_info$domain_label)

  eng_pred <- res_eng$pred_df
  hum_pred <- if (!is.null(res_hum)) res_hum$pred_df else subset(eng_pred, CatchID == humber_info$catch_id)
  tha_pred <- if (!is.null(res_tha)) res_tha$pred_df else subset(eng_pred, CatchID == thames_info$catch_id)
  sev_pred <- if (!is.null(res_sev)) res_sev$pred_df else subset(eng_pred, CatchID == severn_info$catch_id)

  # Fit any user-specified additional basins (--basin=...)
  extra_fits <- lapply(ctx$extra_basins, function(b_info) {
    sub_df <- subset(df_full, CatchID %in% b_info$catch_ids)
    fit_b  <- fit_domain(sub_df, b_info$slug, b_info$domain_label)
    pred_b <- if (!is.null(fit_b)) fit_b$pred_df else subset(eng_pred, CatchID %in% b_info$catch_ids)
    list(info = b_info, fit = fit_b, pred_df = pred_b)
  })

  dom_list <- c(list(res_eng, res_hum, res_tha, res_sev), lapply(extra_fits, function(x) x$fit))
  all_cv     <- do.call(rbind, lapply(dom_list, function(x) if (!is.null(x)) x$cv_summary else NULL))
  all_params <- do.call(rbind, lapply(dom_list, function(x) if (!is.null(x)) x$param_df else NULL))
  det_title_name <- if (period_tag == "all") ctx$det_name else sprintf("%s [%s]", shorten_det_name(ctx$det_name), period_tag)
  all_cv$Source           <- ctx$source_abbr
  all_cv$Source_Desc      <- ctx$source_label
  all_cv$Period           <- period_tag
  all_cv$Determinand_Code <- ctx$det_code
  all_cv$Determinand_Name <- det_title_name
  all_cv$Unit             <- ctx$det_unit
  all_cv$Scale            <- if (ctx$trans_obj$is_transformed) ctx$trans_obj$short_label else "raw"
  all_cv$Transform_Type   <- ctx$trans_obj$type
  all_cv$Transform_Lambda <- ifelse(is.na(ctx$trans_obj$lambda), "", as.character(round(ctx$trans_obj$lambda, 4)))
  all_cv$Transform_Label  <- ctx$trans_obj$label
  all_cv$Network_Abbr     <- ctx$net_meta$abbr
  all_cv$Network_Label    <- ctx$net_meta$short_label
  all_cv$Spatial_Method   <- if (active_meas_error) "Measurement Error Kriging" else "Standard S3N"
  all_cv$Meas_Method      <- if (active_meas_error) "gam" else "standard"
  all_params$Source           <- ctx$source_abbr
  all_params$Source_Desc      <- ctx$source_label
  all_params$Period           <- period_tag
  all_params$Transform_Type   <- ctx$trans_obj$type
  all_params$Transform_Lambda <- ifelse(is.na(ctx$trans_obj$lambda), "", as.character(round(ctx$trans_obj$lambda, 4)))
  all_params$Transform_Label  <- ctx$trans_obj$label
  all_params$Network_Abbr     <- ctx$net_meta$abbr
  all_params$Network_Label    <- ctx$net_meta$short_label
  all_params$Spatial_Method   <- if (active_meas_error) "Measurement Error Kriging" else "Standard S3N"
  all_params$Meas_Method      <- if (active_meas_error) "gam" else "standard"

  # Compute normalised results (Z-score, Robust Z-score, Min-Max, Percentile) using England-wide reference stats
  eng_pred <- add_normalised_columns(eng_pred, ref_df = eng_pred)
  hum_pred <- add_normalised_columns(hum_pred, ref_df = eng_pred)
  tha_pred <- add_normalised_columns(tha_pred, ref_df = eng_pred)
  sev_pred <- add_normalised_columns(sev_pred, ref_df = eng_pred)
  if (length(extra_fits) > 0) {
    for (k in seq_along(extra_fits)) {
      extra_fits[[k]]$pred_df <- add_normalised_columns(extra_fits[[k]]$pred_df, ref_df = eng_pred)
    }
  }

  write.csv(eng_pred,   file.path(out_subdir, sprintf("predictions_%s.csv", ctx$det_code)), row.names = FALSE)
  write.csv(hum_pred,   file.path(out_subdir, sprintf("predictions_humber_%s.csv", ctx$det_code)), row.names = FALSE)
  write.csv(tha_pred,   file.path(out_subdir, sprintf("predictions_thames_%s.csv", ctx$det_code)), row.names = FALSE)
  write.csv(sev_pred,   file.path(out_subdir, sprintf("predictions_severn_%s.csv", ctx$det_code)), row.names = FALSE)
  if (!is.null(res_eng$cv_pts)) {
    write.csv(res_eng$cv_pts, file.path(out_subdir, sprintf("cv_predictions_england_%s.csv", ctx$det_code)), row.names = FALSE)
  }
  if (!is.null(res_hum) && !is.null(res_hum$cv_pts)) {
    write.csv(res_hum$cv_pts, file.path(out_subdir, sprintf("cv_predictions_humber_%s.csv", ctx$det_code)), row.names = FALSE)
  }
  write.csv(all_params, file.path(out_subdir, sprintf("parameters_%s.csv", ctx$det_code)),  row.names = FALSE)
  write.csv(all_cv,     file.path(out_subdir, sprintf("cv_comparison_%s.csv", ctx$det_code)), row.names = FALSE)
  if (length(extra_fits) > 0) {
    for (k in seq_along(extra_fits)) {
      ef <- extra_fits[[k]]
      pred_csv_k <- file.path(out_subdir, sprintf("predictions_basin_%s_%s.csv", ef$info$slug, ctx$det_code))
      write.csv(ef$pred_df, pred_csv_k, row.names = FALSE)
    }
  }

  list(
    out_subdir     = out_subdir,
    period_tag     = period_tag,
    det_title_name = det_title_name,
    eng_pred       = eng_pred,
    hum_pred       = hum_pred,
    tha_pred       = tha_pred,
    sev_pred       = sev_pred,
    res_eng        = res_eng,
    res_hum        = res_hum,
    extra_fits     = extra_fits,
    cv             = all_cv,
    params         = all_params
  )
}

plot_temporal_slice <- function(sl, sc_val, sc_se_eng = NULL, sc_se_hum = NULL, lims_a = NULL, ylim_d = NULL) {
  if (is.null(sl)) return(invisible(NULL))
  out_subdir     <- sl$out_subdir
  det_title_name <- sl$det_title_name
  eng_pred       <- sl$eng_pred
  hum_pred       <- sl$hum_pred
  tha_pred       <- sl$tha_pred
  sev_pred       <- sl$sev_pred
  res_eng        <- sl$res_eng
  res_hum        <- sl$res_hum
  extra_fits     <- sl$extra_fits
  all_cv         <- sl$cv
  all_params     <- sl$params

  n_figs <- 4L + length(extra_fits)
  scale_disp_min <- if (norm_meta$mode == "none") trans_obj$inv(sc_val$vmin) else sc_val$vmin
  scale_disp_max <- if (norm_meta$mode == "none") trans_obj$inv(sc_val$vmax) else sc_val$vmax
  cat(sprintf("\nGenerating %d publication-quality figures in %s (shared scale [%.3f, %.3f]%s) ...\n",
              n_figs, out_subdir, scale_disp_min, scale_disp_max,
              if (norm_meta$mode == "none") "" else paste0(", norm = ", norm_meta$mode)))

  ftag      <- norm_meta$file_tag
  fig1_path <- file.path(out_subdir, sprintf("fig1_england_stream_%s_map%s.png", det_code, ftag))
  fig2_path <- file.path(out_subdir, sprintf("fig2_humber_basin_%s_map%s.png", det_code, ftag))
  fig3_path <- file.path(out_subdir, sprintf("fig3_thames_severn_%s_maps%s.png", det_code, ftag))
  fig4_path <- file.path(out_subdir, sprintf("fig4_cv_and_diagnostics_%s.png", det_code))
  vars_sub  <- format_covariates_subtitle(var_specs = var_specs, param_tab = all_params)
  net_lbl   <- net_meta$short_label

  plot_fig1_england_map(eng_pred, det_code, det_title_name, det_unit, trans_obj, fig1_path,
                        reverse_color = reverse_color, sc_val = sc_val, sc_se = sc_se_eng,
                        vars_label = vars_sub, network_label = net_lbl, source_label = source_label,
                        norm_mode = norm_cli)
  plot_fig2_humber_map(hum_pred, all_cv, sc_val, det_code, det_title_name, det_unit, trans_obj, fig2_path,
                       reverse_color = reverse_color, sc_se = sc_se_hum,
                       vars_label = vars_sub, network_label = net_lbl, source_label = source_label,
                       norm_mode = norm_cli)
  plot_fig3_thames_severn_map(tha_pred, sev_pred, sc_val, det_code, det_title_name, det_unit, trans_obj, fig3_path,
                              reverse_color = reverse_color, vars_label = vars_sub, network_label = net_lbl, source_label = source_label,
                              norm_mode = norm_cli)

  extra_cv_for_fig4 <- lapply(extra_fits, function(ef) {
    list(
      short_title  = ef$info$short_title,
      domain_label = ef$info$domain_label,
      cv_pts       = if (!is.null(ef$fit)) ef$fit$cv_pts else NULL
    )
  })
  plot_fig4_diagnostics(eng_pred, res_eng$cv_pts, if (!is.null(res_hum)) res_hum$cv_pts else NULL,
                        all_cv, all_params, det_code, det_title_name, det_unit, trans_obj, fig4_path,
                        extra_cv_list = extra_cv_for_fig4, lims_a = lims_a, ylim_d = ylim_d,
                        vars_label = vars_sub, network_label = net_lbl, source_label = source_label)

  if (length(extra_fits) > 0) {
    for (k in seq_along(extra_fits)) {
      ef <- extra_fits[[k]]
      fig_k_path <- file.path(out_subdir, sprintf("fig%d_%s_basin_%s_map%s.png", 4L + k, ef$info$slug, det_code, ftag))
      plot_basin_map(ef$pred_df, ef$info$short_title, ef$info$domain_label,
                     all_cv, sc_val, det_code, det_title_name, det_unit, trans_obj, fig_k_path,
                     reverse_color = reverse_color, sc_se = sc_se_hum,
                     vars_label = vars_sub, network_label = net_lbl, source_label = source_label,
                     norm_mode = norm_cli)
      cat(sprintf("  + Saved custom basin map for '%s': %s\n", ef$info$domain_label, fig_k_path))
    }
  }
  cat(sprintf("Completed '%s'! Saved tables and all %d figures to %s/\n", sl$period_tag, n_figs, out_subdir))
}

compute_batch_scales <- function(slice_list) {
  valid_s <- Filter(Negate(is.null), slice_list)
  all_pred_vals <- unlist(lapply(valid_s, function(s) s$eng_pred[[norm_meta$pred_col]]))
  all_se_eng    <- unlist(lapply(valid_s, function(s) s$eng_pred$se_model_sc))
  all_se_hum    <- unlist(lapply(valid_s, function(s) s$hum_pred$se_model_sc))
  all_cv_vals   <- unlist(lapply(valid_s, function(s) c(s$res_eng$cv_pts$obs_y, s$res_eng$cv_pts$cv_comb)))
  all_raw_preds <- unlist(lapply(valid_s, function(s) quantile(s$eng_pred$pred_combined, probs = c(0.01, 0.96), na.rm = TRUE)))

  sc_val <- if (norm_meta$mode == "none") {
    compute_color_scale(all_pred_vals, trans_obj, vmin_raw = vmin_raw, vmax_raw = vmax_raw)
  } else {
    compute_norm_color_scale(all_pred_vals, norm_meta$mode, vmin_raw = vmin_raw, vmax_raw = vmax_raw)
  }
  sc_se_eng <- compute_se_scale(all_se_eng)
  sc_se_hum <- compute_se_scale(all_se_hum)

  lims_a <- range(all_cv_vals, na.rm = TRUE)
  pad_a  <- 0.08 * diff(lims_a)
  lims_a <- c(lims_a[1] - pad_a, lims_a[2] + pad_a)

  yr_d   <- range(all_raw_preds, na.rm = TRUE)
  ylim_d <- c(max(0, yr_d[1] - 0.05 * diff(yr_d)), yr_d[2] + 0.15 * diff(yr_d))

  list(sc_val = sc_val, sc_se_eng = sc_se_eng, sc_se_hum = sc_se_hum, lims_a = lims_a, ylim_d = ylim_d)
}

# Apply any explicit temporal filters (--year, --season, --month) before optional --by= slicing
base_samples <- det_samples
filter_tags  <- character(0)
filter_dirs  <- character(0)
if (!is.na(year_f)) {
  base_samples <- base_samples[!is.na(base_samples$year) & base_samples$year == year_f, ]
  filter_tags  <- c(filter_tags, paste0("Year ", year_f))
  filter_dirs  <- c(filter_dirs, paste0("year_", year_f))
}
if (!is.na(season_f)) {
  base_samples <- base_samples[!is.na(base_samples$season) & base_samples$season == season_f, ]
  filter_tags  <- c(filter_tags, paste0("Season ", season_f))
  filter_dirs  <- c(filter_dirs, paste0("season_", season_f))
}
if (!is.na(month_f)) {
  base_samples <- base_samples[!is.na(base_samples$month) & base_samples$month == month_f, ]
  filter_tags  <- c(filter_tags, sprintf("Month %02d", month_f))
  filter_dirs  <- c(filter_dirs, sprintf("month_%02d", month_f))
}

if (!is.na(by_mode)) {
  slice_root <- if (length(filter_dirs) == 0L) OUT_DIR else file.path(OUT_DIR, paste(filter_dirs, collapse = "_"))
  if (by_mode == "year") {
    periods <- sort(unique(base_samples$year[!is.na(base_samples$year)]))
    batch_res <- lapply(periods, function(yr) {
      sub_s <- base_samples[!is.na(base_samples$year) & base_samples$year == yr, ]
      fit_temporal_slice(sub_s, file.path(slice_root, paste0("year_", yr)), paste0("Year ", yr), slice_ctx)
    })
  } else if (by_mode == "season") {
    periods <- c("winter", "spring", "summer", "autumn")
    batch_res <- lapply(periods, function(seas) {
      sub_s <- base_samples[!is.na(base_samples$season) & base_samples$season == seas, ]
      fit_temporal_slice(sub_s, file.path(slice_root, paste0("season_", seas)), paste0("Season: ", seas), slice_ctx)
    })
  } else if (by_mode == "month") {
    periods <- sort(unique(base_samples$month[!is.na(base_samples$month)]))
    batch_res <- lapply(periods, function(mo) {
      sub_s <- base_samples[!is.na(base_samples$month) & base_samples$month == mo, ]
      fit_temporal_slice(sub_s, file.path(slice_root, sprintf("month_%02d", mo)), sprintf("Month %02d", mo), slice_ctx)
    })
  } else {
    stop(sprintf("Invalid --by='%s'. Choose from: year, season, month.", by_mode))
  }
  shared_sc <- compute_batch_scales(batch_res)
  for (sl in batch_res) {
    plot_temporal_slice(sl, shared_sc$sc_val, shared_sc$sc_se_eng, shared_sc$sc_se_hum,
                        shared_sc$lims_a, shared_sc$ylim_d)
  }
  grid_items <- lapply(Filter(Negate(is.null), batch_res), function(sl) {
    list(label = sl$period_tag, eng_pred = sl$eng_pred, cv_tab = sl$cv)
  })
  fig5_grid_path <- file.path(slice_root, sprintf("fig5_england_by_%s_%s_map%s.png", by_mode, det_code, norm_meta$file_tag))
  cat(sprintf("\nGenerating Figure 5 (%d-panel All England %s comparison grid): %s ...\n",
              length(grid_items), by_mode, fig5_grid_path))
  comb_cv <- do.call(rbind, lapply(batch_res, function(x) if (!is.null(x)) x$cv else NULL))
  comb_pa <- do.call(rbind, lapply(batch_res, function(x) if (!is.null(x)) x$params else NULL))
  vars_sub_grid <- format_covariates_subtitle(var_specs = var_specs, param_tab = comb_pa)
  plot_fig5_temporal_england_grid(grid_items, by_mode, shared_sc$sc_val,
                                  det_code, det_name, det_unit, trans_obj,
                                  fig5_grid_path, reverse_color = reverse_color,
                                  vars_label = vars_sub_grid, network_label = net_meta$short_label,
                                  source_label = source_label, norm_mode = norm_cli)
  write.csv(comb_cv, file.path(slice_root, sprintf("cv_summary_by_%s_%s.csv", by_mode, det_code)), row.names = FALSE)
  write.csv(comb_pa, file.path(slice_root, sprintf("parameters_by_%s_%s.csv", by_mode, det_code)), row.names = FALSE)
} else {
  sub_s  <- base_samples
  p_tag  <- if (length(filter_tags) == 0) "all" else paste(filter_tags, collapse = ", ")
  s_dir  <- if (length(filter_dirs) == 0) OUT_DIR else file.path(OUT_DIR, paste(filter_dirs, collapse = "_"))
  sl_one <- fit_temporal_slice(sub_s, s_dir, p_tag, slice_ctx)
  if (!is.null(sl_one)) {
    # When running a single temporal slice (--year, --season, --month), pool with any existing
    # sibling period predictions in OUT_DIR so scale bars stay consistent across periods
    self_eng_p <- normalizePath(file.path(s_dir, sprintf("predictions_%s.csv", det_code)), mustWork = FALSE)
    sib_files  <- list.files(OUT_DIR, pattern = sprintf("^predictions_%s\\.csv$", det_code),
                             recursive = TRUE, full.names = TRUE)
    sib_files  <- sib_files[normalizePath(sib_files, mustWork = FALSE) != self_eng_p]
    sib_vals <- sl_one$eng_pred[[norm_meta$pred_col]]
    sib_se_e <- sl_one$eng_pred$se_model_sc
    sib_se_h <- sl_one$hum_pred$se_model_sc
    if (length(sib_files) > 0) {
      for (sf in sib_files) {
        sdf <- tryCatch(read.csv(sf, stringsAsFactors = FALSE), error = function(e) NULL)
        if (!is.null(sdf)) {
          if (norm_meta$mode != "none" && !(norm_meta$pred_col %in% names(sdf))) {
            sdf <- add_normalised_columns(sdf)
          }
          sib_vals <- c(sib_vals, sdf[[norm_meta$pred_col]])
          sib_se_e <- c(sib_se_e, sdf$se_model_sc)
        }
        hf <- file.path(dirname(sf), sprintf("predictions_humber_%s.csv", det_code))
        if (file.exists(hf)) {
          hdf <- tryCatch(read.csv(hf, stringsAsFactors = FALSE), error = function(e) NULL)
          if (!is.null(hdf) && "se_model_sc" %in% names(hdf)) {
            sib_se_h <- c(sib_se_h, hdf$se_model_sc)
          }
        }
      }
    }
    sc_val_one <- if (norm_meta$mode == "none") {
      compute_color_scale(sib_vals, trans_obj, vmin_raw = vmin_raw, vmax_raw = vmax_raw)
    } else {
      compute_norm_color_scale(sib_vals, norm_meta$mode, vmin_raw = vmin_raw, vmax_raw = vmax_raw)
    }
    sc_se_e1   <- compute_se_scale(sib_se_e)
    sc_se_h1   <- compute_se_scale(sib_se_h)
    plot_temporal_slice(sl_one, sc_val_one, sc_se_e1, sc_se_h1)
  }
}
