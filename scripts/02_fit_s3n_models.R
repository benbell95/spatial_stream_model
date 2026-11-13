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
source(file.path(BASE_DIR, "scripts/s3n_model_functions.R"))
BASE_DIR <- resolve_base_dir()

args    <- commandArgs(trailingOnly = TRUE)
dir_cfg <- parse_custom_dir_args(args, base_dir = BASE_DIR)
GIS_DIR  <- dir_cfg$gis_dir
DATA_DIR <- dir_cfg$data_dir
OUT_BASE <- dir_cfg$out_dir

meas_error_arg <- any(c("--meas-error", "--gam") %in% args)
n_omp_cli      <- parse_n_omp_arg(args, default = 1L)

src_cli      <- parse_source_arg(args)
source_cli   <- resolve_source_tag(src_cli, base_dir = BASE_DIR)
source_abbr  <- format_source_abbr(src_cli, base_dir = BASE_DIR)
source_label <- format_source_label(src_cli, base_dir = BASE_DIR)
source_tag   <- format_source_folder_tag(source_cli)
range_arg    <- grep("^--range-mode=", args, value = TRUE)
range_cli    <- if (length(range_arg) > 0) tolower(sub("^--range-mode=", "", range_arg[1])) else NA_character_
if (!is.na(range_cli) && !(range_cli %in% c("auto", "none"))) {
  src_tokens <- unlist(strsplit(source_cli, "[+,]+"))
  has_tok <- any(grepl(paste0("_", range_cli, "$"), src_tokens))
  if (!has_tok) {
    if (nzchar(source_tag) && !grepl(paste0("_", range_cli, "$"), source_tag)) {
      source_tag <- paste0(source_tag, "_", range_cli)
    } else if (!nzchar(source_tag)) {
      source_tag <- range_cli
    }
    if (!grepl("[+,]", source_cli) && !grepl(paste0("_", range_cli, "$"), source_cli)) {
      source_cli <- paste0(source_cli, "_", range_cli)
    }
  }
  source_abbr  <- format_source_abbr(source_cli, base_dir = BASE_DIR)
  source_label <- format_source_label(source_cli, base_dir = BASE_DIR)
}

if ("--list-sources" %in% args) {
  print_available_sources()
  quit(save = "no", status = 0)
}

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
net_arg_vec <- grep("^--network=", args, value = TRUE)
network_cli <- if (length(net_arg_vec) > 0) sub("^--network=", "", net_arg_vec[1]) else "orn"
PROC_DIR    <- resolve_network_processed_dir(BASE_DIR, network_cli, data_dir = DATA_DIR, gis_dir = GIS_DIR, source = source_cli)
net_meta    <- load_network_metadata(PROC_DIR)

vars_args <- grep("^--(vars|covariates)=", args, value = TRUE)
var_specs <- if (length(vars_args) > 0) {
  unlist(strsplit(sub("^--(vars|covariates)=", "", vars_args), ",", fixed = TRUE))
} else {
  character(0)
}
var_specs <- trimws(var_specs[nchar(trimws(var_specs)) > 0])

trans_arg  <- grep("^--(transform|trans|transformation)=", args, value = TRUE)
trans_cli  <- if (length(trans_arg) > 0) tolower(sub("^--(transform|trans|transformation)=", "", trans_arg[1])) else NA_character_
lambda_arg <- grep("^--(lambda|transform-lambda)=", args, value = TRUE)
lambda_cli <- if (length(lambda_arg) > 0) as.numeric(sub("^--(lambda|transform-lambda)=", "", lambda_arg[1])) else NULL
scale_arg  <- grep("^--(scale|transform-scale)=", args, value = TRUE)
scale_cli  <- if (length(scale_arg) > 0) as.numeric(sub("^--(scale|transform-scale)=", "", scale_arg[1])) else NULL

if (is.na(trans_cli)) {
  if ("--raw" %in% args) {
    trans_cli <- "none"
  } else if ("--log" %in% args) {
    trans_cli <- "log"
  } else {
    trans_cli <- "log"
  }
}

trans_obj <- build_transformer(trans_cli, lambda = lambda_cli, scale = scale_cli)
use_log   <- trans_obj$is_log

boot_arg   <- grep("^--(boot|n-boot)=", args, value = TRUE)
n_boot_cli <- if (length(boot_arg) > 0) as.integer(sub("^--(boot|n-boot)=", "", boot_arg[1])) else 100L

load_s3n_backend(BASE_DIR)

cat("======================================================================\n")
cat(sprintf("  Scalable Spatial Stream Network (S3N) Nitrate Modeling Pipeline [%s | Source: %s]\n",
            net_meta$label, source_label))
cat("======================================================================\n\n")

t_start <- proc.time()
df_full <- read.csv(file.path(PROC_DIR, "england_rsn_reaches.csv"))
catch_lu <- if (file.exists(file.path(PROC_DIR, "rsn_catchments_lookup.csv"))) {
  read.csv(file.path(PROC_DIR, "rsn_catchments_lookup.csv"), stringsAsFactors = FALSE)
} else {
  NULL
}

# Optimize profile Maximum Likelihood lambda for Box-Cox or Yeo-Johnson BEFORE computing folder tags
if (trans_obj$type %in% c("boxcox", "yeojohnson") && is.na(trans_obj$lambda)) {
  obs_raw <- df_full$nitrate_mean[df_full$has_obs == 1 & !is.na(df_full$nitrate_mean)]
  trans_obj <- build_transformer(trans_cli, y_obs = obs_raw, lambda = lambda_cli, scale = scale_cli)
  cat(sprintf("  Estimated Profile Maximum Likelihood lambda = %.4f for Nitrate (%s)\n",
              trans_obj$lambda, trans_obj$label))
}

trans_tag <- format_transform_folder_tag(trans_obj, default_type = "log")
vars_tag  <- format_covariates_folder_tag(var_specs)
sub_tag <- paste0(net_meta$abbr, "_", source_tag)
if (nzchar(trans_tag)) {
  sub_tag <- paste0(sub_tag, "_", trans_tag)
}
if (nzchar(vars_tag)) {
  sub_tag <- paste0(sub_tag, "_", vars_tag)
}
OUT_ROOT <- rel_path(OUT_BASE, "results", sub_tag)
dir.create(file.path(OUT_ROOT, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "predictions"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "models"), recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Loaded %d stream reaches (%d with RSN Nitrate observations) across %d catchments.\n",
            nrow(df_full), sum(df_full$has_obs == 1), length(unique(df_full$CatchID))))

# Fitting/CV engine lives in scripts/s3n_model_functions.R
analyze_domain <- function(df, domain_id, domain_label, m_neighbors = 15L, n_boot = n_boot_cli, n_omp = n_omp_cli) {
  fit_s3n_domain(
    df = df,
    domain_id = domain_id,
    domain_label = domain_label,
    trans_obj = trans_obj,
    var_specs = var_specs,
    m_neighbors = m_neighbors,
    n_boot = n_boot,
    n_omp = n_omp,
    use_meas_error = isTRUE(meas_error_arg),
    net_meta = net_meta,
    obs_col = "nitrate_mean",
    legacy_nitrate_mode = TRUE,
    out_root = OUT_ROOT
  )
}

humber_info <- resolve_core_basin(df_full, catch_lu, net_meta, "humber")
thames_info <- resolve_core_basin(df_full, catch_lu, net_meta, "thames")
severn_info <- resolve_core_basin(df_full, catch_lu, net_meta, "severn")

eng_domain_label <- switch(source_cli,
  rsn = "All England RSN Network",
  ssn = "All England SSN Network",
  rss = "All England RSS Combined Network",
  sprintf("All England %s Network", source_abbr)
)
res_england <- analyze_domain(df_full, "england", eng_domain_label, m_neighbors = 15L, n_boot = n_boot_cli)
res_humber  <- analyze_domain(subset(df_full, CatchID == humber_info$catch_id), "humber", humber_info$domain_label, m_neighbors = 15L, n_boot = n_boot_cli)
res_thames  <- analyze_domain(subset(df_full, CatchID == thames_info$catch_id), "thames", thames_info$domain_label, m_neighbors = 15L, n_boot = n_boot_cli)
res_severn  <- analyze_domain(subset(df_full, CatchID == severn_info$catch_id), "severn", severn_info$domain_label, m_neighbors = 15L, n_boot = n_boot_cli)

all_cv <- do.call(rbind, Filter(Negate(is.null), list(
  if (!is.null(res_england)) res_england$cv_summary else NULL,
  if (!is.null(res_humber))  res_humber$cv_summary  else NULL,
  if (!is.null(res_thames))  res_thames$cv_summary  else NULL,
  if (!is.null(res_severn))  res_severn$cv_summary  else NULL
)))
all_params <- do.call(rbind, Filter(Negate(is.null), list(
  if (!is.null(res_england)) res_england$param_df else NULL,
  if (!is.null(res_humber))  res_humber$param_df  else NULL,
  if (!is.null(res_thames))  res_thames$param_df  else NULL,
  if (!is.null(res_severn))  res_severn$param_df  else NULL
)))
all_cv$Source               <- source_abbr
all_cv$Source_Desc          <- source_label
all_cv$Transform_Type       <- trans_obj$type
all_cv$Transform_Lambda     <- ifelse(is.na(trans_obj$lambda), "", as.character(round(trans_obj$lambda, 4)))
all_cv$Transform_Label      <- trans_obj$label
all_cv$Scale                <- if (trans_obj$is_transformed) trans_obj$short_label else "raw"
all_cv$Network_Abbr         <- net_meta$abbr
all_cv$Network_Label        <- net_meta$short_label

all_params$Source           <- source_abbr
all_params$Source_Desc      <- source_label
all_params$Transform_Type   <- trans_obj$type
all_params$Transform_Lambda <- ifelse(is.na(trans_obj$lambda), "", as.character(round(trans_obj$lambda, 4)))
all_params$Transform_Label  <- trans_obj$label
all_params$Network_Abbr     <- net_meta$abbr
all_params$Network_Label    <- net_meta$short_label

write.csv(all_cv, file.path(OUT_ROOT, "tables/cv_model_comparison.csv"), row.names = FALSE)
write.csv(all_params, file.path(OUT_ROOT, "tables/s3n_parameter_estimates.csv"), row.names = FALSE)
if (!is.null(res_england) && !is.null(res_england$cv_pts)) {
  write.csv(res_england$cv_pts, file.path(OUT_ROOT, "predictions/england_cv_predictions.csv"), row.names = FALSE)
}
if (!is.null(res_humber) && !is.null(res_humber$cv_pts)) {
  write.csv(res_humber$cv_pts, file.path(OUT_ROOT, "predictions/humber_cv_predictions.csv"), row.names = FALSE)
}

saveRDS(list(england = res_england, humber = res_humber, thames = res_thames, severn = res_severn),
        file = file.path(OUT_ROOT, "models/s3n_model_results.rds"))

cat(sprintf("\n======================================================================\n"))
cat(sprintf("  All S3N models, bootstrap CIs, CV, and predictions finished in %.2f sec.\n",
            (proc.time() - t_start)[3]))
cat(sprintf("======================================================================\n"))

