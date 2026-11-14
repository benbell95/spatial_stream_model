#!/usr/bin/env Rscript
# ==============================================================================
# run_s3n_cokriging.R
#
# Dedicated Multivariate Spatial Stream Network (S3N) Co-Kriging Runner.
# Jointly models P >= 2 determinands (including heterotopic cross-dataset targets
# such as water quality + macroinvertebrates) via Scalable Factor-LMC NNGP
# (with optional fast --method=colocated residual co-kriging) without altering
# the single-determinand modeling workflow.
#
# Usage:
#   Rscript scripts/run_s3n_cokriging.R --help
#   Rscript scripts/run_s3n_cokriging.R --targets=0117,0180
#   Rscript scripts/run_s3n_cokriging.R --targets=inverts:whpt_aspt,wq:0117,wq:0180,wq:9924 --factors=3
# ==============================================================================

args     <- commandArgs(trailingOnly = TRUE)
BASE_DIR <- if (file.exists("scripts/s3n_plot_functions.R")) "." else if (file.exists("s3n_plot_functions.R")) ".." else "."
source(file.path(BASE_DIR, "scripts/s3n_plot_functions.R"))
source(file.path(BASE_DIR, "scripts/s3n_model_functions.R"))
source(file.path(BASE_DIR, "scripts/s3n_multivariate_functions.R"))
BASE_DIR <- resolve_base_dir()

print_cokrig_usage <- function() {
  cat("Multivariate Stream Co-Kriging (Scalable Factor-LMC Joint NNGP)\n\n")
  cat("Usage:\n")
  cat("  Rscript scripts/run_s3n_cokriging.R --targets=<t1>,<t2>[,<t3>,...] [OPTIONS]\n\n")
  cat("Targets:\n")
  cat("  --targets=...          Comma-separated P >= 2 targets (or first positional arg).\n")
  cat("                         Prefix with dataset if needed: wq:0117, inverts:whpt_aspt,\n")
  cat("                         diatoms:tdi4, macrophytes:rmni, rhs:hms, lcms:<cas>\n")
  cat("  --dataset=<ds>         Default dataset when unprefixed (auto-detected if omitted)\n\n")
  cat("Co-Kriging Model Options:\n")
  cat("  --method=lmc|colocated Joint Factor-LMC NNGP (default: lmc) or fast residual co-kriging\n")
  cat("  --factors=<K>          Latent spatial factors K <= P (default: 3; full-rank for P <= 3,\n")
  cat("                         scalable latent factors for P >= 4)\n")
  cat("  --neighbors=<m>        Nearest stream neighbors per Tail-Up/Down DAG (default: 10)\n")
  cat("  --transforms=<t1,...>  Per-target transform (or single broadcast); default: log (wq/lcms),\n")
  cat("                         none (ecology indices)\n")
  cat("  --kernel=<k>           Stream correlation kernel (default: exponential; see --list-kernels)\n")
  cat("  --vars=<v1,v2,...>     Catchment covariates (see --list-vars)\n")
  cat("  --cv-mode=<mode>       Out-of-sample CV holdout mode (default: random; see --list-cv)\n")
  cat("  --threshold=uk|<vals>  UK regulatory EQS auto-lookup or comma-separated numeric thresholds\n")
  cat("  --threshold-dir=<d>    above or below (comma-separated or broadcast; default: above / UK dir)\n")
  cat("  --censor=ros|half_lod  Left-censored ('<') handling (default: ros)\n\n")
  cat("Domain, Temporal & Runtime Options:\n")
  cat("  --source=rsn|ssn|rss   Monitoring network source (default: rsn; see --list-sources)\n")
  cat("  --prefer-source=<s..>  Priority order across co-located reaches (e.g. --prefer-source=rsn)\n")
  cat("  --source-weight=<s:w>  Per-source reliability weights in (0,1] (e.g. rsn:1,waterblitz:0.25)\n")
  cat("  --source-se=<s:sd,...> Per-source observation SE on raw units (e.g. rsn:0.01,waterblitz:0.5)\n")
  cat("  --network=<name>       Stream network (default: orn; see --list-networks)\n")
  cat("  --basin=<CatchID|Name> Additional custom basin(s) to fit & map (see --list-basins)\n")
  cat("  --year=YYYY            Filter all targets to survey year\n")
  cat("  --season=winter|spring|summer|autumn\n")
  cat("  --month=1..12          Filter all targets to month\n")
  cat("  --omp=<n>              OpenMP threads (default: 1; use 4-8 on multi-core workstations)\n")
  cat("  --reverse-color        Invert choropleth palettes on maps\n")
  cat("  --norm=<mode>          Map value normalization (default: none)\n")
  cat("  --fig-format=<fmt>     Output figure format: png (default), pdf, or svg\n")
  cat("  --dpi=<N>              Output PNG resolution in DPI (default: 200; e.g. 300)\n")
  cat("  --fig-size=<pct>       Scale figure dimensions (e.g. 125%, +25%, -20%, or 1.25)\n\n")
  cat("Examples:\n")
  cat("  Rscript scripts/run_s3n_cokriging.R --targets=0117,0180\n")
  cat("  Rscript scripts/run_s3n_cokriging.R --targets=inverts:whpt_aspt,wq:0117,wq:0180,wq:9924 \\\n")
  cat("    --factors=3 --cv-mode=spatial --threshold=uk --omp=4\n")
}

if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
  print_cokrig_usage()
  quit(save = "no", status = 0)
}

if ("--list-sources" %in% args) { print_available_sources(base_dir = BASE_DIR); quit(save = "no", status = 0) }
if ("--list-vars" %in% args || "--list-covariates" %in% args) { print_available_covariates(); quit(save = "no", status = 0) }
if ("--list-transforms" %in% args || "--list-transformations" %in% args) { print_available_transforms(); quit(save = "no", status = 0) }
if ("--list-kernels" %in% args) { print_available_kernels(); quit(save = "no", status = 0) }
if (any(c("--list-validation", "--list-cv") %in% args)) { print_available_validation_options(); quit(save = "no", status = 0) }
if ("--list-thresholds" %in% args) { print_available_thresholds(); quit(save = "no", status = 0) }

dir_cfg  <- parse_custom_dir_args(args, base_dir = BASE_DIR)
GIS_DIR  <- dir_cfg$gis_dir
DATA_DIR <- dir_cfg$data_dir
OUT_BASE <- dir_cfg$out_dir

custom_samples_arg <- grep("^--custom-samples=", args, value = TRUE)
custom_samples_cli <- if (length(custom_samples_arg) > 0) sub("^--custom-samples=", "", custom_samples_arg[1]) else NULL
custom_sites_arg   <- grep("^--custom-sites=", args, value = TRUE)
custom_sites_cli   <- if (length(custom_sites_arg) > 0) sub("^--custom-sites=", "", custom_sites_arg[1]) else NULL
custom_dir_arg     <- grep("^--custom-dir=", args, value = TRUE)
custom_dir_cli     <- if (length(custom_dir_arg) > 0) sub("^--custom-dir=", "", custom_dir_arg[1]) else NULL

src_cli      <- parse_source_arg(args)
source_cli   <- resolve_source_tag(src_cli, base_dir = BASE_DIR)
source_abbr  <- format_source_abbr(src_cli, base_dir = BASE_DIR)
source_label <- format_source_label(src_cli, base_dir = BASE_DIR)
source_tag   <- format_source_folder_tag(source_cli)
range_arg    <- grep("^--range-mode=", args, value = TRUE)
range_cli    <- if (length(range_arg) > 0) tolower(sub("^--range-mode=", "", range_arg[1])) else NA_character_
if (!is.na(range_cli) && !(range_cli %in% c("auto", "none"))) {
  rm_res       <- apply_range_mode_to_source(source_cli, range_cli, base_dir = BASE_DIR)
  source_cli   <- rm_res$source_tag
  source_abbr  <- format_source_abbr(source_cli, base_dir = BASE_DIR)
  source_label <- format_source_label(source_cli, base_dir = BASE_DIR)
  source_tag   <- format_source_folder_tag(source_cli)
}
src_dirs     <- resolve_data_source_dirs(DATA_DIR, source = source_cli, base_dir = BASE_DIR,
                                         custom_dir = custom_dir_cli,
                                         custom_samples = custom_samples_cli,
                                         custom_sites = custom_sites_cli)
RSN_DIR      <- src_dirs$rsn_dir
n_omp_cli    <- parse_n_omp_arg(args, default = 1L)

net_arg_vec <- grep("^--network=", args, value = TRUE)
network_cli <- if (length(net_arg_vec) > 0) sub("^--network=", "", net_arg_vec[1]) else "orn"

if ("--list-networks" %in% args) {
  run_python_script(rel_path(BASE_DIR, "scripts", "01_preprocess_orn.py"),
                    c("--list-networks", paste0("--gis-dir=", shQuote(GIS_DIR))))
  quit(save = "no", status = 0)
}

PROC_DIR       <- resolve_network_processed_dir(BASE_DIR, network_cli, data_dir = DATA_DIR, gis_dir = GIS_DIR,
                                                source = source_cli, processed_dir = dir_cfg$processed_dir,
                                                custom_dir = custom_dir_cli, custom_samples = custom_samples_cli,
                                                custom_sites = custom_sites_cli)
net_meta       <- load_network_metadata(PROC_DIR)
REACHES_CSV    <- rel_path(PROC_DIR, "england_rsn_reaches.csv")
SITES_CSV      <- rel_path(PROC_DIR, "rsn_sites.csv")
CATCHMENTS_CSV <- rel_path(PROC_DIR, "rsn_catchments_lookup.csv")

if ("--list-basins" %in% args) {
  c_lu <- fast_read_csv(CATCHMENTS_CSV)
  rbd_tab <- aggregate(cbind(n_reaches, n_obs_reaches) ~ RBD_ID + RBD_Name, data = c_lu, FUN = sum)
  rbd_tab <- rbd_tab[order(-rbd_tab$n_obs_reaches), ]
  cat("=== River Basin Districts ===\n")
  print(rbd_tab, row.names = FALSE)
  cat("\n=== Top 30 Monitored Catchments ===\n")
  print(head(c_lu[order(-c_lu$n_obs_reaches, -c_lu$n_reaches),
                  c("CatchID", "Main_River", "RBD_Name", "n_reaches", "n_obs_reaches")], 30), row.names = FALSE)
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

# Parse multivariate arguments
targets_arg <- grep("^--targets=", args, value = TRUE)
skip_pos_i  <- integer(0)
src_sp_i    <- which(args == "--source")
if (length(src_sp_i) > 0L) skip_pos_i <- c(skip_pos_i, src_sp_i + 1L)
pos_args    <- args[setdiff(which(!grepl("^--", args)), skip_pos_i)]
target_spec <- if (length(targets_arg) > 0) {
  sub("^--targets=", "", targets_arg[1])
} else if (length(pos_args) > 0) {
  pos_args[1]
} else {
  stop("Specify >= 2 targets via --targets=t1,t2[,...]. Run with --help for details.")
}

ds_arg     <- grep("^--dataset=", args, value = TRUE)
default_ds <- if (length(ds_arg) > 0) tolower(sub("^--dataset=", "", ds_arg[1])) else NA_character_

# Auto-detect dataset per token if --dataset is omitted and token lacks ds: prefix
raw_tokens <- trimws(unlist(strsplit(target_spec, ",", fixed = TRUE)))
raw_tokens <- raw_tokens[nzchar(raw_tokens)]
detect_target_dataset <- function(tok) {
  if (grepl(":", tok, fixed = TRUE)) return(tok)
  if (!is.na(default_ds)) return(paste0(default_ds, ":", tok))
  if (grepl("^[0-9]{4}$", tok)) return(paste0("wq:", tok))
  if (grepl("^[0-9]{2,7}-[0-9]{2}-[0-9]$", tok)) return(paste0("lcms:", tok))
  for (cand_ds in c("inverts", "diatoms", "macrophytes", "rhs")) {
    for (s_name in names(src_dirs$active_dirs)) {
      s_dir  <- src_dirs$active_dirs[[s_name]]
      cand_p <- file.path(s_dir, paste0(toupper(s_name), " ", DATASET_FILES[[cand_ds]]$file))
      if (!file.exists(cand_p)) cand_p <- file.path(s_dir, DATASET_FILES[[cand_ds]]$file)
      if (file.exists(cand_p)) {
        hdr <- colnames(fast_read_csv(cand_p, nrows = 1L, fileEncoding = "UTF-8-BOM"))
        if (tok %in% hdr) return(paste0(cand_ds, ":", tok))
      }
    }
  }
  paste0("wq:", tok)
}
resolved_spec <- paste(vapply(raw_tokens, detect_target_dataset, character(1)), collapse = ",")
targets_parsed <- parse_cokrig_targets(resolved_spec, default_dataset = "wq")
P <- length(targets_parsed)

method_arg <- grep("^--method=", args, value = TRUE)
method_cli <- if (length(method_arg) > 0) tolower(trimws(sub("^--method=", "", method_arg[1]))) else "lmc"
if (!(method_cli %in% c("lmc", "colocated"))) {
  stop(sprintf("Invalid --method='%s'. Choose 'lmc' or 'colocated'.", method_cli))
}

factors_arg <- grep("^--factors=", args, value = TRUE)
k_factors   <- if (length(factors_arg) > 0) max(1L, as.integer(sub("^--factors=", "", factors_arg[1]))) else 3L
k_factors   <- min(P, k_factors)

nb_arg    <- grep("^--neighbors=", args, value = TRUE)
m_neighbors <- if (length(nb_arg) > 0) max(2L, as.integer(sub("^--neighbors=", "", nb_arg[1]))) else 10L

kernel_cli  <- parse_kernel_arg(args, default = "exponential")
cv_mode_cli <- parse_cv_mode_arg(args, default = "random")
reverse_color <- any(args %in% c("--reverse-color", "--reverse-colour", "--rev-color", "--rev-colour", "--reverse-scale"))
norm_cli    <- parse_norm_mode(args, default = "none")
norm_meta   <- get_norm_meta(norm_cli)

censor_arg <- grep("^--censor=", args, value = TRUE)
censor_m   <- if (length(censor_arg) > 0) tolower(sub("^--censor=", "", censor_arg[1])) else "ros"
year_arg   <- grep("^--year=", args, value = TRUE)
year_f     <- if (length(year_arg) > 0) as.integer(sub("^--year=", "", year_arg[1])) else NA_integer_
season_arg <- grep("^--season=", args, value = TRUE)
season_f   <- if (length(season_arg) > 0) tolower(sub("^--season=", "", season_arg[1])) else NA_character_
if (!is.na(season_f) && season_f == "fall") season_f <- "autumn"
month_arg  <- grep("^--month=", args, value = TRUE)
month_f    <- if (length(month_arg) > 0) as.integer(sub("^--month=", "", month_arg[1])) else NA_integer_

prefer_sources_cli <- parse_prefer_source_arg(args)
src_weights_cli    <- parse_source_kv_arg(args, "--source-weight=", min_val = 1e-4, normalize_gt1 = TRUE)
src_se_cli         <- parse_source_kv_arg(args, "--source-se=", min_val = 0.0, normalize_gt1 = FALSE)

basin_args  <- grep("^--basin=", args, value = TRUE)
basin_specs <- if (length(basin_args) > 0) trimws(unlist(strsplit(sub("^--basin=", "", basin_args), ",", fixed = TRUE))) else character(0)
basin_specs <- basin_specs[nzchar(basin_specs)]

vars_args <- grep("^--(vars|covariates)=", args, value = TRUE)
var_specs <- if (length(vars_args) > 0) trimws(unlist(strsplit(sub("^--(vars|covariates)=", "", vars_args), ",", fixed = TRUE))) else character(0)
var_specs <- var_specs[nzchar(var_specs)]

trans_arg <- grep("^--(transforms|transform|trans)=", args, value = TRUE)
trans_tokens <- if (length(trans_arg) > 0) {
  tolower(trimws(unlist(strsplit(sub("^--(transforms|transform|trans)=", "", trans_arg[1]), ",", fixed = TRUE))))
} else {
  character(0)
}
if (length(trans_tokens) == 1L) trans_tokens <- rep(trans_tokens, P)
if (length(trans_tokens) > 0L && length(trans_tokens) != P) {
  stop(sprintf("--transforms length (%d) must be 1 or match number of targets (%d).", length(trans_tokens), P))
}

thresh_raw_arg <- grep("^--threshold=", args, value = TRUE)
thresh_raw_val <- if (length(thresh_raw_arg) > 0) trimws(sub("^--threshold=", "", thresh_raw_arg[1])) else ""
tdir_raw_arg   <- grep("^--threshold-dir=", args, value = TRUE)
tdir_raw_val   <- if (length(tdir_raw_arg) > 0) tolower(trimws(sub("^--threshold-dir=", "", tdir_raw_arg[1]))) else ""

load_s3n_backend(BASE_DIR)

# impute_censored_ros() is provided by scripts/s3n_transforms.R (sourced via s3n_plot_functions.R)
.wq_csv_cache <- new.env(parent = emptyenv())

find_first_col <- function(cols, candidates) {
  cols_low <- tolower(trimws(cols))
  for (cand in candidates) {
    hit <- which(cols_low == tolower(trimws(cand)))
    if (length(hit) > 0) return(cols[hit[1]])
  }
  NULL
}

read_wq_target_file <- function(f_path, source_label_tag, det_code, allow_substr_det = FALSE) {
  if (!file.exists(f_path)) return(NULL)
  cache_k <- normalizePath(f_path, winslash = "/", mustWork = FALSE)
  if (exists(cache_k, envir = .wq_csv_cache, inherits = FALSE)) {
    raw <- get(cache_k, envir = .wq_csv_cache, inherits = FALSE)
  } else {
    raw <- fast_read_csv(f_path, fileEncoding = "UTF-8-BOM")
    assign(cache_k, raw, envir = .wq_csv_cache)
  }
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
    keep_rows <- !is.na(det_vals) & ((det_vals == det_code) | (sub("^0+", "", det_vals) == sub("^0+", "", det_code)) |
      grepl(paste0("^", det_code, "$"), det_vals, ignore.case = TRUE))
    if (allow_substr_det) keep_rows <- keep_rows | (!is.na(det_vals) & grepl(det_code, det_vals, fixed = TRUE))
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

extract_target_samples <- function(tgt, trans_type_override = NULL) {
  ds_mode  <- tgt$dataset
  det_code <- tgt$code
  ds_info  <- DATASET_FILES[[ds_mode]]
  if (ds_mode == "wq") {
    s_dfs <- list()
    for (s_name in names(src_dirs$active_dirs)) {
      s_dir  <- src_dirs$active_dirs[[s_name]]
      f_path <- file.path(s_dir, paste0(toupper(s_name), " Water Quality Samples.csv"))
      if (!file.exists(f_path)) f_path <- file.path(s_dir, "samples.csv")
      if (!file.exists(f_path)) f_path <- file.path(s_dir, "Water Quality Samples.csv")
      if (!file.exists(f_path) && !is.null(custom_samples_cli) && file.exists(custom_samples_cli)) f_path <- custom_samples_cli
      sub_df <- read_wq_target_file(f_path, toupper(s_name), det_code, allow_substr_det = FALSE)
      if (!is.null(sub_df)) s_dfs[[s_name]] <- sub_df
    }
    if (length(s_dfs) == 0) stop(sprintf("No WQ samples found for '%s' in %s.", det_code, source_label))
    det_samples <- do.call(rbind, s_dfs)
    det_name <- det_samples$det_desc[1]
    det_unit <- det_samples$unit_short_desc[1]
    raw_vals <- suppressWarnings(as.numeric(det_samples$meas_result))
    is_cens  <- (det_samples$meas_sign %in% c("<", "yes", "true", "1")) & !grepl("tobit", det_samples$meas_method, ignore.case = TRUE)
  } else if (ds_mode == "lcms") {
    csv_p <- file.path(RSN_DIR, "RSN LCMS GCMS Samples.csv")
    all_samples <- fast_read_csv(csv_p, colClasses = c(cas_number = "character"), fileEncoding = "UTF-8-BOM")
    keep_m <- (!is.na(all_samples$cas_number) & all_samples$cas_number == det_code) |
              (!is.na(all_samples$compound_name) & grepl(det_code, all_samples$compound_name, ignore.case = TRUE))
    det_samples <- all_samples[keep_m, c("network_id", "sample_datetime", "compound_name", "less_than_flag", "meas_result", "unit_short_desc")]
    rm(all_samples)
    if (nrow(det_samples) == 0) stop(sprintf("No LCMS samples found for '%s'.", det_code))
    det_samples$data_source      <- "RSN"
    det_samples$meas_uncertainty <- NA_real_
    det_samples$meas_method      <- ""
    det_name <- det_samples$compound_name[1]
    det_unit <- det_samples$unit_short_desc[1]
    raw_vals <- suppressWarnings(as.numeric(det_samples$meas_result))
    is_cens  <- !is.na(det_samples$less_than_flag) & (tolower(trimws(det_samples$less_than_flag)) == "yes")
  } else {
    s_dfs <- list()
    for (s_name in names(src_dirs$active_dirs)) {
      s_dir  <- src_dirs$active_dirs[[s_name]]
      f_path <- file.path(s_dir, paste0(toupper(s_name), " ", ds_info$file))
      if (!file.exists(f_path)) f_path <- file.path(s_dir, ds_info$file)
      if (!file.exists(f_path)) f_path <- file.path(s_dir, "samples.csv")
      if (file.exists(f_path)) {
        raw <- fast_read_csv(f_path, fileEncoding = "UTF-8-BOM")
        fnames <- colnames(raw)
        id_c <- find_first_col(fnames, c("network_id", "site_id", "station_id", "id"))
        dt_c <- find_first_col(fnames, c(ds_info$dt_col, "sample_datetime", "datetime", "date_time", "date", "survey_date"))
        if (is.null(id_c)) id_c <- fnames[1]
        if (det_code %in% fnames) {
          s_dfs[[s_name]] <- data.frame(
            network_id       = as.character(raw[[id_c]]),
            sample_datetime  = if (!is.null(dt_c)) as.character(raw[[dt_c]]) else rep("2023-06-01", nrow(raw)),
            meas_result      = as.character(raw[[det_code]]),
            meas_uncertainty = rep(NA_real_, nrow(raw)),
            meas_method      = rep("", nrow(raw)),
            data_source      = toupper(s_name),
            stringsAsFactors = FALSE
          )
        }
        rm(raw)
      }
    }
    if (length(s_dfs) == 0) stop(sprintf("Column '%s' not found in %s.", det_code, ds_info$file))
    det_samples <- do.call(rbind, s_dfs)
    m_info   <- if (det_code %in% names(WIDE_VAR_META)) WIDE_VAR_META[[det_code]] else c(short = det_code, desc = det_code, unit = "index")
    det_name <- unname(m_info["short"])
    det_unit <- unname(m_info["unit"])
    raw_vals <- suppressWarnings(as.numeric(det_samples$meas_result))
    is_cens  <- rep(FALSE, length(raw_vals))
  }

  valid_k     <- !is.na(raw_vals) & !is.na(det_samples$sample_datetime) & nchar(trimws(det_samples$sample_datetime)) >= 7
  det_samples <- det_samples[valid_k, , drop = FALSE]
  raw_vals    <- raw_vals[valid_k]
  is_cens     <- is_cens[valid_k]

  det_samples$is_cens <- is_cens
  det_samples$raw_lod <- raw_vals
  det_samples$val <- if (any(is_cens) && tolower(censor_m) %in% c("ros", "em", "tobit")) {
    impute_censored_ros(raw_vals, is_cens)
  } else if (any(is_cens)) {
    ifelse(is_cens, 0.5 * raw_vals, raw_vals)
  } else {
    raw_vals
  }

  det_samples$year   <- suppressWarnings(as.integer(substr(det_samples$sample_datetime, 1, 4)))
  det_samples$month  <- suppressWarnings(as.integer(substr(det_samples$sample_datetime, 6, 7)))
  season_lu          <- c("winter", "winter", "spring", "spring", "spring",
                          "summer", "summer", "summer", "autumn", "autumn", "autumn", "winter")
  det_samples$season <- ifelse(!is.na(det_samples$month) & det_samples$month >= 1L & det_samples$month <= 12L,
                               season_lu[det_samples$month], NA_character_)

  if (!is.na(year_f))   det_samples <- det_samples[!is.na(det_samples$year)   & det_samples$year == year_f, , drop = FALSE]
  if (!is.na(season_f)) det_samples <- det_samples[!is.na(det_samples$season) & det_samples$season == season_f, , drop = FALSE]
  if (!is.na(month_f))  det_samples <- det_samples[!is.na(det_samples$month)  & det_samples$month == month_f, , drop = FALSE]
  if (nrow(det_samples) == 0L) {
    stop(sprintf("No remaining samples for target '%s' after temporal filtering.", det_code))
  }

  t_type <- if (!is.null(trans_type_override) && nzchar(trans_type_override)) {
    trans_type_override
  } else if (ds_mode %in% c("wq", "lcms")) {
    "log"
  } else {
    "none"
  }
  tr_tmp <- build_transformer(t_type)
  valid_t <- tr_tmp$valid_fn(det_samples$val)
  det_samples <- det_samples[valid_t, , drop = FALSE]
  tr_obj <- build_transformer(t_type, y_obs = det_samples$val)

  det_samples$y_sample <- tr_obj$fwd(det_samples$val)
  det_samples$y_sc     <- det_samples$y_sample
  label_short <- shorten_det_name(det_name)
  list(
    samples = det_samples,
    meta    = list(
      code    = det_code,
      dataset = ds_mode,
      slug    = tgt$slug,
      name    = det_name,
      label   = label_short,
      unit    = det_unit
    ),
    trans = tr_obj
  )
}

cat("======================================================================\n")
cat(sprintf("  Multivariate Stream Co-Kriging (%s, P = %d, K <= %d)\n", toupper(method_cli), P, k_factors))
cat("======================================================================\n")

sites_lu   <- fast_read_csv(SITES_CSV)
sites_lu   <- sites_lu[!is.na(sites_lu$included_250m) & sites_lu$included_250m == 1, c("site_id", "RivID", "obs_up_dist_km")]
df_reaches <- fast_read_csv(REACHES_CSV)
catch_lu   <- if (file.exists(CATCHMENTS_CSV)) fast_read_csv(CATCHMENTS_CSV) else NULL

extracted    <- vector("list", P)
targets_meta <- vector("list", P)
trans_list   <- vector("list", P)
reach_tables <- vector("list", P)

for (p in seq_len(P)) {
  t_ov <- if (length(trans_tokens) >= p) trans_tokens[p] else NULL
  res_p <- extract_target_samples(targets_parsed[[p]], trans_type_override = t_ov)
  extracted[[p]]    <- res_p
  targets_meta[[p]] <- res_p$meta
  trans_list[[p]]   <- res_p$trans

  merged_p <- merge(res_p$samples, sites_lu, by.x = "network_id", by.y = "site_id")
  if (nrow(merged_p) == 0) {
    stop(sprintf("No spatially snapped samples found for target '%s'.", res_p$meta$code))
  }
  if (length(prefer_sources_cli) > 0L) {
    pf_p <- filter_samples_by_source_priority(merged_p, prefer_sources_cli)
    merged_p <- pf_p$df
    if (pf_p$n_dropped > 0L) {
      cat(sprintf("    [Source Priority %s] Dropped %d lower-priority sample(s) on %d co-located reach(es).\n",
                  res_p$meta$code, pf_p$n_dropped, pf_p$n_reaches))
    }
  }
  agg_p <- aggregate_samples_by_reach_weighted(
    merged_p, res_p$trans,
    src_weights = src_weights_cli,
    src_se      = src_se_cli
  )
  v_sc_p <- rep(0.0, nrow(agg_p))
  ok_u_p <- !is.na(agg_p$meas_uncertainty) & is.finite(agg_p$meas_uncertainty) & agg_p$meas_uncertainty >= 0
  if (any(ok_u_p)) {
    u_raw_p <- agg_p$meas_uncertainty[ok_u_p]
    if ("meas_method" %in% names(agg_p)) {
      is_log_sd <- tolower(trimws(as.character(agg_p$meas_method[ok_u_p]))) %in% c("tobit_interval", "geometric_midpoint")
      u_raw_p[is_log_sd] <- u_raw_p[is_log_sd] * pmax(1e-6, agg_p$val[ok_u_p][is_log_sd])
    }
    v_sc_p[ok_u_p] <- calc_response_meas_variance(agg_p$val[ok_u_p], u_raw_p, res_p$trans)
  }
  reach_tables[[p]] <- data.frame(
    RivID          = agg_p$RivID,
    y_sc           = agg_p$y_resp,
    v_sc           = v_sc_p,
    source_weight  = agg_p$source_weight,
    obs_up_dist_km = agg_p$obs_up_dist_km,
    stringsAsFactors = FALSE
  )
  cat(sprintf("  [%d/%d] %-24s (%s:%s | trans=%s): %d samples -> %d reaches\n",
              p, P, res_p$meta$label, res_p$meta$dataset, res_p$meta$code,
              res_p$trans$type, nrow(merged_p), nrow(reach_tables[[p]])))
}

# Build union of observed RivIDs across all P targets (supporting heterotopic designs!)
all_obs_rivs <- sort(unique(unlist(lapply(reach_tables, function(tb) tb$RivID))))
obs_reach_rows <- match(all_obs_rivs, df_reaches$RivID)
keep_valid <- !is.na(obs_reach_rows)
all_obs_rivs   <- all_obs_rivs[keep_valid]
obs_reach_rows <- obs_reach_rows[keep_valid]

n_obs_union <- length(all_obs_rivs)
Y_obs_raw_sc  <- matrix(NA_real_, nrow = n_obs_union, ncol = P)
V_obs_phys_sc <- matrix(0.0,      nrow = n_obs_union, ncol = P)
W_obs_mat     <- matrix(1.0,      nrow = n_obs_union, ncol = P)
up_mat        <- matrix(NA_real_, nrow = n_obs_union, ncol = P)
for (p in seq_len(P)) {
  ix <- match(reach_tables[[p]]$RivID, all_obs_rivs)
  Y_obs_raw_sc[ix, p]  <- reach_tables[[p]]$y_sc
  V_obs_phys_sc[ix, p] <- reach_tables[[p]]$v_sc
  W_obs_mat[ix, p]     <- reach_tables[[p]]$source_weight
  up_mat[ix, p]        <- reach_tables[[p]]$obs_up_dist_km
}
obs_up_dist_obs <- apply(up_mat, 1L, function(r) {
  ok <- r[is.finite(r)]
  if (length(ok) > 0L) mean(ok) else NA_real_
})
miss_up <- !is.finite(obs_up_dist_obs)
if (any(miss_up)) obs_up_dist_obs[miss_up] <- df_reaches$up_dist_mid_km[obs_reach_rows[miss_up]]

cat(sprintf("  Union spatial support: %d unique reaches across P = %d targets.\n", n_obs_union, P))

# Resolve exceedance thresholds per target
thresholds_vec     <- rep(NA_real_, P)
threshold_dirs_vec <- rep("above", P)
if (nzchar(thresh_raw_val)) {
  th_toks <- trimws(unlist(strsplit(thresh_raw_val, ",", fixed = TRUE)))
  if (length(th_toks) == 1L) th_toks <- rep(th_toks, P)
  if (length(th_toks) != P) stop("--threshold must be a single value ('uk' or numeric) or match P targets.")
  td_toks <- if (nzchar(tdir_raw_val)) trimws(unlist(strsplit(tdir_raw_val, ",", fixed = TRUE))) else character(0)
  if (length(td_toks) == 1L) td_toks <- rep(td_toks, P)
  for (p in seq_len(P)) {
    dir_in <- if (length(td_toks) >= p) td_toks[p] else ""
    res_t  <- resolve_threshold_spec(th_toks[p], if (nzchar(dir_in)) dir_in else "above", targets_meta[[p]]$code, args)
    thresholds_vec[p]     <- res_t$val
    threshold_dirs_vec[p] <- res_t$dir
  }
} else {
  thresholds_vec <- NULL
}

# Build output directory
slug_joined <- paste(vapply(targets_meta, function(m) gsub("[^a-zA-Z0-9_]+", "_", m$code), character(1)), collapse = "_")
if (nchar(slug_joined) > 64L) slug_joined <- paste0(substr(slug_joined, 1L, 52L), "_p", P)
out_folder <- paste0("cokrig_", slug_joined, "_", net_meta$abbr)
if (nzchar(source_tag)) out_folder <- paste0(out_folder, "_", source_tag)
if (method_cli != "lmc") out_folder <- paste0(out_folder, "_", method_cli)
if (k_factors < P && method_cli == "lmc") out_folder <- paste0(out_folder, "_k", k_factors)
OUT_DIR <- rel_path(OUT_BASE, "results", out_folder)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# resolve_basin_spec() is shared via scripts/s3n_covariates.R
humber_info <- resolve_core_basin(df_reaches, catch_lu, net_meta, "humber")
thames_info <- resolve_core_basin(df_reaches, catch_lu, net_meta, "thames")
severn_info <- resolve_core_basin(df_reaches, catch_lu, net_meta, "severn")

domains <- list(
  list(slug = "england", domain_label = sprintf("All England %s Network", source_abbr), catch_ids = NULL),
  list(slug = "humber",  domain_label = humber_info$domain_label, catch_ids = humber_info$catch_id),
  list(slug = "thames",  domain_label = thames_info$domain_label, catch_ids = thames_info$catch_id),
  list(slug = "severn",  domain_label = severn_info$domain_label, catch_ids = severn_info$catch_id)
)
for (b_sp in basin_specs) {
  b_inf <- resolve_basin_spec(b_sp, catch_lu)
  domains[[length(domains) + 1L]] <- list(slug = b_inf$slug, domain_label = b_inf$domain_label, catch_ids = b_inf$catch_ids)
}

all_cv_list     <- list()
all_cross_list  <- list()
all_cov_list    <- list()
all_exceed_list <- list()
var_names       <- vapply(targets_meta, function(m) m$label, character(1))
vars_sub        <- format_covariates_subtitle(var_specs = var_specs)
method_lbl      <- if (method_cli == "colocated") "Sequential Co-Located Residual Co-Kriging" else sprintf("Factor-LMC Joint NNGP (K=%d)", k_factors)

eng_pred_ref <- NULL
eng_norm_stats_list <- vector("list", P)
for (dom in domains) {
  res_dom <- fit_multivariate_s3n_domain(
    df_reaches      = df_reaches,
    Y_obs_raw_sc    = Y_obs_raw_sc,
    obs_up_dist_obs = obs_up_dist_obs,
    obs_reach_rows  = obs_reach_rows,
    targets_meta    = targets_meta,
    trans_list      = trans_list,
    var_specs       = var_specs,
    domain_id       = dom$slug,
    domain_label    = dom$domain_label,
    catch_filter    = dom$catch_ids,
    m_neighbors     = m_neighbors,
    k_factors       = k_factors,
    kernel          = kernel_cli,
    method          = method_cli,
    cv_mode         = cv_mode_cli,
    thresholds      = thresholds_vec,
    threshold_dirs  = threshold_dirs_vec,
    n_omp           = n_omp_cli,
    V_obs_phys_sc   = V_obs_phys_sc,
    W_obs_mat       = W_obs_mat
  )
  if (is.null(res_dom)) next

  all_cv_list[[length(all_cv_list) + 1L]]       <- res_dom$cv_compare_df
  all_cross_list[[length(all_cross_list) + 1L]] <- res_dom$cross_cov_df
  if (!is.null(res_dom$coverage_df)) {
    all_cov_list[[length(all_cov_list) + 1L]] <- res_dom$coverage_df
  }
  if (!is.null(res_dom$exceed_summary)) {
    all_exceed_list[[length(all_exceed_list) + 1L]] <- res_dom$exceed_summary
  }

  for (p in seq_len(P)) {
    m_p <- targets_meta[[p]]
    c_slug <- gsub("[^a-zA-Z0-9_]+", "_", m_p$code)
    if (dom$slug == "england") {
      eng_norm_stats_list[[p]] <- compute_norm_ref_stats(res_dom$pred_dfs[[p]])
      res_dom$pred_dfs[[p]] <- add_normalised_columns(res_dom$pred_dfs[[p]], ref_df = res_dom$pred_dfs[[p]], ref_stats = eng_norm_stats_list[[p]])
    } else {
      ref_p <- if (is.null(eng_pred_ref)) res_dom$pred_dfs[[p]] else eng_pred_ref[[p]]
      res_dom$pred_dfs[[p]] <- add_normalised_columns(res_dom$pred_dfs[[p]], ref_df = ref_p, ref_stats = eng_norm_stats_list[[p]])
    }
    out_csv <- file.path(OUT_DIR, sprintf("predictions_%s_%s.csv", dom$slug, c_slug))
    fast_write_csv(res_dom$pred_dfs[[p]], out_csv)

    if (dom$slug == "england") {
      ftag <- norm_meta$file_tag
      fig_map <- file.path(OUT_DIR, sprintf("fig_cokrig_england_%s%s.png", c_slug, ftag))
      plot_fig1_england_map(
        res_dom$pred_dfs[[p]], m_p$code, paste0(m_p$name, " [Co-Kriging]"), m_p$unit, trans_list[[p]], fig_map,
        reverse_color = reverse_color, vars_label = vars_sub,
        network_label = net_meta$short_label, source_label = source_label, norm_mode = norm_cli
      )
      if (!is.null(thresholds_vec) && is.finite(thresholds_vec[p])) {
        fig_ex <- file.path(OUT_DIR, sprintf("fig_exceedance_england_%s.png", c_slug))
        plot_fig_exceedance_england(
          res_dom$pred_dfs[[p]], m_p$code, paste0(m_p$name, " [Co-Kriging]"), m_p$unit, fig_ex,
          threshold_val = thresholds_vec[p], threshold_dir = threshold_dirs_vec[p],
          vars_label = vars_sub, network_label = net_meta$short_label, source_label = source_label
        )
      }
    }
  }

  if (dom$slug == "england") {
    eng_pred_ref <- res_dom$pred_dfs
    fig_diag <- file.path(OUT_DIR, "fig_cokrig_cross_correlation.png")
    plot_fig_cokrig_diagnostics(
      cross_cov_df  = res_dom$cross_cov_df,
      cv_compare_df = res_dom$cv_compare_df,
      var_names     = var_names,
      out_file      = fig_diag,
      network_label = net_meta$short_label,
      vars_label    = vars_sub,
      source_label  = source_label,
      method_label  = method_lbl
    )
  }
  rm(res_dom); gc(FALSE)
}

cv_all  <- if (length(all_cv_list) > 0L) do.call(rbind, all_cv_list) else NULL
cc_all  <- if (length(all_cross_list) > 0L) do.call(rbind, all_cross_list) else NULL
cov_all <- if (length(all_cov_list) > 0L) do.call(rbind, all_cov_list) else NULL
ex_all  <- if (length(all_exceed_list) > 0L) do.call(rbind, all_exceed_list) else NULL

if (!is.null(cv_all)) fast_write_csv(cv_all, file.path(OUT_DIR, "cv_univariate_vs_cokriging.csv"))
if (!is.null(cc_all)) fast_write_csv(cc_all, file.path(OUT_DIR, "cokriging_cross_covariance.csv"))
if (!is.null(ex_all)) fast_write_csv(ex_all, file.path(OUT_DIR, "exceedance_summary_cokriging.csv"))

report_path <- file.path(OUT_DIR, "cokriging_evaluation_report.md")
generate_cokriging_report(
  cv_all       = cv_all,
  cc_all       = cc_all,
  coverage_df  = cov_all,
  exceed_df    = ex_all,
  out_file     = report_path,
  method_label = method_lbl,
  cv_mode      = cv_mode_cli
)
cat(sprintf("  + Saved Co-Kriging value-add report: %s\n", report_path))

cat(sprintf("\nMultivariate Co-Kriging completed! Outputs saved to %s/\n", OUT_DIR))
