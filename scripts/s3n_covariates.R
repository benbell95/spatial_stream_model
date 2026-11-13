# ==============================================================================
# s3n_covariates.R
#
# Shared fixed-effect regression covariate catalog, design matrix builder,
# --vars=... parser helpers, and core basin resolution.
# Automatically sourced by scripts/s3n_plot_functions.R.
# ==============================================================================

COVARIATE_CATALOG <- list(
  Elevation_z = list(
    aliases = c("elevation", "elev", "elevation_z"),
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

apply_covariate_specs <- function(specs, base_terms) {
  is_mod <- grepl("^[+-]", specs)
  if (all(is_mod)) {
    terms <- base_terms
  } else {
    terms <- character(0)
    for (sp in specs[!is_mod]) {
      if (tolower(sp) %in% c("intercept", "none")) next
      nm <- match_covariate_token(sp)
      if (!(nm %in% terms)) terms <- c(terms, nm)
    }
    specs <- specs[is_mod]
  }
  for (sp in specs) {
    op <- substr(sp, 1, 1)
    nm <- match_covariate_token(substr(sp, 2, nchar(sp)))
    if (op == "+") {
      if (nm == "Is_Aquifer") {
        terms <- setdiff(terms, c("Is_Chalk_Limestone", "Is_Sandstone"))
      } else if (nm %in% c("Is_Chalk_Limestone", "Is_Sandstone")) {
        terms <- setdiff(terms, "Is_Aquifer")
      }
      if (!(nm %in% terms)) terms <- c(terms, nm)
    } else {
      if (nm == "Is_Aquifer") {
        terms <- setdiff(terms, c("Is_Aquifer", "Is_Chalk_Limestone", "Is_Sandstone"))
      } else {
        terms <- setdiff(terms, nm)
      }
    }
  }
  terms
}

resolve_covariate_terms <- function(df, include_peat = TRUE, var_specs = NULL) {
  sd_peat <- if ("Is_Peat" %in% colnames(df)) suppressWarnings(sd(df$Is_Peat[df$has_obs == 1], na.rm = TRUE)) else NA_real_
  sd_sand <- if ("Is_Sandstone" %in% colnames(df)) suppressWarnings(sd(df$Is_Sandstone[df$has_obs == 1], na.rm = TRUE)) else NA_real_
  has_split_geol <- include_peat && is.finite(sd_peat) && sd_peat > 0 && is.finite(sd_sand) && sd_sand > 0

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

  apply_covariate_specs(specs, default_terms)
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
    if (any(!is.finite(vec))) {
      fin <- vec[is.finite(vec)]
      fill_v <- if (length(fin) > 0) mean(fin) else 0.0
      vec[!is.finite(vec)] <- fill_v
    }
    sd_obs <- suppressWarnings(sd(vec[df$has_obs == 1], na.rm = TRUE))
    if (is.finite(sd_obs) && sd_obs > 0) {
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

format_covariates_subtitle <- function(var_specs = NULL, param_tab = NULL) {
  default_eng <- c("Elevation_z", "Slope_z", "log10_US_Accum_z",
                   "Is_Chalk_Limestone", "Is_Sandstone", "Is_Peat",
                   "Easting_z", "Northing_z")
  default_bas <- c("Elevation_z", "Slope_z", "log10_US_Accum_z",
                   "Is_Aquifer", "Easting_z", "Northing_z")

  has_cli_vars <- !is.null(var_specs) && length(var_specs) > 0 &&
                  !(length(var_specs) == 1 && tolower(trimws(var_specs[1])) == "default")

  if (!is.null(param_tab) && nrow(param_tab) > 0) {
    ref_dom <- if (any(grepl("All England", param_tab$Domain))) {
      param_tab$Domain[grepl("All England", param_tab$Domain)][1]
    } else {
      param_tab$Domain[1]
    }
    eng_p <- unique(as.character(param_tab$Parameter[
      param_tab$Domain == ref_dom &
      !grepl("^(Intercept|sigma\\.sq|tau\\.sq|phi|mean_meas_variance|gam_edf)", param_tab$Parameter)
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
  terms <- apply_covariate_specs(specs, default_eng)
  if (length(terms) == 0) return("None (Intercept-only spatial model)")
  paste(unname(COVARIATE_SHORT_LABEL[terms]), collapse = ", ")
}

resolve_core_basin <- function(df_reaches, catch_lookup = NULL, net_meta = NULL, target = c("humber", "thames", "severn")) {
  target <- match.arg(target)
  avail_cids <- unique(df_reaches$CatchID)

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

  tab_c <- sort(table(df_reaches$CatchID), decreasing = TRUE)
  cid <- as.integer(names(tab_c)[1])
  lbl <- sprintf("%s Basin (CatchID=%d)", paste0(toupper(substr(target, 1, 1)), substr(target, 2, nchar(target))), cid)
  list(catch_id = cid, domain_label = lbl)
}
